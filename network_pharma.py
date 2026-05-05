# -*- coding: utf-8 -*-
"""
网络药理学模块 - SMILES查询 / SwissTargetPrediction靶点预测 / gseapy富集分析
依赖: gseapy, pandas, numpy; PubChem API调用使用标准库 urllib（无需 requests）
"""

import time
import re
import json
from urllib.request import urlopen, Request
from urllib.parse import quote
from urllib.error import URLError, HTTPError
import pandas as pd
import numpy as np

try:
    import streamlit as st
    _ST_AVAILABLE = True
except ImportError:
    _ST_AVAILABLE = False
    st = None

# ===================== 常量 =====================
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_ST_BASE = "https://www.swisstargetprediction.ch/api"
_REQUEST_DELAY = 1.2  # SwissTargetPrediction 请求间隔（秒）
_MAX_RETRIES = 3     # SwissTargetPrediction 重试次数
_PROB_THRESHOLD = 0.01  # 靶点预测概率阈值

# ===================== 通用 HTTP 工具（urllib，不依赖 requests）=====================

def _http_get_json(url: str, timeout: int = 15) -> dict:
    """
    用 urllib GET JSON，替代 requests.get().json()。
    兼容 Python 内置库，不依赖 requests。
    """
    try:
        req = Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (URLError, HTTPError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return {"_error": str(e)}


def _http_post_json(url: str, payload: dict, timeout: int = 30) -> dict:
    """用 urllib POST JSON"""
    try:
        data = json.dumps(payload).encode("utf-8")
        req = Request(url, data=data, headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "Mozilla/5.0",
        })
        with urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (URLError, HTTPError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return {"_error": str(e)}

# ===================== 缓存 =====================

# 简单的内存缓存（替代 st.cache_data，兼容非 Streamlit 环境）
_smiles_cache = {}
_cache_max_size = 2000


def _smiles_cached(key: str, fn, *args):
    if key in _smiles_cache:
        return _smiles_cache[key]
    result = fn(*args)
    if len(_smiles_cache) < _cache_max_size:
        _smiles_cache[key] = result
    return result

# ===================== SMILES 查询 =====================

def _parse_first_name(name: str) -> str:
    """从 compound_name 字段中提取第一个合法名称（处理分号、空格等）"""
    if not name:
        return ''
    name = str(name).strip()
    # 取分号/逗号前第一段
    for sep in (';', ',', '/'):
        if sep in name:
            name = name.split(sep)[0].strip()
    return name


def _name_to_query(name: str) -> str:
    """
    将用户输入的化合物名称转换为适合 PubChem 查询的字符串。
    去掉括号内容、-等常见格式问题。
    """
    if not name:
        return ''
    name = str(name).strip()
    # 去掉括号及其内容，如 "(2R)-Naringenin" -> "Naringenin"
    name = re.sub(r'\([^)]*\)', '', name)
    # 去掉结尾的连字符和多余空格
    name = re.sub(r'\s+-\s*$', '', name).strip()
    # 多个空格变一个
    name = re.sub(r'\s+', ' ', name)
    return name


@st.cache_data(ttl=3600)
def query_smiles_by_name(compound_name: str) -> str:
    """
    通过 PubChem REST API 按化合物名称查询 SMILES。
    使用 urllib（Python 内置），不依赖 requests。

    Args:
        compound_name: 原始化合物名称

    Returns:
        SMILES 字符串，查不到返回 "NOT_FOUND"
    """
    if not compound_name or str(compound_name).strip() in ('-', '', 'nan', 'NaN'):
        return "NOT_FOUND"

    # 解析出第一个合法名称
    first_name = _parse_first_name(compound_name)
    if not first_name:
        return "NOT_FOUND"

    # 查询用名称列表（按优先级）
    names_to_try = []
    q_name = _name_to_query(first_name)
    if q_name and q_name != first_name:
        names_to_try.append(q_name)
    names_to_try.append(first_name)

    for name in names_to_try:
        if not name:
            continue

        # 1. 精确查询
        url = f"{_PUBCHEM_BASE}/compound/name/{quote(name)}/property/IsomericSMILES/JSON"
        result = _http_get_json(url)
        if "_error" not in result:
            try:
                smiles = result['PropertyTable']['Properties'][0].get('IsomericSMILES', '')
                if smiles:
                    return smiles
            except (KeyError, IndexError, TypeError):
                pass

        # 2. 模糊匹配
        url_fuzzy = f"{_PUBCHEM_BASE}/compound/name/{quote(name)}/property/IsomericSMILES/JSON?algorithm=fuzzy"
        result = _http_get_json(url_fuzzy)
        if "_error" not in result:
            try:
                smiles = result['PropertyTable']['Properties'][0].get('IsomericSMILES', '')
                if smiles:
                    return smiles
            except (KeyError, IndexError, TypeError):
                pass

        # 3. CAS 号查询
        if re.match(r'^\d+-\d+-\d+$', name):
            cas_url = f"{_PUBCHEM_BASE}/compound/cas/{name}/property/IsomericSMILES/JSON"
            result = _http_get_json(cas_url)
            if "_error" not in result:
                try:
                    smiles = result['PropertyTable']['Properties'][0].get('IsomericSMILES', '')
                    if smiles:
                        return smiles
                except (KeyError, IndexError, TypeError):
                    pass

    return "NOT_FOUND"


def query_smiles_by_metabolite(metabolite_name: str) -> str:
    """
    直接用代谢物名称（Metabolite列）查询 PubChem SMILES。
    """
    if not metabolite_name or str(metabolite_name).strip() in ('-', '', 'nan'):
        return "NOT_FOUND"

    name = _name_to_query(str(metabolite_name).strip())
    if not name:
        return "NOT_FOUND"

    names_to_try = [name, str(metabolite_name).strip()]

    for query_name in names_to_try:
        # 精确
        url = f"{_PUBCHEM_BASE}/compound/name/{quote(query_name)}/property/IsomericSMILES/JSON"
        result = _http_get_json(url)
        if "_error" not in result:
            try:
                smiles = result['PropertyTable']['Properties'][0].get('IsomericSMILES', '')
                if smiles:
                    return smiles
            except (KeyError, IndexError, TypeError):
                pass

        # 模糊
        url_fuzzy = f"{_PUBCHEM_BASE}/compound/name/{quote(query_name)}/property/IsomericSMILES/JSON?algorithm=fuzzy"
        result = _http_get_json(url_fuzzy)
        if "_error" not in result:
            try:
                smiles = result['PropertyTable']['Properties'][0].get('IsomericSMILES', '')
                if smiles:
                    return smiles
            except (KeyError, IndexError, TypeError):
                pass

    return "NOT_FOUND"


def batch_query_smiles(df: pd.DataFrame, name_col: str = 'compound_name',
                       progress_callback=None) -> pd.DataFrame:
    """
    批量查询 SMILES，支持双列 fallback：
    优先用 compound_name 查询，若 NOT_FOUND 再用 Metabolite 列名称查询。

    Args:
        df: 输入 DataFrame（需含 Metabolite 和 name_col 列）
        name_col: 化合物名列名（默认 compound_name）
        progress_callback: 回调函数，接收 (current, total)

    Returns:
        新增 SMILES 和 SMILES_Source 列的 DataFrame
    """
    result_df = df.copy()
    if 'SMILES' not in result_df.columns:
        result_df['SMILES'] = 'PENDING'
    if 'SMILES_Source' not in result_df.columns:
        result_df['SMILES_Source'] = ''

    invalid_vals = {None, '', '-', 'nan', 'NaN', 'Na', 'NULL'}
    n_total = len(result_df)

    for i, (_, row) in enumerate(result_df.iterrows()):
        smiles_col = row.get('SMILES', '')
        # 已有有效 SMILES 且非占位符则跳过
        if smiles_col and smiles_col not in ('PENDING', 'NOT_IN_DATA', 'NO_NAME', 'NOT_FOUND'):
            if progress_callback:
                progress_callback(i + 1, n_total)
            continue

        metab = str(row.get('Metabolite', '')).strip()
        name = str(row.get(name_col, '')).strip()

        smiles = 'NOT_FOUND'
        source = ''

        # 策略1：compound_name 有效 → 优先查它
        if name and name not in invalid_vals and not pd.isna(row.get(name_col)):
            smiles = query_smiles_by_name(name)
            source = 'compound_name'
            time.sleep(0.3)

        # 策略2：compound_name 无效 或 NOT_FOUND → 用 Metabolite 名称
        if smiles == 'NOT_FOUND' and metab and metab not in invalid_vals:
            smiles = query_smiles_by_metabolite(metab)
            source = 'Metabolite'
            time.sleep(0.3)

        if not smiles or smiles == 'NOT_FOUND':
            # 两者都查不到 → 标记
            if not name or name in invalid_vals:
                smiles = 'NO_NAME'
                source = 'no_input'
            else:
                smiles = 'NOT_FOUND'
                source = source or 'not_found'

        result_df.at[_, 'SMILES'] = smiles
        result_df.at[_, 'SMILES_Source'] = source
        if progress_callback:
            progress_callback(i + 1, n_total)

    return result_df


# ===================== SwissTargetPrediction 靶点预测 =====================

def query_swiss_target_prediction(smiles: str, species: str = "Homo sapiens",
                                  max_retries: int = _MAX_RETRIES) -> tuple:
    """
    调用 SwissTargetPrediction API 获取靶点预测结果。

    Args:
        smiles: 化合物的 SMILES
        species: 物种（默认 "Homo sapiens"）
        max_retries: 429错误重试次数

    Returns:
        (genes_str, target_count, error_msg)
        genes_str: 分号分隔的基因名列表（prob > 0.01），查不到返回 "NOT_FOUND" 或 "ERROR"
        target_count: 靶点数量
        error_msg: 错误信息（无错误为 ""）
    """
    if not smiles or smiles in ('NOT_FOUND', 'PENDING', '', 'nan'):
        return "NOT_FOUND", 0, "No SMILES provided"

    url = f"{_ST_BASE}/search"
    payload = {"smiles": smiles, "species": species}

    for attempt in range(max_retries + 1):
        result = _http_post_json(url, payload, timeout=30)
        err = result.get("_error", "")

        if not err:
            try:
                results = result.get('Results', [])
                filtered = [r for r in results if float(r.get('Probability', 0)) > _PROB_THRESHOLD]
                genes = []
                for r in filtered:
                    for key in ('Gene', 'Target', 'Symbol'):
                        g = r.get(key, '')
                        if g:
                            genes.append(str(g).strip())
                            break
                genes = sorted(set(g for g in genes if g))
                genes_str = ";".join(genes)
                return genes_str, len(genes), ""
            except (KeyError, ValueError, TypeError) as e:
                return "ERROR", 0, f"Parse error: {e}"

        # 处理错误
        if "429" in err or "429" in str(result):
            wait_time = 2 ** attempt + 1
            time.sleep(wait_time)
            continue
        elif "timeout" in err.lower() or "timed out" in err.lower():
            if attempt < max_retries:
                time.sleep(2)
                continue
            return "ERROR", 0, "Timeout"
        else:
            return "ERROR", 0, err

    return "ERROR", 0, "Max retries exceeded"


def batch_swiss_target_prediction(df: pd.DataFrame, smiles_col: str = 'SMILES',
                                   progress_callback=None,
                                   delay: float = _REQUEST_DELAY) -> tuple:
    """
    批量调用 SwissTargetPrediction API。

    Args:
        df: 输入 DataFrame（需含 SMILES 列）
        smiles_col: SMILES 列名
        progress_callback: 回调函数，接收 (current, total)
        delay: 请求间隔（秒）

    Returns:
        (result_df, errors_list)
    """
    result_df = df.copy()
    if 'Predicted_Targets' not in result_df.columns:
        result_df['Predicted_Targets'] = 'PENDING'
    if 'Target_Count' not in result_df.columns:
        result_df['Target_Count'] = 0

    n = len(result_df)
    errors = []
    for i, (_, row) in enumerate(result_df.iterrows()):
        smiles = str(row.get(smiles_col, ''))
        genes_str, count, err = query_swiss_target_prediction(smiles)
        result_df.at[_, 'Predicted_Targets'] = genes_str
        result_df.at[_, 'Target_Count'] = count
        if err:
            errors.append(f"Row {i}: {err}")
        if progress_callback:
            progress_callback(i + 1, n)
        if i < n - 1:
            time.sleep(delay)

    return result_df, errors


# ===================== gseapy 富集分析 =====================

def run_gseapy_enrichment(gene_list: list, organism: str = "human",
                          databases: list = None) -> dict:
    """
    使用 gseapy 做 KEGG + GO 富集分析。

    Args:
        gene_list: 基因名列表（来自 Predicted_Targets 列的分号分隔基因）
        organism: "human"
        databases: 要分析的数据库列表，默认包含 KEGG + 三个 GO 子库

    Returns:
        dict: {db_name: pd.DataFrame} 富集结果 DataFrame 字典
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError("请安装 gseapy: pip install gseapy")

    if databases is None:
        databases = [
            "KEGG_2021_Human",
            "GO_Biological_Process_2021",
            "GO_Molecular_Function_2021",
            "GO_Cellular_Component_2021",
        ]

    results = {}
    for db in databases:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=db,
                organism=organism,
                outdir=None,
                no_plot=True,
            )
            if enr.results is not None and len(enr.results) > 0:
                res_df = enr.results.copy()
                res_df['Database'] = db
                # 标准化列名
                col_map = {}
                for c in res_df.columns:
                    cl = c.lower()
                    if 'term' in cl or 'pathway' in cl:
                        col_map[c] = 'Term'
                    elif 'pvalue' in cl or 'p-value' in cl:
                        col_map[c] = 'P_value'
                    elif 'adjusted' in cl or 'bh' in cl:
                        col_map[c] = 'Adjusted_P_value'
                    elif 'genes' in cl:
                        col_map[c] = 'Genes'
                    elif 'overlap' in cl or 'hits' in cl:
                        col_map[c] = 'Overlap'
                res_df = res_df.rename(columns=col_map)
                results[db] = res_df
        except Exception:
            results[db] = pd.DataFrame()
    return results


def merge_enrichment_results(enr_dict: dict) -> pd.DataFrame:
    """
    合并多个数据库的富集结果。
    """
    dfs = []
    for db, df in enr_dict.items():
        if df is not None and len(df) > 0:
            df = df.copy()
            if 'Database' not in df.columns:
                df['Database'] = db
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    merged = pd.concat(dfs, ignore_index=True)
    # 确保必要列存在
    for col in ['Term', 'Adjusted_P_value', 'Genes', 'Overlap']:
        if col not in merged.columns:
            merged[col] = None
    merged = merged.sort_values('Adjusted_P_value').reset_index(drop=True)
    return merged


# ===================== 可视化 =====================

def plot_enrichment_dotplot(enr_df: pd.DataFrame, top_n: int = 30) -> 'go.Figure':
    """
    绘制 Enrichr 风格的 dotplot（气泡图）。
    X轴 = Count/Overlap，Y轴 = Term，按 P-value 颜色编码。

    Args:
        enr_df: 富集结果 DataFrame，需含 Term, Adjusted_P_value, Genes, Overlap 列
        top_n: 显示前 N 条

    Returns:
        plotly.graph_objects.Figure
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("请安装 plotly: pip install plotly")

    if enr_df is None or len(enr_df) == 0:
        return None

    plot_df = enr_df.head(top_n).copy()

    # 解析 Overlap 列为数值（格式如 "15/200" -> 15）
    def parse_overlap(val):
        if pd.isna(val):
            return 0
        s = str(val)
        if '/' in s:
            try:
                return float(s.split('/')[0])
            except:
                return 0
        try:
            return float(val)
        except:
            return 0

    # 解析 Genes 列获取基因数量
    def gene_count(val):
        if pd.isna(val):
            return 0
        return len([g for g in str(val).split(';') if g.strip()])

    plot_df['Count'] = plot_df.apply(lambda r: gene_count(r.get('Genes', '')), axis=1)
    plot_df['X_val'] = plot_df.apply(lambda r: parse_overlap(r.get('Overlap', 0)), axis=1)

    # Term 名称截断
    plot_df['Term_short'] = plot_df['Term'].apply(lambda x: x[:50] + '...' if len(str(x)) > 50 else str(x))

    # 取 -log10(P) 用于颜色
    p_col = 'Adjusted_P_value' if 'Adjusted_P_value' in plot_df.columns else 'P_value'
    pvals = pd.to_numeric(plot_df[p_col], errors='coerce').fillna(1)
    pvals = pvals.replace(0, 1e-10)
    plot_df['-log10P'] = -np.log10(pvals)

    fig = go.Figure(go.Scatter(
        x=plot_df['X_val'],
        y=plot_df['Term_short'],
        mode='markers',
        marker=dict(
            size=plot_df['Count'] * 3 + 5,
            color=plot_df['-log10P'],
            colorscale='RdYlBu_r',
            showscale=True,
            colorbar=dict(title='-log10(Padj)'),
            sizemin=4,
        ),
        text=plot_df.apply(
            lambda r: f"<b>{r['Term']}</b><br>Count: {r['Count']}<br>AdjP: {r.get('Adjusted_P_value', r.get('P_value', 'N/A'))}<br>Genes: {r.get('Genes', '')}",
            axis=1
        ),
        hoverinfo='text',
    ))
    fig.update_layout(
        title=f'Enrichment Dotplot (Top {top_n} by Adjusted P-value)',
        xaxis_title='Count / Overlap',
        yaxis_title='Term',
        height=max(400, len(plot_df) * 22),
        width=900,
        font=dict(size=10),
        yaxis=dict(autorange='reversed'),
    )
    return fig