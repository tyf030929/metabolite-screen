# -*- coding: utf-8 -*-
"""
网络药理学模块 - SMILES查询 / SwissTargetPrediction靶点预测 / gseapy富集分析
依赖: requests, gseapy, pandas, numpy, streamlit (st.cache_data)
"""

import time
import re
import requests
import pandas as pd
import numpy as np
import streamlit as st

# ===================== 常量 =====================
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_ST_BASE = "https://www.swisstargetprediction.ch/api"
_REQUEST_DELAY = 1.2  # SwissTargetPrediction 请求间隔（秒）
_MAX_RETRIES = 3     # SwissTargetPrediction 重试次数
_PROB_THRESHOLD = 0.01  # 靶点预测概率阈值

# ===================== SMILES 查询 =====================

@st.cache_data(ttl=3600)
def query_smiles_by_name(compound_name: str) -> str:
    """
    通过 PubChem REST API 按化合物名称查询 SMILES。
    使用 @st.cache_data 缓存，避免重复查询。

    Args:
        compound_name: 化合物名称

    Returns:
        SMILES 字符串，查不到返回 "NOT_FOUND"
    """
    if not compound_name or str(compound_name).strip() in ('-', '', 'nan', 'NaN'):
        return "NOT_FOUND"

    name = str(compound_name).strip()

    # 1. 精确查询
    url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES/JSON"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            smiles = data['PropertyTable']['Properties'][0].get('IsomericSMILES', 'NOT_FOUND')
            return smiles if smiles else "NOT_FOUND"
    except Exception:
        pass

    # 2. 模糊匹配
    fuzzy_url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES/JSON?algorithm=fuzzy"
    try:
        resp = requests.get(fuzzy_url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            smiles = data['PropertyTable']['Properties'][0].get('IsomericSMILES', 'NOT_FOUND')
            return smiles if smiles else "NOT_FOUND"
    except Exception:
        pass

    # 3. CAS 号查询
    if re.match(r'^\d+-\d+-\d+$', name):
        cas_url = f"{_PUBCHEM_BASE}/compound/cas/{name}/property/IsomericSMILES/JSON"
        try:
            resp = requests.get(cas_url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                smiles = data['PropertyTable']['Properties'][0].get('IsomericSMILES', 'NOT_FOUND')
                return smiles if smiles else "NOT_FOUND"
        except Exception:
            pass

    return "NOT_FOUND"


def batch_query_smiles(df: pd.DataFrame, name_col: str = 'compound_name',
                       progress_callback=None) -> pd.DataFrame:
    """
    批量查询 SMILES。

    Args:
        df: 输入 DataFrame
        name_col: 化合物名列名
        progress_callback: 回调函数，接收 (current, total)

    Returns:
        新增 SMILES 列的 DataFrame
    """
    result_df = df.copy()
    if 'SMILES' not in result_df.columns:
        result_df['SMILES'] = 'NOT_IN_DATA'

    # 预过滤：跳过无效名称
    invalid_vals = {None, '', '-', 'nan', 'NaN', 'Na', 'NULL'}
    n_total = len(result_df)

    for i, (_, row) in enumerate(result_df.iterrows()):
        name = row.get(name_col, None)
        # 跳过无效名称
        if str(name).strip() in invalid_vals or pd.isna(name):
            result_df.at[_, 'SMILES'] = 'NO_NAME'
            if progress_callback:
                progress_callback(i + 1, n_total)
            continue

        # 已有有效结果则跳过
        cur = result_df.at[_, 'SMILES']
        if cur not in ('PENDING', 'NOT_IN_DATA', 'NO_NAME'):
            if progress_callback:
                progress_callback(i + 1, n_total)
            continue

        smiles = query_smiles_by_name(str(name).strip())
        result_df.at[_, 'SMILES'] = smiles
        if progress_callback:
            progress_callback(i + 1, n_total)
        time.sleep(0.3)  # 避免对 PubChem 请求过快

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
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(max_retries + 1):
        try:
            resp = requests.post(url, json=payload, headers=headers, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get('Results', [])
                filtered = [r for r in results if float(r.get('Probability', 0)) > _PROB_THRESHOLD]
                genes = []
                for r in filtered:
                    # SwissTargetPrediction 返回字段: Target, Gene, Probability
                    for key in ('Gene', 'Target', 'Symbol'):
                        g = r.get(key, '')
                        if g:
                            genes.append(str(g).strip())
                            break
                genes = sorted(set(g for g in genes if g))
                genes_str = ";".join(genes)
                return genes_str, len(genes), ""
            elif resp.status_code == 429:
                wait_time = 2 ** attempt + 1
                time.sleep(wait_time)
                continue
            else:
                return "ERROR", 0, f"HTTP {resp.status_code}"
        except requests.exceptions.Timeout:
            if attempt < max_retries:
                time.sleep(2)
                continue
            return "ERROR", 0, "Timeout"
        except Exception as e:
            return "ERROR", 0, str(e)
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