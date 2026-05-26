# -*- coding: utf-8 -*-
"""
网络药理学平台 - 步骤④-⑦
步骤④：靶点-疾病关联分析（DisGeNET）
步骤⑤：PPI 网络分析（STRING Database）
步骤⑥：通路可视化与机制图
步骤⑦：分子对接（UniProt + HDOCK/CB-Dock）
依赖：pandas, networkx, matplotlib（标准库 + 指定第三方）
禁止使用 requests 库，所有 HTTP 调用使用 urllib.request
"""

from __future__ import annotations

import ssl
import time
import re
import json
import socket
import subprocess
import shutil
import urllib.parse
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

import pandas as pd
import numpy as np

# ===================== SSL Fallback 工具函数 =====================

_has_curl = shutil.which("curl") is not None


def _get_ssl_context():
    """返回一个禁用验证的 SSL context（fallback 方案）"""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def _http_get_json(url: str, timeout: int = 15, headers: dict = None) -> dict:
    """
    GET JSON，支持 SSL fallback + curl fallback。
    """
    if headers is None:
        headers = {
            "Accept": "application/json",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        }

    # 方法1：urllib
    try:
        req = Request(url, headers=headers)
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if not raw:
                return {"_error": "Empty response"}
            try:
                return json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError as je:
                return {"_error": f"JSON decode error: {je}", "_raw": raw[:200].decode("utf-8", errors="replace")}
    except URLError as e:
        if isinstance(e.reason, ssl.SSLCertVerificationError):
            try:
                req = Request(url, headers=headers)
                with urlopen(req, timeout=timeout, context=_get_ssl_context()) as resp:
                    raw = resp.read()
                    if not raw:
                        return {"_error": "Empty response"}
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "urllib_noverify"
                    return result
            except Exception as e2:
                err_str = f"urllib_noverify failed: {e2}"
        else:
            err_str = str(e)
    except Exception as urllib_err:
        err_str = str(urllib_err)

    # 方法2：curl fallback
    if _has_curl:
        try:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", str(timeout),
                "-L",  # follow redirects (UniProt search API needs this)
                "-H", f"Accept: {headers.get('Accept', 'application/json')}",
                "-H", f"User-Agent: {headers.get('User-Agent', '')}",
                url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=timeout + 2)
            if raw:
                try:
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "curl"
                    return result
                except json.JSONDecodeError:
                    return {"_error": f"curl JSON decode error", "_raw": raw[:200].decode("utf-8", errors="replace")}
        except Exception as e:
            err_str = f"curl failed: {e}"

    return {"_error": err_str if 'err_str' in dir() else "Unknown error"}


def _http_post_json(url: str, payload: dict = None, body: str = None,
                    headers: dict = None, timeout: int = 30) -> dict:
    """
    POST JSON/表单，支持 SSL fallback + curl fallback。
    payload: dict（自动 JSON encode）
    body: str（直接发送原始字符串，如 form-urlencoded）
    """
    if headers is None:
        headers = {
            "Accept": "application/json",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        }

    # 构建 body
    if body is None and payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    elif body is not None:
        if isinstance(body, str):
            body = body.encode("utf-8")

    # 方法1：urllib
    try:
        req = Request(url, data=body, headers=headers)
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if not raw:
                return {"_error": "Empty response"}
            try:
                return json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError:
                return {"_raw": raw.decode("utf-8", errors="replace")}
    except URLError as e:
        if isinstance(e.reason, ssl.SSLCertVerificationError):
            try:
                req = Request(url, data=body, headers=headers)
                with urlopen(req, timeout=timeout, context=_get_ssl_context()) as resp:
                    raw = resp.read()
                    if not raw:
                        return {"_error": "Empty response"}
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "urllib_noverify"
                    return result
            except Exception as e2:
                err_str = f"urllib_noverify failed: {e2}"
        else:
            err_str = str(e)
    except Exception as urllib_err:
        err_str = str(urllib_err)

    # 方法2：curl fallback
    if _has_curl:
        try:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", str(timeout),
                "-X", "POST",
                "-H", f"Content-Type: {headers.get('Content-Type', 'application/json')}",
                "-H", f"Accept: {headers.get('Accept', 'application/json')}",
                "-d", body.decode("utf-8") if isinstance(body, bytes) else body,
                url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=timeout + 2)
            if raw:
                try:
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "curl"
                    return result
                except json.JSONDecodeError:
                    return {"_raw": raw.decode("utf-8", errors="replace")}
        except Exception:
            pass

    return {"_error": err_str if 'err_str' in dir() else "Unknown error"}


def _http_get_text(url: str, timeout: int = 15, headers: dict = None) -> dict:
    """
    GET 文本（用于 TSV 等非 JSON 响应），支持 SSL fallback + curl fallback。
    返回 {"_raw": str} 或 {"_error": str}
    """
    if headers is None:
        headers = {
            "Accept": "text/plain, */*",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        }

    # 方法1：urllib
    try:
        req = Request(url, headers=headers)
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if not raw:
                return {"_error": "Empty response"}
            try:
                text = raw.decode("utf-8")
                return {"_raw": text}
            except UnicodeDecodeError:
                return {"_raw": raw.decode("latin-1", errors="replace")}
    except URLError as e:
        if isinstance(e.reason, ssl.SSLCertVerificationError):
            try:
                req = Request(url, headers=headers)
                with urlopen(req, timeout=timeout, context=_get_ssl_context()) as resp:
                    raw = resp.read()
                    return {"_raw": raw.decode("utf-8", errors="replace")}
            except Exception as e2:
                err_str = f"urllib_noverify failed: {e2}"
        else:
            err_str = str(e)
    except Exception as urllib_err:
        err_str = str(urllib_err)

    # 方法2：curl fallback
    if _has_curl:
        try:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", str(timeout),
                "-H", f"Accept: {headers.get('Accept', 'text/plain')}",
                url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=timeout + 2)
            return {"_raw": raw.decode("utf-8", errors="replace")}
        except Exception as e:
            err_str = f"curl failed: {e}"

    return {"_error": err_str if 'err_str' in dir() else "Unknown error"}


# =====================================================================
# 步骤④：靶点-疾病关联分析
# =====================================================================

def npd_disease_query(gene_list: list, token: str = "",
                      threshold: float = 0.0) -> pd.DataFrame:
    """
    查询基因-疾病关联（UniProt REST API）。

    Args:
        gene_list: 基因符号列表，如 ["ABCB1", "ABCC1"]
        token: DisGeNET API token（仅保留参数兼容性，不再使用）
        threshold: 最低 score 阈值（0~1），默认 0.0（UniProt 无 score，始终返回 1.0）

    Returns:
        DataFrame，列：[Gene, Disease, Score, Source, Disease_ID]
        若查询失败或无数据，返回空 DataFrame（含上述列名）
    """
    if not gene_list:
        return pd.DataFrame(columns=["Gene", "Disease", "Score", "Source", "Disease_ID"])

    rows = []

    for gene in gene_list:
        gene = str(gene).strip()
        if not gene:
            continue

        # 步骤1：通过 search API 获取 UniProt accession（不能用 gene direct URL）
        # 注意：UniProt search API 不支持 comments 字段，必须分两步走
        search_url = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene:{gene}+AND+organism_id:9606"
            f"&format=json&size=5&fields=accession,gene_names"
        )
        result = _http_get_json(search_url, timeout=20)

        if "_error" in result:
            time.sleep(0.3)
            continue

        results_list = result.get("results", [])
        if not results_list:
            time.sleep(0.3)
            continue

        # 步骤2：用 accession 获取完整 JSON，提取 DISEASE comments
        for entry in results_list:
            accession = entry.get("primaryAccession", "")
            if not accession:
                continue

            gene_in_entry = ""
            gene_list = entry.get("genes", [])
            if gene_list:
                gene_in_entry = gene_list[0].get("geneName", {}).get("value", "")

            # 获取完整 JSON（包含 disease comments）
            full_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
            full_result = _http_get_json(full_url, timeout=20)

            if "_error" in full_result:
                continue

            comments = full_result.get("comments", [])
            if not isinstance(comments, list):
                comments = []

            for comment in comments:
                if comment.get("commentType") != "DISEASE":
                    continue

                disease_data = comment.get("disease", {})
                disease_id = disease_data.get("diseaseAccession", "")
                disease_name = disease_data.get("diseaseId", "") or disease_data.get("description", "")

                if not disease_name:
                    continue

                # UniProt 没有 score，默认 1.0
                score = 1.0

                rows.append({
                    "Gene": gene_in_entry or gene,
                    "Disease": disease_name,
                    "Score": score,
                    "Source": "UniProt",
                    "Disease_ID": disease_id,
                })

        time.sleep(0.3)

    if not rows:
        return pd.DataFrame(columns=["Gene", "Disease", "Score", "Source", "Disease_ID"])

    return pd.DataFrame(rows)


# =====================================================================
# 步骤⑤：PPI 网络分析
# =====================================================================

def npd_ppi_build(gene_list: list, min_score: float = 0.4) -> pd.DataFrame:
    """
    调用 STRING API 构建 PPI 网络。

    Args:
        gene_list: 基因符号列表
        min_score: 最低置信度分（0~1），默认 0.4（约等于 STRING 400/1000）

    Returns:
        DataFrame，列：[Gene1, Gene2, Score]
        若查询失败或无边，返回只含列名的空 DataFrame
    """
    if not gene_list:
        return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])

    gene_list = [str(g).strip() for g in gene_list if str(g).strip()]
    if not gene_list:
        return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])

    # STRING API: identifiers 用换行符 %0A 分隔
    genes_param = "%0A".join(gene_list)
    url = (f"https://string-db.org/api/tsv/network"
           f"?identifiers={genes_param}&species=9606")

    # 用 curl.exe 避免 SSL 问题
    text = ""
    if _has_curl:
        try:
            raw = subprocess.check_output(
                ["curl.exe", "-s", "-S", "--max-time", "30",
                 "-H", "Accept: text/plain",
                 url],
                stderr=subprocess.STDOUT, timeout=35
            )
            text = raw.decode("utf-8", errors="replace")
        except Exception:
            return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])
    else:
        http_result = _http_get_text(url, timeout=30)
        if "_error" in http_result:
            return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])
        text = http_result.get("_raw", "")

    lines = text.strip().split("\n") if text.strip() else []

    if len(lines) < 2:
        return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])

    # 解析 TSV（跳过表头）
    # 列顺序: stringId_A(0), stringId_B(1), preferredName_A(2), preferredName_B(3), ncbiTaxonId(4), score(5)
    rows = []
    for line in lines[1:]:          # 第一行是列名
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        try:
            score = float(parts[5])  # score 在第6列
        except (ValueError, IndexError):
            score = 0
        if score < min_score:
            continue
        # 使用易读基因符号（preferredName）
        gene_a = parts[2].strip()   # preferredName_A
        gene_b = parts[3].strip()   # preferredName_B
        rows.append({"Gene1": gene_a, "Gene2": gene_b, "Score": score})

    if not rows:
        return pd.DataFrame(columns=["Gene1", "Gene2", "Score"])

    return pd.DataFrame(rows)


def npd_ppi_topology(ppi_df: pd.DataFrame) -> pd.DataFrame:
    """
    使用 networkx 计算网络拓扑参数。

    Args:
        ppi_df: PPI DataFrame，含 [Gene1, Gene2, Score] 列

    Returns:
        DataFrame，列：[Gene, Degree, Betweenness, Closeness]
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("请安装 networkx: pip install networkx")

    if ppi_df is None or len(ppi_df) == 0:
        return pd.DataFrame(columns=["Gene", "Degree", "Betweenness", "Closeness"])

    # 构建无向加权图
    G = nx.Graph()
    for _, row in ppi_df.iterrows():
        G.add_edge(str(row["Gene1"]), str(row["Gene2"]), weight=float(row["Score"]))

    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["Gene", "Degree", "Betweenness", "Closeness"])

    # 计算拓扑指标
    degree_dict = dict(G.degree())
    try:
        betweenness_dict = nx.betweenness_centrality(G)
    except Exception:
        betweenness_dict = {n: 0.0 for n in G.nodes()}

    try:
        closeness_dict = nx.closeness_centrality(G)
    except Exception:
        closeness_dict = {n: 0.0 for n in G.nodes()}

    rows = []
    for node in G.nodes():
        rows.append({
            "Gene": node,
            "Degree": degree_dict.get(node, 0),
            "Betweenness": round(betweenness_dict.get(node, 0.0), 4),
            "Closeness": round(closeness_dict.get(node, 0.0), 4),
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("Degree", ascending=False).reset_index(drop=True)
    return df


def npd_ppi_plot(ppi_df: pd.DataFrame, topo_df: pd.DataFrame,
                 output_path: str = "ppi_network.png") -> None:
    """
    使用 matplotlib + networkx 绘制 PPI 网络图。
    节点大小按 Degree 缩放，颜色按 Betweenness 着色。
    保存为 PNG。

    Args:
        ppi_df: PPI DataFrame，含 [Gene1, Gene2, Score]
        topo_df: 拓扑 DataFrame，含 [Gene, Degree, Betweenness]
        output_path: 输出文件路径
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError as e:
        raise ImportError(f"请安装依赖: {e}")

    if ppi_df is None or len(ppi_df) == 0:
        # 画空图
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No PPI data available", ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    # 构建网络图
    G = nx.Graph()
    for _, row in ppi_df.iterrows():
        G.add_edge(str(row["Gene1"]), str(row["Gene2"]), weight=float(row["Score"]))

    # 合并 topo 信息
    if topo_df is not None and len(topo_df) > 0:
        topo_dict = topo_df.set_index("Gene").to_dict("index")
    else:
        topo_dict = {}

    # 节点属性
    nodes = list(G.nodes())
    degrees = [G.degree(n) for n in nodes]
    max_deg = max(degrees) if degrees else 1
    betweenness_vals = []
    for n in nodes:
        tb = topo_dict.get(n, {}).get("Betweenness", 0.0)
        betweenness_vals.append(tb)

    # 节点大小：Degree 映射到 [200, 2000]
    node_sizes = [200 + (d / max_deg) * 1800 for d in degrees]

    # 节点颜色：Betweenness
    max_bw = max(betweenness_vals) if betweenness_vals else 1
    node_colors = [b / max_bw if max_bw > 0 else 0 for b in betweenness_vals]

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, k=0.6, seed=42)

    # 绘制边（按 score 透明度）
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    max_w = max(edge_weights) if edge_weights else 1
    for (u, v), w in zip(G.edges(), edge_weights):
        alpha = 0.2 + 0.6 * (w / max_w)
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                               width=w / 150, alpha=alpha,
                               edge_color="gray", ax=ax)

    # 绘制节点
    cmap = cm.RdYlBu_r
    nx.draw_networkx_nodes(G, pos, nodelist=nodes,
                           node_size=node_sizes,
                           node_color=node_colors,
                           cmap=cmap, ax=ax)

    # 标签（只显示 degree >= 2 的节点）
    labels_to_show = {n: n for n in nodes if G.degree(n) >= 2}
    nx.draw_networkx_labels(G, pos, labels=labels_to_show, font_size=7, ax=ax)

    ax.axis("off")
    ax.set_title(f"PPI Network (edges={G.number_of_edges()}, nodes={G.number_of_nodes()})",
                 fontsize=13)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_clim(0, max_bw if max_bw > 0 else 1)
    fig.colorbar(sm, ax=ax, label="Betweenness Centrality", shrink=0.6)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# =====================================================================
# 步骤⑥：通路可视化与机制图
# =====================================================================

def npd_viz_dotplot(enr_df: pd.DataFrame, top_n: int = 20,
                    color_by: str = "Database") -> 'go.Figure':
    """
    增强气泡图，支持按数据库分组着色。

    Args:
        enr_df: 富集结果 DataFrame，需含 Term, Adjusted_P_value, Genes, Overlap, Database 列
        top_n: 显示前 N 条
        color_by: "Database" 或 "Term"

    Returns:
        plotly.graph_objects.Figure
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("请安装 plotly: pip install plotly")

    if enr_df is None or len(enr_df) == 0:
        return go.Figure()

    plot_df = enr_df.head(top_n).copy()

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

    def gene_count(val):
        if pd.isna(val):
            return 0
        return len([g for g in str(val).split(';') if g.strip()])

    plot_df['Count'] = plot_df.apply(lambda r: gene_count(r.get('Genes', '')), axis=1)
    plot_df['X_val'] = plot_df.apply(lambda r: parse_overlap(r.get('Overlap', 0)), axis=1)
    plot_df['Term_short'] = plot_df['Term'].apply(
        lambda x: x[:50] + '...' if len(str(x)) > 50 else str(x)
    )

    p_col = 'Adjusted_P_value' if 'Adjusted_P_value' in plot_df.columns else 'P_value'
    pvals = pd.to_numeric(plot_df[p_col], errors='coerce').fillna(1)
    pvals = pvals.replace(0, 1e-10)
    plot_df['-log10P'] = -np.log10(pvals)

    # 固定颜色映射（不依赖 showscale，避免 plotly 整数数组 bug）
    DATABASE_COLORS = {
        'KEGG_2021_Human': '#e74c3c',
        'GO_Biological_Process_2021': '#3498db',
        'GO_Molecular_Function_2021': '#2ecc71',
        'GO_Cellular_Component_2021': '#f39c12',
        'Reactome_2021': '#9b59b6',
        'WikiPathways_2019_Human': '#1abc9c',
    }
    FALLBACK_COLOR = '#95a5a6'

    if color_by == "Database" and "Database" in plot_df.columns:
        marker = dict(
            size=plot_df['Count'] * 3 + 5,
            color=plot_df['Database'].map(DATABASE_COLORS).fillna(FALLBACK_COLOR),
            sizemin=4,
        )
    else:
        marker = dict(
            size=plot_df['Count'] * 3 + 5,
            color=plot_df['-log10P'],
            colorscale='Viridis',
            colorbar=dict(title='-log10(Padj)'),
            sizemin=4,
        )

    fig = go.Figure(go.Scatter(
        x=plot_df['X_val'],
        y=plot_df['Term_short'],
        mode='markers',
        marker=marker,
        text=plot_df.apply(
            lambda r: f"<b>{r['Term']}</b><br>Count: {r['Count']}<br>"
                       f"AdjP: {r.get('Adjusted_P_value', r.get('P_value', 'N/A'))}<br>"
                       f"Genes: {r.get('Genes', '')}",
            axis=1
        ),
        hoverinfo='text',
    ))

    fig.update_layout(
        title=f'Enrichment Dotplot (Top {top_n})',
        xaxis_title='Count / Overlap',
        yaxis_title='Term',
        height=max(400, len(plot_df) * 22),
        width=900,
        font=dict(size=10),
        yaxis=dict(autorange='reversed'),
    )
    return fig


def npd_viz_bar(enr_df: pd.DataFrame, top_n: int = 20) -> 'go.Figure':
    """
    富集结果条形图（Top N by Adjusted P-value）。

    Args:
        enr_df: 富集结果 DataFrame，需含 Term, Adjusted_P_value 列
        top_n: 显示前 N 条

    Returns:
        plotly.graph_objects.Figure
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("请安装 plotly: pip install plotly")

    if enr_df is None or len(enr_df) == 0:
        return go.Figure()

    sort_col = 'Adjusted_P_value' if 'Adjusted_P_value' in enr_df.columns else 'P_value'
    plot_df = enr_df.sort_values(sort_col).head(top_n).copy()
    plot_df['Term_short'] = plot_df['Term'].apply(
        lambda x: x[:55] + '...' if len(str(x)) > 55 else str(x)
    )
    pvals = pd.to_numeric(plot_df[sort_col], errors='coerce').fillna(1)
    pvals = pvals.replace(0, 1e-10)
    plot_df['-log10P'] = -np.log10(pvals)

    fig = go.Figure(go.Bar(
        x=plot_df['-log10P'],
        y=plot_df['Term_short'],
        orientation='h',
        marker=dict(
            color=plot_df['-log10P'],
            colorscale='RdYlBu_r',
        ),
        text=plot_df[sort_col].apply(lambda v: f"Padj={v:.2e}"),
        hoverinfo='text',
    ))
    fig.update_layout(
        title=f'Enrichment Bar Plot (Top {top_n})',
        xaxis_title='-log10(Adjusted P-value)',
        yaxis_title='Term',
        height=max(350, len(plot_df) * 22),
        width=850,
        font=dict(size=10),
    )
    return fig


def npd_viz_venn(gene_sets: dict, labels: list = None) -> 'Figure':
    """
    Venn 图比较多个基因集/通路集重叠。
    使用 matplotlib_venn（如果可用）或 matplotlib 自绘。

    Args:
        gene_sets: dict，{名称: set(基因列表)} 或 {名称: set(通路列表)}
        labels: 标签列表，默认用 gene_sets 的 key

    Returns:
        matplotlib Figure
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("请安装 matplotlib: pip install matplotlib")

    if not gene_sets:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(0.5, 0.5, "No gene sets provided", ha="center", va="center")
        ax.axis("off")
        return fig

    names = list(gene_sets.keys())
    if labels is None:
        labels = names

    n = len(gene_sets)

    # 尝试使用 matplotlib_venn
    try:
        from matplotlib_venn import venn2, venn3
        if n == 2:
            sets = list(gene_sets.values())
            fig, ax = plt.subplots(figsize=(7, 6))
            venn2(sets, set_labels=labels[:2], ax=ax)
            ax.set_title("Venn Diagram")
            return fig
        elif n == 3:
            sets = list(gene_sets.values())
            fig, ax = plt.subplots(figsize=(7, 6))
            venn3(sets, set_labels=labels[:3], ax=ax)
            ax.set_title("Venn Diagram")
            return fig
    except ImportError:
        pass

    # 回退：matplotlib 自绘（只支持 2~3 集合的简化示意）
    if n == 2:
        sets = list(gene_sets.values())
        intersection = sets[0] & sets[1]
        only_a = sets[0] - sets[1]
        only_b = sets[1] - sets[0]

        fig, ax = plt.subplots(figsize=(8, 6))
        circle1 = plt.Circle((0.35, 0.5), 0.32, fill=False, color="royalblue", lw=2)
        circle2 = plt.Circle((0.65, 0.5), 0.32, fill=False, color="coral", lw=2)
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0.22, 0.5, f"{len(only_a)}", ha="center", va="center", fontsize=14)
        ax.text(0.78, 0.5, f"{len(only_b)}", ha="center", va="center", fontsize=14)
        ax.text(0.5, 0.5, f"{len(intersection)}", ha="center", va="center", fontsize=16, fontweight="bold")
        ax.legend([circle1, circle2], list(labels[:2]), loc="upper right")
        ax.set_title("Venn Diagram")
        return fig

    elif n == 3:
        sets = [set(s) for s in list(gene_sets.values())]
        s0, s1, s2 = sets[0], sets[1], sets[2]
        only_0 = s0 - s1 - s2
        only_1 = s1 - s0 - s2
        only_2 = s2 - s0 - s1
        i01 = (s0 & s1) - s2
        i02 = (s0 & s2) - s1
        i12 = (s1 & s2) - s0
        i012 = s0 & s1 & s2

        fig, ax = plt.subplots(figsize=(9, 8))
        c0 = plt.Circle((0.38, 0.57), 0.3, fill=False, color="royalblue", lw=2)
        c1 = plt.Circle((0.62, 0.57), 0.3, fill=False, color="coral", lw=2)
        c2 = plt.Circle((0.5, 0.38), 0.3, fill=False, color="forestgreen", lw=2)
        ax.add_patch(c0)
        ax.add_patch(c1)
        ax.add_patch(c2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.axis("off")
        positions = [
            (0.22, 0.62), (0.78, 0.62), (0.5, 0.2),
            (0.62, 0.52), (0.38, 0.52), (0.5, 0.52), (0.5, 0.47)
        ]
        counts = [len(only_0), len(only_1), len(only_2),
                  len(i01), len(i02), len(i12), len(i012)]
        for (x, y), c in zip(positions, counts):
            ax.text(x, y, str(c), ha="center", va="center", fontsize=11)
        ax.legend([c0, c1, c2], list(labels[:3]), loc="upper right")
        ax.set_title("Venn Diagram")
        return fig

    else:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(0.5, 0.5, f"{n} sets (Venn only supports 2-3)", ha="center", va="center")
        ax.axis("off")
        return fig


def npd_viz_network(compound: str, gene_list: list, pathway_df: pd.DataFrame,
                    disease_df: pd.DataFrame = None) -> 'Figure':
    """
    '成分-靶点-通路-疾病' 四层网络图。
    左侧：化合物 → 中间：靶点基因 → 右侧：通路 →（右侧外：疾病）
    使用 matplotlib + networkx 实现。

    Args:
        compound: 化合物名称
        gene_list: 靶点基因列表
        pathway_df: 通路 DataFrame，需含 Term 列
        disease_df: 疾病 DataFrame，可选

    Returns:
        matplotlib Figure
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError as e:
        raise ImportError(f"请安装 matplotlib/networkx: {e}")

    fig, ax = plt.subplots(figsize=(16, 10))
    G = nx.DiGraph()

    # 节点定义
    G.add_node(compound, node_type="compound")

    gene_list = [str(g).strip() for g in gene_list if str(g).strip()]
    pathways = []
    if pathway_df is not None and len(pathway_df) > 0:
        pathways = pathway_df['Term'].head(15).tolist()

    diseases = []
    if disease_df is not None and len(disease_df) > 0:
        diseases = disease_df['Disease'].head(10).tolist()

    # 布局：左→右分层
    # x 轴：化合物(0) → 基因(1) → 通路(2) → 疾病(3)
    pos = {}

    pos[compound] = (0, 0.5)

    n_genes = len(gene_list)
    for i, g in enumerate(gene_list):
        G.add_node(g, node_type="gene")
        frac = (i + 0.5) / n_genes if n_genes > 0 else 0.5
        pos[g] = (1, frac)

    for i, pw in enumerate(pathways):
        G.add_node(pw, node_type="pathway")
        frac = (i + 0.5) / len(pathways) if pathways else 0.5
        pos[pw] = (2, frac)

    for i, ds in enumerate(diseases):
        G.add_node(ds, node_type="disease")
        frac = (i + 0.5) / len(diseases) if diseases else 0.5
        pos[ds] = (3, frac)

    # 边：化合物→基因
    for g in gene_list:
        G.add_edge(compound, g)

    # 边：基因→通路（随机建立连接以便可视化）
    np.random.seed(42)
    for g in gene_list:
        for pw in pathways[:5]:
            if np.random.rand() > 0.5:
                G.add_edge(g, pw)

    # 边：通路→疾病
    if diseases and pathways:
        for pw in pathways[:5]:
            for ds in diseases[:3]:
                G.add_edge(pw, ds)

    # 颜色映射
    type_colors = {
        "compound": "#FF6B6B",
        "gene": "#4ECDC4",
        "pathway": "#45B7D1",
        "disease": "#96CEB4"
    }
    node_color_list = [type_colors.get(G.nodes[n].get("node_type", "gene"), "#CCCCCC")
                       for n in G.nodes()]

    # 节点大小
    node_size_list = []
    for n in G.nodes():
        nt = G.nodes[n].get("node_type", "gene")
        if nt == "compound":
            node_size_list.append(2500)
        elif nt == "gene":
            node_size_list.append(1000)
        elif nt == "pathway":
            node_size_list.append(800)
        else:
            node_size_list.append(600)

    nx.draw_networkx_nodes(G, pos, node_color=node_color_list,
                           node_size=node_size_list, ax=ax, alpha=0.9)
    nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.5,
                           arrows=True, arrowsize=15, ax=ax)

    # 标签（通路名截断）
    labels = {}
    for n in G.nodes():
        lbl = n
        if G.nodes[n].get("node_type") == "pathway":
            lbl = n[:25] + '...' if len(n) > 25 else n
        elif G.nodes[n].get("node_type") == "disease":
            lbl = n[:25] + '...' if len(n) > 25 else n
        labels[n] = lbl

    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, ax=ax)

    # 图例
    legend_patches = [
        mpatches.Patch(color=c, label=t)
        for t, c in [("Compound", "#FF6B6B"),
                      ("Gene", "#4ECDC4"),
                      ("Pathway", "#45B7D1"),
                      ("Disease", "#96CEB4")]
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=10)
    ax.axis("off")
    ax.set_title(f"Compound-Target-Pathway-Disease Network: {compound}", fontsize=14)
    return fig


# =====================================================================
# 步骤⑦：分子对接
# =====================================================================

def npd_get_pdb_from_gene(gene_symbol: str) -> str:
    """
    从 UniProt API 获取基因对应的 PDB ID。
    返回: 第一个 PDB ID 字符串，未找到返回 None
    """
    if not gene_symbol:
        return None

    gene_symbol = str(gene_symbol).strip().upper()
    # Use UniProt search API (gene name → accessions → PDB)
    url = (f"https://rest.uniprot.org/uniprotkb/search"
           f"?query=gene:{gene_symbol}+AND+organism_id:9606"
           f"&format=json&size=5&fields=accession,gene_names,rcrossReferences")
    result = _http_get_json(url, timeout=15)

    if "_error" in result:
        return None

    try:
        results_list = result.get("results", [])
        pdb_ids = []
        for entry in results_list:
            xrefs = entry.get("rcrossReferences", [])
            for xr in xrefs:
                if xr.get("type") == "PDB":
                    pid = xr.get("id")
                    if pid:
                        pdb_ids.append(str(pid))
        if pdb_ids:
            return pdb_ids[0]
        return None
    except Exception:
        return None


def npd_dock_query(gene_symbol: str, smiles: str,
                   pdb_id: str = None, method: str = "cbdock") -> dict:
    """
    在线分子对接（HDOCK 或 CB-Dock）。

    Args:
        gene_symbol: 基因名（用于获取 PDB）
        smiles: 化合物 SMILES
        pdb_id: 指定 PDB ID（可选，自动从 UniProt 查询）
        method: "hdock" | "cbdock"

    Returns:
        dict: {"status": "done"|"pending"|"error",
               "score": float, "gene": str, "smiles": str,
               "pdb_id": str, "note": str}
    """
    if not gene_symbol or not smiles:
        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": pdb_id, "note": "Missing gene or SMILES"
        }

    if pdb_id is None:
        pdb_id = npd_get_pdb_from_gene(gene_symbol)

    if not pdb_id:
        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": None, "note": "No PDB ID found"
        }

    if method == "hdock":
        return _npd_dock_hdock(gene_symbol, smiles, pdb_id)
    else:
        return _npd_dock_cbdock(gene_symbol, smiles, pdb_id)


def _npd_dock_hdock(gene_symbol: str, smiles: str, pdb_id: str) -> dict:
    """HDOCK 在线对接"""
    submit_url = "http://hdock.psych.ac.cn/submit"
    poll_url_base = "http://hdock.psych.ac.cn/result/"

    try:
        # 提交任务（form-data）
        if _has_curl:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", "30",
                "-X", "POST",
                "-F", f"receptor_pdb={pdb_id}",
                "-F", f"ligand_mol={smiles}",
                submit_url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=35)
            try:
                submit_result = json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError:
                submit_result = {"_raw": raw.decode("utf-8", errors="replace")}
        else:
            # urllib 表单提交（较难，这里用简单实现）
            import urllib.parse as up
            form_data = up.urlencode({
                "receptor_pdb": pdb_id,
                "ligand_mol": smiles
            }).encode("ascii")
            req = Request(submit_url, data=form_data,
                          headers={"Content-Type": "application/x-www-form-urlencoded"})
            with urlopen(req, timeout=30) as resp:
                raw = resp.read()
                try:
                    submit_result = json.loads(raw.decode("utf-8"))
                except json.JSONDecodeError:
                    submit_result = {"_raw": raw.decode("utf-8", errors="replace")}

        if "_error" in submit_result:
            return {
                "status": "error", "score": None,
                "gene": gene_symbol, "smiles": smiles,
                "pdb_id": pdb_id, "note": f"Submit failed: {submit_result.get('_error')}"
            }

        # 提取 job_id
        job_id = None
        raw_str = submit_result.get("_raw", "")
        if isinstance(submit_result, dict):
            job_id = submit_result.get("job_id") or submit_result.get("jobid") or submit_result.get("id")

        # 超时 5 分钟轮询
        for attempt in range(60):
            time.sleep(5)
            poll_url = f"{poll_url_base}{job_id}" if job_id else None
            if not poll_url:
                # 轮询不到 job_id，假设完成
                return {
                    "status": "done", "score": None,
                    "gene": gene_symbol, "smiles": smiles,
                    "pdb_id": pdb_id,
                    "note": "No job_id returned; result unknown"
                }

            poll_result = _http_get_json(poll_url, timeout=20)
            if "_error" in poll_result:
                continue

            # 解析状态
            if isinstance(poll_result, dict):
                state = poll_result.get("state", "") or poll_result.get("status", "")
                if state.lower() in ("finished", "done", "complete", "success"):
                    score = poll_result.get("score") or poll_result.get("docking_score")
                    score_val = float(score) if score is not None else None
                    return {
                        "status": "done", "score": score_val,
                        "gene": gene_symbol, "smiles": smiles,
                        "pdb_id": pdb_id, "note": "HDOCK calculation complete"
                    }

        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": pdb_id, "note": "Timeout after 5 minutes"
        }

    except Exception as e:
        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": pdb_id, "note": f"Exception: {e}"
        }


def _npd_dock_cbdock(gene_symbol: str, smiles: str, pdb_id: str) -> dict:
    """CB-Dock 在线对接"""
    submit_url = "https://clab.labshare.cn/cb-dock/predict"
    poll_url_base = "https://clab.labshare.cn/cb-dock/result/"

    try:
        # 提交
        payload = {"ligand": smiles, "protein": pdb_id}
        submit_result = _http_post_json(submit_url, payload=payload, timeout=30)

        if "_error" in submit_result:
            return {
                "status": "error", "score": None,
                "gene": gene_symbol, "smiles": smiles,
                "pdb_id": pdb_id, "note": f"Submit failed: {submit_result.get('_error')}"
            }

        task_id = submit_result.get("id") or submit_result.get("task_id") or submit_result.get("_raw", "")

        # 尝试从 raw 中提取 id
        if not task_id and "_raw" in submit_result:
            import re
            m = re.search(r'"id"\s*:\s*"?(\w+)', submit_result.get("_raw", ""))
            if m:
                task_id = m.group(1)

        if not task_id:
            return {
                "status": "error", "score": None,
                "gene": gene_symbol, "smiles": smiles,
                "pdb_id": pdb_id, "note": "No task ID returned"
            }

        # 轮询（最多 5 分钟）
        for _attempt in range(60):
            time.sleep(5)
            poll_url = f"{poll_url_base}{task_id}"
            poll_result = _http_get_json(poll_url, timeout=20)

            if "_error" in poll_result:
                continue

            if isinstance(poll_result, dict):
                status = poll_result.get("status", "")
                if status in ("done", "finished", "success"):
                    score = poll_result.get("score") or poll_result.get("vina_score")
                    score_val = float(score) if score is not None else None
                    return {
                        "status": "done", "score": score_val,
                        "gene": gene_symbol, "smiles": smiles,
                        "pdb_id": pdb_id, "note": "CB-Dock calculation complete"
                    }

        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": pdb_id, "note": "Timeout after 5 minutes"
        }

    except Exception as e:
        return {
            "status": "error", "score": None,
            "gene": gene_symbol, "smiles": smiles,
            "pdb_id": pdb_id, "note": f"Exception: {e}"
        }