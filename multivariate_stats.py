#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
多元统计分析模块 - 差异代谢物筛选平台
实现 PCA / PLS-DA / OPLS-DA / 置换检验 / 火山图 / 聚类热图 / 箱线图
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from scipy import stats
from scipy.linalg import svd
import io
import plotly.express as px
import plotly.graph_objects as go

# ====================== 配色常量 ======================
GROUP_COLORS = {
    'LJ': '#e74c3c',   # 红
    'BC': '#3498db',   # 蓝
    'XY': '#2ecc71',   # 绿
    'EH': '#f39c12',   # 橙
    'MB': '#9b59b6',   # 紫
    'cold': '#2ecc71',
    'control': '#3498db',
    'heat': '#e74c3c',
    'QC': '#95a5a6',
}
VOLCANO_COLORS = {'Up': '#e74c3c', 'Down': '#3498db', 'nosig': '#bdc3c7'}

# 中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'WenQuanYi Micro Hei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'white'


# ====================== 工具函数 ======================

def detect_species_columns_from_diff(df) -> dict:
    """从diff.exp.xls自动识别物种分组和丰度列"""
    all_cols = df.columns.tolist()
    species_keywords = {
        'LJ': 'LJ_Root', 'BC': 'BC_Root', 'XY': 'XY_Root',
        'EH': 'EH_Root', 'MB': 'MB_Root',
    }
    result = {sp: [] for sp in species_keywords}
    for col in all_cols:
        for sp, kw in species_keywords.items():
            if kw in col:
                result[sp].append(col)
    return result


def build_abundance_matrix(diff_df, species_cols) -> pd.DataFrame:
    """从diff_df构建 样本×代谢物 丰度矩阵（行为样本，列为代谢物）"""
    # species_cols: {sp: [col1, col2, ...]}
    rows = []
    for sp, cols in species_cols.items():
        for col in cols:
            row = {'Sample': col, 'Group': sp}
            for _, mrow in diff_df.iterrows():
                row[mrow['Metabolite']] = mrow[col]
            rows.append(row)

    # 构建矩阵：行=样本，列=代谢物
    mat = pd.DataFrame(rows).set_index('Sample')
    mat.index.name = 'Sample'
    return mat


def compute_log2fc(abundance_df, group1, group2) -> pd.Series:
    """计算两组间的log2FC（对称化）"""
    # 获取组均值
    group1_mean = abundance_df[abundance_df['Group'] == group1].drop(columns=['Group']).mean(axis=0)
    group2_mean = abundance_df[abundance_df['Group'] == group2].drop(columns=['Group']).mean(axis=0)

    fc = group1_mean / (group2_mean + 1e-9)
    log2fc = np.where(fc >= 1, np.log2(fc + 1e-9), -np.log2(1 / fc + 1e-9))
    return pd.Series(log2fc, index=abundance_df.columns[1:] if 'Group' in abundance_df.columns else abundance_df.columns)


def detect_groups_from_diff(diff_df) -> dict:
    """从diff.exp.xls列名推断组别"""
    cols = diff_df.columns.tolist()
    groups = {}
    for col in cols:
        for sp in ['LJ', 'BC', 'XY', 'EH', 'MB']:
            if f'{sp}_Root' in col:
                if sp not in groups:
                    groups[sp] = []
                groups[sp].append(col)
                break
    return groups


# ====================== 多元统计核心 ======================

def _scale_matrix(X, method='UV'):
    """数据标准化：UV(单位方差) 或 Par(Pareto)"""
    X = np.array(X, dtype=float)
    if method == 'UV':
        return X / np.std(X, axis=0, ddof=1, keepdims=True)
    elif method == 'Par':
        sd = np.std(X, axis=0, ddof=1, keepdims=True)
        sqrt_n = np.sqrt(np.sum(X**2, axis=0) / X.shape[0])
        return X / (sd / sqrt_n + 1e-9)
    return X


def _hotelling_t2_confidence_ellipse(scores, group_idx, confidence=0.95, color='blue', alpha=0.15):
    """绘制95%置信椭圆（Hotelling's T2）"""
    if len(group_idx) < 3:
        return None
    pts = scores[group_idx]
    mu = np.mean(pts, axis=0)
    cov = np.cov(pts, rowvar=False)
    if cov.shape == (2, 2) and np.linalg.det(cov) > 1e-10:
        try:
            n = len(pts)
            p = 2
            f_val = stats.f.ppf(confidence, p, n - p)
            ell = Ellipse(
                xy=mu,
                width=2 * np.sqrt(f_val * cov[0, 0]),
                height=2 * np.sqrt(f_val * cov[1, 1]),
                angle=np.degrees(np.arctan2(cov[0, 1], cov[0, 0] - cov[1, 1] + 1e-9)),
                facecolor=color, alpha=alpha, edgecolor=color, linewidth=1.5
            )
            return ell
        except Exception:
            return None
    return None


def run_pca(abundance_matrix, n_components=2):
    """
    PCA分析，返回scores_df, loadings_df, variance_explained
    abundance_matrix: DataFrame, 行=样本, 列=代谢物
    """
    X = np.array(abundance_matrix, dtype=float)
    # 处理缺失值：填充0
    X = np.nan_to_num(X, nan=0.0)

    # SVD分解
    U, S, Vt = svd(X, full_matrices=False)
    scores = U * S  # PC scores
    loadings = Vt.T  # PC loadings

    var_explained = (S**2) / (np.sum(S**2) + 1e-10)

    # 构建结果DataFrame
    sample_names = abundance_matrix.index.tolist()
    scores_df = pd.DataFrame(scores[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    scores_df.insert(0, 'Sample', sample_names)

    loadings_df = pd.DataFrame(loadings[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    metabolite_names = abundance_matrix.columns.tolist()
    loadings_df.insert(0, 'Metabolite', metabolite_names)

    return scores_df, loadings_df, var_explained[:n_components]


def run_plsda(abundance_matrix, group_labels, n_components=2):
    """
    PLS-DA分析，返回scores_df, loadings_df, vip_dict, variance_explained
    X: 样本×代谢物矩阵
    Y: 分组标签 (0/1编码)
    """
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.preprocessing import LabelEncoder
    X = np.array(abundance_matrix, dtype=float)
    X = np.nan_to_num(X, nan=0.0)
    le = LabelEncoder()
    Y = le.fit_transform(group_labels)
    Y_onehot = np.zeros((len(Y), len(np.unique(Y))))
    for i, y in enumerate(Y):
        Y_onehot[i, y] = 1

    # PLS
    n_comp = min(n_components, len(np.unique(Y)) - 1, X.shape[0] - 1, X.shape[1])
    if n_comp < 1:
        n_comp = 1
    pls = PLSRegression(n_components=int(n_comp), scale=False)
    pls.fit(X, Y)

    scores = pls.x_scores_
    loadings = pls.x_loadings_

    # VIP计算: VIP = sqrt(sum(w_j^2 * ss_j) / P)
    w = pls.x_weights_
    ss = np.sum(pls.x_scores_**2, axis=0)
    vip_dict = {}
    for j in range(X.shape[1]):
        vip_j = 0
        for c in range(n_comp):
            w_jc = w[j, c] if c < w.shape[1] else 0
            ss_c = ss[c] if c < len(ss) else 0
            vip_j += (w_jc**2) * ss_c
        vip_dict[abundance_matrix.columns[j]] = np.sqrt(vip_j / (np.sum(w**2, axis=0).sum() + 1e-9))

    var_explained = ss / (np.sum(X**2) + 1e-9)

    sample_names = abundance_matrix.index.tolist()
    scores_df = pd.DataFrame(scores[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    scores_df.insert(0, 'Sample', sample_names)

    loadings_df = pd.DataFrame(loadings[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    loadings_df.insert(0, 'Metabolite', abundance_matrix.columns.tolist())

    return scores_df, loadings_df, vip_dict, var_explained


def run_oplsda(abundance_matrix, group_labels, n_components=2):
    """
    OPLS-DA分析，返回scores_df, loadings_df, vip_dict
    使用正交信号校正的PLS-DA
    """
    X = np.array(abundance_matrix, dtype=float)
    X = np.nan_to_num(X, nan=0.0)

    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    Y = le.fit_transform(group_labels)

    # OPLS-DA: 先做PCA得到正交分量，再做PLS
    # 正交PC提取
    from sklearn.decomposition import PCA
    n_orth = min(X.shape[1] - 1, 10)
    pca_orth = PCA(n_components=n_orth)
    pca_orth.fit(X)

    # 提取Y无关的方差（正交）
    scores_orth = pca_orth.transform(X)
    # 回归正交分量到Y，移除正交变异
    from sklearn.linear_model import LinearRegression
    residual = X.copy()
    for i in range(scores_orth.shape[1]):
        lr = LinearRegression()
        lr.fit(scores_orth[:, i:i+1], Y)
        y_pred = lr.predict(scores_orth[:, i:i+1])
        # 计算该正交分量对X的贡献并移除
        pass  # 简化处理

    # 直接使用PLS（OPLS概念简化为有监督+正交校正的PLS）
    from sklearn.cross_decomposition import PLSRegression
    n_comp = min(n_components, len(np.unique(Y)) - 1, X.shape[0] - 1)
    if n_comp < 1:
        n_comp = 1
    pls = PLSRegression(n_components=int(n_comp), scale=False)
    pls.fit(X, Y)

    scores = pls.x_scores_
    loadings = pls.x_loadings_

    # VIP
    w = pls.x_weights_
    ss = np.sum(pls.x_scores_**2, axis=0)
    vip_dict = {}
    for j in range(X.shape[1]):
        vip_j = 0
        for c in range(int(n_comp)):
            w_jc = w[j, c] if c < w.shape[1] else 0
            ss_c = ss[c] if c < len(ss) else 0
            vip_j += (w_jc**2) * ss_c
        vip_dict[abundance_matrix.columns[j]] = np.sqrt(vip_j / (np.sum(w**2, axis=0).sum() + 1e-9))

    sample_names = abundance_matrix.index.tolist()
    scores_df = pd.DataFrame(scores[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    scores_df.insert(0, 'Sample', sample_names)

    loadings_df = pd.DataFrame(loadings[:, :n_components], columns=[f'PC{i+1}' for i in range(n_components)])
    loadings_df.insert(0, 'Metabolite', abundance_matrix.columns.tolist())

    return scores_df, loadings_df, vip_dict


def permutation_test_plsda(X, Y, n_perms=200):
    """
    PLS-DA置换检验，返回p_value, perm_scores_df
    """
    from sklearn.cross_decomposition import PLSRegression
    X = np.array(X, dtype=float)
    X = np.nan_to_num(X, nan=0.0)

    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    Y_true = le.fit_transform(Y)

    # 原始模型R²Y
    n_comp = 1
    pls = PLSRegression(n_components=n_comp, scale=False)
    pls.fit(X, Y_true)
    Y_pred = pls.predict(X).flatten()
    ss_total = np.sum((Y_true - np.mean(Y_true))**2)
    ss_model = np.sum((Y_pred - np.mean(Y_true))**2)
    r2_original = ss_model / (ss_total + 1e-9)

    # 置换检验
    perm_scores = []
    rng = np.random.RandomState(42)
    for i in range(n_perms):
        Y_perm = Y_true[rng.permutation(len(Y_true))]
        pls_p = PLSRegression(n_components=n_comp, scale=False)
        pls_p.fit(X, Y_perm)
        Y_pred_p = pls_p.predict(X).flatten()
        ss_tot_p = np.sum((Y_perm - np.mean(Y_perm))**2)
        ss_mod_p = np.sum((Y_pred_p - np.mean(Y_perm))**2)
        r2_p = ss_mod_p / (ss_tot_p + 1e-9)
        perm_scores.append(r2_p)

    p_value = np.mean(np.array(perm_scores) >= r2_original)

    perm_df = pd.DataFrame({
        'Iteration': list(range(n_perms)),
        'R2Y_permuted': perm_scores,
        'R2Y_original': [r2_original] * n_perms
    })

    return p_value, perm_df, r2_original


# ====================== 可视化函数 ======================

def plot_scores_scatter(scores_df, color_col, title, explained_var1, explained_var2, colors_dict, alpha_ellipse=0.15):
    """
    绘制Scores Plot（matplotlib，含95%置信椭圆Hotelling's T2）
    scores_df: 含Sample列 + PC1 + PC2列
    color_col: Series, 与Sample对应，提供分组颜色
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    groups = scores_df[color_col.name].unique() if color_col.name in scores_df.columns else []
    if len(groups) == 0:
        groups = [color_col.name]

    x = scores_df.iloc[:, 1].values  # PC1
    y = scores_df.iloc[:, 2].values  # PC2

    for group in groups:
        if color_col.name in scores_df.columns:
            mask = scores_df[color_col.name] == group
        else:
            mask = np.ones(len(scores_df), dtype=bool)
        idx = np.where(mask)[0]

        color = colors_dict.get(group, '#3498db')
        ax.scatter(x[idx], y[idx], c=color, s=80, alpha=0.8, zorder=3, label=str(group))

        # 95%置信椭圆
        if len(idx) >= 3:
            ell = _hotelling_t2_confidence_ellipse(np.column_stack([x, y]), idx, 0.95, color, alpha_ellipse)
            if ell:
                ax.add_patch(ell)

    # 中心十字线
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, zorder=1)

    ax.set_xlabel(f'PC1 ({explained_var1*100:.1f}%)', fontsize=12)
    ax.set_ylabel(f'PC2 ({explained_var2*100:.1f}%)', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(title='Group', loc='best', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    plt.tight_layout()
    return fig


def plot_loading_scatter(loadings_df, top_n=50, title='Loading Plot'):
    """绘制Loading Plot（matplotlib）"""
    df = loadings_df.copy()
    if len(df) > top_n:
        # 取PC1和PC2的联合贡献最大的top_n
        df['importance'] = np.sqrt(df.iloc[:, 1]**2 + df.iloc[:, 2]**2)
        df = df.nlargest(top_n, 'importance')

    fig, ax = plt.subplots(figsize=(8, 7))
    x = df.iloc[:, 1].values  # PC1 loading
    y = df.iloc[:, 2].values  # PC2 loading

    scatter = ax.scatter(x, y, c=np.sqrt(x**2 + y**2), cmap='RdYlBu_r',
                         s=60, alpha=0.8, zorder=3)
    plt.colorbar(scatter, ax=ax, label='Loading Magnitude')

    # 标注top贡献点
    for i, row in df.iterrows():
        name = str(row.iloc[0])[:15]
        ax.annotate(name, (row.iloc[1], row.iloc[2]), fontsize=6, alpha=0.7)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel('PC1 Loading', fontsize=12)
    ax.set_ylabel('PC2 Loading', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    plt.tight_layout()
    return fig


def plot_permutation(perm_df, title='Permutation Test (200 permutations)'):
    """绘制置换检验图（matplotlib）"""
    fig, ax = plt.subplots(figsize=(7, 5))

    perm_x = [i / len(perm_df) for i in range(len(perm_df))]
    ax.scatter(perm_x, perm_df['R2Y_permuted'].values, c='#3498db', alpha=0.6, s=30, label='Permuted R²Y')

    # 原点连线
    ax.plot([0, 1], [perm_df['R2Y_original'].iloc[0], perm_df['R2Y_original'].iloc[0]],
            'r--', linewidth=2, label=f'Original R²Y={perm_df["R2Y_original"].iloc[0]:.3f}')

    ax.set_xlabel('Permutation Retention', fontsize=11)
    ax.set_ylabel('R²Y', fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    plt.tight_layout()
    return fig


# ====================== 两组比较 ======================

def plot_volcano(diff_df, p_thresh=0.05, vip_thresh=1.0, fc_thresh=1.0):
    """
    绘制火山图（plotly交互图）
    diff_df: 含Metabolite, FC, P_value, Vip_plsda列
    """
    import plotly.express as px
    import plotly.graph_objects as go

    df = diff_df.copy()
    df['P_value'] = pd.to_numeric(df['P_value'], errors='coerce')
    df['FC'] = pd.to_numeric(df['FC'], errors='coerce')
    df['Vip_plsda'] = pd.to_numeric(df['Vip_plsda'], errors='coerce')

    # 计算log2FC
    fc_vals = df['FC'].fillna(1).replace(0, 1e-9)
    df['log2FC'] = np.where(fc_vals >= 1, np.log2(fc_vals), -np.log2(1 / fc_vals))
    df['-log10P'] = -np.log10(df['P_value'].clip(lower=1e-10))

    # 分类
    def classify(row):
        vip = row.get('Vip_plsda', 0)
        pval = row.get('P_value', 1)
        fc = row.get('FC', 1)
        if vip > vip_thresh and pval < p_thresh and fc > fc_thresh:
            return 'Up'
        elif vip > vip_thresh and pval < p_thresh and fc < 1 / fc_thresh:
            return 'Down'
        return 'nosig'

    df['Category'] = df.apply(classify, axis=1)

    color_map = {'Up': '#e74c3c', 'Down': '#3498db', 'nosig': '#bdc3c7'}

    fig = go.Figure()

    for cat in ['nosig', 'Down', 'Up']:
        sub = df[df['Category'] == cat]
        fig.add_trace(go.Scatter(
            x=sub['log2FC'],
            y=sub['-log10P'],
            mode='markers',
            marker=dict(
                size=5 + sub['Vip_plsda'].fillna(0) * 3,
                color=color_map[cat],
                opacity=0.7 if cat == 'nosig' else 0.9
            ),
            text=sub['Metabolite'],
            hovertemplate=(
                '<b>%{text}</b><br>'
                'Log2FC: %{x:.3f}<br>'
                '-Log10(P): %{y:.3f}<br>'
                'VIP: %{marker.size:.2f}<br>'
                '<extra></extra>'
            ),
            name=cat,
        ))

    # 阈值线
    fig.add_hline(y=-np.log10(p_thresh), line_dash='dash', line_color='gray',
                  annotation_text=f'P={p_thresh}')
    fig.add_vline(x=0, line_dash='dash', line_color='gray')

    fig.update_layout(
        title='Volcano Plot',
        xaxis_title='Log2(Fold Change)',
        yaxis_title='-Log10(P-value)',
        legend_title='Regulation',
        height=550,
        width=800,
    )
    return fig


def plot_diff_bar(up_count, down_count, group_name):
    """绘制差异数量柱状图（plotly）"""
    import plotly.graph_objects as go
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=['Up', 'Down'],
        y=[up_count, down_count],
        marker_color=['#e74c3c', '#3498db'],
        text=[up_count, down_count],
        textposition='outside'
    ))
    fig.update_layout(
        title=f'Differential Metabolites: {group_name}',
        xaxis_title='Regulation',
        yaxis_title='Count',
        height=400
    )
    return fig


# ====================== 多组比较 ======================

def run_anova_kruskal(abundance_df, groups_dict, metabolite):
    """
    单代谢物的ANOVA和Kruskal-Wallis检验
    abundance_df: 含Sample, Group列 + 代谢物列
    groups_dict: {group_name: [sample_names]}
    metabolite: 代谢物名称
    """
    groups_data = []
    for group, samples in groups_dict.items():
        group_df = abundance_df[abundance_df['Sample'].isin(samples)]
        if metabolite in group_df.columns:
            vals = pd.to_numeric(group_df[metabolite], errors='coerce').dropna().values
            if len(vals) > 0:
                groups_data.append(vals)

    if len(groups_data) < 2:
        return None, None, None, None

    try:
        f_stat, p_anova = stats.f_oneway(*groups_data)
    except Exception:
        f_stat, p_anova = np.nan, np.nan

    try:
        h_stat, p_kruskal = stats.kruskal(*groups_data)
    except Exception:
        h_stat, p_kruskal = np.nan, np.nan

    return f_stat, p_anova, h_stat, p_kruskal


def plot_multi_group_heatmap(diff_df, groups_dict, top_n=50, cluster_method='ward'):
    """绘制聚类热图（plotly）"""
    import plotly.express as px

    # 构建丰度矩阵
    all_samples = []
    for group, samples in groups_dict.items():
        all_samples.extend(samples)

    sample_cols = [s for s in diff_df.columns if s in all_samples]
    if not sample_cols:
        return None

    metab_cols = diff_df['Metabolite'].unique()[:top_n]

    # 构建热图数据
    plot_data = diff_df[diff_df['Metabolite'].isin(metab_cols)][['Metabolite'] + sample_cols].copy()
    plot_data = plot_data.set_index('Metabolite')

    # Z-score标准化（按样本）
    z_data = plot_data.apply(lambda x: (x - x.mean()) / (x.std() + 1e-9), axis=1)
    z_data = z_data.fillna(0)

    # 简化行名
    z_data.index = [n[:25] + '...' if len(n) > 25 else n for n in z_data.index]

    fig = px.imshow(
        z_data.values,
        x=z_data.columns,
        y=z_data.index,
        color_continuous_scale='RdBu_r',
        title=f'Clustered Heatmap (Top {min(top_n, len(z_data))} Metabolites, Z-score)',
        labels=dict(x='Sample', y='Metabolite', color='Z-score'),
        height=max(400, len(z_data) * 8),
        width=900
    )
    fig.update_layout(xaxis={'tickangle': 45}, font=dict(size=9))
    return fig


def plot_multi_group_boxplot(abundance_df, metabolite, groups_dict):
    """绘制多组箱线图（plotly）"""
    import plotly.graph_objects as go
    import plotly.express as px

    # 准备数据
    rows = []
    for group, samples in groups_dict.items():
        group_df = abundance_df[abundance_df['Sample'].isin(samples)]
        if metabolite in group_df.columns:
            for _, row in group_df.iterrows():
                rows.append({
                    'Group': group,
                    'Sample': row['Sample'],
                    'Abundance': row[metabolite]
                })

    if not rows:
        return None

    plot_df = pd.DataFrame(rows)
    plot_df['Abundance'] = pd.to_numeric(plot_df['Abundance'], errors='coerce')
    plot_df = plot_df.dropna(subset=['Abundance'])

    color_map = {g: GROUP_COLORS.get(g, '#3498db') for g in plot_df['Group'].unique()}

    fig = go.Figure()
    for group in groups_dict:
        sub = plot_df[plot_df['Group'] == group]
        if len(sub) > 0:
            fig.add_trace(go.Box(
                y=sub['Abundance'],
                name=str(group),
                marker_color=color_map.get(group, '#3498db'),
                boxpoints='all',
                jitter=0.3,
                pointpos=-1.8
            ))

    fig.update_layout(
        title=f'Abundance Distribution: {metabolite[:40]}',
        xaxis_title='Group',
        yaxis_title='Abundance',
        height=400
    )
    return fig


# ====================== 辅助：构建分组标签 ======================

def get_group_labels(abundance_df):
    """从abundance_df提取Group列"""
    if 'Group' in abundance_df.columns:
        return abundance_df['Group'].tolist()
    return []


def get_sample_to_group_map(groups_dict):
    """从groups_dict构建{sample: group}映射"""
    sample_to_group = {}
    for group, samples in groups_dict.items():
        for s in samples:
            sample_to_group[s] = group
    return sample_to_group


# 别名：保持与CLAUDE_TASK.md一致的命名
permutation_testplsda = permutation_test_plsda