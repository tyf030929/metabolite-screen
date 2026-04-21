# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    content = f.read()

start_marker = '# ====== MiniMax风格美化CSS ======'
end_marker = '# 中文字体支持'

start_idx = content.find(start_marker)
end_idx = content.find(end_marker)

if start_idx == -1 or end_idx == -1:
    print(f'NOT FOUND: start={start_idx}, end={end_idx}')
    exit(1)

print(f'Block found: start={start_idx}, end={end_idx}, len={end_idx-start_idx}')

new_css = '''# ====== MetaboLab Aura Design System CSS ======
# 参考：stitch_high_end_website_prototype - Tailwind+Material Design 3 原型设计
# 精确复刻配色、圆角、阴影、字体、动效全部细节
st.markdown("""
<style>
/* ====== Google Fonts: Inter ====== */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800;900&family=Noto+Sans+SC:wght@400;500;700&display=swap');

/* ====== 全局基础 ====== */
html, body {
    font-family: 'Inter', 'Noto Sans SC', -apple-system, sans-serif !important;
    -webkit-font-smoothing: antialiased;
    letter-spacing: -0.01em;
}

h1, h2, h3, h4 {
    letter-spacing: -0.02em !important;
}

/* ====== Streamlit 根背景 - 薰衣草白 ====== */
.stApp {
    background-color: #fef7ff !important;
}

/* ====== 顶部固定导航栏 - 玻璃拟态 ====== */
[data-testid="stHeader"] {
    background: rgba(254, 247, 255, 0.82) !important;
    backdrop-filter: blur(24px) !important;
    -webkit-backdrop-filter: blur(24px) !important;
    border-bottom: none !important;
    box-shadow: 0 40px 60px -15px rgba(99, 14, 212, 0.05) !important;
}

/* ====== 侧边栏 - 薰衣草紫背景 ====== */
[data-testid="stSidebar"] {
    background-color: #f3edf7 !important;
    border-right: none !important;
}

/* 侧边栏内部所有文字 */
[data-testid="stSidebar"] .stMarkdown p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] label {
    color: #1d1a24 !important;
    font-family: 'Inter', sans-serif !important;
}

/* 侧边栏标题 - 品牌紫 */
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 {
    color: #630ed4 !important;
    font-weight: 700 !important;
}

/* 侧边栏分隔线 */
[data-testid="stSidebar"] hr {
    border: none !important;
    border-top: 1px solid rgba(204, 195, 216, 0.5) !important;
}

/* 侧边栏 Selectbox */
[data-testid="stSidebar"] [data-baseweb="select"] > div,
[data-testid="stSidebar"] [data-baseweb="multiselect"] > div {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99, 14, 212, 0.1) !important;
}

/* 侧边栏标签 */
[data-testid="stSidebar"] .stSelectbox > label,
[data-testid="stSidebar"] .stMultiSelect > label {
    color: #4a4455 !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
}

/* 侧边栏 Slider 轨道 */
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-track {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    border: none !important;
    height: 4px !important;
    border-radius: 9999px !important;
}
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-thumb {
    background-color: #630ed4 !important;
    width: 16px !important;
    height: 16px !important;
    box-shadow: 0 2px 8px rgba(99, 14, 212, 0.4) !important;
}
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-rail {
    background-color: #e8dfee !important;
    height: 4px !important;
    border-radius: 9999px !important;
    opacity: 1 !important;
}

/* ====== Tab 标签导航 - 胶囊式 ====== */
.stTabs [data-baseweb="tab-list"] {
    background-color: #f3ebfa !important;
    border-radius: 1.5rem !important;
    padding: 6px !important;
    gap: 4px !important;
    border: none !important;
    box-shadow: none !important;
}

.stTabs [data-baseweb="tab"] {
    border-radius: 1rem !important;
    font-weight: 500 !important;
    font-size: 13px !important;
    color: #4a4455 !important;
    background: transparent !important;
    border: none !important;
    padding: 8px 18px !important;
    transition: all 0.25s ease !important;
    letter-spacing: 0 !important;
}

.stTabs [data-baseweb="tab"]:hover {
    background-color: rgba(99, 14, 212, 0.08) !important;
    color: #630ed4 !important;
}

.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important;
    font-weight: 600 !important;
    box-shadow: 0 4px 20px rgba(99, 14, 212, 0.4) !important;
    border: none !important;
}

/* ====== 主按钮 - 紫色渐变胶囊 ====== */
.stButton > button[kind="primary"],
div[data-testid="stMainBlockContainer"] button[kind="primary"] {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 9999px !important;
    padding: 0.6rem 2rem !important;
    font-weight: 600 !important;
    font-size: 14px !important;
    letter-spacing: 0 !important;
    box-shadow: 0 8px 16px -4px rgba(99, 14, 212, 0.35) !important;
    transition: all 0.2s ease !important;
    font-family: 'Inter', sans-serif !important;
}

.stButton > button[kind="primary"]:hover {
    background: linear-gradient(135deg, #5509c6 0%, #6d28d9 100%) !important;
    box-shadow: 0 12px 24px -4px rgba(99, 14, 212, 0.45) !important;
    transform: translateY(-1px) !important;
    opacity: 0.95 !important;
}

.stButton > button[kind="primary"]:active {
    transform: scale(0.98) translateY(0) !important;
}

/* ====== 普通按钮 ====== */
.stButton > button {
    border-radius: 1rem !important;
    font-weight: 500 !important;
    font-size: 13px !important;
    border: none !important;
    background-color: #f3ebfa !important;
    color: #1d1a24 !important;
    padding: 0.5rem 1.2rem !important;
    transition: all 0.2s ease !important;
    font-family: 'Inter', sans-serif !important;
    box-shadow: none !important;
}

.stButton > button:hover {
    background-color: #ede5f4 !important;
    box-shadow: 0 2px 12px rgba(99, 14, 212, 0.12) !important;
}

/* ====== 下载按钮 ====== */
.stDownloadButton > button {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 9999px !important;
    font-weight: 600 !important;
    font-size: 13px !important;
    padding: 0.5rem 1.5rem !important;
    box-shadow: 0 4px 16px rgba(99, 14, 212, 0.3) !important;
    font-family: 'Inter', sans-serif !important;
}
.stDownloadButton > button:hover {
    box-shadow: 0 6px 24px rgba(99, 14, 212, 0.4) !important;
    transform: translateY(-1px) !important;
}

/* ====== Metric 大数字 ====== */
[data-testid="stMetricValue"] {
    color: #630ed4 !important;
    font-weight: 800 !important;
    font-size: 36px !important;
    letter-spacing: -0.04em !important;
    font-family: 'Inter', sans-serif !important;
}
[data-testid="stMetricLabel"] {
    color: #4a4455 !important;
    font-weight: 500 !important;
    font-size: 12px !important;
    letter-spacing: 0.03em !important;
    text-transform: uppercase !important;
}

/* ====== DataFrame 表格 ====== */
.stDataFrame {
    border-radius: 1.5rem !important;
    overflow: hidden !important;
    border: none !important;
    box-shadow: 0 20px 40px -15px rgba(99, 14, 212, 0.06) !important;
    font-family: 'Inter', sans-serif !important;
}

.stDataFrame thead tr th {
    background-color: #f3ebfa !important;
    color: #4a4455 !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    border-bottom: none !important;
    padding: 14px 20px !important;
}

.stDataFrame tbody tr {
    background-color: #ffffff !important;
    transition: background-color 0.15s ease !important;
}
.stDataFrame tbody tr:hover {
    background-color: rgba(249, 241, 255, 0.6) !important;
}
.stDataFrame tbody tr:nth-child(even) {
    background-color: rgba(237, 229, 244, 0.15) !important;
}
.stDataFrame tbody tr:nth-child(even):hover {
    background-color: rgba(249, 241, 255, 0.6) !important;
}

.stDataFrame tbody td {
    color: #1d1a24 !important;
    font-size: 13px !important;
    padding: 12px 20px !important;
    border-bottom: none !important;
    font-family: 'Inter', sans-serif !important;
}

/* 表头圆角 */
.stDataFrame thead th:first-child {
    border-radius: 1.5rem 0 0 0 !important;
}
.stDataFrame thead th:last-child {
    border-radius: 0 1.5rem 0 0 !important;
}

/* ====== 展开面板 ====== */
details {
    background-color: #ffffff !important;
    border-radius: 1.5rem !important;
    border: none !important;
    box-shadow: 0 4px 20px rgba(99, 14, 212, 0.07) !important;
    overflow: hidden !important;
    font-family: 'Inter', sans-serif !important;
}
details summary {
    background-color: #f9f1ff !important;
    color: #1d1a24 !important;
    font-weight: 600 !important;
    font-size: 13px !important;
    border-radius: 1.5rem !important;
    padding: 14px 20px !important;
    border-bottom: none !important;
    letter-spacing: 0 !important;
}
details[open] summary {
    border-radius: 1.5rem 1.5rem 0 0 !important;
}
details > div {
    background-color: #ffffff !important;
    padding: 16px 20px !important;
}

/* ====== Alert ====== */
.stAlert {
    border-radius: 1.5rem !important;
    border: none !important;
    box-shadow: 0 4px 20px rgba(99, 14, 212, 0.08) !important;
    font-family: 'Inter', sans-serif !important;
}

/* ====== Selectbox / Multiselect ====== */
.stSelectbox [data-baseweb="select"] > div,
.stMultiSelect [data-baseweb="select"] > div {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99, 14, 212, 0.1) !important;
}
.stSelectbox label,
.stMultiSelect label {
    color: #4a4455 !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
}

/* ====== Slider 轨道 ====== */
.stSlider [data-baseweb="slider"] .MuiSlider-track {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    border: none !important;
    height: 4px !important;
    border-radius: 9999px !important;
}
.stSlider [data-baseweb="slider"] .MuiSlider-thumb {
    background-color: #630ed4 !important;
    width: 18px !important;
    height: 18px !important;
    box-shadow: 0 2px 10px rgba(99, 14, 212, 0.45) !important;
}
.stSlider [data-baseweb="slider"] .MuiSlider-rail {
    background-color: #e8dfee !important;
    height: 4px !important;
    border-radius: 9999px !important;
    opacity: 1 !important;
}

/* ====== Text Input / Number Input ====== */
.stTextInput input,
.stTextArea textarea,
.stNumberInput input {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99, 14, 212, 0.1) !important;
    font-size: 13px !important;
    color: #1d1a24 !important;
    font-family: 'Inter', sans-serif !important;
    padding: 8px 14px !important;
}
.stTextInput input:focus,
.stTextArea textarea:focus,
.stNumberInput input:focus {
    outline: 2px solid #7c3aed !important;
    outline-offset: 0 !important;
    box-shadow: 0 2px 12px rgba(99, 14, 212, 0.15) !important;
}

/* ====== Spinner ====== */
.stSpinner > div {
    border-color: rgba(99, 14, 212, 0.2) !important;
    border-top-color: #630ed4 !important;
}

/* ====== Progress Bar ====== */
.stProgress > div > div > div {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    border-radius: 9999px !important;
    height: 6px !important;
}

/* ====== Plotly 图表 ====== */
.js-plotly-plot .plotly {
    border-radius: 1.5rem !important;
    overflow: hidden !important;
    box-shadow: 0 4px 24px rgba(99, 14, 212, 0.1) !important;
}
.js-plotly-plot .plotly .modebar {
    background: rgba(254, 247, 255, 0.9) !important;
    border-radius: 1rem !important;
    padding: 4px 8px !important;
}

/* ====== 分隔线隐藏 ====== */
hr { border: none !important; }

/* ====== 全局滚动条 ====== */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: #ccc3d8; border-radius: 9999px; }
::-webkit-scrollbar-thumb:hover { background: #7b7487; }

/* ====== Streamlit block 卡片 ====== */
[data-testid="stVerticalBlock"] [data-testid="stHorizontalBlock"] > div {
    background-color: #ffffff !important;
    border-radius: 1.5rem !important;
    box-shadow: 0 2px 16px rgba(99, 14, 212, 0.06) !important;
    padding: 20px !important;
}

/* ====== 文件上传器 ====== */
[data-testid="stFileUploader"] > div {
    border-radius: 1.5rem !important;
    border: 2px dashed rgba(99, 14, 212, 0.25) !important;
    background-color: rgba(249, 241, 255, 0.5) !important;
    backdrop-filter: blur(8px) !important;
    transition: all 0.3s ease !important;
}
[data-testid="stFileUploader"] > div:hover {
    border-color: rgba(99, 14, 212, 0.5) !important;
    background-color: rgba(249, 241, 255, 0.8) !important;
}
</style>
""", unsafe_allow_html=True)
'''

new_content = content[:start_idx] + new_css + '\n' + content[end_idx:]

with open(r'E:\基因组组装\网络药理学\11\app.py', 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f'Replacement done. New file length: {len(new_content)} chars')
print(f'Replaced {end_idx - start_idx} chars with {len(new_css)} chars')
