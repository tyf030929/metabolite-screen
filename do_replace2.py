# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    content = f.read()

# Find the st.markdown CSS block precisely
# It starts at: # ====== comment line, then st.markdown("""<style>
# It ends at: """, unsafe_allow_html=True)

# Step 1: Find start of CSS comment block (the line starting with # ===)
lines = content.split('\n')
css_comment_line_idx = -1
css_st_markdown_start = -1
css_end_literal = -1

for i, line in enumerate(lines):
    if line.startswith('# ===') and 'CSS' in line and i > 300 and i < 600:
        css_comment_line_idx = i
        print(f'CSS comment line {i+1}: {line[:60]}')
    if css_comment_line_idx >= 0 and css_st_markdown_start < 0:
        if 'st.markdown("""' in line and i == css_comment_line_idx + 1:
            css_st_markdown_start = i
            print(f'st.markdown start line {i+1}')
    if css_end_literal < 0 and 'unsafe_allow_html=True' in line and i > css_comment_line_idx:
        # This is the end of the CSS block
        css_end_literal = i
        print(f'CSS end literal line {i+1}: {line[:60]}')
        break

# Build char positions
start_pos = sum(len(l)+1 for l in lines[:css_comment_line_idx])
end_pos = sum(len(l)+1 for l in lines[:css_end_literal+1])

print(f'\nChar range: {start_pos} to {end_pos}, removing {end_pos-start_pos} chars')

old_block = '\n'.join(lines[css_comment_line_idx:cms_end_literal+1])
print(f'Old block ({len(old_block)} chars): {repr(old_block[:200])}')

new_css = '''# ====== MetaboLab Aura Design System CSS ======
# 参考：stitch_high_end_website_prototype - Tailwind MD3 原型设计精确复刻
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

h1, h2, h3, h4 { letter-spacing: -0.02em !important; }

/* ====== 根背景 - 薰衣草白 ====== */
.stApp { background-color: #fef7ff !important; }

/* ====== 顶部导航栏 - 玻璃拟态 ====== */
[data-testid="stHeader"] {
    background: rgba(254, 247, 255, 0.82) !important;
    backdrop-filter: blur(24px) !important;
    -webkit-backdrop-filter: blur(24px) !important;
    border-bottom: none !important;
    box-shadow: 0 40px 60px -15px rgba(99, 14, 212, 0.05) !important;
}

/* ====== 侧边栏 - 薰衣草紫 ====== */
[data-testid="stSidebar"] {
    background-color: #f3edf7 !important;
    border-right: none !important;
}
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 { color: #630ed4 !important; font-weight: 700 !important; }
[data-testid="stSidebar"] hr { border: none !important; border-top: 1px solid rgba(204,195,216,0.4) !important; }
[data-testid="stSidebar"] .stMarkdown p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] label { color: #1d1a24 !important; }

/* 侧边栏 Selectbox */
[data-testid="stSidebar"] [data-baseweb="select"] > div,
[data-testid="stSidebar"] [data-baseweb="multiselect"] > div {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important;
}
[data-testid="stSidebar"] .stSelectbox > label,
[data-testid="stSidebar"] .stMultiSelect > label,
[data-testid="stSidebar"] label {
    color: #4a4455 !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
}

/* 侧边栏 Slider */
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-track {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    height: 4px !important;
    border-radius: 9999px !important;
    border: none !important;
}
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-thumb {
    background-color: #630ed4 !important;
    width: 16px !important;
    height: 16px !important;
    box-shadow: 0 2px 8px rgba(99,14,212,0.4) !important;
}
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-rail {
    background-color: #e8dfee !important;
    height: 4px !important;
    border-radius: 9999px !important;
    opacity: 1 !important;
}

/* ====== Tab 标签 - 胶囊式 ====== */
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
}
.stTabs [data-baseweb="tab"]:hover {
    background-color: rgba(99,14,212,0.08) !important;
    color: #630ed4 !important;
}
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important;
    font-weight: 600 !important;
    box-shadow: 0 4px 20px rgba(99,14,212,0.4) !important;
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
    box-shadow: 0 8px 16px -4px rgba(99,14,212,0.35) !important;
    transition: all 0.2s ease !important;
    font-family: 'Inter', sans-serif !important;
}
.stButton > button[kind="primary"]:hover {
    background: linear-gradient(135deg, #5509c6 0%, #6d28d9 100%) !important;
    box-shadow: 0 12px 24px -4px rgba(99,14,212,0.45) !important;
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
    box-shadow: none !important;
}
.stButton > button:hover {
    background-color: #ede5f4 !important;
    box-shadow: 0 2px 12px rgba(99,14,212,0.12) !important;
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
    box-shadow: 0 4px 16px rgba(99,14,212,0.3) !important;
}
.stDownloadButton > button:hover {
    box-shadow: 0 6px 24px rgba(99,14,212,0.4) !important;
    transform: translateY(-1px) !important;
}

/* ====== Metric ====== */
[data-testid="stMetricValue"] {
    color: #630ed4 !important;
    font-weight: 800 !important;
    font-size: 36px !important;
    letter-spacing: -0.04em !important;
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
    box-shadow: 0 20px 40px -15px rgba(99,14,212,0.06) !important;
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
.stDataFrame thead th:first-child { border-radius: 1.5rem 0 0 0 !important; }
.stDataFrame thead th:last-child { border-radius: 0 1.5rem 0 0 !important; }
.stDataFrame tbody tr {
    background-color: #ffffff !important;
    transition: background-color 0.15s ease !important;
}
.stDataFrame tbody tr:hover { background-color: rgba(249,241,255,0.6) !important; }
.stDataFrame tbody tr:nth-child(even) { background-color: rgba(237,229,244,0.15) !important; }
.stDataFrame tbody tr:nth-child(even):hover { background-color: rgba(249,241,255,0.6) !important; }
.stDataFrame tbody td {
    color: #1d1a24 !important;
    font-size: 13px !important;
    padding: 12px 20px !important;
    border-bottom: none !important;
}

/* ====== 展开面板 ====== */
details {
    background-color: #ffffff !important;
    border-radius: 1.5rem !important;
    border: none !important;
    box-shadow: 0 4px 20px rgba(99,14,212,0.07) !important;
    overflow: hidden !important;
}
details summary {
    background-color: #f9f1ff !important;
    color: #1d1a24 !important;
    font-weight: 600 !important;
    font-size: 13px !important;
    border-radius: 1.5rem !important;
    padding: 14px 20px !important;
}
details[open] summary { border-radius: 1.5rem 1.5rem 0 0 !important; }
details > div { background-color: #ffffff !important; padding: 16px 20px !important; }

/* ====== Alert ====== */
.stAlert {
    border-radius: 1.5rem !important;
    border: none !important;
    box-shadow: 0 4px 20px rgba(99,14,212,0.08) !important;
}

/* ====== Selectbox / Multiselect ====== */
.stSelectbox [data-baseweb="select"] > div,
.stMultiSelect [data-baseweb="select"] > div {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important;
}
.stSelectbox label, .stMultiSelect label {
    color: #4a4455 !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
}

/* ====== Slider ====== */
.stSlider [data-baseweb="slider"] .MuiSlider-track {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    height: 4px !important;
    border-radius: 9999px !important;
    border: none !important;
}
.stSlider [data-baseweb="slider"] .MuiSlider-thumb {
    background-color: #630ed4 !important;
    width: 18px !important;
    height: 18px !important;
    box-shadow: 0 2px 10px rgba(99,14,212,0.45) !important;
}
.stSlider [data-baseweb="slider"] .MuiSlider-rail {
    background-color: #e8dfee !important;
    height: 4px !important;
    border-radius: 9999px !important;
    opacity: 1 !important;
}

/* ====== Text / Number Input ====== */
.stTextInput input, .stTextArea textarea, .stNumberInput input {
    background-color: #ffffff !important;
    border: none !important;
    border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important;
    font-size: 13px !important;
    color: #1d1a24 !important;
    padding: 8px 14px !important;
}
.stTextInput input:focus, .stTextArea textarea:focus, .stNumberInput input:focus {
    outline: 2px solid #7c3aed !important;
    box-shadow: 0 2px 12px rgba(99,14,212,0.15) !important;
}

/* ====== Spinner ====== */
.stSpinner > div { border-color: rgba(99,14,212,0.2) !important; border-top-color: #630ed4 !important; }

/* ====== Progress ====== */
.stProgress > div > div > div {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important;
    border-radius: 9999px !important;
    height: 6px !important;
}

/* ====== Plotly 图表 ====== */
.js-plotly-plot .plotly {
    border-radius: 1.5rem !important;
    overflow: hidden !important;
    box-shadow: 0 4px 24px rgba(99,14,212,0.1) !important;
}
.js-plotly-plot .plotly .modebar {
    background: rgba(254,247,255,0.9) !important;
    border-radius: 1rem !important;
    padding: 4px 8px !important;
}

/* ====== 全局滚动条 ====== */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: #ccc3d8; border-radius: 9999px; }
::-webkit-scrollbar-thumb:hover { background: #7b7487; }

/* ====== Streamlit container 卡片 ====== */
[data-testid="stVerticalBlock"] [data-testid="stHorizontalBlock"] > div {
    background-color: #ffffff !important;
    border-radius: 1.5rem !important;
    box-shadow: 0 2px 16px rgba(99,14,212,0.06) !important;
    padding: 20px !important;
}

/* ====== 文件上传器 ====== */
[data-testid="stFileUploader"] > div {
    border-radius: 1.5rem !important;
    border: 2px dashed rgba(99,14,212,0.25) !important;
    background-color: rgba(249,241,255,0.5) !important;
    backdrop-filter: blur(8px) !important;
    transition: all 0.3s ease !important;
}
[data-testid="stFileUploader"] > div:hover {
    border-color: rgba(99,14,212,0.5) !important;
    background-color: rgba(249,241,255,0.8) !important;
}
</style>
""", unsafe_allow_html=True)

'''

new_content = content[:start_pos] + new_css + '\n' + content[end_pos:]

with open(r'E:\基因组组装\网络药理学\11\app.py', 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f'\nDone. New file: {len(new_content)} chars')
print(f'Replaced {end_pos - start_pos} chars with {len(new_css)} chars')
