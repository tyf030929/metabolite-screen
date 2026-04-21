# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    lines = f.readlines()

print(f'Current lines: {len(lines)}')
# Current structure:
# idx 408: '# ====== MetaboLab Aura Design System CSS ======\n'
# idx 409: '# 参考：...'
# idx 410: '\n'
# idx 411-519: old st.markdown CSS block
# idx 520: '"", unsafe_allow_html=True)\n'
# idx 521: '\n'
# idx 522: plt.rcParams['font.sans-serif']...
# idx 523: plt.rcParams['axes.unicode_minus']...

# We want:
# lines[0:408] + NEW CSS block + lines[522:]

new_css = [
    '# ====== MetaboLab Aura Design System CSS ======\n',
    '# 参考：stitch_high_end_website_prototype - Tailwind MD3 原型设计精确复刻\n',
    '\n',
    'st.markdown("""\n',
    '<style>\n',
    '/* ====== Google Fonts: Inter ====== */\n',
    "@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800;900&family=Noto+Sans+SC:wght@400;500;700&display=swap');\n",
    '\n',
    '/* ====== 全局基础 ====== */\n',
    'html, body { font-family: "Inter", "Noto Sans SC", -apple-system, sans-serif !important; -webkit-font-smoothing: antialiased; letter-spacing: -0.01em; }\n',
    'h1, h2, h3, h4 { letter-spacing: -0.02em !important; }\n',
    '\n',
    '/* ====== 根背景 - 薰衣草白 ====== */\n',
    '.stApp { background-color: #fef7ff !important; }\n',
    '\n',
    '/* ====== 顶部导航栏 - 玻璃拟态 ====== */\n',
    '[data-testid="stHeader"] { background: rgba(254, 247, 255, 0.82) !important; backdrop-filter: blur(24px) !important; -webkit-backdrop-filter: blur(24px) !important; border-bottom: none !important; box-shadow: 0 40px 60px -15px rgba(99, 14, 212, 0.05) !important; }\n',
    '\n',
    '/* ====== 侧边栏 - 薰衣草紫 ====== */\n',
    '[data-testid="stSidebar"] { background-color: #f3ebfa !important; border-right: none !important; }\n',
    '[data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 { color: #630ed4 !important; font-weight: 700 !important; }\n',
    '[data-testid="stSidebar"] hr { border: none !important; border-top: 1px solid rgba(204,195,216,0.4) !important; }\n',
    '[data-testid="stSidebar"] .stMarkdown p, [data-testid="stSidebar"] span, [data-testid="stSidebar"] label { color: #1d1a24 !important; }\n',
    '\n',
    '/* ====== 侧边栏 Selectbox ====== */\n',
    '[data-testid="stSidebar"] [data-baseweb="select"] > div, [data-testid="stSidebar"] [data-baseweb="multiselect"] > div { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; }\n',
    '[data-testid="stSidebar"] .stSelectbox > label, [data-testid="stSidebar"] .stMultiSelect > label, [data-testid="stSidebar"] label { color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.06em !important; text-transform: uppercase !important; }\n',
    '\n',
    '/* ====== 侧边栏 Slider ====== */\n',
    '[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-track { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; height: 4px !important; border-radius: 9999px !important; border: none !important; }\n',
    '[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-thumb { background-color: #630ed4 !important; width: 16px !important; height: 16px !important; box-shadow: 0 2px 8px rgba(99,14,212,0.4) !important; }\n',
    '[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-rail { background-color: #e8dfee !important; height: 4px !important; border-radius: 9999px !important; opacity: 1 !important; }\n',
    '\n',
    '/* ====== Tab 胶囊式 ====== */\n',
    '.stTabs [data-baseweb="tab-list"] { background-color: #f3ebfa !important; border-radius: 1.5rem !important; padding: 6px !important; gap: 4px !important; border: none !important; }\n',
    '.stTabs [data-baseweb="tab"] { border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important; color: #4a4455 !important; background: transparent !important; border: none !important; padding: 8px 18px !important; transition: all 0.25s ease !important; }\n',
    '.stTabs [data-baseweb="tab"]:hover { background-color: rgba(99,14,212,0.08) !important; color: #630ed4 !important; }\n',
    '.stTabs [aria-selected="true"] { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; font-weight: 600 !important; box-shadow: 0 4px 20px rgba(99,14,212,0.4) !important; border: none !important; }\n',
    '\n',
    '/* ====== 主按钮 ====== */\n',
    '.stButton > button[kind="primary"], div[data-testid="stMainBlockContainer"] button[kind="primary"] { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; border: none !important; border-radius: 9999px !important; padding: 0.6rem 2rem !important; font-weight: 600 !important; font-size: 14px !important; box-shadow: 0 8px 16px -4px rgba(99,14,212,0.35) !important; transition: all 0.2s ease !important; }\n',
    '.stButton > button[kind="primary"]:hover { background: linear-gradient(135deg, #5509c6 0%, #6d28d9 100%) !important; box-shadow: 0 12px 24px -4px rgba(99,14,212,0.45) !important; transform: translateY(-1px) !important; opacity: 0.95 !important; }\n',
    '.stButton > button[kind="primary"]:active { transform: scale(0.98) translateY(0) !important; }\n',
    '\n',
    '/* ====== 普通按钮 ====== */\n',
    '.stButton > button { border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important; border: none !important; background-color: #f3ebfa !important; color: #1d1a24 !important; padding: 0.5rem 1.2rem !important; transition: all 0.2s ease !important; }\n',
    '.stButton > button:hover { background-color: #ede5f4 !important; box-shadow: 0 2px 12px rgba(99,14,212,0.12) !important; }\n',
    '\n',
    '/* ====== 下载按钮 ====== */\n',
    '.stDownloadButton > button { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; border: none !important; border-radius: 9999px !important; font-weight: 600 !important; font-size: 13px !important; padding: 0.5rem 1.5rem !important; box-shadow: 0 4px 16px rgba(99,14,212,0.3) !important; }\n',
    '.stDownloadButton > button:hover { box-shadow: 0 6px 24px rgba(99,14,212,0.4) !important; transform: translateY(-1px) !important; }\n',
    '\n',
    '/* ====== Metric ====== */\n',
    '[data-testid="stMetricValue"] { color: #630ed4 !important; font-weight: 800 !important; font-size: 36px !important; letter-spacing: -0.04em !important; }\n',
    '[data-testid="stMetricLabel"] { color: #4a4455 !important; font-weight: 500 !important; font-size: 12px !important; letter-spacing: 0.03em !important; text-transform: uppercase !important; }\n',
    '\n',
    '/* ====== DataFrame ====== */\n',
    '.stDataFrame { border-radius: 1.5rem !important; overflow: hidden !important; border: none !important; box-shadow: 0 20px 40px -15px rgba(99,14,212,0.06) !important; }\n',
    '.stDataFrame thead tr th { background-color: #f3ebfa !important; color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.08em !important; text-transform: uppercase !important; border-bottom: none !important; padding: 14px 20px !important; }\n',
    '.stDataFrame thead th:first-child { border-radius: 1.5rem 0 0 0 !important; }\n',
    '.stDataFrame thead th:last-child { border-radius: 0 1.5rem 0 0 !important; }\n',
    '.stDataFrame tbody tr { background-color: #ffffff !important; transition: background-color 0.15s ease !important; }\n',
    '.stDataFrame tbody tr:hover { background-color: rgba(249,241,255,0.6) !important; }\n',
    '.stDataFrame tbody tr:nth-child(even) { background-color: rgba(237,229,244,0.15) !important; }\n',
    '.stDataFrame tbody tr:nth-child(even):hover { background-color: rgba(249,241,255,0.6) !important; }\n',
    '.stDataFrame tbody td { color: #1d1a24 !important; font-size: 13px !important; padding: 12px 20px !important; border-bottom: none !important; }\n',
    '\n',
    '/* ====== 展开面板 ====== */\n',
    'details { background-color: #ffffff !important; border-radius: 1.5rem !important; border: none !important; box-shadow: 0 4px 20px rgba(99,14,212,0.07) !important; overflow: hidden !important; }\n',
    'details summary { background-color: #f9f1ff !important; color: #1d1a24 !important; font-weight: 600 !important; font-size: 13px !important; border-radius: 1.5rem !important; padding: 14px 20px !important; }\n',
    'details[open] summary { border-radius: 1.5rem 1.5rem 0 0 !important; }\n',
    'details > div { background-color: #ffffff !important; padding: 16px 20px !important; }\n',
    '\n',
    '/* ====== Alert ====== */\n',
    '.stAlert { border-radius: 1.5rem !important; border: none !important; box-shadow: 0 4px 20px rgba(99,14,212,0.08) !important; }\n',
    '\n',
    '/* ====== Selectbox / Multiselect ====== */\n',
    '.stSelectbox [data-baseweb="select"] > div, .stMultiSelect [data-baseweb="select"] > div { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; }\n',
    '.stSelectbox label, .stMultiSelect label { color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.06em !important; text-transform: uppercase !important; }\n',
    '\n',
    '/* ====== Slider ====== */\n',
    '.stSlider [data-baseweb="slider"] .MuiSlider-track { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; height: 4px !important; border-radius: 9999px !important; border: none !important; }\n',
    '.stSlider [data-baseweb="slider"] .MuiSlider-thumb { background-color: #630ed4 !important; width: 18px !important; height: 18px !important; box-shadow: 0 2px 10px rgba(99,14,212,0.45) !important; }\n',
    '.stSlider [data-baseweb="slider"] .MuiSlider-rail { background-color: #e8dfee !important; height: 4px !important; border-radius: 9999px !important; opacity: 1 !important; }\n',
    '\n',
    '/* ====== Text / Number Input ====== */\n',
    '.stTextInput input, .stTextArea textarea, .stNumberInput input { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; font-size: 13px !important; color: #1d1a24 !important; padding: 8px 14px !important; }\n',
    '.stTextInput input:focus, .stTextArea textarea:focus, .stNumberInput input:focus { outline: 2px solid #7c3aed !important; box-shadow: 0 2px 12px rgba(99,14,212,0.15) !important; }\n',
    '\n',
    '/* ====== Spinner ====== */\n',
    '.stSpinner > div { border-color: rgba(99,14,212,0.2) !important; border-top-color: #630ed4 !important; }\n',
    '\n',
    '/* ====== Progress ====== */\n',
    '.stProgress > div > div > div { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; border-radius: 9999px !important; height: 6px !important; }\n',
    '\n',
    '/* ====== Plotly ====== */\n',
    '.js-plotly-plot .plotly { border-radius: 1.5rem !important; overflow: hidden !important; box-shadow: 0 4px 24px rgba(99,14,212,0.1) !important; }\n',
    '.js-plotly-plot .plotly .modebar { background: rgba(254,247,255,0.9) !important; border-radius: 1rem !important; padding: 4px 8px !important; }\n',
    '\n',
    '/* ====== 全局滚动条 ====== */\n',
    '::-webkit-scrollbar { width: 6px; height: 6px; }\n',
    '::-webkit-scrollbar-track { background: transparent; }\n',
    '::-webkit-scrollbar-thumb { background: #ccc3d8; border-radius: 9999px; }\n',
    '::-webkit-scrollbar-thumb:hover { background: #7b7487; }\n',
    '\n',
    '/* ====== Container 卡片 ====== */\n',
    '[data-testid="stVerticalBlock"] [data-testid="stHorizontalBlock"] > div { background-color: #ffffff !important; border-radius: 1.5rem !important; box-shadow: 0 2px 16px rgba(99,14,212,0.06) !important; padding: 20px !important; }\n',
    '\n',
    '/* ====== 文件上传器 ====== */\n',
    '[data-testid="stFileUploader"] > div { border-radius: 1.5rem !important; border: 2px dashed rgba(99,14,212,0.25) !important; background-color: rgba(249,241,255,0.5) !important; backdrop-filter: blur(8px) !important; transition: all 0.3s ease !important; }\n',
    '[data-testid="stFileUploader"] > div:hover { border-color: rgba(99,14,212,0.5) !important; background-color: rgba(249,241,255,0.8) !important; }\n',
    '</style>\n',
    '""", unsafe_allow_html=True)\n',
    '\n',
]

new_lines = list(lines[0:408]) + new_css + list(lines[522:])
result = ''.join(new_lines)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(result)

print(f'Done. New file: {len(result)} chars, {len(new_lines)} lines')

import ast
try:
    ast.parse(result)
    print('Syntax: OK')
except SyntaxError as e:
    print(f'Syntax Error at line {e.lineno}: {e.msg}')
    rl = result.splitlines()
    for i in range(max(0,e.lineno-3), min(len(rl),e.lineno+2)):
        print(f'  L{i+1}: {rl[i][:80]}')

print(f'Has #630ed4: {"630ed4" in result}')
print(f'Has plt.rcParams: {"plt.rcParams" in result}')
print(f'Has def main: {"def main" in result}')
