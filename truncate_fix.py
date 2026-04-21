# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    lines = f.readlines()

print(f'Total lines in current file: {len(lines)}')

# The file structure:
# Lines ~0-906: good code (ends with "return fig" from radar chart)
# Lines ~907-910: blank + "# ====== Streamlit UI ======" + "def main():"
# Lines ~911+: DUPLICATE corrupt content (brand banner + duplicate CSS)
# We want to keep lines 0-908 (up to "return fig"), then add the NEW correct def main()

# Find where "return fig" (the last good one) is
for i in range(len(lines)-1, -1, -1):
    if 'return fig' in lines[i] and i > 900:
        print(f'Last "return fig" at line {i+1}: {repr(lines[i][:50])}')
        last_return_fig = i
        break

# Find where "# ====== Streamlit UI ======" appears after the pharma rules (second occurrence)
ui_marker_count = 0
streamlit_ui_line = -1
for i, line in enumerate(lines):
    if '# ====== Streamlit UI ======' in line:
        ui_marker_count += 1
        print(f'Streamlit UI marker #{ui_marker_count} at line {i+1}')
        if ui_marker_count == 2:
            streamlit_ui_line = i
            break

print(f'Streamlit UI marker 2 at line: {streamlit_ui_line+1}')

# New clean def main() with correct CSS
new_main = '''
# ====================== Streamlit UI ======================

def main():
    # ===== 品牌头部横幅 =====
    st.markdown("""
    <div style="background: linear-gradient(135deg, #630ed4 0%, #7c3aed 50%, #9f67f5 100%);
                padding: 24px 32px; border-radius: 1.5rem; margin-bottom: 24px;
                box-shadow: 0 8px 32px rgba(99, 14, 212, 0.3);">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div style="display: flex; align-items: center; gap: 16px;">
                <div style="font-size: 36px; filter: drop-shadow(0 2px 8px rgba(0,0,0,0.2));">&#129704;</div>
                <div>
                    <h1 style="color: #ffffff; margin: 0; font-size: 22px; font-weight: 700; letter-spacing: -0.03em; text-shadow: 0 1px 4px rgba(0,0,0,0.15);">
                        &#24046;&#24322;&#20195;&#35874;&#29289;&#33647;&#29992;&#31579;&#36873;&#24179;&#21488;
                    </h1>
                    <p style="color: rgba(255,255,255,0.82); margin: 6px 0 0 0; font-size: 13px; letter-spacing: 0.01em;">
                        Multivariate Stats &middot; Volcano Plot &middot; KEGG Enrichment &middot; Network Pharmacology
                    </p>
                </div>
            </div>
            <div style="background: rgba(255,255,255,0.15); border-radius: 1rem; padding: 8px 16px; backdrop-filter: blur(8px);">
                <span style="color: rgba(255,255,255,0.9); font-size: 12px; font-weight: 600; letter-spacing: 0.05em; text-transform: uppercase;">MetaboLab</span>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

'''

# Build new file: lines 0 to last_return_fig (inclusive) + new_main
new_lines = lines[:last_return_fig + 1]
new_lines.append(new_main)

result = ''.join(new_lines)
with open(filepath, 'w', encoding='utf-8') as f:
    f.write(result)

print(f'\nDone. New file: {len(result)} chars, {len(new_lines)} lines')

# Verify
with open(filepath, encoding='utf-8') as f:
    verify = f.read()
import ast
try:
    ast.parse(verify)
    print('Syntax: OK')
except SyntaxError as e:
    print(f'Syntax Error at line {e.lineno}: {e.msg}')
    # Show context
    vl = verify.splitlines()
    for i in range(max(0, e.lineno-3), min(len(vl), e.lineno+2)):
        print(f'  L{i+1}: {vl[i][:80]}')
