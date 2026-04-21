# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    content = f.read()

# Find plt.rcParams (the FIRST occurrence, which comes after our CSS)
plt_pos = content.find("plt.rcParams['font.sans-serif']")
print(f'plt.rcParams at: {plt_pos}')

# Find the '# ====== MetaboLab Aura' comment line that starts our CSS block
meta_start = content.find('# ====== MetaboLab Aura Design System CSS')
print(f'MetaboLab Aura at: {meta_start}')

# The new CSS block content we want (from the current file's MetaboLab block to unsafe_allow_html)
# But we need to find where the st.markdown STARTS (after the comment)
# It starts with 'st.markdown("""\n'
sm_start = content.find('st.markdown("""\n', meta_start)
print(f'st.markdown for our CSS at: {sm_start}')

# Find end of our CSS block
ua_end = content.find('\nunsafe_allow_html=True)', sm_start)
print(f'unsafe_allow_html at: {ua_end}')

# The block we want to replace: from '# ====== MetaboLab Aura' to 'unsafe_allow_html=True)\n'
old_css_start = meta_start
old_css_end = ua_end + len('\nunsafe_allow_html=True)\n')
print(f'Old CSS block: {old_css_start} to {old_css_end}, length={old_css_end-old_css_start}')

# The OLD Aura CSS block also needs to be removed
# Find '# ====================== Aura Scientific' (the OLD CSS block header)
aura_old_start = content.find('# ====================== Aura Scientific')
print(f'Old Aura CSS at: {aura_old_start}')

# Find where this old block ENDS (unsafe_allow_html=True) - need to find the FIRST one after aura_old_start
# But there might be multiple. Let's find ALL st.markdown blocks
import re
all_starts = [m.start() for m in re.finditer(r'st\.markdown\("""\s*<style>', content)]
print(f'All st.markdown<style> starts: {all_starts}')

# The first one is the old CSS block
old_sm_start = all_starts[0]
# The second one is the new CSS block
new_sm_start = all_starts[1] if len(all_starts) > 1 else -1
print(f'Old st.markdown at: {old_sm_start}, New st.markdown at: {new_sm_start}')

# Find end of old CSS block
old_ua = content.find('\nunsafe_allow_html=True)', old_sm_start)
print(f'Old unsafe_allow_html at: {old_ua}')
old_css_end = old_ua + len('\nunsafe_allow_html=True)\n')

# Find where plt.rcParams is (comes after both CSS blocks in the clean file)
plt_pos = content.find("plt.rcParams['font.sans-serif']")
print(f'plt.rcParams at: {plt_pos}')

# New def main() function to append
new_def_main = '''

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

# Build new content:
# Keep: content[0:plt_pos] (everything up to plt.rcParams['font.sans-serif'])
# Add: plt.rcParams + '\n\n' + new_def_main + '\n\n' + content[old_sm_start:old_css_end] (new CSS block)
# Wait, actually new_def_main ALREADY has the brand banner but we need to include the NEW correct CSS block
# The new CSS block is: content[old_sm_start:old_css_end] (from the SECOND st.markdown)
# And we want: plt.rcParams + new CSS + new_def_main

# But wait - new_def_main has the brand banner STYLED with inline HTML
# We also need the st.markdown CSS from lines 408-521
# That's content[meta_start:old_css_end]

# Correct plan:
# 1. content[0:plt_pos] = up to plt.rcParams['font.sans-serif']
# 2. plt.rcParams code
# 3. blank
# 4. New CSS block (the SECOND st.markdown, from old_sm_start to old_css_end)
# 5. blank
# 6. New def main() with brand banner HTML

# Find plt.rcParams['font.sans-serif'] line start
plt_line_start = content.rfind('\n', 0, plt_pos) + 1
print(f'plt.rcParams line starts at: {plt_line_start}')

# NEW CSS block: content[old_sm_start:old_css_end]
new_css_block = content[old_sm_start:old_css_end]

# New content
new_content = (
    content[:plt_line_start] +  # up to plt.rcParams line
    '\n\n' +
    new_css_block +
    new_def_main
)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f'\nDone. New file: {len(new_content)} chars')

# Verify
with open(filepath, encoding='utf-8') as f:
    verify = f.read()
import ast
try:
    ast.parse(verify)
    print('Syntax: OK')
except SyntaxError as e:
    print(f'Syntax Error at line {e.lineno}: {e.msg}')
    vl = verify.splitlines()
    for i in range(max(0,e.lineno-3), min(len(vl),e.lineno+2)):
        print(f'  L{i+1}: {vl[i][:80]}')
print(f'File size: {len(verify)} chars')
print(f'Lines: {len(verify.splitlines())}')
