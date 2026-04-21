# -*- coding: utf-8 -*-
import sys, io, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    content = f.read()

# Find ALL st.markdown("<style> blocks
pattern = r'st\.markdown\("""\s*<style>.*?</style>""",\s*unsafe_allow_html=True\)'
matches = list(re.finditer(pattern, content, re.DOTALL))
print(f'Found {len(matches)} CSS blocks')
for i, m in enumerate(matches):
    print(f'  Block {i+1}: chars {m.start()}-{m.end()}, len={m.end()-m.start()}')
    print(f'  Preview: {repr(m.group()[:80])}')

if len(matches) < 2:
    print('ERROR: need at least 2 blocks')
    exit(1)

# Block 1 = OLD CSS (to remove)
# Block 2 = NEW CSS (to keep)
old_css = matches[0].group()
old_css_start = matches[0].start()
old_css_end = matches[0].end()

# Find plt.rcParams position
plt_pos = content.find("plt.rcParams['font.sans-serif']")
if plt_pos < 0:
    print('ERROR: plt.rcParams not found')
    exit(1)
print(f'plt.rcParams at: {plt_pos}')

# Find def main() - first occurrence
main_match = re.search(r'\ndef main\(\):', content)
if main_match:
    main_pos = main_match.start()
    print(f'def main at: {main_pos}')
else:
    print('ERROR: def main not found')
    exit(1)

# New brand banner HTML for def main()
brand_banner = '''
    # ===== 品牌头部横幅 =====
    st.markdown("""
    <div style="background: linear-gradient(135deg, #630ed4 0%, #7c3aed 50%, #9f67f5 100%);
                padding: 24px 32px; border-radius: 1.5rem; margin-bottom: 24px;
                box-shadow: 0 8px 32px rgba(99, 14, 212, 0.3);">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div style="display: flex; align-items: center; gap: 16px;">
                <div style="font-size: 36px;">&#129704;</div>
                <div>
                    <h1 style="color: #ffffff; margin: 0; font-size: 22px; font-weight: 700; letter-spacing: -0.03em;">
                        &#24046;&#24322;&#20195;&#35874;&#29289;&#33647;&#29992;&#31579;&#36873;&#24179;&#21488;
                    </h1>
                    <p style="color: rgba(255,255,255,0.82); margin: 6px 0 0 0; font-size: 13px;">
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
# [0 : old_css_start]  - everything before OLD CSS
# [old_css_end : plt_pos] - plt.rcParams area (old CSS block removed!)
# plt.rcParams code
# [plt_pos] onwards - everything from plt.rcParams to end
# BUT we need to replace def main() brand banner part

# Strategy: 
# Part 1: content[0 : old_css_start]  (before old CSS)
# Part 2: content[plt_pos : main_pos]  (plt.rcParams area + everything up to def main)
# Part 3: def main() with NEW brand banner
# Part 4: content[main_pos_end : ]     (rest after "def main():" including the NEW CSS block if needed)

# Actually simpler:
# 1. Remove OLD CSS block completely
# 2. The rest (plt.rcParams + def main + NEW CSS block) stays as-is
# 3. Update def main()'s brand banner

# Find end of "def main():" in the new content
# Look for the brand banner line in def main()
banner_in_main = content.find('    # ===== 品牌头部横幅 =====', main_pos)
print(f'Brand banner in main at: {banner_in_main}')

# Find end of def main() brand banner (the st.markdown that follows)
# It starts with '    st.markdown("""'
banner_sm_start = content.find('    st.markdown("""', banner_in_main)
banner_sm_end = content.find('"", unsafe_allow_html=True)', banner_sm_start)
banner_sm_end = banner_sm_end + len('"", unsafe_allow_html=True)')
print(f'Brand banner st.markdown: {banner_sm_start} to {banner_sm_end}')

# Now build new content:
# Part 1: content[0 : old_css_start] 
# Part 2: content[old_css_end : banner_sm_start]  (everything after old CSS up to old brand banner)
# Part 3: brand_banner (new banner)
# Part 4: content[banner_sm_end : ] (rest of file)

new_content = (
    content[0 : old_css_start] +
    content[old_css_end : banner_sm_start] +
    brand_banner +
    content[banner_sm_end : ]
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
print(f'File: {len(verify)} chars, {len(verify.splitlines())} lines')
print(f'Old CSS removed: {"#630ed4" not in verify[:old_css_start]}')
