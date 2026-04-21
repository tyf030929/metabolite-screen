# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    app_content = f.read()
with open(r'E:\基因组组装\网络药理学\11\fix_css.py', encoding='utf-8') as f:
    new_css = f.read()

# Find the old CSS block: starts at 'st.markdown("""\n<style>', ends before '# 中文字体支持'
css_start = app_content.find('st.markdown("""\n<style>')
if css_start < 0:
    print('ERROR: old CSS block not found')
    exit(1)

# Find the comment line just before the CSS block
comment_start = app_content.rfind('\n# ===', 0, css_start)
if comment_start < 0:
    comment_start = app_content.rfind('\n# =====', 0, css_start)
print(f'CSS block start (comment): char {comment_start}')

# Find end: unsafe_allow_html=True followed by the marker
end_marker = '# \u4e2d\u6587\u5b57\u4f53\u652f\u6301'
end_pos = app_content.find('\n' + end_marker, css_start)
if end_pos < 0:
    print('ERROR: end marker not found')
    exit(1)
end_pos = end_pos + 1  # include the newline before marker
print(f'CSS block end: char {end_pos}, removing {end_pos - comment_start} chars')

# Build new content
new_content = app_content[:comment_start] + '\n' + new_css + '\n' + app_content[end_pos:]

with open(r'E:\基因组组装\网络药理学\11\app.py', 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f'Done. New file: {len(new_content)} chars')
# Quick verify
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    v = f.read()
count = v.count('# ====== MetaboLab Aura')
print(f'MetaboLab Aura markers: {count}')
print(f'Has # 中文字体支持: {"# 中文字体支持" in v}')
