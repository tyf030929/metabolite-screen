# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    content = f.read()

# Find the CSS block: starts at st.markdown("""\n<style> and ends before '# 中文字体支持'
# Use a line-based approach
lines = content.split('\n')
new_lines = []
in_css_block = False
css_block_start_found = False

for i, line in enumerate(lines):
    if not in_css_block and 'st.markdown("""' in line and i > 300:
        # Start of CSS block - skip until end marker
        in_css_block = True
        css_block_start_found = True
        # Add the new CSS comment line as first line
        new_lines.append('# ====== MetaboLab Aura Design System CSS ======')
        new_lines.append('# 参考：stitch_high_end_website_prototype - Tailwind MD3 原型设计')
        new_lines.append('')
        continue
    if in_css_block:
        # Continue skipping until we find the end marker
        if '# 中文字体支持' in line:
            in_css_block = False
            new_lines.append(line)  # Add the end marker line
            continue
        else:
            continue  # Skip lines in CSS block
    new_lines.append(line)

result = '\n'.join(new_lines)
with open(r'E:\基因组组装\网络药理学\11\app.py', 'w', encoding='utf-8') as f:
    f.write(result)

print(f'Done. New file: {len(result)} chars, {len(result.split(chr(10)))} lines')
# Verify st.markdown CSS is gone
if 'st.markdown("""' in result:
    print('WARNING: st.markdown still in file')
else:
    print('st.markdown block removed OK')
if '# 中文字体支持' in result:
    print('# 中文字体支持 still present OK')
