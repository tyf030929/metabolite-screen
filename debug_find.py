# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    content = f.read()

# Aura Scientific starts at 17852, let's look around it
aura = 'Aura Scientific'
idx = content.find(aura)
print(f'Aura Scientific at: {idx}')

# Find start of CSS comment block
# Go back to find a line starting with # =====
search_start = max(0, idx - 500)
hash_line = content.rfind('\n# =====', search_start, idx)
print(f'CSS comment line starts at: {hash_line}')

# The comment line itself
# Find the start of the line
line_start = content.rfind('\n', 0, hash_line) + 1
print(f'CSS block full start at: {line_start}')
print(f'Content at start: {repr(content[line_start:line_start+80])}')

# Find end marker '# 中文字体支持'
end_marker = '# 中文字体支持'
end_idx = content.find(end_marker)
print(f'End marker at: {end_idx}')

# The end of the old CSS block - just before '# 中文字体支持'
old_block = content[line_start:end_idx]
print(f'Old block length: {len(old_block)} chars')
print(f'Old block end: {repr(old_block[-100:])}')

# Now build new CSS block  
new_css = content  # placeholder
