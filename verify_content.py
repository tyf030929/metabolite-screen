# -*- coding: utf-8 -*-
import sys, io, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    content = f.read()

# Find st.markdown occurrences
idx = 0
count = 0
while True:
    idx = content.find('st.markdown', idx)
    if idx < 0:
        break
    print(f'st.markdown #{count+1} at {idx}: {repr(content[idx:idx+30])}')
    idx += 1
    count += 1
    if count > 5:
        break

print()

# Try simple pattern
matches = list(re.finditer(r'stm\.markdown', content))
print(f'stm.matches: {len(matches)}')

# Try the triple-quote pattern
m = re.search(r'st\.markdow[n]\("""', content)
print(f'With [n]: {m}')

# Check what 'st.markdown' actually looks like in the file
for test in [r'st.markdown', r'st\.markdown', 'st.markdown']:
    idx = content.find(test)
    if idx >= 0:
        print(f'Found {test!r} at {idx}: {repr(content[idx:idx+40])}')
        break
