# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    lines = f.readlines()

print(f'Total lines: {len(lines)}')
# Show lines 406-420
for i in range(406, 420):
    print(f'L{i+1}: {repr(lines[i][:80])}')
print('...')
for i in range(775, 785):
    print(f'L{i+1}: {repr(lines[i][:80])}')
