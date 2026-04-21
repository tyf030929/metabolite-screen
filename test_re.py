# -*- coding: utf-8 -*-
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

filepath = r'E:\基因组组装\网络药理学\11\app.py'
with open(filepath, encoding='utf-8') as f:
    content = f.read()

# Find all st.markdown blocks with simple string search
idx = 0
while True:
    idx = content.find('st.markdown', idx)
    if idx < 0:
        break
    snippet = content[idx:idx+50]
    snippet_clean = snippet.encode('ascii', errors='replace').decode('ascii')
    print(f'st.markdown at {idx}: {repr(snippet_clean)}')
    idx += 1

print()
# Try actual regex
try:
    matches = list(re.finditer(r'st\.markdown\("""\s*<style>.*?</style>""",\s*unsafe_allow_html=True\)', content, re.DOTALL))
    print(f'Regex matches: {len(matches)}')
except Exception as e:
    print(f'Re regex error: {e}')

# Try simpler regex
try:
    matches2 = list(re.finditer(r'st\.markdown\(\"\"\"', content))
    print(f'Simple regex matches: {len(matches2)}')
    for m in matches2:
        print(f'  at {m.start()}: {repr(content[m.start():m.start()+40])}')
except Exception as e:
    print(f'Simple regex error: {e}')
