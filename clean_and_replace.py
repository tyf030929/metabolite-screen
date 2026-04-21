# -*- coding: utf-8 -*-
with open(r'E:\基因组组装\网络药理学\11\app.py', encoding='utf-8') as f:
    content = f.read()

# The file has:
# Line ~409: old Aura Scientific CSS comment (STAYS)
# Line ~412: '# ====== MetaboLab Aura Design...' (REMOVE - duplicate inserted by do_replace.py)
# Line ~413: '# 参考...' (REMOVE)
# Line ~415: '# 中文字体支持' (REMOVE)
# After that: plt.rcParams (STAYS)
# But ALSO the st.markdown block from line ~409 through ~778 is STILL THERE

# We need to:
# 1. Remove the duplicate lines 412-417 (the 'MetaboLab Aura' block)
# 2. Replace the actual old st.markdown CSS block with new CSS

# Find the old st.markdown block
sm_start = content.find('st.markdown(' + chr(34) + chr(34) + chr(34) + '\n<style>')
print(f'st.markdown start: {sm_start}')
if sm_start < 0:
    # Try alternative
    sm_start = content.find('st.markdown(' + '"\\""\\"')
    print(f'alt st.markdown start: {sm_start}')
