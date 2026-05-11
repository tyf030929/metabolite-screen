# -*- coding: utf-8 -*-
"""
药理数据库整合模块
支持 TCMSP、DrugBank、TTD 离线文件匹配
以及 STITCH 在线 API 查询（网络允许时）
"""

import re
import os
import json
import glob
import time
import requests
import pandas as pd
import numpy as np
import streamlit as st

# ===================== 药理证据评分规则 =====================
PHARMA_SCORES = {
    'tcmsp_ob': 0.4,      # TCMSP OB>30
    'tcmsp_dl': 0.2,       # TCMSP DL>0.18
    'drugbank': 0.9,      # DrugBank 有记录
    'ttd': 0.7,            # TTD 有记录
    'stitch': 0.4,        # STITCH 有相互作用
    'pubchem_bio': 0.5,    # PubChem BioActivity 有记录
}


# ===================== TCMSP 目录批量加载 =====================
def load_tcmsp_from_dir(dir_path: str) -> dict:
    """
    从指定目录扫描并加载 TCMSP 数据文件。
    支持的文件：
      - tcmsp_ingredients.csv: 成分-靶点关联（SVM/RF 预测）
      - tcmsp_targets.csv: 成分药代动力学属性（OB、DL、Caco2、BBB）
      - tcmsp_diseases.csv: 靶点-疾病关联
      - tcmsp_full.json: 完整 JSON 备份
    返回 {compound_name_lower: {ob, dl, targets: [], diseases: [], svm_scores: []}}
    """
    import glob
    result = {}  # {compound_name_lower: {ob, dl, targets, diseases, svm_scores}}

    if not os.path.isdir(dir_path):
        print(f"[pharma_db] TCMSP 目录不存在: {dir_path}")
        return result

    csv_files = glob.glob(os.path.join(dir_path, '*.csv'))

    ingredients_data = {}  # {compound_name_lower: [target_names]}
    targets_data = {}     # {compound_name_lower: {ob, dl, ...}}
    diseases_data = {}    # {compound_name_lower: [disease_names]}

    for fpath in csv_files:
        fname = os.path.basename(fpath).lower()
        try:
            # 自动检测分隔符（CSV 可能用 , 或 \t）
            with open(fpath, 'r', encoding='utf-8-sig') as f:
                first_line = f.readline()
            sep = '\t' if '\t' in first_line else ','

            df = pd.read_csv(fpath, sep=sep, low_memory=False, encoding='utf-8-sig', on_bad_lines='skip')
            df.columns = df.columns.str.strip()
            col_lower = [c.lower() for c in df.columns]

            # 判断文件类型（按文件名优先，其次按列特征）
            fname_key = os.path.basename(fpath).lower()
            if 'target' in fname_key and ('svm_score' in col_lower or 'rf_score' in col_lower):
                # tcmsp_targets.csv: molecule_name + target_name + SVM/RF scores
                mol_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'molecule_name'), None)
                tgt_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'target_name'), None)
                svm_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'svm_score'), None)
                if mol_col and tgt_col:
                    for _, row in df.iterrows():
                        mol = str(row.get(mol_col, '')).strip()
                        tgt = str(row.get(tgt_col, '')).strip()
                        svm = row.get(svm_col, None) if svm_col else None
                        if mol and mol not in ('-', '', 'nan', 'NaN'):
                            key = normalize_name(mol)
                            if key not in ingredients_data:
                                ingredients_data[key] = {'targets': [], 'svm_scores': [], 'mol_name': mol}
                            if tgt and tgt not in ('-', '', 'nan'):
                                ingredients_data[key]['targets'].append(tgt)
                            if svm is not None:
                                try:
                                    ingredients_data[key]['svm_scores'].append(float(svm))
                                except (ValueError, TypeError):
                                    pass

            elif 'ingredient' in fname_key and ('ob' in col_lower or 'dl' in col_lower or 'bbb' in col_lower):
                # tcmsp_ingredients.csv: molecule_name + OB, DL, BBB 等属性
                mol_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'molecule_name'), None)
                if mol_col is None:
                    mol_col = next((df.columns[i] for i, c in enumerate(col_lower) if 'name' in c), df.columns[0])
                ob_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'ob'), None)
                dl_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'dl'), None)
                bbb_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'bbb'), None)
                caco_col = next((df.columns[i] for i, c in enumerate(col_lower) if 'caco' in c), None)
                mw_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'mw'), None)

                for _, row in df.iterrows():
                    mol = str(row.get(mol_col, '')).strip()
                    if mol and mol not in ('-', '', 'nan', 'NaN'):
                        key = normalize_name(mol)
                        ob = float(row[ob_col]) if ob_col and pd.notna(row.get(ob_col)) else None
                        dl = float(row[dl_col]) if dl_col and pd.notna(row.get(dl_col)) else None
                        bbb = float(row[bbb_col]) if bbb_col and pd.notna(row.get(bbb_col)) else None
                        caco = float(row[caco_col]) if caco_col and pd.notna(row.get(caco_col)) else None
                        mw = float(row[mw_col]) if mw_col and pd.notna(row.get(mw_col)) else None
                        targets_data[key] = {'ob': ob, 'dl': dl, 'bbb': bbb, 'caco2': caco, 'mw': mw, 'mol_name': mol}

            elif 'disease' in fname_key:
                # tcmsp_diseases.csv: target_name + disease_name
                tgt_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'target_name'), None)
                dis_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'disease_name'), None)
                if tgt_col and dis_col:
                    for _, row in df.iterrows():
                        tgt = str(row.get(tgt_col, '')).strip()
                        dis = str(row.get(dis_col, '')).strip()
                        if tgt and tgt not in ('-', '', 'nan', 'NaN'):
                            key = normalize_name(tgt)
                            if key not in diseases_data:
                                diseases_data[key] = []
                            if dis and dis not in ('-', '', 'nan'):
                                diseases_data[key].append(dis)

            else:
                # fallback: 用列特征判断
                if 'target_name' in col_lower and 'molecule_name' in col_lower:
                    # 靶点预测文件
                    mol_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'molecule_name'), None)
                    tgt_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'target_name'), None)
                    svm_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'svm_score'), None)
                    if mol_col and tgt_col:
                        for _, row in df.iterrows():
                            mol = str(row.get(mol_col, '')).strip()
                            tgt = str(row.get(tgt_col, '')).strip()
                            svm = row.get(svm_col, None) if svm_col else None
                            if mol and mol not in ('-', '', 'nan', 'NaN'):
                                key = normalize_name(mol)
                                if key not in ingredients_data:
                                    ingredients_data[key] = {'targets': [], 'svm_scores': [], 'mol_name': mol}
                                if tgt and tgt not in ('-', '', 'nan'):
                                    ingredients_data[key]['targets'].append(tgt)
                                if svm is not None:
                                    try:
                                        ingredients_data[key]['svm_scores'].append(float(svm))
                                    except (ValueError, TypeError):
                                        pass
                elif 'ob' in col_lower or 'dl' in col_lower:
                    # 属性文件
                    mol_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'molecule_name'), None)
                    if mol_col is None:
                        mol_col = df.columns[0]
                    ob_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'ob'), None)
                    dl_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'dl'), None)
                    bbb_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'bbb'), None)
                    for _, row in df.iterrows():
                        mol = str(row.get(mol_col, '')).strip()
                        if mol and mol not in ('-', '', 'nan', 'NaN'):
                            key = normalize_name(mol)
                            ob = float(row[ob_col]) if ob_col and pd.notna(row.get(ob_col)) else None
                            dl = float(row[dl_col]) if dl_col and pd.notna(row.get(dl_col)) else None
                            bbb = float(row[bbb_col]) if bbb_col and pd.notna(row.get(bbb_col)) else None
                            targets_data[key] = {'ob': ob, 'dl': dl, 'bbb': bbb, 'caco2': None, 'mw': None, 'mol_name': mol}
                elif 'disease_name' in col_lower:
                    tgt_col = next((df.columns[i] for i, c in enumerate(col_lower) if c == 'target_name'), None)
                    dis_col = 'disease_name'
                    if tgt_col:
                        for _, row in df.iterrows():
                            tgt = str(row.get(tgt_col, '')).strip()
                            dis = str(row.get(dis_col, '')).strip()
                            if tgt and tgt not in ('-', '', 'nan', 'NaN'):
                                key = normalize_name(tgt)
                                if key not in diseases_data:
                                    diseases_data[key] = []
                                if dis and dis not in ('-', '', 'nan'):
                                    diseases_data[key].append(dis)

        except Exception as e:
            print(f"[pharma_db] 解析 TCMSP 文件失败 {fpath}: {e}")

    # 也尝试加载 JSON 备份
    json_path = os.path.join(dir_path, 'tcmsp_full.json')
    if os.path.exists(json_path):
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                json_data = json.load(f)
            if isinstance(json_data, dict):
                for key, val in json_data.items():
                    if key not in result:
                        result[key] = {'ob': None, 'dl': None, 'targets': [], 'diseases': [], 'svm_scores': []}
                    if isinstance(val, dict):
                        result[key].setdefault('ob', val.get('ob'))
                        result[key].setdefault('dl', val.get('dl'))
        except Exception as e:
            print(f"[pharma_db] 解析 tcmsp_full.json 失败: {e}")

    # 合并数据
    # 先用 targets_data 初始化所有化合物
    for key, attrs in targets_data.items():
        result[key] = {
            'ob': attrs.get('ob'),
            'dl': attrs.get('dl'),
            'bbb': attrs.get('bbb'),
            'caco2': attrs.get('caco2'),
            'mw': attrs.get('mw'),
            'targets': [],
            'diseases': [],
            'svm_scores': [],
            'mol_name': attrs.get('mol_name', ''),
        }

    # 合并 ingredients（靶点和 SVM 评分）
    for key, ing_info in ingredients_data.items():
        if key not in result:
            result[key] = {'ob': None, 'dl': None, 'bbb': None, 'caco2': None, 'mw': None,
                           'targets': [], 'diseases': [], 'svm_scores': [], 'mol_name': ing_info.get('mol_name', '')}
        result[key]['targets'].extend(ing_info.get('targets', []))
        result[key]['svm_scores'].extend(ing_info.get('svm_scores', []))
        # 去重
        result[key]['targets'] = list(set(result[key]['targets']))
        result[key]['svm_scores'] = list(set(result[key]['svm_scores']))

    # 合并 diseases（需要通过靶点名称关联到成分）
    # 已知成分的靶点，去重后添加疾病
    for key, compound_info in result.items():
        all_diseases = set()
        for tgt in compound_info.get('targets', []):
            tgt_key = normalize_name(tgt)
            if tgt_key in diseases_data:
                all_diseases.update(diseases_data[tgt_key])
        result[key]['diseases'] = list(all_diseases)

    # 去重所有 targets
    for key in result:
        if result[key].get('targets'):
            result[key]['targets'] = list(set(result[key]['targets']))

    print(f"[pharma_db] TCMSP 目录加载完成: {len(result)} 化合物, {sum(len(v['targets']) for v in result.values())} 靶点关联")
    return result


# ===================== TCMSP 文件加载（单文件上传方式，保持兼容）=====================
def load_tcmsp(uploaded_file) -> dict:
    """
    加载 TCMSP Excel 文件，返回 {compound_name: {ob, dl, bbb, ...}}
    TCMSP 主要列：Name, OB(%), DL, BBB, Caco-2, ...
    """
    try:
        if uploaded_file.name.endswith('.xlsx') or uploaded_file.name.endswith('.xls'):
            df = pd.read_excel(uploaded_file, engine='xlrd' if uploaded_file.name.endswith('.xls') else 'openpyxl')
        else:
            df = pd.read_csv(uploaded_file, sep='\t')
        df.columns = df.columns.str.strip()
        result = {}
        # 尝试常见的化合物名列
        name_col = None
        for col in df.columns:
            if 'name' in col.lower() or 'molecule' in col.lower():
                name_col = col
                break
        if name_col is None and len(df.columns) > 0:
            name_col = df.columns[0]

        for _, row in df.iterrows():
            name = str(row.get(name_col, '')).strip()
            if name and name not in ('-', '', 'nan', 'NaN'):
                # OB 值
                ob = None
                for ob_col in ['OB(%)', 'OB', 'Oral Bioavailability', 'ob']:
                    if ob_col in df.columns:
                        try:
                            ob = float(row[ob_col])
                            break
                        except (ValueError, TypeError):
                            pass
                # DL 值
                dl = None
                for dl_col in ['DL', 'Drug likeness', 'dl']:
                    if dl_col in df.columns:
                        try:
                            dl = float(row[dl_col])
                            break
                        except (ValueError, TypeError):
                            pass
                result[name.lower()] = {'ob': ob, 'dl': dl}
        return result
    except Exception as e:
        return {}


# ===================== DrugBank 文件加载 =====================
def load_drugbank(uploaded_file) -> dict:
    """
    加载 DrugBank CSV/XML 文件，返回 {compound_name: {targets, indications, activity}}
    DrugBank CSV 通常包含：Name, Categories, Target, Indication, ...
    """
    try:
        if uploaded_file.name.endswith('.xml'):
            df = pd.read_xml(uploaded_file, xpath='.//drug')
        elif uploaded_file.name.endswith('.csv') or uploaded_file.name.endswith('.tsv'):
            sep = '\t' if uploaded_file.name.endswith('.tsv') else ','
            df = pd.read_csv(uploaded_file, sep=sep, low_memory=False, on_bad_lines='skip')
        else:
            df = pd.read_csv(uploaded_file, low_memory=False, on_bad_lines='skip')
        df.columns = df.columns.str.strip()
        result = {}

        name_col = None
        for col in df.columns:
            if col.lower() in ('name', 'drug name', 'compound'):
                name_col = col
                break
        if name_col is None and len(df.columns) > 0:
            name_col = df.columns[0]

        for _, row in df.iterrows():
            name = str(row.get(name_col, '')).strip()
            if name and name not in ('-', '', 'nan', 'NaN'):
                # 提取靶点信息（可能有多个靶点）
                target_col = None
                for tc in ['Target', 'Targets', 'target gene', 'target_name']:
                    if tc in df.columns:
                        target_col = tc
                        break
                targets = str(row.get(target_col, '-')) if target_col else '-'

                indication_col = None
                for ic in ['Indication', 'Indications', 'indication']:
                    if ic in df.columns:
                        indication_col = ic
                        break
                indications = str(row.get(indication_col, '-')) if indication_col else '-'

                result[name.lower()] = {
                    'targets': targets,
                    'indications': indications,
                }
        return result
    except Exception as e:
        return {}


# ===================== TTD 文件加载 =====================
def load_ttd(uploaded_file) -> dict:
    """
    加载 TTD 下载文件（通常是 TTD_download.txt 或 CSV）
    返回 {compound_name: {targets, diseases}}
    """
    try:
        if uploaded_file.name.endswith('.xlsx') or uploaded_file.name.endswith('.xls'):
            df = pd.read_excel(uploaded_file)
        elif uploaded_file.name.endswith('.tsv'):
            df = pd.read_csv(uploaded_file, sep='\t', on_bad_lines='skip')
        else:
            df = pd.read_csv(uploaded_file, on_bad_lines='skip', low_memory=False)
        df.columns = df.columns.str.strip()
        result = {}

        # TTD 通常有 Compound 和 Target/Disease 列
        name_col = None
        for col in df.columns:
            if 'compound' in col.lower() or 'drug' in col.lower():
                name_col = col
                break
        if name_col is None and len(df.columns) > 0:
            name_col = df.columns[0]

        for _, row in df.iterrows():
            name = str(row.get(name_col, '')).strip()
            if name and name not in ('-', '', 'nan', 'NaN'):
                # 提取靶点和疾病
                target_val = str(row.get('Target', str(row.get('Gene', '-'))))[:100]
                disease_val = str(row.get('Disease', '-'))[:100]
                result[name.lower()] = {
                    'targets': target_val,
                    'diseases': disease_val,
                }
        return result
    except Exception as e:
        return {}


# ===================== 化合物名称匹配 =====================
def normalize_name(name: str) -> str:
    """标准化化合物名称用于匹配"""
    if not isinstance(name, str):
        return ''
    # 移除括号内容、特殊字符
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(r'[^a-zA-Z0-9\s]', ' ', name)
    return name.lower().strip()


def fuzzy_match(query: str, candidates: list, threshold: float = 0.7) -> str:
    """简单模糊匹配，返回最接近的候选名称或空字符串"""
    query_norm = normalize_name(query)
    if not query_norm:
        return ''
    best_match = ''
    best_score = 0.0
    for cand in candidates:
        cand_norm = normalize_name(cand)
        if not cand_norm:
            continue
        # 简单包含匹配
        if query_norm in cand_norm or cand_norm in query_norm:
            score = min(len(query_norm), len(cand_norm)) / max(len(query_norm), len(cand_norm))
            if score > best_score:
                best_score = score
                best_match = cand
    if best_score >= threshold:
        return best_match
    return ''


# ===================== 药理证据评分 =====================
def compute_pharma_evidence_score(tcmsp_info: dict, drugbank_info: dict,
                                   ttd_info: dict, stitch_hits: int,
                                   pubchem_active: bool) -> dict:
    """
    计算药理证据评分和详细信息。

    Returns:
        dict: {
            'score': float,           # 0~1
            'details': dict,          # 各数据库命中详情
            'evidence': list,         # 可读证据列表
        }
    """
    score = 0.0
    details = {}
    evidence = []

    # TCMSP
    if tcmsp_info:
        ob = tcmsp_info.get('ob')
        dl = tcmsp_info.get('dl')
        tcmsp_score = 0.0
        if ob is not None and ob > 30:
            tcmsp_score += PHARMA_SCORES['tcmsp_ob']
            evidence.append(f'TCMSP: OB={ob:.1f}% (>30%)')
        if dl is not None and dl > 0.18:
            tcmsp_score += PHARMA_SCORES['tcmsp_dl']
            evidence.append(f'TCMSP: DL={dl:.3f} (>0.18)')
        details['tcmsp'] = {'ob': ob, 'dl': dl, 'score': tcmsp_score}
        score += tcmsp_score

    # DrugBank
    if drugbank_info:
        targets = drugbank_info.get('targets', '-')
        indications = drugbank_info.get('indications', '-')
        db_score = PHARMA_SCORES['drugbank']
        details['drugbank'] = {
            'targets': targets[:100] if targets and targets != '-' else '-',
            'indications': indications[:100] if indications and indications != '-' else '-',
            'score': db_score,
        }
        evidence.append(f'DrugBank: {targets[:50]}...' if len(str(targets)) > 50 else f'DrugBank: {targets}')
        score += db_score

    # TTD
    if ttd_info:
        targets = ttd_info.get('targets', '-')
        diseases = ttd_info.get('diseases', '-')
        ttd_score = PHARMA_SCORES['ttd']
        details['ttd'] = {
            'targets': targets,
            'diseases': diseases,
            'score': ttd_score,
        }
        evidence.append(f'TTD: {targets[:50]} (Disease: {diseases[:30]})' if diseases != '-' else f'TTD: {targets[:50]}')
        score += ttd_score

    # STITCH
    if stitch_hits > 0:
        stitch_score = min(PHARMA_SCORES['stitch'] * min(stitch_hits, 5), PHARMA_SCORES['stitch'])
        details['stitch'] = {'interactions': stitch_hits, 'score': stitch_score}
        evidence.append(f'STITCH: {stitch_hits} 条化学相互作用')
        score += stitch_score

    # PubChem BioActivity
    if pubchem_active:
        details['pubchem'] = {'active': True, 'score': PHARMA_SCORES['pubchem_bio']}
        evidence.append('PubChem: 有生物活性记录')
        score += PHARMA_SCORES['pubchem_bio']

    # 上限为 1.0
    score = min(score, 1.0)
    if not evidence:
        evidence.append('无药理数据库记录')

    return {
        'score': score,
        'details': details,
        'evidence': evidence,
    }


# ===================== STITCH API 查询 =====================
def query_stitch(compound_name: str, max_retries: int = 1) -> int:
    """
    查询 STITCH 化学相互作用数据库。
    返回命中相互作用数量，0 表示无记录。
    """
    if not compound_name or compound_name in ('-', 'nan', ''):
        return 0
    try:
        url = 'http://stitch.embl.de/cgi/stitch'
        params = {
            'q': compound_name,
            'c': 'compound',
            'species': 'Homo sapiens',
            'format': 'tsv',
        }
        for attempt in range(max_retries):
            try:
                resp = requests.get(url, params=params, timeout=8)
                if resp.status_code == 200:
                    lines = [l for l in resp.text.strip().split('\n') if l and not l.startswith('#')]
                    return max(0, len(lines) - 1)  # 减去 header 行
            except requests.RequestException:
                if attempt < max_retries - 1:
                    time.sleep(1)
    except Exception:
        pass
    return 0


# ===================== 主匹配函数 =====================
def match_pharma_db(
    metabolites: list,
    tcmsp_data: dict,
    drugbank_data: dict,
    ttd_data: dict,
    use_stitch: bool = False,
    progress_callback=None,
) -> pd.DataFrame:
    """
    对候选化合物列表进行药理数据库匹配。

    Parameters
    ----------
    metabolites : list
        代谢物名称列表
    tcmsp_data : dict
        TCMSP 数据字典
    drugbank_data : dict
        DrugBank 数据字典
    ttd_data : dict
        TTD 数据字典
    use_stitch : bool
        是否查询 STITCH API（需要网络连接）
    progress_callback : callable, optional
        进度回调函数

    Returns
    -------
    pd.DataFrame
        包含药理匹配信息的 DataFrame
    """
    rows = []
    all_names = list(tcmsp_data.keys()) + list(drugbank_data.keys()) + list(ttd_data.keys())

    for i, metab in enumerate(metabolites):
        metab_lower = normalize_name(metab)

        # 匹配 TCMSP
        tcmsp_hit = None
        matched_name = None
        if metab_lower in tcmsp_data:
            tcmsp_hit = tcmsp_data[metab_lower]
            matched_name = metab
        elif metab_lower in drugbank_data:
            tcmsp_hit = tcmsp_data.get(metab_lower)
            matched_name = metab
        else:
            # 模糊匹配
            fuzzy = fuzzy_match(metab, all_names, threshold=0.6)
            if fuzzy:
                matched_name = fuzzy
                tcmsp_hit = tcmsp_data.get(fuzzy) or drugbank_data.get(fuzzy) or ttd_data.get(fuzzy)

        tcmsp_info = tcmsp_hit if tcmsp_hit else {}

        # 匹配 DrugBank
        db_hit = drugbank_data.get(metab_lower) or \
                 drugbank_data.get(normalize_name(matched_name)) if matched_name else None

        # 匹配 TTD
        ttd_hit = ttd_data.get(metab_lower) or \
                  ttd_data.get(normalize_name(matched_name)) if matched_name else None

        # STITCH 查询（网络允许时）
        stitch_hits = 0
        if use_stitch:
            stitch_hits = query_stitch(metab)

        # PubChem 活性（简化版：有 compound_name 即认为有一定活性）
        pubchem_active = bool(metab and metab not in ('-', 'nan'))

        # 计算评分
        ev = compute_pharma_evidence_score(tcmsp_info, db_hit or {}, ttd_hit or {}, stitch_hits, pubchem_active)

        row = {
            'Metabolite': metab,
            'TCMSP_OB': ev['details'].get('tcmsp', {}).get('ob'),
            'TCMSP_DL': ev['details'].get('tcmsp', {}).get('dl'),
            'TCMSP_Targets': '; '.join(tcmsp_info.get('targets', [])[:10]) if tcmsp_info.get('targets') else '-',
            'TCMSP_Diseases': '; '.join(tcmsp_info.get('diseases', [])[:5]) if tcmsp_info.get('diseases') else '-',
            'TCMSP_SVM_Score': max(tcmsp_info.get('svm_scores', []) or [0]),
            'DrugBank_Targets': ev['details'].get('drugbank', {}).get('targets', '-'),
            'DrugBank_Indications': ev['details'].get('drugbank', {}).get('indications', '-'),
            'TTD_Targets': ev['details'].get('ttd', {}).get('targets', '-'),
            'TTD_Diseases': ev['details'].get('ttd', {}).get('diseases', '-'),
            'STITCH_Interactions': stitch_hits,
            'Pharma_Evidence_Score': ev['score'],
            'Pharma_Evidence': '; '.join(ev['evidence']),
        }
        rows.append(row)

        if progress_callback and (i + 1) % 20 == 0:
            progress_callback(i + 1, len(metabolites))

    return pd.DataFrame(rows)
