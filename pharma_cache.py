# -*- coding: utf-8 -*-
"""
药理数据库自动下载与缓存模块
自动从 DrugCentral 下载药物-靶点互作数据，并集成 PubChem PUG REST API
在首次查询时自动下载，仅下载一次，之后从内存缓存读取
"""

import io
import os
import re
import time
import gzip
import hashlib
import requests
import pandas as pd
import numpy as np

# ===================== 全局缓存（单例） =====================
_drugcentral_loaded = False
_drugcentral_targets = {}   # {compound_name_lower: {targets, indications, diseases}}
_drugcentral_structures = {}  # {compound_name_lower: {smiles, name}}

_pubchem_cache = {}          # {compound_name_lower: {properties, bioactivity}}

# ===================== DrugCentral 自动下载 =====================
_DRUGCENTRAL_INTERACTION_URL = (
    "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
)
_DRUGCENTRAL_STRUCTURES_URL = (
    "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv"
)


def _normalize_name(name: str) -> str:
    """标准化化合物名称用于匹配"""
    if not isinstance(name, str):
        return ''
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(r'[^a-zA-Z0-9\s]', ' ', name)
    return name.lower().strip()


def _load_drugcentral(force_download=False) -> bool:
    """
    加载 DrugCentral 数据到内存。
    首次调用时自动下载（~800KB），之后直接返回缓存。
    返回 True 表示加载成功，False 表示失败。
    """
    global _drugcentral_loaded, _drugcentral_targets, _drugcentral_structures

    if _drugcentral_loaded and not force_download:
        return True

    print("[pharma_cache] 首次使用，正在下载 DrugCentral 数据库（仅下载一次）...")

    # 下载 1: drug-target interaction
    interaction_data = _download_with_retry(_DRUGCENTRAL_INTERACTION_URL)
    if interaction_data:
        _parse_interaction_data(interaction_data)
    else:
        print("[pharma_cache] DrugCentral interaction 下载失败，将仅使用 PubChem 数据")

    # 下载 2: structures
    structure_data = _download_with_retry(_DRUGCENTRAL_STRUCTURES_URL)
    if structure_data:
        _parse_structure_data(structure_data)

    _drugcentral_loaded = True
    print(f"[pharma_cache] DrugCentral 加载完成: {len(_drugcentral_targets)} 靶点记录, {len(_drugcentral_structures)} 结构记录")
    return True


def _download_with_retry(url: str, max_retries: int = 2) -> bytes:
    """带重试的下载"""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=60, stream=True)
            if r.status_code == 200:
                return r.content
        except Exception as e:
            print(f"[pharma_cache] 下载失败（尝试 {attempt+1}）: {e}")
            if attempt < max_retries - 1:
                time.sleep(2)
    return b''


def _parse_interaction_data(data: bytes) -> None:
    """解析 DrugCentral drug-target interaction TSV"""
    global _drugcentral_targets
    try:
        # 如果是 gz 压缩
        if data[:2] == b'\x1f\x8b':
            f = gzip.GzipFile(fileobj=io.BytesIO(data))
            content = f.read()
        else:
            content = data

        df = pd.read_csv(io.BytesIO(content), sep='\t', low_memory=False,
                         on_bad_lines='skip')
        df.columns = df.columns.str.strip()

        # 找关键列
        name_col = None
        target_col = None
        disease_col = None
        act_col = None

        for col in df.columns:
            cl = col.lower()
            if cl in ('drug_name', 'compound_name', 'name', 'ingredient', 'drug'):
                name_col = col
            elif cl in ('target', 'target_name', 'gene', 'uniprot'):
                target_col = col
            elif cl in ('disease', 'indication', 'do_name'):
                disease_col = col
            elif cl in ('action', 'activity', 'relationship'):
                act_col = col

        if name_col is None:
            name_col = df.columns[0]
        if target_col is None:
            for col in df.columns:
                if col != name_col:
                    target_col = col
                    break

        for _, row in df.iterrows():
            name = str(row.get(name_col, '')).strip()
            if not name or name in ('-', '', 'nan'):
                continue
            name_key = _normalize_name(name)
            target = str(row.get(target_col, '-'))[:200] if target_col else '-'
            disease = str(row.get(disease_col, '-'))[:200] if disease_col else '-'
            action = str(row.get(act_col, '-'))[:100] if act_col else '-'

            if name_key not in _drugcentral_targets:
                _drugcentral_targets[name_key] = {
                    'targets': [],
                    'diseases': [],
                    'actions': [],
                }
            _drugcentral_targets[name_key]['targets'].append(target)
            _drugcentral_targets[name_key]['diseases'].append(disease)
            _drugcentral_targets[name_key]['actions'].append(action)

        # 去重
        for name_key in _drugcentral_targets:
            info = _drugcentral_targets[name_key]
            info['targets'] = list(set(info['targets']))
            info['diseases'] = list(set(info['diseases']))
            info['actions'] = list(set(info['actions']))

    except Exception as e:
        print(f"[pharma_cache] 解析 DrugCentral interaction 失败: {e}")


def _parse_structure_data(data: bytes) -> None:
    """解析 DrugCentral structures SMILES TSV"""
    global _drugcentral_structures
    try:
        df = pd.read_csv(io.BytesIO(data), sep='\t', low_memory=False,
                         on_bad_lines='skip')
        df.columns = df.columns.str.strip()

        name_col = None
        for col in df.columns:
            if col.lower() in ('name', 'drug_name', 'compound'):
                name_col = col
                break
        if name_col is None:
            name_col = df.columns[0]

        smiles_col = None
        for col in df.columns:
            if col.lower() in ('smiles', 'structure'):
                smiles_col = col
                break

        for _, row in df.iterrows():
            name = str(row.get(name_col, '')).strip()
            if not name or name in ('-', '', 'nan'):
                continue
            name_key = _normalize_name(name)
            smiles = str(row.get(smiles_col, '-'))[:500] if smiles_col else '-'
            _drugcentral_structures[name_key] = {'smiles': smiles, 'name': name}

    except Exception as e:
        print(f"[pharma_cache] 解析 DrugCentral structures 失败: {e}")


# ===================== PubChem PUG REST API =====================
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_RATE_LIMIT_DELAY = 0.25  # 每次请求间隔（秒）


def _pubchem_rate_limit():
    """简单的速率限制"""
    time.sleep(_RATE_LIMIT_DELAY)


def query_pubchem_compound(compound_name: str, max_retries: int = 2) -> dict:
    """
    通过 PubChem PUG REST API 查询化合物属性。
    返回包含 CIDs、SMILES、XLogP、分子量等信息。
    """
    global _pubchem_cache
    name_key = _normalize_name(compound_name)

    # 检查内存缓存
    if name_key in _pubchem_cache:
        return _pubchem_cache[name_key]

    result = {
        'cid': None, 'smiles': None, 'xlogp': None,
        'molecular_weight': None, 'complexity': None,
        'bioactivity_count': 0, 'active': False,
    }

    for attempt in range(max_retries):
        try:
            _pubchem_rate_limit()

            # Step 1: 通过名称获取 CID
            url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
            r = requests.get(url, timeout=10)
            if r.status_code != 200:
                break
            data = r.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            if not cids:
                break
            cid = cids[0]
            result['cid'] = cid

            # Step 2: 获取化合物属性
            _pubchem_rate_limit()
            prop_url = (
                f"{_PUBCHEM_BASE}/compound/cid/{cid}/property"
                f"/CanonicalSMILES,XLogP,MolecularWeight,Complexity/JSON"
            )
            r2 = requests.get(prop_url, timeout=10)
            if r2.status_code == 200:
                props = r2.json().get('PropertyTable', {}).get('Properties', [{}])[0]
                result['smiles'] = props.get('CanonicalSMILES')
                result['xlogp'] = props.get('XLogP')
                result['molecular_weight'] = props.get('MolecularWeight')
                result['complexity'] = props.get('Complexity')

            # Step 3: 获取 BioActivity 简要信息
            _pubchem_rate_limit()
            bio_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON"
            r3 = requests.get(bio_url, timeout=10)
            if r3.status_code == 200:
                bio_data = r3.json()
                count = bio_data.get('InformationList', {}).get('Information', [{}])
                if count and isinstance(count, list):
                    result['bioactivity_count'] = int(count[0].get('TotalResults', 0))
                    result['active'] = result['bioactivity_count'] > 0

            break
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                pass

    _pubchem_cache[name_key] = result
    return result


def query_pubchem_bioactivity_summary(compound_name: str, max_results: int = 5) -> list:
    """
    获取化合物的 BioActivity 记录摘要。
    返回活性类型描述列表。
    """
    name_key = _normalize_name(compound_name)
    if name_key in _pubchem_cache:
        cached = _pubchem_cache[name_key]
        if not cached.get('active'):
            return []

    try:
        _pubchem_rate_limit()
        cid = None
        url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            cids = r.json().get('IdentifierList', {}).get('CID', [])
            cid = cids[0] if cids else None

        if not cid:
            return []

        _pubchem_rate_limit()
        bio_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON"
        r2 = requests.get(bio_url, timeout=15)
        if r2.status_code != 200:
            return []

        data = r2.json()
        assays = data.get('AssaySummaries', {}).get('AssaySummary', [])
        descriptions = []
        seen = set()
        for assay in assays:
            desc = str(assay.get('Description', ''))[:100]
            if desc and desc not in seen:
                seen.add(desc)
                descriptions.append(desc)
                if len(descriptions) >= max_results:
                    break
        return descriptions

    except Exception:
        return []


def query_pubchem_targets(compound_name: str, max_targets: int = 10) -> list:
    """
    通过 PubChem PUG REST 获取化合物的靶点信息（如果有）。
    """
    name_key = _normalize_name(compound_name)
    if name_key in _pubchem_cache:
        cached = _pubchem_cache[name_key]
        if not cached.get('bioactivity_count', 0) > 0:
            return []

    targets = []
    try:
        _pubchem_rate_limit()
        url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return []
        cids = r.json().get('IdentifierList', {}).get('CID', [])
        if not cids:
            return []
        cid = cids[0]

        # 通过 PUG REST 获取靶点
        _pubchem_rate_limit()
        target_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/targets/JSON"
        r2 = requests.get(target_url, timeout=15)
        if r2.status_code == 200:
            data = r2.json()
            target_list = data.get('TargetList', [])
            for t in target_list[:max_targets]:
                tid = t.get('TargetID', '')
                tname = t.get('TargetName', '')
                if tname and tid:
                    targets.append(f"{tname} (ID:{tid})")

    except Exception:
        pass
    return targets


# ===================== 统一查询接口 =====================

def query_pharma_info(compound_name: str) -> dict:
    """
    统一查询接口：同时查询 DrugCentral + PubChem。
    返回完整的药理信息字典。
    """
    name_key = _normalize_name(compound_name)

    # 确保 DrugCentral 已加载
    _load_drugcentral()

    info = {
        'source': [],
        'drugcentral_targets': [],
        'drugcentral_diseases': [],
        'drugcentral_actions': [],
        'pubchem_cid': None,
        'pubchem_smiles': None,
        'pubchem_xlogp': None,
        'pubchem_bioactivity_count': 0,
        'pubchem_bioactive': False,
        'pubchem_targets': [],
        'pubchem_activities': [],
    }

    # DrugCentral
    dc = _drugcentral_targets.get(name_key)
    if dc:
        info['source'].append('DrugCentral')
        info['drugcentral_targets'] = dc.get('targets', [])[:10]
        info['drugcentral_diseases'] = dc.get('diseases', [])[:5]
        info['drugcentral_actions'] = dc.get('actions', [])[:5]

    # PubChem
    pc = query_pubchem_compound(compound_name)
    if pc.get('cid'):
        info['source'].append('PubChem')
        info['pubchem_cid'] = pc['cid']
        info['pubchem_smiles'] = pc['smiles']
        info['pubchem_xlogp'] = pc['xlogp']
        info['pubchem_bioactivity_count'] = pc.get('bioactivity_count', 0)
        info['pubchem_bioactive'] = pc.get('active', False)

    # PubChem 靶点（单独请求）
    if info['pubchem_bioactive']:
        targets = query_pubchem_targets(compound_name)
        info['pubchem_targets'] = targets
        activities = query_pubchem_bioactivity_summary(compound_name)
        info['pubchem_activities'] = activities

    return info


def compute_pharma_evidence_score(info: dict) -> tuple:
    """
    根据药理信息计算证据评分。
    返回 (score: float, evidence: list[str])
    """
    score = 0.0
    evidence = []

    # DrugCentral 有靶点记录
    if info['drugcentral_targets']:
        score += 0.9
        n_targets = len(info['drugcentral_targets'])
        evidence.append(
            f"DrugCentral: {n_targets}靶点 "
            f"({'; '.join(info['drugcentral_targets'][:3])})"
        )
        if info['drugcentral_diseases']:
            evidence.append(f"  适应症: {info['drugcentral_diseases'][0]}")

    # PubChem BioActivity 有活性
    if info['pubchem_bioactive']:
        bio_score = min(info['pubchem_bioactivity_count'] / 50.0, 1.0) * 0.5
        score += bio_score
        evidence.append(
            f"PubChem: {info['pubchem_bioactivity_count']}条生物活性记录"
        )
        if info['pubchem_targets']:
            evidence.append(f"  靶点: {'; '.join(info['pubchem_targets'][:3])}")
        if info['pubchem_activities']:
            evidence.append(f"  活性: {info['pubchem_activities'][0][:60]}")

    # DrugCentral + PubChem 均无记录
    if not info['source']:
        evidence.append('无数据库记录（需实验验证）')

    # 上限
    score = min(score, 1.0)
    return score, evidence


def match_pharma_online(metabolites: list, progress_callback=None) -> pd.DataFrame:
    """
    对候选化合物列表进行在线药理数据库查询。

    Parameters
    ----------
    metabolites : list
        代谢物名称列表
    progress_callback : callable, optional
        进度回调函数，接受 (current, total) 两个参数

    Returns
    -------
    pd.DataFrame
        包含药理匹配信息的 DataFrame
    """
    rows = []
    total = len(metabolites)

    for i, metab in enumerate(metabolites):
        info = query_pharma_info(metab)
        score, evidence = compute_pharma_evidence_score(info)

        row = {
            'Metabolite': metab,
            'Data_Sources': ', '.join(info['source']),
            'DrugCentral_Targets': '; '.join(info['drugcentral_targets'][:5]),
            'DrugCentral_Diseases': '; '.join(info['drugcentral_diseases'][:3]),
            'DrugCentral_Actions': '; '.join(info['drugcentral_actions'][:3]),
            'PubChem_CID': info['pubchem_cid'],
            'PubChem_XLogP': info['pubchem_xlogp'],
            'PubChem_BioActivity_Count': info['pubchem_bioactivity_count'],
            'PubChem_Targets': '; '.join(info['pubchem_targets'][:5]),
            'PubChem_Activities': '; '.join(info['pubchem_activities'][:3]),
            'Pharma_Evidence_Score': score,
            'Pharma_Evidence': ' | '.join(evidence),
        }
        rows.append(row)

        if progress_callback and (i + 1) % 10 == 0:
            progress_callback(i + 1, total)

    return pd.DataFrame(rows)


def is_drugcentral_available() -> bool:
    """检查 DrugCentral 数据是否已加载"""
    return _drugcentral_loaded and len(_drugcentral_targets) > 0


def get_cache_status() -> dict:
    """获取缓存状态"""
    return {
        'drugcentral_loaded': _drugcentral_loaded,
        'drugcentral_targets_count': len(_drugcentral_targets),
        'drugcentral_structures_count': len(_drugcentral_structures),
        'pubchem_cache_count': len(_pubchem_cache),
    }
