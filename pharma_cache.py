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
import json
import hashlib
import requests
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ===================== DrugCentral 磁盘缓存路径 =====================
_CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.pharma_cache')
_DRUGCENTRAL_CACHE_FILE = os.path.join(_CACHE_DIR, 'drugcentral_data.json')

# ===================== 全局缓存（单例） =====================
_drugcentral_loaded = False
_drugcentral_targets = {}
_drugcentral_structures = {}
_pubchem_cache = {}

# ===================== DrugCentral 自动下载 =====================
_DRUGCENTRAL_INTERACTION_URL = (
    "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
)
_DRUGCENTRAL_STRUCTURES_URL = (
    "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv"
)

# ===================== 连接池（全局复用） =====================
_session = None

def _get_session() -> requests.Session:
    global _session
    if _session is None:
        _session = requests.Session()
        adapter = HTTPAdapter(
            pool_connections=10, pool_maxsize=10,
            max_retries=Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
        )
        _session.mount('http://', adapter)
        _session.mount('https://', adapter)
    return _session

def _normalize_name(name: str) -> str:
    if not isinstance(name, str): return ''
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(r'[^a-zA-Z0-9\s]', ' ', name)
    return name.lower().strip()

def _ensure_cache_dir():
    if not os.path.exists(_CACHE_DIR): os.makedirs(_CACHE_DIR, exist_ok=True)

def _load_drugcentral_from_disk() -> bool:
    if not os.path.exists(_DRUGCENTRAL_CACHE_FILE): return False
    try:
        with open(_DRUGCENTRAL_CACHE_FILE, 'r', encoding='utf-8') as f: data = json.load(f)
        global _drugcentral_targets, _drugcentral_structures
        _drugcentral_targets = data.get('targets', {})
        _drugcentral_structures = data.get('structures', {})
        print(f"[pharma_cache] DrugCentral 已从磁盘缓存加载: {len(_drugcentral_targets)} 靶点记录")
        return True
    except Exception as e:
        print(f"[pharma_cache] 磁盘缓存读取失败: {e}"); return False

def _save_drugcentral_to_disk():
    _ensure_cache_dir()
    try:
        with open(_DRUGCENTRAL_CACHE_FILE, 'w', encoding='utf-8') as f:
            json.dump({'targets': _drugcentral_targets, 'structures': _drugcentral_structures}, f, ensure_ascii=False)
        print(f"[pharma_cache] DrugCentral 已保存到磁盘缓存")
    except Exception as e:
        print(f"[pharma_cache] 磁盘缓存保存失败: {e}")

def _download_with_retry(session: requests.Session, url: str, max_retries: int = 3, timeout: int = 120) -> bytes:
    for attempt in range(max_retries):
        try:
            r = session.get(url, timeout=timeout, stream=True)
            if r.status_code == 200:
                content = r.content
                print(f"[pharma_cache] 下载完成: {url.split('/')[-1]} ({len(content)/(1024*1024):.1f} MB)")
                return content
            elif r.status_code == 404:
                print(f"[pharma_cache] 文件不存在: {url}"); return b''
            else:
                print(f"[pharma_cache] 下载失败（状态码 {r.status_code}），重试 {attempt+1}/{max_retries}")
        except Exception as e:
            print(f"[pharma_cache] 下载异常（尝试 {attempt+1}/{max_retries}）: {e}")
            if attempt < max_retries - 1: time.sleep(2 ** attempt)
    return b''

def _load_drugcentral(force_download=False) -> bool:
    global _drugcentral_loaded, _drugcentral_targets, _drugcentral_structures
    if _drugcentral_loaded and not force_download: return True
    if not force_download and _load_drugcentral_from_disk():
        _drugcentral_loaded = True; return True
    print("[pharma_cache] 首次使用，正在下载 DrugCentral 数据库（并发下载中，仅下载一次）...")
    session = _get_session()
    with ThreadPoolExecutor(max_workers=2) as executor:
        f1 = executor.submit(_download_with_retry, session, _DRUGCENTRAL_INTERACTION_URL)
        f2 = executor.submit(_download_with_retry, session, _DRUGCENTRAL_STRUCTURES_URL)
        try: interaction_data = f1.result(timeout=300)
        except Exception as e: print(f"[pharma_cache] interaction 文件下载超时或失败: {e}"); interaction_data = b''
        try: structure_data = f2.result(timeout=300)
        except Exception as e: print(f"[pharma_cache] structures 文件下载超时或失败: {e}"); structure_data = b''
    if interaction_data: _parse_interaction_data(interaction_data)
    else: print("[pharma_cache] DrugCentral interaction 下载失败，将仅使用 PubChem 数据")
    if structure_data: _parse_structure_data(structure_data)
    _save_drugcentral_to_disk()
    _drugcentral_loaded = True
    print(f"[pharma_cache] DrugCentral 加载完成: {len(_drugcentral_targets)} 靶点记录, {len(_drugcentral_structures)} 结构记录")
    return True

def _parse_interaction_data(data: bytes) -> None:
    global _drugcentral_targets
    try:
        content = gzip.GzipFile(fileobj=io.BytesIO(data)).read() if data[:2]==b'\x1f\x8b' else data
        df = pd.read_csv(io.BytesIO(content), sep='\t', low_memory=False, on_bad_lines='skip')
        df.columns = df.columns.str.strip()
        name_col = next((c for c in df.columns if c.lower() in ('drug_name','compound_name','name','ingredient','drug')), df.columns[0])
        target_col = next((c for c in df.columns if c.lower() in ('target','target_name','gene','uniprot') and c!=name_col), None)
        disease_col = next((c for c in df.columns if c.lower() in ('disease','indication','do_name')), None)
        act_col = next((c for c in df.columns if c.lower() in ('action','activity','relationship')), None)
        if target_col is None:
            for c in df.columns:
                if c != name_col: target_col = c; break
        for _, row in df.iterrows():
            name = str(row.get(name_col,'')).strip()
            if not name or name in ('-','','nan'): continue
            nk = _normalize_name(name)
            if nk not in _drugcentral_targets: _drugcentral_targets[nk] = {'targets':[],'diseases':[],'actions':[]}
            _drugcentral_targets[nk]['targets'].append(str(row.get(target_col,'-'))[:200] if target_col else '-')
            _drugcentral_targets[nk]['diseases'].append(str(row.get(disease_col,'-'))[:200] if disease_col else '-')
            _drugcentral_targets[nk]['actions'].append(str(row.get(act_col,'-'))[:100] if act_col else '-')
        for nk in _drugcentral_targets:
            d=_drugcentral_targets[nk]; d['targets']=list(set(d['targets'])); d['diseases']=list(set(d['diseases'])); d['actions']=list(set(d['actions']))
        print(f"[pharma_cache] 解析 interaction 数据完成: {len(_drugcentral_targets)} 条记录")
    except Exception as e: print(f"[pharma_cache] 解析 DrugCentral interaction 失败: {e}")

def _parse_structure_data(data: bytes) -> None:
    global _drugcentral_structures
    try:
        df = pd.read_csv(io.BytesIO(data), sep='\t', low_memory=False, on_bad_lines='skip')
        df.columns = df.columns.str.strip()
        name_col = next((c for c in df.columns if c.lower() in ('name','drug_name','compound')), df.columns[0])
        smiles_col = next((c for c in df.columns if c.lower() in ('smiles','structure')), None)
        for _, row in df.iterrows():
            name = str(row.get(name_col,'')).strip()
            if not name or name in ('-','','nan'): continue
            nk = _normalize_name(name)
            _drugcentral_structures[nk] = {'smiles': str(row.get(smiles_col,'-'))[:500] if smiles_col else '-', 'name': name}
        print(f"[pharma_cache] 解析 structures 数据完成: {len(_drugcentral_structures)} 条记录")
    except Exception as e: print(f"[pharma_cache] 解析 DrugCentral structures 失败: {e}")

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_RATE_LIMIT_DELAY = 0.25
_pubchem_session = None

def _get_pubchem_session() -> requests.Session:
    global _pubchem_session
    if _pubchem_session is None:
        _pubchem_session = requests.Session()
        adapter = HTTPAdapter(pool_connections=10, pool_maxsize=10, max_retries=Retry(total=2, backoff_factor=0.3))
        _pubchem_session.mount('http://', adapter); _pubchem_session.mount('https://', adapter)
    return _pubchem_session

def _pubchem_rate_limit(): time.sleep(_RATE_LIMIT_DELAY)

def query_pubchem_compound(compound_name: str, max_retries: int = 2) -> dict:
    global _pubchem_cache
    nk = _normalize_name(compound_name)
    if nk in _pubchem_cache: return _pubchem_cache[nk]
    result = {'cid':None,'smiles':None,'xlogp':None,'molecular_weight':None,'complexity':None,'bioactivity_count':0,'active':False}
    session = _get_pubchem_session()
    for attempt in range(max_retries):
        try:
            _pubchem_rate_limit()
            url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
            r = session.get(url, timeout=10)
            if r.status_code!=200: break
            cids = r.json().get('IdentifierList',{}).get('CID',[])
            if not cids: break
            cid = cids[0]; result['cid'] = cid
            _pubchem_rate_limit()
            prop_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/property/CanonicalSMILES,XLogP,MolecularWeight,Complexity/JSON"
            r2 = session.get(prop_url, timeout=10)
            if r2.status_code==200:
                props = r2.json().get('PropertyTable',{}).get('Properties',[{}])[0]
                result.update({'smiles':props.get('CanonicalSMILES'),'xlogp':props.get('XLogP'),'molecular_weight':props.get('MolecularWeight'),'complexity':props.get('Complexity')})
            _pubchem_rate_limit()
            bio_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON"
            r3 = session.get(bio_url, timeout=10)
            if r3.status_code==200:
                bio_data = r3.json()
                count = bio_data.get('InformationList',{}).get('Information',[{}])
                if count and isinstance(count,list):
                    result['bioactivity_count']=int(count[0].get('TotalResults',0)); result['active']=result['bioactivity_count']>0
            break
        except Exception as e:
            if attempt < max_retries-1: time.sleep(1)
    _pubchem_cache[nk] = result
    return result

def query_pubchem_bioactivity_summary(compound_name: str, max_results: int = 5) -> list:
    nk = _normalize_name(compound_name)
    if nk in _pubchem_cache and not _pubchem_cache[nk].get('active'): return []
    session = _get_pubchem_session()
    try:
        _pubchem_rate_limit()
        url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
        r = session.get(url, timeout=10)
        if r.status_code!=200: return []
        cids = r.json().get('IdentifierList',{}).get('CID',[]); cid=cids[0] if cids else None
        if not cid: return []
        _pubchem_rate_limit()
        bio_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON"
        r2 = session.get(bio_url, timeout=15)
        if r2.status_code!=200: return []
        assays = r2.json().get('AssaySummaries',{}).get('AssaySummary',[])
        descriptions=[]; seen=set()
        for assay in assays:
            desc = str(assay.get('Description',''))[:100]
            if desc and desc not in seen:
                seen.add(desc); descriptions.append(desc)
                if len(descriptions)>=max_results: break
        return descriptions
    except: return []

def query_pubchem_targe
...(truncated)...
