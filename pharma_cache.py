# -*- coding: utf-8 -*-
"""
药理数据库自动下载与缓存模块 - 修复版
去除悬空代码，添加并发下载优化
已修复：DrugCentral 磁盘缓存持久化，避免每次重新下载
"""
import io, os, re, time, gzip, json, hashlib, requests, pandas as pd, numpy as np
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ===================== 磁盘缓存路径 =====================
_CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.pharma_cache')
_DRUGCENTRAL_CACHE_FILE = os.path.join(_CACHE_DIR, 'drugcentral_data.json')

# ===================== 全局缓存 =====================
_drugcentral_loaded = False
_drugcentral_targets = {}
_drugcentral_structures = {}
_pubchem_cache = {}

# ===================== 下载地址 =====================
_DRUGCENTRAL_INTERACTION_URL = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
_DRUGCENTRAL_STRUCTURES_URL = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv"

# ===================== 连接池 =====================
_session = None
def _get_session():
    global _session
    if _session is None:
        _session = requests.Session()
        adapter = HTTPAdapter(pool_connections=10, pool_maxsize=10, max_retries=Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504]))
        _session.mount('http://', adapter); _session.mount('https://', adapter)
    return _session

def _normalize_name(name):
    if not isinstance(name, str): return ''
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(r'[^a-zA-Z0-9\s]', ' ', name)
    return name.lower().strip()

def _ensure_cache_dir():
    if not os.path.exists(_CACHE_DIR): os.makedirs(_CACHE_DIR, exist_ok=True)

def _validate_disk_cache():
    """
    验证磁盘缓存是否完整有效（不依赖内存中的 _drugcentral_loaded 标志）。
    解决 Streamlit 重载时全局变量重置导致缓存失效的问题。
    """
    cache_path = _DRUGCENTRAL_CACHE_FILE
    if not os.path.exists(cache_path):
        return False
    # 文件太小说明下载不完整（实际 drugcentral_data.json 应 > 5MB）
    if os.path.getsize(cache_path) < 1024 * 1024:  # < 1MB
        return False
    # 尝试解析验证结构完整性
    try:
        with open(cache_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if not isinstance(data, dict):
            return False
        targets = data.get('targets', {})
        structures = data.get('structures', {})
        # 至少要有实质数据（DrugCentral 解析后应有数万条记录）
        if not isinstance(targets, dict) or len(targets) < 100:
            return False
        # structures 允许为空（structures 文件可能下载失败）
        if not isinstance(structures, dict):
            return False
        return True
    except Exception:
        return False

def _load_drugcentral_from_disk():
    if not os.path.exists(_DRUGCENTRAL_CACHE_FILE): return False
    try:
        with open(_DRUGCENTRAL_CACHE_FILE, 'r', encoding='utf-8') as f: data = json.load(f)
        global _drugcentral_targets, _drugcentral_structures
        _drugcentral_targets = data.get('targets', {})
        _drugcentral_structures = data.get('structures', {})
        print(f"[pharma_cache] DrugCentral loaded from disk: {len(_drugcentral_targets)} targets")
        return True
    except Exception as e:
        print(f"[pharma_cache] Disk cache read failed: {e}"); return False

def _save_drugcentral_to_disk():
    _ensure_cache_dir()
    try:
        with open(_DRUGCENTRAL_CACHE_FILE, 'w', encoding='utf-8') as f:
            json.dump({'targets': _drugcentral_targets, 'structures': _drugcentral_structures}, f, ensure_ascii=False)
        print(f"[pharma_cache] DrugCentral saved to disk")
    except Exception as e:
        print(f"[pharma_cache] Disk cache write failed: {e}")

def _download_with_retry(session, url, max_retries=3, timeout=120):
    for attempt in range(max_retries):
        try:
            r = session.get(url, timeout=timeout, stream=True)
            if r.status_code == 200:
                content = r.content
                print(f"[pharma_cache] Downloaded: {url.split('/')[-1]} ({len(content)/(1024*1024):.1f} MB)")
                return content
            elif r.status_code == 404:
                print(f"[pharma_cache] File not found: {url}"); return b''
            else:
                print(f"[pharma_cache] Download failed ({r.status_code}), retry {attempt+1}/{max_retries}")
        except Exception as e:
            print(f"[pharma_cache] Download attempt {attempt+1} failed: {e}")
            if attempt < max_retries-1: time.sleep(2 ** attempt)
    return b''

def _load_drugcentral(force_download=False):
    global _drugcentral_loaded, _drugcentral_targets, _drugcentral_structures
    # 如果已加载且非强制下载，直接返回（内存命中）
    if _drugcentral_loaded and not force_download and len(_drugcentral_targets) > 0:
        return True
    # 优先检查磁盘缓存（不依赖 _drugcentral_loaded 标志，解决 Streamlit 重载问题）
    if not force_download and _validate_disk_cache():
        if _load_drugcentral_from_disk():
            _drugcentral_loaded = True
            return True
    # 以下为下载逻辑...
    session = _get_session()
    with ThreadPoolExecutor(max_workers=2) as executor:
        f1 = executor.submit(_download_with_retry, session, _DRUGCENTRAL_INTERACTION_URL)
        f2 = executor.submit(_download_with_retry, session, _DRUGCENTRAL_STRUCTURES_URL)
        try: interaction_data = f1.result(timeout=300)
        except Exception as e: print(f"[pharma_cache] interaction download failed: {e}"); interaction_data = b''
        try: structure_data = f2.result(timeout=300)
        except Exception as e: print(f"[pharma_cache] structures download failed: {e}"); structure_data = b''
    if interaction_data: _parse_interaction_data(interaction_data)
    else: print("[pharma_cache] Interaction download failed, using PubChem only")
    if structure_data: _parse_structure_data(structure_data)
    _save_drugcentral_to_disk()
    _drugcentral_loaded = True
    print(f"[pharma_cache] DrugCentral loaded: {len(_drugcentral_targets)} targets, {len(_drugcentral_structures)} structures")
    return True

def _parse_interaction_data(data):
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
        print(f"[pharma_cache] Parsed {len(_drugcentral_targets)} interaction records")
    except Exception as e: print(f"[pharma_cache] Interaction parse failed: {e}")

def _parse_structure_data(data):
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
        print(f"[pharma_cache] Parsed {len(_drugcentral_structures)} structure records")
    except Exception as e: print(f"[pharma_cache] Structure parse failed: {e}")

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_RATE_LIMIT_DELAY = 0.25
_pubchem_session = None

def _get_pubchem_session():
    global _pubchem_session
    if _pubchem_session is None:
        _pubchem_session = requests.Session()
        adapter = HTTPAdapter(pool_connections=10, pool_maxsize=10, max_retries=Retry(total=2, backoff_factor=0.3))
        _pubchem_session.mount('http://', adapter); _pubchem_session.mount('https://', adapter)
    return _pubchem_session

def _pubchem_rate_limit(): time.sleep(_RATE_LIMIT_DELAY)

def query_pubchem_compound(compound_name, max_retries=2):
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

def query_pubchem_bioactivity_summary(compound_name, max_results=5):
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

def query_pubchem_targets(compound_name, max_targets=10):
    nk = _normalize_name(compound_name)
    if nk in _pubchem_cache and not _pubchem_cache[nk].get('bioactivity_count',0)>0: return []
    session = _get_pubchem_session()
    targets=[]
    try:
        _pubchem_rate_limit()
        url = f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
        r = session.get(url, timeout=10)
        if r.status_code!=200: return []
        cids = r.json().get('IdentifierList',{}).get('CID',[]); cid=cids[0] if cids else None
        if not cid: return []
        _pubchem_rate_limit()
        target_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/targets/JSON"
        r2 = session.get(target_url, timeout=15)
        if r2.status_code==200:
            for t in r2.json().get('TargetList',[])[:max_targets]:
                tid=t.get('TargetID',''); tname=t.get('TargetName','')
                if tname and tid: targets.append(f"{tname} (ID:{tid})")
    except: pass
    return targets

def query_pharma_info(compound_name):
    nk = _normalize_name(compound_name)
    _load_drugcentral()
    info = {'source':[],'drugcentral_targets':[],'drugcentral_diseases':[],'drugcentral_actions':[],
            'pubchem_cid':None,'pubchem_smiles':None,'pubchem_xlogp':None,
            'pubchem_bioactivity_count':0,'pubchem_bioactive':False,'pubchem_targets':[],'pubchem_activities':[]}
    dc = _drugcentral_targets.get(nk)
    if dc:
        info['source'].append('DrugCentral')
        info['drugcentral_targets']=dc.get('targets',[])[:10]
        info['drugcentral_diseases']=dc.get('diseases',[])[:5]
        info['drugcentral_actions']=dc.get('actions',[])[:5]
    pc = query_pubchem_compound(compound_name)
    if pc.get('cid'):
        info['source'].append('PubChem')
        info.update({'pubchem_cid':pc['cid'],'pubchem_smiles':pc['smiles'],'pubchem_xlogp':pc['xlogp'],
                      'pubchem_bioactivity_count':pc.get('bioactivity_count',0),'pubchem_bioactive':pc.get('active',False)})
    if info['pubchem_bioactive']:
        info['pubchem_targets']=query_pubchem_targets(compound_name)
        info['pubchem_activities']=query_pubchem_bioactivity_summary(compound_name)
    return info

def compute_pharma_evidence_score(info):
    score=0.0; evidence=[]
    if info['drugcentral_targets']:
        score+=0.9; n=len(info['drugcentral_targets'])
        tgt='; '.join(info['drugcentral_targets'][:3])
        evidence.append(f"DrugCentral: {n} targets ({tgt})")
        if info['drugcentral_diseases']: evidence.append(f"  Disease: {info['drugcentral_diseases'][0]}")
    if info['pubchem_bioactive']:
        bio_score=min(info['pubchem_bioactivity_count']/50.0,1.0)*0.5; score+=bio_score
        evidence.append(f"PubChem: {info['pubchem_bioactivity_count']} bioactivity records")
        if info['pubchem_targets']: evidence.append(f"  Targets: {'; '.join(info['pubchem_targets'][:3])}")
        if info['pubchem_activities']: evidence.append(f"  Activity: {info['pubchem_activities'][0][:60]}")
    if not info['source']: evidence.append('No database records')
    return min(score,1.0), evidence

def match_pharma_online(metabolites, progress_callback=None):
    rows=[]; total=len(metabolites)
    for i, metab in enumerate(metabolites):
        info=query_pharma_info(metab); score, evidence=compute_pharma_evidence_score(info)
        rows.append({'Metabolite':metab,'Data_Sources':', '.join(info['source']),
            'DrugCentral_Targets':'; '.join(info['drugcentral_targets'][:5]),
            'DrugCentral_Diseases':'; '.join(info['drugcentral_diseases'][:3]),
            'DrugCentral_Actions':'; '.join(info['drugcentral_actions'][:3]),
            'PubChem_CID':info['pubchem_cid'],'PubChem_XLogP':info['pubchem_xlogp'],
            'PubChem_BioActivity_Count':info['pubchem_bioactivity_count'],
            'PubChem_Targets':'; '.join(info['pubchem_targets'][:5]),
            'PubChem_Activities':'; '.join(info['pubchem_activities'][:3]),
            'Pharma_Evidence_Score':score,'Pharma_Evidence':' | '.join(evidence)})
        if progress_callback and (i+1)%10==0: progress_callback(i+1,total)
    return pd.DataFrame(rows)

def is_drugcentral_available(): return _drugcentral_loaded and len(_drugcentral_targets)>0

def get_cache_status():
    """
    返回缓存状态，包含磁盘缓存信息。
    UI 层可据此判断是否显示"使用本地缓存"而非"下载"。
    """
    disk_cache_valid = _validate_disk_cache()
    disk_count = 0
    if disk_cache_valid:
        try:
            with open(_DRUGCENTRAL_CACHE_FILE, 'r', encoding='utf-8') as f:
                data = json.load(f)
            disk_count = len(data.get('targets', {}))
        except Exception:
            pass
    return {
        'drugcentral_loaded': _drugcentral_loaded and len(_drugcentral_targets) > 0,
        'drugcentral_targets_count': len(_drugcentral_targets),
        'drugcentral_structures_count': len(_drugcentral_structures),
        'disk_cache_valid': disk_cache_valid,
        'disk_cache_path': _DRUGCENTRAL_CACHE_FILE if disk_cache_valid else None,
        'disk_cache_targets_count': disk_count,
        'pubchem_cache_count': len(_pubchem_cache)
    }

