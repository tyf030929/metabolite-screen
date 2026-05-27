# -*- coding: utf-8 -*-
"""
网络药理学模块 - SMILES查询 / SwissTargetPrediction靶点预测 / gseapy富集分析
依赖: gseapy, pandas, numpy; PubChem API调用使用标准库 urllib（无需 requests）
"""

import ssl
import time
import re
import json
import socket
import subprocess
import shutil
from urllib.request import urlopen, Request
from urllib.parse import quote
from urllib.error import URLError, HTTPError
import pandas as pd
import numpy as np

try:
    import streamlit as st
    _ST_AVAILABLE = True
except ImportError:
    _ST_AVAILABLE = False
    st = None

# ===================== 常量 =====================
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_ST_BASE = "https://www.swisstargetprediction.ch/api"
_REQUEST_DELAY = 1.2  # SwissTargetPrediction 请求间隔（秒）
_MAX_RETRIES = 3     # SwissTargetPrediction 重试次数
_PROB_THRESHOLD = 0.01  # 靶点预测概率阈值

# ===================== curl 探测 =====================
_has_curl = shutil.which("curl") is not None


def _get_ssl_context():
    """返回一个禁用验证的 SSL context（仅用于 fallback，应对 Windows Python 3.13 证书问题）"""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def _http_get_json(url: str, timeout: int = 15) -> dict:
    """
    GET JSON，支持 SSL fallback + curl fallback。
    优先用 urllib → SSL 失败则禁用验证重试 → 最后 curl。
    """
    headers = {
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }

    # 方法1：urllib（正常 SSL）
    try:
        req = Request(url, headers=headers)
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if not raw:
                return {"_error": "Empty response"}
            try:
                return json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError as je:
                return {"_error": f"JSON decode error: {je}", "_raw": raw[:200].decode("utf-8", errors="replace")}
    except URLError as e:
        # 检查是否是 SSL 证书验证失败（被 URLError 包裹）
        if isinstance(e.reason, ssl.SSLCertVerificationError):
            # 方法1b：SSL 验证失败，禁用验证重试
            try:
                req = Request(url, headers=headers)
                with urlopen(req, timeout=timeout, context=_get_ssl_context()) as resp:
                    raw = resp.read()
                    if not raw:
                        return {"_error": "Empty response"}
                    try:
                        result = json.loads(raw.decode("utf-8"))
                        result["_via"] = "urllib_noverify"
                        return result
                    except json.JSONDecodeError as je:
                        return {"_error": f"JSON decode error (noverify): {je}", "_raw": raw[:200].decode("utf-8", errors="replace")}
            except Exception as e:
                err_str = f"urllib_noverify failed: {e}"
        else:
            err_str = str(e)
    except Exception as urllib_err:
        err_str = str(urllib_err)

    # 方法2：curl fallback
    if _has_curl:
        try:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", str(timeout),
                "-H", "Accept: application/json",
                "-H", "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
                url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=timeout + 2)
            if raw:
                try:
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "curl"
                    return result
                except json.JSONDecodeError as je:
                    return {"_error": f"curl JSON decode error: {je}", "_raw": raw[:200].decode("utf-8", errors="replace")}
        except Exception as e:
            err_str = f"curl failed: {e}"

    return {"_error": err_str}


def _http_post_json(url: str, payload: dict, timeout: int = 30) -> dict:
    """POST JSON，支持 SSL fallback + curl fallback"""
    body = json.dumps(payload).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }

    # 方法1：urllib（正常 SSL）
    try:
        req = Request(url, data=body, headers=headers)
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if not raw:
                return {"_error": "Empty response"}
            return json.loads(raw.decode("utf-8"))
    except URLError as e:
        # 检查是否是 SSL 证书验证失败（被 URLError 包裹）
        if isinstance(e.reason, ssl.SSLCertVerificationError):
            # 方法1b：SSL 验证失败，禁用验证重试
            try:
                req = Request(url, data=body, headers=headers)
                with urlopen(req, timeout=timeout, context=_get_ssl_context()) as resp:
                    raw = resp.read()
                    if not raw:
                        return {"_error": "Empty response"}
                    result = json.loads(raw.decode("utf-8"))
                    result["_via"] = "urllib_noverify"
                    return result
            except Exception as e:
                err_str = f"urllib_noverify failed: {e}"
        else:
            err_str = str(e)
    except Exception as urllib_err:
        err_str = str(urllib_err)

    # 方法2：curl fallback
    if _has_curl:
        try:
            curl_cmd = [
                "curl.exe", "-s", "-S",
                "--max-time", str(timeout),
                "-X", "POST",
                "-H", "Content-Type: application/json",
                "-H", "Accept: application/json",
                "-d", body.decode("utf-8"),
                url
            ]
            raw = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT, timeout=timeout + 2)
            if raw:
                result = json.loads(raw.decode("utf-8"))
                result["_via"] = "curl"
                return result
        except Exception:
            pass

    return {"_error": err_str}


def check_network_connectivity() -> dict:
    """
    返回网络连通性诊断信息。
    """
    results = {}

    # DNS 测试
    try:
        socket.setdefaulttimeout(5)
        socket.gethostbyname("pubchem.ncbi.nlm.nih.gov")
        results["dns"] = "OK"
    except Exception as e:
        results["dns"] = f"FAIL: {e}"

    # HTTP 测试（不用 urllib，用 curl if available）
    test_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Quercetin/property/IsomericSMILES/JSON"
    try:
        if _has_curl:
            out = subprocess.check_output(
                ["curl.exe", "-s", "-S", "--max-time", "10", test_url],
                stderr=subprocess.STDOUT
            )
            data = json.loads(out)
            props = data.get("PropertyTable", {}).get("Properties", [])
            sm = (props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES") or props[0].get("SMILES", "")) if props else ""
            results["http"] = f"OK → SMILES: {sm}" if sm else f"OK but no SMILES: {json.dumps(data)[:200]}"
        else:
            req = Request(test_url, headers={"User-Agent": "Mozilla/5.0"})
            with urlopen(req, timeout=10, context=_get_ssl_context()) as r:
                data = json.loads(r.read())
                props = data.get("PropertyTable", {}).get("Properties", [])
                sm = (props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES") or props[0].get("SMILES", "")) if props else ""
                results["http"] = f"OK → SMILES: {sm}" if sm else f"OK but no SMILES: {json.dumps(data)[:200]}"
    except Exception as e:
        results["http"] = f"FAIL: {e}"

    # SwissTargetPrediction 连通性
    try:
        st_url = "https://www.swisstargetprediction.ch/api/search"
        if _has_curl:
            out = subprocess.check_output(
                ["curl.exe", "-s", "-S", "--max-time", "10", "-X", "POST",
                 "-H", "Content-Type: application/json",
                 "-d", '{"smiles":"CC","species":"Homo sapiens"}',
                 st_url],
                stderr=subprocess.STDOUT
            )
            results["swiss"] = f"OK: {out.decode('utf-8')[:100]}"
        else:
            req = Request(st_url,
                          data=json.dumps({"smiles": "CC", "species": "Homo sapiens"}).encode(),
                          headers={"Content-Type": "application/json"})
            with urlopen(req, timeout=10, context=_get_ssl_context()) as r:
                results["swiss"] = f"OK: {r.read().decode('utf-8')[:100]}"
    except Exception as e:
        results["swiss"] = f"FAIL: {e}"

    return results

# ===================== 本地 SQLite SMILES 数据库 =====================
# 使用 Python 内置 sqlite3，完全不依赖第三方库
# 查询顺序：本地数据库 → PubChem API（找到后自动缓存到本地）

import sqlite3
import os

_DB_DIR = os.path.dirname(os.path.abspath(__file__))
_DB_PATH = os.path.join(_DB_DIR, "smiles_cache.db")


def _get_db() -> sqlite3.Connection:
    """获取数据库连接，自动建表"""
    conn = sqlite3.connect(_DB_PATH, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS compound_smiles (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            name_lower      TEXT    NOT NULL UNIQUE,
            name_raw        TEXT,
            smiles          TEXT    NOT NULL,
            source          TEXT    DEFAULT 'pubchem',
            cid             INTEGER,
            cached_at       TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS metabolite_smiles (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            metab_lower     TEXT    NOT NULL UNIQUE,
            metab_raw       TEXT,
            smiles          TEXT    NOT NULL,
            source          TEXT    DEFAULT 'pubchem',
            cached_at       TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS api_log (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            compound_name   TEXT,
            metabolite_name TEXT,
            smiles          TEXT,
            status          TEXT,
            via             TEXT,
            error_msg       TEXT,
            queried_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    conn.commit()
    return conn


def db_lookup_name(name: str) -> tuple:
    """
    在本地数据库按名称查找 SMILES。
    Returns: (smiles, source) 或 (None, None)
    """
    try:
        conn = _get_db()
        cur = conn.execute(
            "SELECT smiles, source FROM compound_smiles WHERE name_lower = ?",
            (str(name).strip().lower(),)
        )
        row = cur.fetchone()
        conn.close()
        if row:
            return row[0], row[1]
    except Exception:
        pass
    return None, None


def db_lookup_metabolite(metab: str) -> tuple:
    """
    在本地数据库按代谢物名称查找 SMILES。
    Returns: (smiles, source) 或 (None, None)
    """
    try:
        conn = _get_db()
        cur = conn.execute(
            "SELECT smiles, source FROM metabolite_smiles WHERE metab_lower = ?",
            (str(metab).strip().lower(),)
        )
        row = cur.fetchone()
        conn.close()
        if row:
            return row[0], row[1]
    except Exception:
        pass
    return None, None


def db_insert_name(name: str, smiles: str, source: str = "pubchem", cid: int = None) -> bool:
    """插入/更新名称→SMILES 记录"""
    try:
        conn = _get_db()
        conn.execute("""
            INSERT INTO compound_smiles (name_lower, name_raw, smiles, source, cid)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(name_lower) DO UPDATE SET
                smiles = excluded.smiles,
                source = excluded.source,
                cid = excluded.cid,
                cached_at = CURRENT_TIMESTAMP
        """, (str(name).strip().lower(), str(name).strip(), smiles, source, cid))
        conn.commit()
        conn.close()
        return True
    except Exception:
        return False


def db_insert_metabolite(metab: str, smiles: str, source: str = "pubchem") -> bool:
    """插入/更新代谢物名→SMILES 记录"""
    try:
        conn = _get_db()
        conn.execute("""
            INSERT INTO metabolite_smiles (metab_lower, metab_raw, smiles, source)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(metab_lower) DO UPDATE SET
                smiles = excluded.smiles,
                source = excluded.source,
                cached_at = CURRENT_TIMESTAMP
        """, (str(metab).strip().lower(), str(metab).strip(), smiles, source))
        conn.commit()
        conn.close()
        return True
    except Exception:
        return False


def db_log(compound_name: str, metabolite_name: str, smiles: str,
           status: str, via: str = "", error_msg: str = ""):
    """记录 API 查询日志"""
    try:
        conn = _get_db()
        conn.execute("""
            INSERT INTO api_log (compound_name, metabolite_name, smiles, status, via, error_msg)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (compound_name, metabolite_name, smiles, status, via, error_msg))
        conn.commit()
        conn.close()
    except Exception:
        pass


def db_get_stats() -> dict:
    """返回数据库统计信息"""
    try:
        conn = _get_db()
        c1 = conn.execute("SELECT COUNT(*) FROM compound_smiles").fetchone()[0]
        c2 = conn.execute("SELECT COUNT(*) FROM metabolite_smiles").fetchone()[0]
        c3 = conn.execute("SELECT COUNT(*) FROM api_log").fetchone()[0]
        total_size = os.path.getsize(_DB_PATH) if os.path.exists(_DB_PATH) else 0
        conn.close()
        return {
            "compound_name_records": c1,
            "metabolite_name_records": c2,
            "api_log_entries": c3,
            "db_file_size_mb": round(total_size / 1024 / 1024, 2),
            "db_path": _DB_PATH,
        }
    except Exception as e:
        return {"error": str(e)}


def import_smiles_from_csv(csv_path: str, name_col: str = "name",
                            smiles_col: str = "smiles") -> int:
    """
    从 CSV/Excel 文件批量导入 SMILES 数据到本地数据库。
    返回成功导入的记录数。
    """
    try:
        if csv_path.endswith('.csv') or csv_path.endswith('.tsv'):
            sep = '\t' if csv_path.endswith('.tsv') else ','
            df = pd.read_csv(csv_path, sep=sep)
        else:
            df = pd.read_excel(csv_path)
        count = 0
        for _, row in df.iterrows():
            n = str(row.get(name_col, '')).strip()
            s = str(row.get(smiles_col, '')).strip()
            if n and s and s not in ('NOT_FOUND', 'nan', ''):
                db_insert_name(n, s, source="user_upload")
                count += 1
        return count
    except Exception as e:
        raise RuntimeError(f"导入失败: {e}")


def _pubchem_name_to_smiles(name: str, timeout: int = 15) -> tuple:
    """
    多策略 PubChem 查询，按优先级尝试：
    1. 精确查询（PropertyTable）
    2. CID 中转查询（更稳定）
    3. 模糊匹配
    4. CAS 号查询
    5. 相似性搜索（threshold=85）
    Returns: (smiles, strategy, api_error_detail)
    """
    if not name:
        return "", "skipped", ""

    # 先做名称清理
    query_name = _name_to_query(name)
    if not query_name:
        return "", "invalid_name", ""

    # 尝试列表（按优先级）
    strategies = []

    # 策略1：直接精确查询
    strategies.append(("exact", query_name))

    # 策略2：如果清理后名称与原始不同，额外尝试原始名
    raw_name = str(name).strip()
    if query_name.lower() != raw_name.lower():
        strategies.append(("exact", raw_name))

    # 策略3：去除 "3-" "7-" 等数字前缀（如 Pelargonidin 3-Glucoside → Pelargonidin）
    cleaned = re.sub(r'^\d+[-\s]', '', query_name).strip()
    if cleaned and cleaned.lower() != query_name.lower():
        strategies.append(("exact", cleaned))

    # 策略4：CAS 号检测
    if re.match(r'^\d{2,7}-\d{2}-\d$', query_name):
        strategies.append(("cas", query_name))

    # 策略5：相似性搜索（兜底）
    strategies.append(("similarity", query_name))

    last_error = ""
    for strategy, q in strategies:
        try:
            if strategy == "exact":
                url = f"{_PUBCHEM_BASE}/compound/name/{quote(q)}/property/IsomericSMILES/JSON"
                result = _http_get_json(url, timeout=timeout)
                if "_error" in result:
                    last_error = result["_error"]
                    continue
                # 检查是否是 Fault 结构
                if "Fault" in result:
                    fault_msg = result.get('Fault', {})
                    last_error = f"Fault: {fault_msg.get('Message', fault_msg.get('message', str(fault_msg)))}"
                    continue
                props = result.get("PropertyTable", {}).get("Properties", [])
                if not props:
                    last_error = f"No properties returned for '{q}'"
                    continue
                # 尝试 IsomericSMILES，如果不存在则尝试 CanonicalSMILES，最后尝试 SMILES
                smiles = props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES") or props[0].get("SMILES", "")
                if smiles:
                    return smiles, "exact", ""
                last_error = f"Properties exist but no SMILES field. Keys: {list(props[0].keys())}"
                continue

            elif strategy == "cas":
                url = f"{_PUBCHEM_BASE}/compound/cas/{q}/property/IsomericSMILES/JSON"
                result = _http_get_json(url, timeout=timeout)
                if "_error" in result:
                    last_error = result["_error"]
                    continue
                if "Fault" in result:
                    last_error = f"CAS query fault: {result.get('Fault', {}).get('Message', '')}"
                    continue
                props = result.get("PropertyTable", {}).get("Properties", [])
                if props:
                    smiles = props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES") or props[0].get("SMILES", "")
                    if smiles:
                        return smiles, "cas", ""

            elif strategy == "similarity":
                # 先查 CID，再用 CID 查 SMILES（相似性搜索只返回 CID 列表）
                cid_url = f"{_PUBCHEM_BASE}/compound/name/{quote(q)}/cids/JSON?algorithm=similarity&threshold=85"
                cid_result = _http_get_json(cid_url, timeout=timeout)
                if "_error" in cid_result:
                    last_error = f"Similarity CID error: {cid_result['_error']}"
                    continue
                if "IdentifierList" not in cid_result:
                    last_error = f"Similarity: no IdentifierList in response"
                    continue
                cids = cid_result["IdentifierList"].get("CID", [])
                if not cids:
                    last_error = "Similarity: no CIDs found"
                    continue
                # 取第一个 CID 查 SMILES
                sm_url = f"{_PUBCHEM_BASE}/compound/cid/{cids[0]}/property/IsomericSMILES/JSON"
                sm_result = _http_get_json(sm_url, timeout=timeout)
                if "_error" in sm_result:
                    last_error = f"Similarity SMILES error: {sm_result['_error']}"
                    continue
                props = sm_result.get("PropertyTable", {}).get("Properties", [])
                if props:
                    smiles = props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES") or props[0].get("SMILES", "")
                    if smiles:
                        return smiles, "similarity", ""
                last_error = "Similarity: no SMILES found"
                continue

        except Exception as e:
            last_error = f"{strategy} exception: {e}"
            continue

    return "", "all_failed", last_error


def db_query_compound_with_fallback(name: str, metabolite: str = "") -> tuple:
    """
    带 fallback 的查询：
    1. 查本地 compound_name 数据库
    2. 查本地 metabolite 数据库
    3. 尝试 API（compound_name 优先）
    4. API 失败则尝试 metabolite 名称 API
    5. 所有结果自动写入本地数据库

    Returns: (smiles, source, via)
      source: 'local_db' | 'pubchem_api' | 'not_found'
      via: 'compound_name' | 'metabolite' | ''
    """
    invalid_vals = {None, '', '-', 'nan', 'NaN'}

    # 1. 本地查找 compound_name
    if name and str(name).strip() not in invalid_vals:
        sm, src = db_lookup_name(str(name).strip())
        if sm:
            db_log(name, metabolite, sm, "HIT_local_db", via="compound_name")
            return sm, "local_db", "compound_name"

    # 2. 本地查找 metabolite
    if metabolite and str(metabolite).strip() not in invalid_vals:
        sm, src = db_lookup_metabolite(str(metabolite).strip())
        if sm:
            db_log(name, metabolite, sm, "HIT_local_db", via="metabolite")
            return sm, "local_db", "metabolite"

    # 3. API 查询（compound_name 优先）- 多策略
    via = ""
    if name and str(name).strip() not in invalid_vals:
        via = "compound_name"
        smiles, strategy, api_err = _pubchem_name_to_smiles(str(name).strip())
        if smiles:
            db_insert_name(name, smiles, source=f"pubchem_{strategy}")
            db_log(name, metabolite, smiles, "HIT_api", via=f"compound_name:{strategy}")
            return smiles, "pubchem_api", "compound_name"
        # 未查到，记录真实错误原因
        db_log(name, metabolite, "", f"MISS_api", via="compound_name", error_msg=api_err or "no_result")

    # 4. API fallback：用 metabolite 名称查（同样多策略）
    if metabolite and str(metabolite).strip() not in invalid_vals:
        via = "metabolite"
        smiles, strategy, api_err = _pubchem_name_to_smiles(str(metabolite).strip())
        if smiles:
            db_insert_metabolite(metabolite, smiles, source=f"pubchem_{strategy}")
            db_log(name, metabolite, smiles, "HIT_api", via=f"metabolite:{strategy}")
            return smiles, "pubchem_api", "metabolite"
        db_log(name, metabolite, "", f"MISS_api", via="metabolite", error_msg=api_err or "no_result")

    # 5. 两者都查不到
    db_log(name, metabolite, "", "NOT_FOUND", via=via or "none")
    return "NOT_FOUND", "not_found", ""


# ===================== 缓存 =====================

# 简单的内存缓存
_smiles_cache = {}
_cache_max_size = 2000


def _smiles_cached(key: str, fn, *args):
    if key in _smiles_cache:
        return _smiles_cache[key]
    result = fn(*args)
    if len(_smiles_cache) < _cache_max_size:
        _smiles_cache[key] = result
    return result

# ===================== SMILES 查询 =====================

def _parse_first_name(name: str) -> str:
    """从 compound_name 字段中提取第一个合法名称（处理分号、空格等）"""
    if not name:
        return ''
    name = str(name).strip()
    # 取分号/逗号前第一段
    for sep in (';', ',', '/'):
        if sep in name:
            name = name.split(sep)[0].strip()
    return name


def _name_to_query(name: str) -> str:
    """
    将用户输入的化合物名称转换为适合 PubChem 查询的字符串。
    去掉括号内容、-等常见格式问题。
    """
    if not name:
        return ''
    name = str(name).strip()
    # 去掉括号及其内容，如 "(2R)-Naringenin" -> "Naringenin"
    name = re.sub(r'\([^)]*\)', '', name)
    # 去掉结尾的连字符和多余空格
    name = re.sub(r'\s+-\s*$', '', name).strip()
    # 多个空格变一个
    name = re.sub(r'\s+', ' ', name)
    return name


@st.cache_data(ttl=3600)
def query_smiles_by_name(compound_name: str) -> str:
    """
    通过 PubChem REST API 按化合物名称查询 SMILES（多策略 fallback）。
    使用 urllib（Python 内置），不依赖 requests。

    Args:
        compound_name: 原始化合物名称

    Returns:
        SMILES 字符串，查不到返回 "NOT_FOUND"
    """
    if not compound_name or str(compound_name).strip() in ('-', '', 'nan', 'NaN'):
        return "NOT_FOUND"

    first_name = _parse_first_name(compound_name)
    if not first_name:
        return "NOT_FOUND"

    # 多策略查询（compound_name 优先）
    smiles, strategy, api_err = _pubchem_name_to_smiles(first_name)
    return smiles if smiles else "NOT_FOUND"


def query_smiles_by_metabolite(metabolite_name: str) -> str:
    """
    通过 PubChem REST API 按代谢物名称查询 SMILES（多策略 fallback）。
    """
    if not metabolite_name or str(metabolite_name).strip() in ('-', '', 'nan'):
        return "NOT_FOUND"

    name = _name_to_query(str(metabolite_name).strip())
    if not name:
        return "NOT_FOUND"

    smiles, strategy, api_err = _pubchem_name_to_smiles(name)
    return smiles if smiles else "NOT_FOUND"


def batch_query_smiles(df: pd.DataFrame, name_col: str = 'compound_name',
                       progress_callback=None) -> pd.DataFrame:
    """批量查询 SMILES：本地 SQLite 缓存优先，找不到再调 API，结果自动落库。"""
    result_df = df.copy()
    if 'SMILES' not in result_df.columns:
        result_df['SMILES'] = 'PENDING'
    if 'SMILES_Source' not in result_df.columns:
        result_df['SMILES_Source'] = ''
    n_total = len(result_df)
    for i, (_, row) in enumerate(result_df.iterrows()):
        cur_smiles = row.get('SMILES', '')
        if cur_smiles and cur_smiles not in ('PENDING', 'NOT_IN_DATA', 'NO_NAME', 'NOT_FOUND'):
            if progress_callback:
                progress_callback(i + 1, n_total)
            continue
        metab = str(row.get('Metabolite', '')).strip()
        name = str(row.get(name_col, '')).strip()
        smiles, source, via = db_query_compound_with_fallback(name, metab)
        if not smiles or smiles in ('NOT_FOUND', '') or smiles == 'NOT_FOUND':
            if not name.strip() or name.strip() in ('-', 'nan', '', 'NaN'):
                smiles = 'NO_NAME'
                source = 'no_input'
            else:
                smiles = 'NOT_FOUND'
                source = 'not_found'
        result_df.at[_, 'SMILES'] = smiles
        result_df.at[_, 'SMILES_Source'] = source
        if progress_callback:
            progress_callback(i + 1, n_total)
    return result_df

# ===================== SwissTargetPrediction 靶点预测 =====================

def query_swiss_target_prediction(smiles: str, species: str = "Homo sapiens",
                                  max_retries: int = _MAX_RETRIES) -> tuple:
    """
    调用 SwissTargetPrediction 获取靶点预测结果（2026-05 改版后：表单提交 + HTML 解析）。

    Args:
        smiles: 化合物的 SMILES
        species: 物种（默认 "Homo sapiens"）
        max_retries: 重试次数

    Returns:
        (genes_str, target_count, error_msg)
        genes_str: 分号分隔的基因名列表（prob > 0.01），查不到返回 "NOT_FOUND" 或 "ERROR"
        target_count: 靶点数量
        error_msg: 错误信息（无错误为 ""）
    """
    if not smiles or smiles in ('NOT_FOUND', 'PENDING', '', 'nan'):
        return "NOT_FOUND", 0, "No SMILES provided"

    # 物种映射
    species_map = {
        "Homo sapiens": "Homo_sapiens",
        "Mus musculus": "Mus_musculus",
        "Rattus norvegicus": "Rattus_norvegicus",
    }
    org = species_map.get(species, "Homo_sapiens")

    import http.cookiejar
    from urllib.parse import urlencode
    from urllib.request import build_opener, HTTPCookieProcessor

    for attempt in range(max_retries + 1):
        try:
            # 使用 cookie jar 自动管理 session
            cookie_jar = http.cookiejar.CookieJar()
            ssl_ctx = _get_ssl_context()
            opener = build_opener(HTTPCookieProcessor(cookie_jar))

            # Step 0: GET 主页获取 session cookie
            main_req = Request("https://www.swisstargetprediction.ch/", headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            opener.open(main_req, timeout=15)

            # Step 1: 提交任务
            form_data = urlencode({
                'smiles': smiles,
                'organism': org,
                'ioi': '2'
            }).encode('ascii')

            submit_req = Request("https://www.swisstargetprediction.ch/predict.php", data=form_data, headers={
                'Content-Type': 'application/x-www-form-urlencoded',
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
                'Referer': 'https://www.swisstargetprediction.ch/',
                'Origin': 'https://www.swisstargetprediction.ch'
            })

            resp = opener.open(submit_req, timeout=60)
            submit_html = resp.read().decode('utf-8')

            # Step 2: 提取 job ID
            job_match = re.search(r'result\.php\?job=(\d+)', submit_html)
            if not job_match:
                # 可能是同步返回了结果页面
                if 'Target' in submit_html and 'Probability' in submit_html:
                    return _parse_stp_result_html(submit_html)
                return "ERROR", 0, "No job ID in response"

            job_id = job_match.group(1)

            # Step 3: 等待计算完成（最多等 90 秒）
            for wait in range(18):
                time.sleep(5)
                check_url = "https://www.swisstargetprediction.ch/result.php?job=%s&organism=%s" % (job_id, org)
                check_req = Request(check_url, headers={
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
                })
                resp = opener.open(check_req, timeout=30)
                check_html = resp.read().decode('utf-8')
                if 'Target' in check_html and 'Probability' in check_html:
                    return _parse_stp_result_html(check_html)
                if 'finished' in check_html.lower():
                    break

            return "ERROR", 0, "Calculation timeout (>90s)"

        except Exception as e:
            if attempt < max_retries:
                time.sleep(2)
                continue
            return "ERROR", 0, str(e)

    return "ERROR", 0, "Max retries exceeded"


def _parse_stp_result_html(html: str) -> tuple:
    """
    解析 SwissTargetPrediction 结果页面的 HTML 表格。

    Returns: (genes_str, target_count, error_msg)
    """
    try:
        # 找到表格
        table_match = re.search(r'<table[^>]*>(.*?)</table>', html, re.DOTALL)
        if not table_match:
            return "ERROR", 0, "No result table found"

        table = table_match.group(1)
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', table, re.DOTALL)

        if len(rows) < 2:
            return "ERROR", 0, "Table has no data rows"

        # 解析表头
        header_row = rows[0]
        headers = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', header_row, re.DOTALL)
        headers_clean = [re.sub(r'<[^>]+>', '', h).strip() for h in headers]

        # 找到 Common name 和 Probability 列的索引
        gene_idx = None
        prob_idx = None
        for i, h in enumerate(headers_clean):
            hl = h.lower()
            if 'common name' in hl or 'gene' in hl:
                gene_idx = i
            elif 'probability' in hl:
                prob_idx = i

        if gene_idx is None:
            gene_idx = 1  # 默认第二列是 Common name
        if prob_idx is None:
            prob_idx = len(headers_clean) - 2  # 倒数第二列

        # 解析数据行
        genes = []
        for row in rows[1:]:
            cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL)
            cleaned = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]

            if len(cleaned) > max(gene_idx, prob_idx):
                gene = cleaned[gene_idx].strip()
                try:
                    prob = float(cleaned[prob_idx])
                except (ValueError, IndexError):
                    prob = 0

                if gene and prob > _PROB_THRESHOLD:
                    genes.append(gene)

        genes = sorted(set(genes))
        genes_str = ";".join(genes)
        return genes_str, len(genes), ""

    except Exception as e:
        return "ERROR", 0, f"Parse error: {e}"


def batch_swiss_target_prediction(df: pd.DataFrame, smiles_col: str = 'SMILES',
                                   progress_callback=None,
                                   delay: float = _REQUEST_DELAY) -> tuple:
    """
    批量调用 SwissTargetPrediction API。

    Args:
        df: 输入 DataFrame（需含 SMILES 列）
        smiles_col: SMILES 列名
        progress_callback: 回调函数，接收 (current, total)
        delay: 请求间隔（秒）

    Returns:
        (result_df, errors_list)
    """
    result_df = df.copy()
    if 'Predicted_Targets' not in result_df.columns:
        result_df['Predicted_Targets'] = 'PENDING'
    if 'Target_Count' not in result_df.columns:
        result_df['Target_Count'] = 0

    n = len(result_df)
    errors = []
    for i, (_, row) in enumerate(result_df.iterrows()):
        smiles = str(row.get(smiles_col, ''))
        genes_str, count, err = query_swiss_target_prediction(smiles)
        result_df.at[_, 'Predicted_Targets'] = genes_str
        result_df.at[_, 'Target_Count'] = count
        if err:
            errors.append(f"Row {i}: {err}")
        if progress_callback:
            progress_callback(i + 1, n)
        if i < n - 1:
            time.sleep(delay)

    return result_df, errors


# ===================== gseapy 富集分析 =====================

def run_gseapy_enrichment(gene_list: list, organism: str = "human",
                          databases: list = None) -> dict:
    """
    使用 gseapy 做 KEGG + GO 富集分析。

    Args:
        gene_list: 基因名列表（来自 Predicted_Targets 列的分号分隔基因）
        organism: "human"
        databases: 要分析的数据库列表，默认包含 KEGG + 三个 GO 子库

    Returns:
        dict: {db_name: pd.DataFrame} 富集结果 DataFrame 字典
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError("请安装 gseapy: pip install gseapy")

    if databases is None:
        databases = [
            "KEGG_2021_Human",
            "GO_Biological_Process_2021",
            "GO_Molecular_Function_2021",
            "GO_Cellular_Component_2021",
        ]

    results = {}
    for db in databases:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=db,
                organism=organism,
                outdir=None,
                no_plot=True,
            )
            if enr.results is not None and len(enr.results) > 0:
                res_df = enr.results.copy()
                res_df['Database'] = db
                # 去除重复列名（保留第一个）- 只在合并前处理
                results[db] = res_df
        except Exception:
            results[db] = pd.DataFrame()
    return results


def merge_enrichment_results(enr_dict: dict) -> pd.DataFrame:
    """
    合并多个数据库的富集结果。
    """
    # 标准化列名映射（从 Enrichr 原始列名 → 标准列名）
    COL_MAP = {
        'term': 'Term', 'pathway': 'Term',
        'p-value': 'P_value', 'pvalue': 'P_value', 'p value': 'P_value',
        'old p-value': 'P_value',
        'adjusted p-value': 'Adjusted_P_value', 'adjusted pvalue': 'Adjusted_P_value',
        'adjusted p-value (bh)': 'Adjusted_P_value', 'fdr': 'Adjusted_P_value',
        'genes': 'Genes', 'gene': 'Genes',
        'overlap': 'Overlap', 'overlap/hits': 'Overlap',
    }
    dfs = []
    for db, df in enr_dict.items():
        if df is not None and len(df) > 0:
            df = df.copy()
            if 'Database' not in df.columns:
                df['Database'] = db
            # 标准化列名
            col_map = {}
            for c in df.columns:
                cl = c.lower().strip()
                if cl in COL_MAP:
                    col_map[c] = COL_MAP[cl]
            df = df.rename(columns=col_map)
            # 去重（保留第一个出现的）
            df = df.loc[:, ~df.columns.duplicated()]
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    merged = pd.concat(dfs, ignore_index=True)
    # 确保必要列存在
    for col in ['Term', 'Adjusted_P_value', 'Genes', 'Overlap']:
        if col not in merged.columns:
            merged[col] = None
    merged = merged.sort_values('Adjusted_P_value').reset_index(drop=True)
    return merged


# ===================== 可视化 =====================

def plot_enrichment_dotplot(enr_df: pd.DataFrame, top_n: int = 30) -> 'go.Figure':
    """
    绘制 Enrichr 风格的 dotplot（气泡图）。
    X轴 = Count/Overlap，Y轴 = Term，按 P-value 颜色编码。

    Args:
        enr_df: 富集结果 DataFrame，需含 Term, Adjusted_P_value, Genes, Overlap 列
        top_n: 显示前 N 条

    Returns:
        plotly.graph_objects.Figure
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("请安装 plotly: pip install plotly")

    if enr_df is None or len(enr_df) == 0:
        return None

    plot_df = enr_df.head(top_n).copy()

    # 解析 Overlap 列为数值（格式如 "15/200" -> 15）
    def parse_overlap(val):
        if pd.isna(val):
            return 0
        s = str(val)
        if '/' in s:
            try:
                return float(s.split('/')[0])
            except:
                return 0
        try:
            return float(val)
        except:
            return 0

    # 解析 Genes 列获取基因数量
    def gene_count(val):
        if pd.isna(val):
            return 0
        return len([g for g in str(val).split(';') if g.strip()])

    plot_df['Count'] = plot_df.apply(lambda r: gene_count(r.get('Genes', '')), axis=1)
    plot_df['X_val'] = plot_df.apply(lambda r: parse_overlap(r.get('Overlap', 0)), axis=1)

    # Term 名称截断
    plot_df['Term_short'] = plot_df['Term'].apply(lambda x: x[:50] + '...' if len(str(x)) > 50 else str(x))

    # 取 -log10(P) 用于颜色
    p_col = 'Adjusted_P_value' if 'Adjusted_P_value' in plot_df.columns else 'P_value'
    pvals = pd.to_numeric(plot_df[p_col], errors='coerce').fillna(1)
    pvals = pvals.replace(0, 1e-10)
    plot_df['-log10P'] = -np.log10(pvals)

    fig = go.Figure(go.Scatter(
        x=plot_df['X_val'],
        y=plot_df['Term_short'],
        mode='markers',
        marker=dict(
            size=plot_df['Count'] * 3 + 5,
            color=plot_df['-log10P'],
            colorscale='RdYlBu_r',
            showscale=True,
            colorbar=dict(title='-log10(Padj)'),
            sizemin=4,
        ),
        text=plot_df.apply(
            lambda r: f"<b>{r['Term']}</b><br>Count: {r['Count']}<br>AdjP: {r.get('Adjusted_P_value', r.get('P_value', 'N/A'))}<br>Genes: {r.get('Genes', '')}",
            axis=1
        ),
        hoverinfo='text',
    ))
    fig.update_layout(
        title=f'Enrichment Dotplot (Top {top_n} by Adjusted P-value)',
        xaxis_title='Count / Overlap',
        yaxis_title='Term',
        height=max(400, len(plot_df) * 22),
        width=900,
        font=dict(size=10),
        yaxis=dict(autorange='reversed'),
    )
