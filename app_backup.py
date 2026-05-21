"""
app_backup.py — Streamlit 数据分析应用的备份/恢复/下载/上传模块

功能：
  - save_backup(session_state_dict)  → 将所有 DataFrame 保存到 SQLite
  - load_backup()                     → 从 SQLite 恢复所有 DataFrame
  - get_backup_info()                 → 返回备份元信息
  - export_all_results(session_state_dict) → 导出所有结果为 Excel (多 sheet)
  - import_results(uploaded_file)     → 从上传的 Excel/SQLite 恢复

技术约束：仅依赖 Python 标准库 + pandas
"""

import io
import os
import sqlite3
import json
import traceback
from datetime import datetime, timezone, timedelta

import pandas as pd

# ── 需要备份的 session_state 键 ──────────────────────────────────
BACKUP_KEYS = [
    "diff_df", "diff_filtered", "annotated_df", "pharma_df",
    "abundance_df", "star_df", "star_df_v2",
    "smiles_df", "smiles_all_df", "target_df",
    "enr_df", "enrichment_df", "pharma_match_df",
    "pca_results", "plsda_results", "oplsda_results",
    "multi_group_results", "volcano_result", "_mat_data_groups",
]

# 有些值可能是 dict / list / None 而非 DataFrame，需要特殊处理
_DICT_LIKE_KEYS = {
    "pca_results", "plsda_results", "oplsda_results",
    "multi_group_results", "volcano_result", "_mat_data_groups",
}

BACKUP_VERSION = "1.0.0"
TZ_CST = timezone(timedelta(hours=8))


# ═══════════════════════════════════════════════════════════════════
# 内部工具函数
# ═══════════════════════════════════════════════════════════════════

def _is_dataframe(obj) -> bool:
    """判断对象是否为 DataFrame"""
    return isinstance(obj, pd.DataFrame) and not obj.empty


def _ensure_str_columns(df: pd.DataFrame) -> pd.DataFrame:
    """将 DataFrame 中所有列名转为 str（SQLite 不接受纯数字列名）"""
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    return df


def _df_to_bytes_for_sqlite(df: pd.DataFrame) -> str:
    """将 DataFrame 序列化为 JSON 字符串，存入 _blob 表"""
    # 需要处理 datetime 列
    df = df.copy()
    for col in df.columns:
        if pd.api.types.is_datetime64_any_dtype(df[col]):
            df[col] = df[col].astype(str)
        elif pd.api.types.is_integer_dtype(df[col]):
            df[col] = df[col].astype(object).where(df[col].notna(), None)
    # NaN → None
    df = df.where(df.notna(), None)
    return df.to_json(orient="records", force_ascii=False, date_format="iso")


def _bytes_to_df(json_str: str) -> pd.DataFrame:
    """从 JSON 字符串还原 DataFrame"""
    records = json.loads(json_str)
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def _serialize_value(value) -> str:
    """将 dict/list/其他类型序列化为 JSON 字符串"""
    if value is None:
        return ""
    return json.dumps(value, ensure_ascii=False, default=str, sort_keys=False)


def _deserialize_value(json_str: str):
    """从 JSON 字符串反序列化"""
    if not json_str:
        return None
    try:
        return json.loads(json_str)
    except (json.JSONDecodeError, TypeError):
        return json_str


# ═══════════════════════════════════════════════════════════════════
# 1. save_backup — 保存到 SQLite
# ═══════════════════════════════════════════════════════════════════

def save_backup(session_state_dict: dict, db_path: str = "analysis_backup.db") -> dict:
    """
    将 session_state 中的 DataFrame 保存到 SQLite 数据库。

    参数:
        session_state_dict: streamlit.session_state 的副本或 dict
        db_path: SQLite 数据库路径

    返回:
        dict: {"success": bool, "message": str, "tables": dict}
    """
    tables_info = {}
    conn = None
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # ── 清理旧表 ──
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        old_tables = [row[0] for row in cursor.fetchall()]
        for tbl in old_tables:
            cursor.execute(f'DROP TABLE IF EXISTS "{tbl}"')

        # ── 写入 DataFrame 到 _blob 表（JSON 序列化）──
        for key in BACKUP_KEYS:
            value = session_state_dict.get(key)
            if value is None:
                continue

            if key in _DICT_LIKE_KEYS:
                # dict / list 类型 → JSON 序列化
                json_str = _serialize_value(value)
                if json_str:
                    cursor.execute(
                        f'CREATE TABLE IF NOT EXISTS "{key}" (json_data TEXT)'
                    )
                    cursor.execute(f'INSERT INTO "{key}" VALUES (?)', (json_str,))
                    tables_info[key] = {"type": "json", "rows": 1}
            elif _is_dataframe(value):
                # DataFrame → JSON records
                df = _ensure_str_columns(value)
                # 转换 datetime 列
                for col in df.columns:
                    if pd.api.types.is_datetime64_any_dtype(df[col]):
                        df[col] = df[col].astype(str)
                # NaN → None
                df = df.where(df.notna(), None)
                records = df.to_dict(orient="records")

                if not records:
                    tables_info[key] = {"type": "dataframe", "rows": 0}
                    continue

                # 提取所有列名，创建表
                columns = list(df.columns)
                col_defs = ", ".join(f'"{col}" TEXT' for col in columns)
                cursor.execute(f'CREATE TABLE IF NOT EXISTS "{key}" ({col_defs})')

                # 批量插入
                placeholders = ", ".join(["?"] * len(columns))
                insert_sql = f'INSERT INTO "{key}" VALUES ({placeholders})'
                rows = []
                for rec in records:
                    row = [str(rec.get(c)) if rec.get(c) is not None else None for c in columns]
                    rows.append(row)
                cursor.executemany(insert_sql, rows)

                tables_info[key] = {"type": "dataframe", "rows": len(rows), "columns": columns}

        # ── 写入元信息 ──
        meta = {
            "version": BACKUP_VERSION,
            "backup_time": datetime.now(TZ_CST).isoformat(),
            "tables": tables_info,
            "total_tables": len(tables_info),
        }
        cursor.execute(
            'CREATE TABLE IF NOT EXISTS "_backup_meta" (meta_json TEXT)'
        )
        cursor.execute('INSERT INTO "_backup_meta" VALUES (?)',
                       (_serialize_value(meta),))

        conn.commit()

        return {
            "success": True,
            "message": f"备份成功，包含 {len(tables_info)} 张表",
            "tables": tables_info,
        }

    except Exception as e:
        tb = traceback.format_exc()
        return {"success": False, "message": f"备份失败: {e}\n{tb}", "tables": {}}
    finally:
        if conn:
            conn.close()


# ═══════════════════════════════════════════════════════════════════
# 2. load_backup — 从 SQLite 恢复
# ═══════════════════════════════════════════════════════════════════

def load_backup(db_path: str = "analysis_backup.db") -> dict:
    """
    从 SQLite 数据库恢复所有 DataFrame。

    返回:
        dict: {"success": bool, "message": str, "data": dict, "meta": dict}
    """
    result = {"success": False, "message": "", "data": {}, "meta": {}}

    if not os.path.exists(db_path):
        result["message"] = "备份文件不存在"
        return result

    conn = None
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 读取元信息
        cursor.execute('SELECT meta_json FROM "_backup_meta" LIMIT 1')
        row = cursor.fetchone()
        if row:
            result["meta"] = _deserialize_value(row[0]) or {}

        # 读取每个键的表
        for key in BACKUP_KEYS:
            try:
                cursor.execute(f'SELECT * FROM "{key}"')
                rows = cursor.fetchall()
                cols = [description[0] for description in cursor.description]

                if not rows:
                    continue

                # 判断是 JSON 序列化还是 DataFrame 表
                if cols == ["json_data"] and len(rows) == 1:
                    # dict/list 类型
                    result["data"][key] = _deserialize_value(rows[0][0])
                else:
                    # DataFrame
                    df = pd.DataFrame(rows, columns=cols)
                    # 转换类型
                    for col in df.columns:
                        # 尝试将数值字符串转为数字
                        try:
                            numeric_vals = pd.to_numeric(df[col], errors="coerce")
                            if numeric_vals.notna().sum() / max(len(df), 1) > 0.5:
                                df[col] = numeric_vals
                                continue
                        except Exception:
                            pass
                        # 尝试将日期字符串转为 datetime
                        try:
                            date_vals = pd.to_datetime(df[col], errors="coerce", utc=True, format="mixed")
                            if date_vals.notna().sum() / max(len(df), 1) > 0.5:
                                df[col] = date_vals.dt.tz_convert("Asia/Shanghai")
                                continue
                        except Exception:
                            pass
                    result["data"][key] = df

            except sqlite3.OperationalError:
                # 表不存在，跳过
                continue

        result["success"] = True
        result["message"] = f"恢复成功，共 {len(result['data'])} 个数据集"

    except Exception as e:
        tb = traceback.format_exc()
        result["message"] = f"恢复失败: {e}\n{tb}"
    finally:
        if conn:
            conn.close()

    return result


# ═══════════════════════════════════════════════════════════════════
# 3. get_backup_info — 获取备份元信息
# ═══════════════════════════════════════════════════════════════════

def get_backup_info(db_path: str = "analysis_backup.db") -> dict:
    """
    返回备份文件的元信息，包括时间、包含哪些表、行数等。

    返回:
        dict: {"exists": bool, "meta": dict, "file_size_kb": float}
    """
    info = {"exists": False, "meta": {}, "file_size_kb": 0.0}

    if not os.path.exists(db_path):
        return info

    info["exists"] = True
    info["file_size_kb"] = round(os.path.getsize(db_path) / 1024, 2)

    conn = None
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 读取元信息
        try:
            cursor.execute('SELECT meta_json FROM "_backup_meta" LIMIT 1')
            row = cursor.fetchone()
            if row:
                info["meta"] = _deserialize_value(row[0]) or {}
        except sqlite3.OperationalError:
            pass

        # 补充实际表信息
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE '_%'")
        tables = [row[0] for row in cursor.fetchall()]
        table_info = {}
        for tbl in tables:
            cursor.execute(f'SELECT COUNT(*) FROM "{tbl}"')
            count = cursor.fetchone()[0]
            table_info[tbl] = {"rows": count}

        info["tables"] = table_info

    except Exception:
        pass
    finally:
        if conn:
            conn.close()

    return info


# ═══════════════════════════════════════════════════════════════════
# 4. export_all_results — 导出为 Excel
# ═══════════════════════════════════════════════════════════════════

def export_all_results(session_state_dict: dict) -> bytes:
    """
    将所有 DataFrame 导出为 Excel 文件（多 sheet），返回 bytes。

    返回:
        bytes: Excel 文件内容
    """
    buffer = io.BytesIO()

    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        sheet_count = 0

        # 写入元信息 sheet
        meta_info = {
            "备份版本": BACKUP_VERSION,
            "导出时间": datetime.now(TZ_CST).strftime("%Y-%m-%d %H:%M:%S"),
            "包含的表": [],
        }

        for key in BACKUP_KEYS:
            value = session_state_dict.get(key)
            if value is None:
                continue

            # 安全化 sheet 名称（Excel 限制 31 字符）
            safe_name = key[:31]

            if _is_dataframe(value):
                df = _ensure_str_columns(value)
                # 转换 datetime 列
                for col in df.columns:
                    if pd.api.types.is_datetime64_any_dtype(df[col]):
                        df[col] = df[col].astype(str)
                df.to_excel(writer, sheet_name=safe_name, index=False)
                sheet_count += 1
                meta_info["包含的表"].append(f"{key} ({len(df)} 行)")

            elif key in _DICT_LIKE_KEYS and value is not None:
                # dict/list → 保存为文本 sheet
                text = json.dumps(value, ensure_ascii=False, indent=2, default=str)
                df_text = pd.DataFrame({"内容": [text]})
                df_text.to_excel(writer, sheet_name=safe_name, index=False)
                sheet_count += 1
                meta_info["包含的表"].append(f"{key} (JSON)")

        # 写入元信息
        meta_df = pd.DataFrame(list(meta_info.items()), columns=["项目", "值"])
        meta_df.to_excel(writer, sheet_name="_meta", index=False)

    buffer.seek(0)
    return buffer.getvalue()


# ═══════════════════════════════════════════════════════════════════
# 5. import_results — 从上传的 Excel/SQLite 恢复
# ═══════════════════════════════════════════════════════════════════

def import_results(uploaded_file, save_db: bool = True,
                   db_path: str = "analysis_backup.db") -> dict:
    """
    从上传的文件（Excel 或 SQLite）恢复数据。

    参数:
        uploaded_file: file-like 对象（如 st.file_uploader 返回的对象）
        save_db: 是否同时保存到 SQLite
        db_path: SQLite 保存路径

    返回:
        dict: {"success": bool, "message": str, "data": dict}
    """
    result = {"success": False, "message": "", "data": {}}

    if uploaded_file is None:
        result["message"] = "未提供文件"
        return result

    try:
        filename = getattr(uploaded_file, "name", "")
        file_ext = filename.lower().rsplit(".", 1)[-1] if "." in filename else ""

        if file_ext == "db":
            # SQLite 文件
            content = uploaded_file.read()
            buffer = io.BytesIO(content)
            # sqlite3 需要文件路径，写临时文件
            import tempfile
            with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
                tmp.write(content)
                tmp_path = tmp.name
            try:
                load_result = load_backup(tmp_path)
                result.update(load_result)
                if save_db and load_result["success"]:
                    # 同时保存到本地
                    with open(db_path, "wb") as f:
                        f.write(content)
            finally:
                os.unlink(tmp_path)

        elif file_ext in ("xlsx", "xls"):
            # Excel 文件
            content = uploaded_file.read()
            buffer = io.BytesIO(content)
            xl = pd.ExcelFile(buffer)
            data = {}

            for sheet_name in xl.sheet_names:
                if sheet_name == "_meta":
                    continue  # 跳过元信息 sheet

                # 映射回原始 key
                matched_key = None
                for key in BACKUP_KEYS:
                    if key[:31] == sheet_name or key == sheet_name:
                        matched_key = key
                        break

                if matched_key is None:
                    # 尝试匹配（sheet 名可能被截断）
                    matched_key = sheet_name

                df = xl.parse(sheet_name)

                # 检查是否是 JSON 文本 sheet
                if len(df) == 1 and "内容" in df.columns:
                    json_str = df["内容"].iloc[0]
                    data[matched_key] = _deserialize_value(json_str)
                else:
                    # 尝试还原类型
                    for col in df.columns:
                        try:
                            numeric_vals = pd.to_numeric(df[col], errors="coerce")
                            if numeric_vals.notna().sum() / max(len(df), 1) > 0.5:
                                df[col] = numeric_vals
                                continue
                        except Exception:
                            pass
                        try:
                            date_vals = pd.to_datetime(df[col], errors="coerce", utc=True, format="mixed")
                            if date_vals.notna().sum() / max(len(df), 1) > 0.5:
                                df[col] = date_vals
                                continue
                        except Exception:
                            pass
                    data[matched_key] = df

            result["data"] = data
            result["success"] = True
            result["message"] = f"从 Excel 恢复成功，共 {len(data)} 个数据集"

        else:
            result["message"] = f"不支持的文件格式: {file_ext}"

    except Exception as e:
        tb = traceback.format_exc()
        result["message"] = f"导入失败: {e}\n{tb}"

    return result


# ═══════════════════════════════════════════════════════════════════
# 便捷函数：用于 Streamlit 侧边栏集成
# ═══════════════════════════════════════════════════════════════════

def get_backup_db_bytes(db_path: str = "analysis_backup.db") -> bytes | None:
    """读取本地 SQLite 文件并返回 bytes，用于 Streamlit 下载按钮"""
    if not os.path.exists(db_path):
        return None
    with open(db_path, "rb") as f:
        return f.read()


def get_backup_summary(session_state_dict: dict) -> list[dict]:
    """返回当前 session_state 中可备份数据的摘要列表"""
    summary = []
    for key in BACKUP_KEYS:
        value = session_state_dict.get(key)
        if value is None:
            continue

        if _is_dataframe(value):
            summary.append({
                "key": key,
                "type": "DataFrame",
                "rows": len(value),
                "columns": len(value.columns),
            })
        elif key in _DICT_LIKE_KEYS and value is not None:
            if isinstance(value, dict):
                summary.append({
                    "key": key,
                    "type": "dict",
                    "keys": list(value.keys())[:5],
                })
            elif isinstance(value, list):
                summary.append({
                    "key": key,
                    "type": "list",
                    "length": len(value),
                })
            else:
                summary.append({
                    "key": key,
                    "type": type(value).__name__,
                })
    return summary


# ═══════════════════════════════════════════════════════════════════
# 测试入口
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # 简单测试
    print("=== app_backup.py 模块测试 ===")

    # 模拟 session_state
    test_state = {
        "diff_df": pd.DataFrame({"gene": ["TP53", "BRCA1"], "pvalue": [0.01, 0.05]}),
        "enrichment_df": pd.DataFrame({"term": ["GO:001"], "pvalue": [0.001]}),
        "pca_results": {"components": [[1, 2], [3, 4]], "variance": [0.5, 0.3]},
    }

    # 测试 save_backup
    r = save_backup(test_state, "test_backup.db")
    print(f"Save: {r}")

    # 测试 load_backup
    r = load_backup("test_backup.db")
    print(f"Load: {r['message']}")
    for k, v in r["data"].items():
        print(f"  {k}: {type(v).__name__} → {v if not isinstance(v, pd.DataFrame) else v.shape}")

    # 测试 get_backup_info
    info = get_backup_info("test_backup.db")
    print(f"Info: {info}")

    # 清理
    os.remove("test_backup.db")
    print("\nDone.")
