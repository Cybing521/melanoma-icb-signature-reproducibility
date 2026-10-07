#!/usr/bin/env python3
"""GEO 数据门禁：拉取 series_matrix 头部元数据，核对样本—患者映射与临床标签字段。

只下载 series_matrix 的元数据部分（KB 级），不下载表达矩阵、不下载原始测序数据。
输出：
  meta/<GSE>_samples.tsv   每样本一行，含 GSM、标题、来源、全部 characteristics 键值
  meta/data_gate_summary.json
  logs/data_gate_report.md

用法：
  python3 scripts/00_fetch_geo_metadata.py --outdir /path/to/project
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

# 方案第五节列出的队列，附带方案中登记的样本量，用于交叉核对
SERIES = [
    ("GSE308434", "纵向 ICB 单核转录组", 42),
    ("GSE294273", "淋巴结 scRNA-seq", 12),
    ("GSE294272", "淋巴结 bulk RNA-seq", 29),
    ("GSE78220", "治疗前 bulk ICB", 28),
    ("GSE91061", "ICB bulk（含两瘤种）", 109),
    ("GSE244982", "进展后 bulk ICB", 41),
]

RESPONSE_KEY_RE = re.compile(
    r"respons|recist|clinical[_ ]?benefit|best[_ ]?response|treatment[_ ]?outcome|"
    r"patient[_ ]?response|therapy[_ ]?response",
    re.I,
)
PATIENT_KEY_RE = re.compile(
    r"patient|subject|case[_ ]?id|individual|donor|specimen|patient_id", re.I
)
TIMEPOINT_KEY_RE = re.compile(
    r"timepoint|time[_ ]?point|visit|pre|post|on[_ ]?treatment|baseline|"
    r"treatment|week|month|day|stage|collection",
    re.I,
)


def series_matrix_url(gse: str) -> str:
    bucket = f"{gse[:-3]}nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/matrix/"
        f"{gse}_series_matrix.txt.gz"
    )


def fetch(url: str, retries: int = 3, timeout: int = 90) -> bytes:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-metadata-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, OSError) as exc:  # 网络类错误重试
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"下载失败 {url}: {last}")


def parse_header(raw: bytes) -> tuple[dict[str, str], dict[str, list[list[str]]], int]:
    """拆分 series_matrix 头部。

    series_matrix 每行以制表符分隔，第 1 列是字段名，其余列依次对应各样本；
    characteristics 等字段会出现多行，每行仍与样本列一一对齐。
    返回 (series 字段, sample 字段 -> 行列表, 样本数)。
    """
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    series: dict[str, str] = {}
    sample_rows: dict[str, list[list[str]]] = {}
    n_samples = 0
    for line in text.splitlines():
        if not line.startswith("!"):
            continue  # 数据矩阵本体不在本次处理范围
        parts = line.rstrip("\r").split("\t")
        key = parts[0][1:]
        values = parts[1:]
        if key.startswith("Series_"):
            series[key] = values[0] if values else ""
        elif key.startswith("Sample_"):
            sample_rows.setdefault(key, []).append(values)
            n_samples = max(n_samples, len(values))
    return series, sample_rows, n_samples


def split_characteristics(body: str) -> dict[str, str]:
    """把 "key: value; key: value" 拆成键值对，键统一小写。"""
    out: dict[str, str] = {}
    for chunk in body.split(";"):
        chunk = chunk.strip()
        if not chunk or ":" not in chunk:
            continue
        key, _, value = chunk.partition(":")
        out[key.strip().lower()] = value.strip()
    return out


def build_sample_table(
    sample_rows: dict[str, list[list[str]]], n_samples: int, titles: list[str]
) -> list[dict[str, str]]:
    """按样本列位置对齐；多行字段逐行合并，characteristics 拆成独立键。"""
    records: list[dict[str, str]] = [{"char": ""} for _ in range(n_samples)]
    chars: list[dict[str, str]] = [{} for _ in range(n_samples)]

    def cell_text(cell: str) -> str:
        """GEO 部分字段取值带 "样本标题: " 前缀与引号，按已知标题剥离并去引号。"""
        for title in titles:
            if title and cell.startswith(f"{title}:"):
                return cell[len(title) + 1 :].strip().strip('"').strip()
        return cell.strip().strip('"').strip()

    for key, rows in sample_rows.items():
        for row in rows:
            for idx in range(min(n_samples, len(row))):
                text = cell_text(row[idx])
                if key == "Sample_characteristics_ch1":
                    chars[idx].update(split_characteristics(text))
                elif key == "Sample_geo_accession":
                    records[idx]["gsm"] = text
                elif key == "Sample_title":
                    records[idx]["title"] = text
                elif key == "Sample_source_name_ch1":
                    records[idx]["source"] = text
                else:
                    records[idx][key] = text

    rows_out = []
    for rec, kv in zip(records, chars):
        rec = dict(rec)
        rec["characteristics"] = "; ".join(f"{k}={v}" for k, v in kv.items())
        rec.update({f"char::{k}": v for k, v in kv.items()})
        rows_out.append(rec)
    return rows_out


def audit_series(gse: str, role: str, expected_n: int, outdir: Path) -> dict:
    record: dict = {
        "gse": gse,
        "role": role,
        "registered_n": expected_n,
        "matrix_n": None,
        "status": "pending",
        "error": None,
        "response_keys": [],
        "patient_keys": [],
        "timepoint_keys": [],
        "sample_value_preview": {},
    }
    try:
        raw = fetch(series_matrix_url(gse))
        series, sample_rows, n_samples = parse_header(raw)
        title_row = sample_rows.get("Sample_title", [[]])[0] if sample_rows.get("Sample_title") else []
        rows = build_sample_table(sample_rows, n_samples, list(title_row))
    except Exception as exc:  # noqa: BLE001 - 记录后继续其他队列
        record.update(status="fetch_failed", error=str(exc))
        return record

    record["matrix_n"] = len(rows)
    record["series_title"] = series.get("Series_title", "")
    record["unique_titles"] = len({r["title"] for r in rows})

    all_keys: set[str] = set()
    values_by_key: dict[str, set[str]] = {}
    for row in rows:
        for key, value in row.items():
            if key.startswith("char::"):
                all_keys.add(key[6:])
                values_by_key.setdefault(key[6:], set()).add(value)
    record["characteristic_keys"] = sorted(all_keys)
    record["response_keys"] = sorted(k for k in all_keys if RESPONSE_KEY_RE.search(k))
    record["patient_keys"] = sorted(k for k in all_keys if PATIENT_KEY_RE.search(k))
    record["timepoint_keys"] = sorted(k for k in all_keys if TIMEPOINT_KEY_RE.search(k))
    record["sample_value_preview"] = {
        k: sorted(v)[:12] for k, v in values_by_key.items() if RESPONSE_KEY_RE.search(k)
    }
    record["status"] = "ok" if len(rows) == expected_n else "n_mismatch"

    out = outdir / "meta" / f"{gse}_samples.tsv"
    char_cols = [f"char::{k}" for k in sorted(all_keys)]
    cols = ["gsm", "title", "source", "characteristics"] + char_cols
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append(
            "\t".join(str(row.get(c, "") or "").replace("\t", " ") for c in cols)
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return record


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="项目根目录（含 meta/ logs/）")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    (outdir / "meta").mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)

    records = []
    for gse, role, expected in SERIES:
        rec = audit_series(gse, role, expected, outdir)
        records.append(rec)
        print(
            f"[{rec['status']:>12}] {gse} {role}: 登记 {expected} / 实得 {rec['matrix_n']}",
            flush=True,
        )
        if rec["error"]:
            print(f"             错误: {rec['error']}", flush=True)

    (outdir / "meta" / "data_gate_summary.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    md = ["# GEO 数据门禁报告", ""]
    md.append(f"生成时间：服务器执行，核查 {len(records)} 个队列。\n")
    md.append("| 队列 | 角色 | 方案登记样本量 | series_matrix 实得 | 状态 | 响应标签字段 | 患者标识字段 |")
    md.append("|---|---|---:|---:|---|---|---|")
    for r in records:
        md.append(
            f"| {r['gse']} | {r['role']} | {r['registered_n']} | "
            f"{r['matrix_n'] if r['matrix_n'] is not None else 'NA'} | {r['status']} | "
            f"{', '.join(r['response_keys']) or '—'} | {', '.join(r['patient_keys']) or '—'} |"
        )
    md.append("")
    for r in records:
        if r["sample_value_preview"]:
            md.append(f"## {r['gse']} 标签取值")
            for key, vals in r["sample_value_preview"].items():
                md.append(f"- `{key}`: {', '.join(vals)}")
            md.append("")
    (outdir / "logs" / "data_gate_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"\n报告已写入 {outdir / 'logs' / 'data_gate_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())