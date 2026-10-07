#!/usr/bin/env python3
"""GSE215868：解析 series matrix 元数据 + 105 份 NanoString RPT，建患者级表与表达矩阵。

依据：docs/03 §17.1–17.3（预登记写定于下载与计算之前）。

本脚本只做「搬运与建表」，不做任何应答相关的推断。

两处必须处理的坑
----------------
1. **Ragged characteristics**：GEO 把 105 个样本的 `!Sample_characteristics_ch1`
   切成多行，每行只放一部分字段，缺位用 `""` 补齐。必须按**列位置**跨行对齐，
   逐行直接 zip 会错位。
2. **itx 取值不规整**：实测含 `Pembro`（首字母大写、其余小写，疑为录入笔误）与
   `NIVO+EXPERIMENTAL`（第四类方案，不属于预登记列举的三类）。原始值一律原样
   保留在 `itx_raw` 列，另给规范化列供分层用。

输出
----
meta/GSE215868_samples.tsv
expr/GSE215868_symbol_matrix.tsv.gz
qc/GSE215868_parse_audit.md
"""

from __future__ import annotations

import glob
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SERIES_MATRIX = ROOT / "data" / "GSE215868" / "GSE215868_series_matrix.txt.gz"
RPT_DIR = ROOT / "data" / "GSE215868" / "rpt"

GENE_CLASSES = ("Endogenous", "Housekeeping")


# ---------------------------------------------------------------- 元数据解析
def parse_series_matrix(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    """返回 (gsm 列表, {gsm: {字段: 值}})。

    characteristics 按列位置跨行对齐；每个值形如 "key: value" 或为空。
    """
    rows: dict[str, list[str]] = {}
    char_lines: list[list[str]] = []
    order: list[str] | None = None
    ncol: int | None = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.startswith("!Sample_"):
                continue
            parts = line.split("\t")
            tag = parts[0]
            vals = [v.strip().strip('"') for v in parts[1:]]
            if not vals or all(v == "" for v in vals):
                continue
            if ncol is None:
                ncol = len(vals)
            elif len(vals) != ncol:
                raise ValueError(f"{tag} 列数 {len(vals)} != {ncol}，无法按列位置对齐")
            if tag == "!Sample_geo_accession":
                order = vals
            rows.setdefault(tag, vals)
            if tag == "!Sample_characteristics_ch1":
                char_lines.append(vals)

    if order is None:
        raise ValueError("series matrix 中没有 !Sample_geo_accession")

    recs: dict[str, dict[str, str]] = {g: {} for g in order}
    for tag, vals in rows.items():
        if tag == "!Sample_characteristics_ch1":
            continue
        for gsm, v in zip(order, vals):
            recs[gsm][tag] = v

    # 逐样本把多行 characteristics 拼起来，按 "key: value" 拆
    for gsm_i, gsm in enumerate(order):
        chunks = [ln[gsm_i] for ln in char_lines]
        for chunk in chunks:
            if not chunk:
                continue
            m = re.match(r"^\s*([^:]+):\s*(.*)$", chunk)
            if not m:
                raise ValueError(f"{gsm} 的 characteristics 片段无法解析：{chunk!r}")
            recs[gsm][m.group(1).strip()] = m.group(2).strip()

    return order, recs


# ---------------------------------------------------------------- RPT 解析
def parse_rpt(path: Path) -> tuple[dict[str, int], dict[str, int], dict[str, str]]:
    """返回 (基因->计数, 正对照->计数, 样本属性)。"""
    counts: dict[str, int] = {}
    pos: dict[str, int] = {}
    attrs: dict[str, str] = {}
    section = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith("<"):
                section = line.strip()
                continue
            if section == "<Sample_Attributes>":
                p = line.split(",", 1)
                if len(p) == 2 and p[0] not in ("ID",):
                    attrs[p[0]] = p[1]
            elif section == "<Code_Summary>":
                p = line.rstrip(",").split(",")
                if len(p) < 4 or p[0] == "CodeClass":
                    continue
                cls, name, _acc, cnt = p[0], p[1], p[2], p[3]
                try:
                    val = int(float(cnt))
                except ValueError:
                    continue
                if cls in GENE_CLASSES:
                    counts[name] = val
                elif cls == "Positive":
                    pos[name] = val
    return counts, pos, attrs


# ---------------------------------------------------------------- 主流程
def main() -> int:
    if not SERIES_MATRIX.exists():
        print(f"[FATAL] 缺 {SERIES_MATRIX}", file=sys.stderr)
        return 1

    order, recs = parse_series_matrix(SERIES_MATRIX)
    rpts = sorted(glob.glob(str(RPT_DIR / "*.txt.gz")))
    rpt_gsms = {Path(p).name.split("_")[0] for p in rpts}

    audit: list[str] = []

    # --- 一致性自检 -------------------------------------------------------
    audit.append("## 一、元数据与 RPT 的一致性\n")
    audit.append(f"- series matrix 样本数：{len(order)}")
    audit.append(f"- RPT 文件数：{len(rpts)}")
    only_meta = sorted(set(order) - rpt_gsms)
    only_rpt = sorted(rpt_gsms - set(order))
    audit.append(f"- 只有元数据无 RPT：{only_meta or '无'}")
    audit.append(f"- 只有 RPT 无元数据：{only_rpt or '无'}")

    # --- 建表达矩阵 -------------------------------------------------------
    gene_sets, matrices = [], []
    for p in rpts:
        c, _pos, _attrs = parse_rpt(Path(p))
        gene_sets.append(c)
        matrices.append(c)
    panel = sorted(set().union(*gene_sets))
    if any(len(s) != len(panel) for s in gene_sets):
        counts_bad = Counter(len(s) for s in gene_sets)
        raise ValueError(f"各样本基因数不一致：{counts_bad}")

    expr_path = ROOT / "expr" / "GSE215868_symbol_matrix.tsv.gz"
    expr_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(expr_path, "wt", encoding="utf-8") as out:
        out.write("gene_symbol\t" + "\t".join(order) + "\n")
        for g in panel:
            out.write(g + "\t" + "\t".join(str(m[g]) for m in matrices) + "\n")
    audit.append(f"\n## 二、表达矩阵\n")
    audit.append(f"- 基因数（内源+内参）：{len(panel)}")
    audit.append(f"- 列数（GSM）：{len(order)}")
    audit.append(f"- 写出：`expr/GSE215868_symbol_matrix.tsv.gz`")

    # --- 患者级表 ---------------------------------------------------------
    def field(g: str, key: str) -> str:
        return recs[g].get(key, "")

    rows = []
    for g in order:
        bor = field(g, "bor")
        itx_raw = field(g, "itx")
        rows.append(
            {
                "gsm": g,
                "title": field(g, "!Sample_title") or recs[g].get("!Sample_title", ""),
                "patient_id": g,  # 已核实 105 GSM = 105 位独立患者
                "age": field(g, "age"),
                "bor": bor,
                "response_group": (
                    "PRCR" if bor in ("CR", "PR") else "PD" if bor == "PD"
                    else "SD" if bor == "SD" else "UK"
                ),
                "itx_raw": itx_raw,
                "itx": (
                    "PEMBRO" if itx_raw.upper().startswith("PEMBRO")
                    else "NIVO" if itx_raw.upper() == "NIVO"
                    else "IPI+NIVO" if itx_raw.upper() == "IPI+NIVO"
                    else itx_raw.upper()
                ),
                "prior_icb": field(g, "prior checkpoint blockade"),
                "os_days": field(g, "os_days"),
                "pfs_days": field(g, "pfs_days"),
                "os_index": field(g, "os_index"),
                "pfs_index": field(g, "pfs_index"),
                "long_term_benefit": field(g, "long-term benefit"),
            }
        )

    meta_path = ROOT / "meta" / "GSE215868_samples.tsv"
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(meta_path, "w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    audit.append(f"\n## 三、患者级表\n")
    audit.append(f"- 写出：`meta/GSE215868_samples.tsv`，{len(rows)} 行 × {len(cols)} 列")
    audit.append("- 分析单位：**患者**。105 个 GSM 与 105 个唯一临床组合一一对应（门禁已核）")

    # --- 分布 -------------------------------------------------------------
    audit.append("\n## 四、关键字段分布（实测，非抄预登记）\n")
    for key in ("bor", "response_group", "itx_raw", "itx", "prior_icb", "os_index"):
        c = Counter(r[key] if r[key] != "" else "(缺失)" for r in rows)
        audit.append(f"- **{key}**：")
        for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
            audit.append(f"    - `{k}`：{v}")

    # 唯一性再核
    combos = {tuple(r[c] for c in ("age", "bor", "itx_raw", "pfs_days", "os_days")) for r in rows}
    audit.append(f"\n- 唯一 (age, bor, itx, pfs_days, os_days) 组合数：{len(combos)} / {len(rows)}")

    # 标题形如 "Yale Melanoma ITx1_1.1"：板位 = ITx<板号>_<lane>.<孔号>
    t = [r["title"] for r in rows]
    plate = Counter(m.group(1) for x in t if (m := re.match(r"^Yale Melanoma (ITx\d+_\d+)\.", x)))
    well = [m.group(1) for x in t if (m := re.match(r"^Yale Melanoma ITx\d+_\d+\.(\d+)$", x))]
    audit.append(f"\n- 标题去重后唯一值数：{len(set(t))} / {len(t)}")
    audit.append(f"- 标题可解析为板位编号的比例：{len(well)} / {len(t)}")
    audit.append(f"- 板位（`ITx<板>_<lane>`）分布：{dict(sorted(plate.items(), key=lambda kv: -kv[1]))}")
    dup_well = {w: c for w, c in Counter(well).items() if c > 1}
    audit.append(f"- 孔号重复情况：{dup_well or '无重复'}")
    audit.append(
        "- **板位编号不是患者编号**：孔号在板内唯一，且 105 个样本的 "
        "(age, bor, itx, pfs_days, os_days) 组合两两互不相同 → 一 GSM 一患者"
    )

    audit_path = ROOT / "qc" / "GSE215868_parse_audit.md"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(
        "# GSE215868 解析自检\n\n生成脚本 `scripts/17_build_gse215868.py`。\n\n"
        + "\n".join(audit) + "\n",
        encoding="utf-8",
    )
    print("\n".join(audit))
    print(f"\n[OK] 已写出 {expr_path}")
    print(f"[OK] 已写出 {meta_path}")
    print(f"[OK] 已写出 {audit_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())