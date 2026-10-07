#!/usr/bin/env python3
"""把四个队列的原始表达文件统一整理成**基因符号 × 样本**矩阵。

要解决的坑
----------
1. **GSE91061 的三个矩阵（raw / fpkm / rld）首列不是行号，而是 Entrez Gene ID**
   （实测 min=1、max=101559451、n=22187），但列名留空——GEO 补充文件没随附注释表。
   本脚本用 NCBI `Homo_sapiens.gene_info` 做 Entrez→HGNC 官方符号映射并记录命中率。
2. 其余三个队列本就是符号，但**列名体系各不相同**（患者名 / 患者名.时间点 / GSM），
   不统一就无法和 `meta/patient_level.tsv` 对上。

列名统一为 **GSM 号**
-------------------
患者级纳排表 `meta/patient_level.tsv` 用 `gsm_list` 关联表达矩阵，
因此矩阵列一律用 GSM，原始列名另存 `meta/<cohort>_expr_columns.tsv` 备查。
这样任何一个模型脚本都只需按 GSM join，不必再处理 `Pt27A`、`Pt1.baseline` 这类特例。

表达值口径
----------
* GSE91061 用 **raw 未标准化计数**（fpkm/rld 已被 Illumina 流程标准化，跨样本方差被压缩，
  用于「样本内基因 z-score 后取签名均值」会引入处理痕迹）；fpkm 仅作敏感性分析备用。
* GSE294272 用 featureCounts **原始计数**。
* GSE78220 是 FPKM、GSE244982 是作者给的标准化值，只有这两种可用——已在
  `docs/03` 记录为跨队列比较的已知限制。

输出
----
expr/<cohort>_symbol_matrix.tsv.gz
meta/<cohort>_expr_columns.tsv
meta/GSE91061_symbol_map.tsv
logs/build_expr_<cohort>.log
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
import tarfile
from pathlib import Path

# 主分析使用的状态打分签名（与 scripts/09_state_scores.py 保持一致）
STATE_SETS_MAIN = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",   # 抗原呈递的代理，MSigDB v7 起原集下架
    "HALLMARK_E2F_TARGETS",           # 增殖
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]


def load_entrez_to_symbol(gene_info: Path) -> dict[str, str]:
    """Entrez Gene ID → HGNC 官方符号；官方符号为空时退回 NCBI 常用符号。

    注意：gene_info 的表头行本身就以 `#tax_id` 开头，若直接用
    `DictReader(ln for ln in fh if not ln.startswith('#'))`，表头会被过滤掉，
    结果把第一条数据行当表头、映射表全空。必须先单独读表头。
    """
    mapping: dict[str, str] = {}
    with gzip.open(gene_info, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        fieldnames = next(csv.reader([header], delimiter="\t"))
        rows = csv.DictReader(
            (ln for ln in fh if not ln.startswith("#")), fieldnames=fieldnames, delimiter="\t"
        )
        for r in rows:
            gid = (r.get("GeneID") or "").strip()
            if not gid:
                continue
            official = (r.get("Symbol_from_nomenclature_authority") or "").strip()
            fallback = (r.get("Symbol") or "").strip()
            s = official or fallback
            if s and s != "-":
                mapping.setdefault(gid, s)
    return mapping


def load_title_to_gsm(samples_tsv: Path) -> tuple[dict[str, str], dict[str, str]]:
    """从 meta/<cohort>_samples.tsv 取 gsm↔title 双向映射。"""
    t2g: dict[str, str] = {}
    g2t: dict[str, str] = {}
    with samples_tsv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gsm = (r.get("gsm") or "").strip()
            title = (r.get("title") or "").strip()
            if gsm and title:
                t2g.setdefault(title, gsm)
                g2t.setdefault(gsm, title)
    return t2g, g2t


def write_matrix(out: Path, samples: list[str], rows: list[tuple[str, list[str]]]) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["symbol"] + samples)
        for sym, vals in rows:
            w.writerow([sym] + vals)


def write_columns(out: Path, rows: list[list[str]]) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gsm", "original_column", "mapped"])
        w.writerows(rows)


def check_hallmark(symbols: set[str], gmt: Path) -> list[str]:
    """统计主分析签名在各矩阵上的覆盖情况——覆盖率不过半就不能用来打分。"""
    if not gmt.exists():
        return []
    with gmt.open(encoding="utf-8") as fh:
        sets = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    out = []
    for name in STATE_SETS_MAIN:
        genes = sets.get(name)
        if not genes:
            continue
        hit = len(genes & symbols)
        out.append(f"    {name:<46} {hit:>3}/{len(genes):<3} {hit / len(genes) * 100:5.1f}%")
    return out


# --------------------------------------------------------------------------- #
# GSE91061：Entrez 原始计数
# --------------------------------------------------------------------------- #
def build_gse91061(outdir: Path, data_dir: Path, gene_info: Path) -> None:
    src = data_dir / "GSE91061" / "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz"
    log: list[str] = []

    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        cols = next(csv.reader([header]))
        orig_samples = [c.strip() for c in cols[1:]]
        raw_rows = []
        for line in fh:
            if not line.strip():
                continue
            parts = next(csv.reader([line.rstrip("\n")]))
            raw_rows.append((parts[0].strip().strip('"'), parts[1:]))

    mapping = load_entrez_to_symbol(gene_info)
    log.append(f"GSE91061 raw 矩阵：{len(raw_rows)} 行 × {len(orig_samples)} 样本")

    kept: list[tuple[str, list[str]]] = []
    map_rows: list[list[str]] = []
    seen: set[str] = set()
    n_miss = 0
    for gid, vals in raw_rows:
        sym = mapping.get(gid)
        if not sym:
            n_miss += 1
            map_rows.append([gid, "", "no_symbol_in_gene_info"])
            continue
        if sym in seen:
            map_rows.append([gid, sym, "duplicate_symbol_dropped"])
            continue
        seen.add(sym)
        map_rows.append([gid, sym, "ok"])
        kept.append((sym, vals))
    log.append(f"Entrez→Symbol 命中 {len(raw_rows) - n_miss}，未命中 {n_miss}，"
               f"映射率 {(len(raw_rows) - n_miss) / len(raw_rows) * 100:.2f}%")
    log.append(f"去重后唯一符号 {len(seen)} 个")

    t2g, _ = load_title_to_gsm(outdir / "meta" / "GSE91061_samples.tsv")
    col_rows, samples = [], []
    for title in orig_samples:
        gsm = t2g.get(title, "")
        samples.append(gsm or title)
        col_rows.append([gsm, title, "ok" if gsm else "title_not_in_meta"])

    write_matrix(outdir / "expr" / "GSE91061_symbol_matrix.tsv.gz", samples, sorted(kept))
    write_columns(outdir / "meta" / "GSE91061_expr_columns.tsv", col_rows)
    with (outdir / "meta" / "GSE91061_symbol_map.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["entrez_id", "symbol", "status"])
        w.writerows(map_rows)

    n_gsm = sum(1 for c in col_rows if c[2] == "ok")
    log.append(f"样本列 {len(samples)} 个，其中映射到 GSM {n_gsm} 个")
    log.append("主分析签名覆盖：")
    log.extend(check_hallmark({s for s, _ in kept}, outdir / "resources" / "hallmark_2024.1.Hs.gmt"))

    (outdir / "logs" / "build_expr_GSE91061.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


# --------------------------------------------------------------------------- #
# GSE78220：xlsx，列名 PtN.timepoint
# --------------------------------------------------------------------------- #
def build_gse78220(outdir: Path, data_dir: Path) -> None:
    import openpyxl

    log: list[str] = []
    src = data_dir / "GSE78220" / "GSE78220_PatientFPKM.xlsx"
    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    it = ws.iter_rows(values_only=True)
    head = next(it)
    orig_samples = [str(c).strip() for c in head[1:] if c is not None]

    rows: dict[str, list[str]] = {}
    for r in it:
        if not r or r[0] in (None, ""):
            continue
        rows.setdefault(str(r[0]).strip(), [("" if v is None else str(v)) for v in r[1:]])
    wb.close()
    log.append(f"GSE78220 FPKM 矩阵：{len(rows)} 个符号 × {len(orig_samples)} 样本")

    t2g, _ = load_title_to_gsm(outdir / "meta" / "GSE78220_samples.tsv")
    samples, col_rows = [], []
    for col in orig_samples:
        base = col.split(".")[0]
        gsm = t2g.get(base, "")
        samples.append(gsm or col)
        col_rows.append([gsm, col, "ok" if gsm else "patient_not_in_meta"])

    n = len(samples)
    kept = [(sym, vals[:n]) for sym, vals in rows.items()]
    write_matrix(outdir / "expr" / "GSE78220_symbol_matrix.tsv.gz", samples, sorted(kept))
    write_columns(outdir / "meta" / "GSE78220_expr_columns.tsv", col_rows)
    log.append(f"样本列 {n} 个，映射到 GSM {sum(1 for c in col_rows if c[2] == 'ok')} 个")
    log.append("主分析签名覆盖：")
    log.extend(check_hallmark(set(rows), outdir / "resources" / "hallmark_2024.1.Hs.gmt"))

    (outdir / "logs" / "build_expr_GSE78220.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


# --------------------------------------------------------------------------- #
# GSE294272：RAW.tar 内逐样本 featureCounts，文件名即 GSM
# --------------------------------------------------------------------------- #
def build_gse294272(outdir: Path, data_dir: Path) -> None:
    log: list[str] = []
    tar_path = data_dir / "GSE294272" / "GSE294272_RAW.tar"
    counts: dict[str, dict[str, int]] = {}
    samples: list[str] = []

    # 该 RAW.tar 实际是未压缩 tar（magic 为 "ustar"），不能用 "r:gz"，交给 * 自动探测
    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf:
            if not member.isfile():
                continue
            gsm = member.name.split("_", 1)[0]
            fh = tf.extractfile(member)
            if fh is None:
                continue
            with gzip.open(fh, "rt", encoding="utf-8", errors="replace") as gz:
                head = gz.readline()  # # Program:featureCounts ...
                cols = gz.readline().rstrip("\n").split("\t")
                gidx = cols.index("Geneid")
                cidx = len(cols) - 1  # 最后一列是计数
                per_sample: dict[str, int] = {}
                for line in gz:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    per_sample[parts[gidx]] = int(float(parts[cidx]))
            counts[gsm] = per_sample
            samples.append(gsm)
            del head

    all_genes = sorted(set().union(*(set(v) for v in counts.values())))
    samples = sorted(samples)
    log.append(f"GSE294272 featureCounts：{len(samples)} 个样本，合并后 {len(all_genes)} 个基因")
    log.append(f"每样本基因数中位 {sorted(len(v) for v in counts.values())[len(samples) // 2]}")

    rows = [(g, [str(counts[s].get(g, 0)) for s in samples]) for g in all_genes]
    write_matrix(outdir / "expr" / "GSE294272_symbol_matrix.tsv.gz", samples, rows)
    write_columns(outdir / "meta" / "GSE294272_expr_columns.tsv",
                  [[s, f"{s}_*_featurecount.txt.gz", "ok"] for s in samples])
    log.append("主分析签名覆盖：")
    log.extend(check_hallmark(set(all_genes), outdir / "resources" / "hallmark_2024.1.Hs.gmt"))

    (outdir / "logs" / "build_expr_GSE294272.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


# --------------------------------------------------------------------------- #
# GSE244982：作者标准化值，列名即患者号
# --------------------------------------------------------------------------- #
def build_gse244982(outdir: Path, data_dir: Path) -> None:
    log: list[str] = []
    src = data_dir / "GSE244982" / "GSE244982_ProcessedData_bulkRNAseq.txt.gz"
    rows: dict[str, list[str]] = {}
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        orig_samples = [c.strip() for c in head[1:] if c.strip()]
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.setdefault(parts[0].strip(), parts[1:])

    log.append(f"GSE244982 标准化值矩阵：{len(rows)} 个符号 × {len(orig_samples)} 样本")
    _, g2t = load_title_to_gsm(outdir / "meta" / "GSE244982_samples.tsv")
    t2g, _ = load_title_to_gsm(outdir / "meta" / "GSE244982_samples.tsv")
    samples = [t2g.get(c, "") or c for c in orig_samples]
    col_rows = [[t2g.get(c, ""), c, "ok" if t2g.get(c) else "patient_not_in_meta"] for c in orig_samples]

    n = len(orig_samples)
    kept = [(sym, vals[:n]) for sym, vals in rows.items()]
    write_matrix(outdir / "expr" / "GSE244982_symbol_matrix.tsv.gz", samples, sorted(kept))
    write_columns(outdir / "meta" / "GSE244982_expr_columns.tsv", col_rows)
    log.append(f"样本列 {n} 个，映射到 GSM {sum(1 for c in col_rows if c[2] == 'ok')} 个")
    log.append("主分析签名覆盖：")
    log.extend(check_hallmark(set(rows), outdir / "resources" / "hallmark_2024.1.Hs.gmt"))

    (outdir / "logs" / "build_expr_GSE244982.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


BUILDERS = {
    "GSE91061": lambda o, d, g: build_gse91061(o, d, g),
    "GSE78220": lambda o, d, g: build_gse78220(o, d),
    "GSE294272": lambda o, d, g: build_gse294272(o, d),
    "GSE244982": lambda o, d, g: build_gse244982(o, d),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--gene-info", default="resources/Homo_sapiens.gene_info.gz")
    ap.add_argument("--cohort", default="all", choices=["all", *BUILDERS])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    data_dir = Path(args.data_dir)
    gene_info = Path(args.gene_info)
    (outdir / "expr").mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)

    todo = list(BUILDERS) if args.cohort == "all" else [args.cohort]
    for c in todo:
        if not gene_info.exists():
            print(f"[跳过 {c}] 缺少 {gene_info}", file=sys.stderr)
            continue
        BUILDERS[c](outdir, data_dir, gene_info)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
