#!/usr/bin/env python3
"""模块 2b：GSE308433（lpWGS）拷贝数剂量 vs 肿瘤内在状态 / TCR 克隆性。

口径全部照《终点口径与分析预登记》第十四节执行，实现前已锁定，未查看任何关联结果。

为什么只做到染色体臂级
----------------------
低深度 WGS 的段跨度远大于单个基因，逐 bin 的 `.cna.seg`（单样本 2486 行）看似精细，
实际是同一事件的重复采样；`segments` 文件才是 FACETS 的共识段，用它做臂级汇总。
**禁止基因级拷贝数判定**——lpWGS 分辨率不支持，强行做等于把段边界当成生物学事实。

为什么用 copy.number 而不是事件名
----------------------------------
⚠️ 曾踩的坑：初稿把分子写成 `(GAIN bp - LOSS bp)`，用 FACETS 的事件字符串做加减。
但本数据集 **39/39 样本的 `LOSS` 事件出现 0 次**——FACETS 把拷贝数 1 的段记作
`HETD`、拷贝数 0 记作 `HOMDEL`，负向事件根本不叫 `LOSS`。按事件名取负值会得到
"只扩增、无缺失"的假象（实测 0/39 样本有任何臂缺失），而同一批 FACETS 判定本身正常
（33/39 众数 CN=2，中性段占 60%）。CN 是 FACETS 的原始判定，事件名只是符号约定。

为什么不用 logR 幅度
-------------------
`logR` 幅度几乎正比于肿瘤纯度，而补充文件里没有 FACETS 的 purity 输出
（tar 内 78 条 = 39 样本 × 2 文件，无 `.purity.txt`）。直接混用 logR 幅度，
会把"样本纯度低"误读成"拷贝数变化小"。

统计量
------
  cn_arm(arm)        = Σ (CN-2)×段长 / Σ 段长        相对二倍体，值域约 [-2, +2]
  cn_dosage_X        = median{ cn_arm(arm) : arm ∈ 签名 X 基因所占的不同臂 }
  cna_burden         = |CN-2| ≥ 1 的段覆盖 bp / 全常染色体覆盖 bp

`cn_dosage` 按**不同臂**取中位数而非按基因：同一臂上的多个基因不是独立观测，
按基因取会把一个臂的信息重复计入、虚增有效样本量。

输出
----
qc/GSE308433_arm_cnv.tsv          每文库 × 每臂
qc/GSE308433_dosage_scores.tsv    每文库的 6 个剂量分 + cna_burden + 众数 CN
qc/module2b_cnv_associations.tsv  关联结果（主比较 3 对 + 探索性 BH-FDR）
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
import tarfile
from pathlib import Path

# GRCh37 标准细胞遗传学臂界（公开参考常量，不做主观划分）
ARMS: dict[tuple[str, str], tuple[int, int]] = {
    ("1", "p"): (1, 123_400_000),          ("1", "q"): (123_400_001, 249_250_621),
    ("2", "p"): (1, 93_300_000),           ("2", "q"): (93_300_001, 243_199_373),
    ("3", "p"): (1, 90_900_000),           ("3", "q"): (90_900_001, 198_022_430),
    ("4", "p"): (1, 49_600_000),           ("4", "q"): (49_600_001, 191_154_276),
    ("5", "p"): (1, 48_500_000),           ("5", "q"): (48_500_001, 180_915_260),
    ("6", "p"): (1, 58_900_000),           ("6", "q"): (58_900_001, 170_805_979),
    ("7", "p"): (1, 60_100_000),           ("7", "q"): (60_100_001, 159_138_663),
    ("8", "p"): (1, 43_200_000),           ("8", "q"): (43_200_001, 145_138_636),
    ("9", "p"): (1, 42_200_000),           ("9", "q"): (42_200_001, 138_394_717),
    ("10", "p"): (1, 39_200_000),          ("10", "q"): (39_200_001, 135_534_747),
    ("11", "p"): (1, 48_300_000),          ("11", "q"): (48_300_001, 135_086_622),
    ("12", "p"): (1, 33_200_000),          ("12", "q"): (33_200_001, 133_851_895),
    ("13", "p"): (1, 16_500_000),          ("13", "q"): (16_500_001, 114_364_328),
    ("14", "p"): (1, 16_100_000),          ("14", "q"): (16_100_001, 107_043_718),
    ("15", "p"): (1, 17_500_000),          ("15", "q"): (17_500_001, 101_991_189),
    ("16", "p"): (1, 35_300_000),          ("16", "q"): (35_300_001, 90_338_345),
    ("17", "p"): (1, 22_700_000),          ("17", "q"): (22_700_001, 81_195_210),
    ("18", "p"): (1, 15_400_000),          ("18", "q"): (15_400_001, 78_077_248),
    ("19", "p"): (1, 14_200_000),          ("19", "q"): (14_200_001, 59_128_983),
    ("20", "p"): (1, 25_700_000),          ("20", "q"): (25_700_001, 63_025_520),
    ("21", "p"): (1, 10_900_000),          ("21", "q"): (10_900_001, 48_129_895),
    ("22", "p"): (1, 13_700_000),          ("22", "q"): (13_700_001, 51_304_566),
}
ARM_LABELS = [f"{c}{a}" for (c, a) in ARMS]
AUTOSOMES = {c for c, _ in ARMS}

SIGNATURES = {
    "ifng_response": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "antigen_presentation": "HALLMARK_ALLOGRAFT_REJECTION",
    "proliferation": "HALLMARK_E2F_TARGETS",
    "emt": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "hypoxia": "HALLMARK_HYPOXIA",
    "tnf_nfkb": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
}
# 事先指定的主比较 3 对，不得事后挑选
PRIMARY = ["ifng_response", "antigen_presentation", "emt"]
CLONALITY_METRICS = ["clonality", "top1_frac", "top10_frac", "clones_for_50pct"]
MIN_TCR_CELLS = 200

GSM_PREFIX = re.compile(r"^GSM\d+_")
SUFFIX = "_lpwgs_tumor.seg.gz"


# --------------------------------------------------------------------------- #
def read_gmt(p: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                out[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return out


def read_tsv(p: Path) -> list[dict]:
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_symbol_to_arm(gene_info: Path) -> dict[str, str]:
    """基因符号 → 染色体臂（如 BRAF→7q）。只保留常染色体 1–22，X/Y/MT 一律排除。"""
    with gzip.open(gene_info, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        fields = next(csv.reader([header], delimiter="\t"))
        rows = csv.DictReader(
            (ln for ln in fh if not ln.startswith("#")), fieldnames=fields, delimiter="\t"
        )
        out: dict[str, str] = {}
        for r in rows:
            sym = (r.get("Symbol_from_nomenclature_authority") or r.get("Symbol") or "").strip()
            chrom = (r.get("chromosome") or "").strip()
            loc = (r.get("map_location") or "").strip()
            if not sym or not loc:
                continue
            m = re.match(r"^(\d+)([pq])", loc)
            if not m:
                continue
            key = (m.group(1), m.group(2))
            if key not in ARMS or chrom != m.group(1):
                continue
            out.setdefault(sym.upper(), f"{key[0]}{key[1]}")
        return out


def rankdata(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    rx, ry = rankdata(x), rankdata(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return float("nan"), float("nan")
    rho = num / (dx * dy)
    z = 0.5 * math.log((1 + rho) / (1 - rho)) * math.sqrt(max(n - 3, 1))
    return rho, math.erfc(abs(z) / math.sqrt(2))


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), 1):
        k = m - rank + 1
        val = min(prev, pvals[idx] * m / k)
        adj[idx] = val
        prev = val
    return adj


# --------------------------------------------------------------------------- #
def parse_segments(tar_path: Path) -> dict[str, list[dict]]:
    """从 RAW.tar 里读出每文库的 FACETS 共识段。"""
    per_sample: dict[str, list[dict]] = {}
    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(SUFFIX):
                continue
            sample = member.name[: -len(SUFFIX)]
            fh = tf.extractfile(member)
            if fh is None:
                continue
            segs = []
            with gzip.open(fh, "rt", encoding="utf-8", errors="replace") as gz:
                gz.readline()  # sample chr start end event copy.number bins median
                for line in gz:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 8:
                        continue
                    segs.append({
                        "chrom": p[1], "start": int(p[2]), "end": int(p[3]),
                        "event": p[4], "cn": int(p[5]), "bins": int(p[6]),
                    })
            per_sample[sample] = segs
    return per_sample


def arm_stats(segs: list[dict]) -> tuple[dict[str, float], float, float, dict[str, float], int]:
    """返回 (每臂 cn_arm, cna_burden, 众数 CN, 事件 bp 占比, 有效段数)。

    cn_arm = Σ (CN-2)×bp / Σ bp，相对二倍体；段可能跨臂，按比例切分。
    """
    num: dict[str, float] = {a: 0.0 for a in ARM_LABELS}
    cov: dict[str, float] = {a: 0.0 for a in ARM_LABELS}
    cn_bp: dict[str, int] = {}
    ev_bp: dict[str, float] = {}
    n_useful = 0

    for s in segs:
        chrom = s["chrom"].replace("chr", "")
        if chrom not in AUTOSOMES:
            continue
        ov = {a: 0 for a in ("p", "q")}
        for arm in ("p", "q"):
            lo, hi = ARMS[(chrom, arm)]
            x = min(s["end"], hi) - max(s["start"], lo) + 1
            if x > 0:
                ov[arm] = x
        if ov["p"] + ov["q"] == 0:
            continue
        n_useful += 1
        span = ov["p"] + ov["q"]
        cn_bp[s["cn"]] = cn_bp.get(s["cn"], 0) + span
        ev_bp[s["event"]] = ev_bp.get(s["event"], 0.0) + span / 1_000_000.0
        for arm in ("p", "q"):
            if ov[arm]:
                label = f"{chrom}{arm}"
                cov[label] += ov[arm]
                num[label] += (s["cn"] - 2) * ov[arm]

    cn_arm = {a: num[a] / cov[a] if cov[a] > 0 else float("nan") for a in ARM_LABELS}
    tot_cov = sum(cov.values())
    altered = sum(v for k, v in cn_bp.items() if k != 2)
    burden = altered / tot_cov if tot_cov else float("nan")
    modal_cn = max(cn_bp.items(), key=lambda x: x[1])[0] if cn_bp else float("nan")
    ev_total = sum(ev_bp.values()) or 1.0
    return cn_arm, burden, float(modal_cn), {k: v / ev_total for k, v in ev_bp.items()}, n_useful


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--data-dir", default="data/GSE308433")
    ap.add_argument("--gmt", default="resources/hallmark_2024.1.Hs.gmt")
    ap.add_argument("--gene-info", default="resources/Homo_sapiens.gene_info.gz")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    qc = outdir / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    tar_path = Path(args.data_dir) / "GSE308433_RAW.tar"

    gmt = read_gmt(Path(args.gmt))
    sym2arm = load_symbol_to_arm(Path(args.gene_info))

    # 签名基因 → 所占的不同臂（只算在常染色体且 gene_info 里有定位的基因）
    sig_arms: dict[str, list[str]] = {}
    sig_cov: dict[str, tuple[int, int]] = {}
    for short, full in SIGNATURES.items():
        genes = gmt[full]
        arms = sorted({sym2arm[g] for g in genes if g in sym2arm})
        mapped = len([g for g in genes if g in sym2arm])
        sig_arms[short] = arms
        sig_cov[short] = (mapped, len(genes))
        print(f"{short:<22} 基因 {mapped}/{len(genes)} 有常染色体定位，涉及 {len(arms)} 条臂")

    segs_by_sample = parse_segments(tar_path)
    print(f"\n读入 {len(segs_by_sample)} 个文库的拷贝数分段")

    # ---- 臂级矩阵 ----
    with (qc / "GSE308433_arm_cnv.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["library", "sample"] + ARM_LABELS)
        for sample in sorted(segs_by_sample):
            lib = GSM_PREFIX.sub("", sample)
            cn_arm, *_ = arm_stats(segs_by_sample[sample])
            w.writerow([lib, sample] + [f"{cn_arm[a]:.4f}" if cn_arm[a] == cn_arm[a] else "nan"
                                        for a in ARM_LABELS])

    # ---- 剂量分 ----
    scores: dict[str, dict] = {}
    with (qc / "GSE308433_dosage_scores.tsv").open("w", encoding="utf-8", newline="") as fh:
        cols = (["library", "sample", "modal_cn", "cna_burden", "n_segments",
                 "n_arms_altered"] + [f"cn_dosage_{s}" for s in SIGNATURES])
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for sample in sorted(segs_by_sample):
            lib = GSM_PREFIX.sub("", sample)
            cn_arm, burden, modal_cn, ev_frac, nseg = arm_stats(segs_by_sample[sample])
            row: dict[str, object] = {"library": lib, "sample": sample,
                                      "modal_cn": int(modal_cn),
                                      "cna_burden": round(burden, 4), "n_segments": nseg}
            n_alt = sum(1 for a in ARM_LABELS
                        if cn_arm[a] == cn_arm[a] and abs(cn_arm[a]) >= 0.25)
            row["n_arms_altered"] = n_alt
            for short, arms in sig_arms.items():
                vals = [cn_arm[a] for a in arms if cn_arm[a] == cn_arm[a]]
                row[f"cn_dosage_{short}"] = round(sorted(vals)[len(vals) // 2], 4) if vals else "nan"
            scores[lib] = row
            w.writerow(row)

    # ---- 事件构成质控（保留事件名，仅作质控，不进统计量）----
    ev_keys = sorted({k for s in segs_by_sample.values() for k in arm_stats(s)[3]})
    with (qc / "GSE308433_event_fractions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["library", "sample"] + [f"frac_{k}" for k in ev_keys])
        for sample in sorted(segs_by_sample):
            ev_frac = arm_stats(segs_by_sample[sample])[3]
            w.writerow([GSM_PREFIX.sub("", sample), sample]
                       + [f"{ev_frac.get(k, 0.0):.4f}" for k in ev_keys])

    # ---- 关联分析 ----
    states: dict[str, dict] = {}
    sf = qc / "GSE308434_state_scores.tsv"
    if sf.exists():
        for r in read_tsv(sf):
            states[GSM_PREFIX.sub("", r["sample"])] = r
    tcr: dict[str, dict] = {}
    tf_ = qc / "GSE308435_tcr_clonality.tsv"
    if tf_.exists():
        for r in read_tsv(tf_):
            tcr[r.get("library") or GSM_PREFIX.sub("", r["sample"])] = r

    results = []

    def test(name_side: str, dose_side: str, pairs: list[tuple[str, dict]], primary: bool) -> None:
        xs = [float(a) for a, _ in pairs]
        ys = [float(b) for _, b in pairs]
        if len(xs) < 5:
            return
        rho, p = spearman(xs, ys)
        results.append({
            "state_side": name_side, "dose_side": dose_side, "n": len(xs),
            "spearman_rho": round(rho, 4), "p_approx": f"{p:.4g}",
            "is_primary": "yes" if primary else "exploratory",
        })

    # 主比较：状态分 × 同名签名剂量分
    for short in PRIMARY:
        pairs = []
        for lib, st in states.items():
            if lib not in scores:
                continue
            a, b = st.get(f"{short}_sd", ""), scores[lib].get(f"cn_dosage_{short}", "")
            if a in ("", "nan") or b in ("", "nan"):
                continue
            pairs.append((a, b))
        test(f"{short}_sd", f"cn_dosage_{short}", pairs, primary=True)

    # 探索性：另 3 个签名
    for short in SIGNATURES:
        if short in PRIMARY:
            continue
        pairs = []
        for lib, st in states.items():
            if lib not in scores:
                continue
            a, b = st.get(f"{short}_sd", ""), scores[lib].get(f"cn_dosage_{short}", "")
            if a in ("", "nan") or b in ("", "nan"):
                continue
            pairs.append((a, b))
        test(f"{short}_sd", f"cn_dosage_{short}", pairs, primary=False)

    # 探索性：cna_burden × TCR 克隆性（仅 TCR ≥200 子集，与模块 2 主分析同一批）
    tcr_ok = {k: v for k, v in tcr.items() if int(v["total_cells"]) >= MIN_TCR_CELLS}
    for cm in CLONALITY_METRICS:
        pairs = []
        for lib, tc in tcr_ok.items():
            if lib not in scores:
                continue
            a, b = scores[lib]["cna_burden"], tc.get(cm, "")
            if a in ("", "nan") or b in ("", "nan"):
                continue
            pairs.append((a, b))
        test("cna_burden", cm, pairs, primary=False)

    expl = [i for i, r in enumerate(results) if r["is_primary"] == "exploratory"]
    if expl:
        adj = bh([float(results[i]["p_approx"]) for i in expl])
        for i, a in zip(expl, adj):
            results[i]["p_fdr_bh"] = f"{a:.4g}"
    for r in results:
        r.setdefault("p_fdr_bh", "")

    order = {"yes": 0, "exploratory": 1}
    results.sort(key=lambda r: (order[r["is_primary"]], float(r["p_approx"])))
    with (qc / "module2b_cnv_associations.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["state_side", "dose_side", "n", "spearman_rho",
                                           "p_approx", "is_primary", "p_fdr_bh"], delimiter="\t")
        w.writeheader()
        w.writerows(results)

    print(f"\n文库数：CNV {len(scores)}、状态分 {len(states)}、TCR 全部 {len(tcr)}、"
          f"TCR≥{MIN_TCR_CELLS} {len(tcr_ok)}")
    print(f"\n主比较（事先指定）：")
    for r in results:
        if r["is_primary"] == "yes":
            print(f"  {r['state_side']:<22}× {r['dose_side']:<26} n={r['n']:<3} "
                  f"rho={r['spearman_rho']:>7} p≈{r['p_approx']}")
    print(f"\n探索性 {len(expl)} 项已做 BH-FDR → qc/module2b_cnv_associations.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
