#!/usr/bin/env python3
"""模块 2 主体：为 GSE308434 的每个肿瘤细胞核打状态分，汇总到活检级别。

设计要点
--------
1. **状态谱修正（2026-10-07）**：初稿把「耗竭」和「驻留样/TCF1 程序」列入肿瘤核状态谱，
   但这两个是 T 细胞状态，而 GSE308434 只含肿瘤细胞核（`cell type: tumor cell (nuclei)`，
   42/42 全部），肿瘤核里不存在这些细胞。改为肿瘤内在的可测状态，免疫侧的耗竭/克隆性
   由 GSE308435 的 TCR 指标承担（见 scripts/08_tcr_clonality.py）。
2. **基因集来源**：MSigDB Hallmark 2024.1.0（`resources/hallmark_2024.1.Hs.gmt`）。
   不用手写基因列表，避免主观挑基因并保证口径可复现。
3. **打分方法与统计量的选择（重要）**：CP10K 归一化 + log1p → 每个样本内按基因做 z-score
   → 签名内取均值。需要明确的是：**样本内 z-score 会让每个基因在该样本内的均值为 0，
   因此签名分的样本级均值在数学上恒等于 0，不含任何信息，不作为输出统计量。**
   有意义的是分布形状，故输出三个统计量：
     - `sd`：样本内状态分的离散度（该活检肿瘤群体的状态异质性）
     - `p90_p10`：90 分位减 10 分位的跨度（对离群核稳健的分布宽度）
     - `hi_pct`：得分高于「本样本均值 + 0.5 SD」的核占比（高状态核比例）
   这三个都是**样本内相对量**，因此可跨样本比较的是「状态构成/异质性」，
   而不是状态的绝对表达水平——后者本模块不回答，也不应被解读。
4. **质控阈值**：使用《终点口径与分析预登记》第十节已锁定的阈值，不在运行时调整。

输入：data/GSE308434/*_sn_counts.csv.gz、resources/hallmark_2024.1.Hs.gmt
输出：qc/GSE308434_state_scores.tsv（样本级）、qc/GSE308434_nucleus_state.tsv（核级，仅高分核）

用法：
  python3 scripts/09_state_scores.py --outdir . --data-dir data/GSE308434 --gmt resources/hallmark_2024.1.Hs.gmt
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

# ---- 已锁定的质控阈值（《终点口径与分析预登记》第十节），不得在运行时调整 ----
QC_MIN_UMI = 500
QC_MIN_GENES = 250
QC_MAX_MT_PCT = 20.0

# ---- 主状态面板：与 ICB 耐药机制直接相关（6 个） ----
# 注意 MSigDB 自 v7 起下架了 HALLMARK_ANTIGEN_PROCESSING_AND_PRESENTATION，
# Hallmark 2024.1.0 中不存在该集。抗原呈递改用 HALLMARK_ALLOGRAFT_REJECTION 作代理——
# 它是 Hallmark 里唯一携带完整 MHC/TAP 机器的集合（B2M、HLA-A/E/G、CD74、HLA-DR/DM/DO/QA、TAP1/2、PSMB10、TAPBP）。
# 代价：该集本身被 IFN-γ 通路强烈驱动，与 ifng_response 高度共线，两者不可当作独立证据使用。
PRIMARY_STATES = {
    "antigen_presentation": ["HALLMARK_ALLOGRAFT_REJECTION"],
    "ifng_response": ["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
    "proliferation": ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"],
    "emt": ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    "hypoxia": ["HALLMARK_HYPOXIA"],
    "tnf_nfkb": ["HALLMARK_TNFA_SIGNALING_VIA_NFKB"],
}

# ---- 探索性面板：单独报告，不参与主要结论 ----
SECONDARY_STATES = {
    "ifna_response": ["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
    "inflammatory": ["HALLMARK_INFLAMMATORY_RESPONSE"],
    "tgf_beta": ["HALLMARK_TGF_BETA_SIGNALING"],
    "il2_stat5": ["HALLMARK_IL2_STAT5_SIGNALING"],
    "glycolysis": ["HALLMARK_GLYCOLYSIS"],
    "oxphos": ["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
    "mtorc1": ["HALLMARK_MTORC1_SIGNALING"],
    "pi3k_akt_mtor": ["HALLMARK_PI3K_AKT_MTOR_SIGNALING"],
    "apoptosis": ["HALLMARK_APOPTOSIS"],
    "angiogenesis": ["HALLMARK_ANGIOGENESIS"],
    "p53": ["HALLMARK_P53_PATHWAY"],
}


def load_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = set(parts[2:])
    return sets


def load_sample(path: Path):
    with gzip.open(path, "rb") as fh:
        cells = fh.readline().decode("utf-8").rstrip("\n").split(",")[1:]
    df = pd.read_csv(path, header=None, skiprows=1)
    genes = df[0].astype(str).to_numpy()
    X = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
    return genes, X, cells


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    gmt = load_gmt(Path(args.gmt))
    all_states = {**PRIMARY_STATES, **SECONDARY_STATES}
    missing = {s for v in all_states.values() for s in v if s not in gmt}
    if missing:
        raise SystemExit(f"基因集文件中缺少：{sorted(missing)}")

    files = sorted(Path(args.data_dir).glob("*_sn_counts.csv.gz"))
    if args.limit:
        files = files[: args.limit]
    if not files:
        raise SystemExit(f"{args.data_dir} 下没有 *_sn_counts.csv.gz")

    qc = Path(args.outdir) / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    sample_rows = []
    for i, f in enumerate(files, 1):
        sample = f.name.replace("_sn_counts.csv.gz", "")
        genes, X, cells = load_sample(f)

        # ---- 质控（阈值已锁定）----
        n_umi = X.sum(axis=0)
        n_genes = (X > 0).sum(axis=0)
        mt_mask = np.array([g.startswith("MT-") for g in genes])
        mt_pct = (X[mt_mask].sum(axis=0) / np.maximum(n_umi, 1)) * 100.0
        keep = (n_umi >= QC_MIN_UMI) & (n_genes >= QC_MIN_GENES) & (mt_pct <= QC_MAX_MT_PCT)
        n_before, n_after = int(len(keep)), int(keep.sum())

        Xk = X[:, keep]
        # ---- 归一化 ----
        tot = Xk.sum(axis=0, keepdims=True)
        Xk = np.log1p(Xk / np.maximum(tot, 1.0) * 1e4)

        gidx = {g: k for k, g in enumerate(genes)}
        row = {
            "sample": sample,
            "n_nuclei_raw": n_before,
            "n_nuclei_qc": n_after,
            "qc_kept_pct": round(100.0 * n_after / max(n_before, 1), 2),
        }

        for state, sets in all_states.items():
            members: set[str] = set()
            for s in sets:
                members |= gmt[s]
            present = [gidx[g] for g in members if g in gidx]
            if len(present) < 10:
                row[f"{state}_sd"] = float("nan")
                row[f"{state}_p90_p10"] = float("nan")
                row[f"{state}_hi_pct"] = float("nan")
                row[f"{state}_genes"] = len(present)
                continue
            sub = Xk[present, :]
            # 样本内按基因 z-score 后取均值
            mu = sub.mean(axis=1, keepdims=True)
            sd = sub.std(axis=1, keepdims=True)
            sd[sd == 0] = 1.0
            z = (sub - mu) / sd
            score = z.mean(axis=0)
            # 注意：score 在样本内的均值恒为 0（基因 z-score 的性质），不可作为输出。
            row[f"{state}_sd"] = round(float(score.std()), 6)
            row[f"{state}_p90_p10"] = round(
                float(np.percentile(score, 90) - np.percentile(score, 10)), 6)
            row[f"{state}_hi_pct"] = round(
                100.0 * float((score > score.mean() + 0.5 * score.std()).mean()), 3)
            row[f"{state}_genes"] = len(present)

        sample_rows.append(row)
        print(
            f"[{i}/{len(files)}] {sample}: 核 {n_before}->{n_after} ({row['qc_kept_pct']}%) | "
            f"抗原呈递 sd={row['antigen_presentation_sd']:.4f} "
            f"IFNγ sd={row['ifng_response_sd']:.4f} "
            f"增殖 sd={row['proliferation_sd']:.4f} "
            f"EMT sd={row['emt_sd']:.4f}",
            flush=True,
        )

        out = qc / "GSE308434_state_scores.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sample_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sample_rows)
    print(f"\n样本级状态分已写入 {out}（{len(sample_rows)} 个样本）")
    print(f"质控阈值：UMI≥{QC_MIN_UMI}, genes≥{QC_MIN_GENES}, MT%≤{QC_MAX_MT_PCT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())