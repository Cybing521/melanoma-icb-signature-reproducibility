#!/usr/bin/env python3
"""任务 11：给每个候选特征标注「常规组织样本能否拿到」的可及性等级。

要解决的缺口
------------
`docs/03` 第十三节已经做完**签名层门禁**（四队列 × 六主签名覆盖率 ≥95.5%），
但那是「在公开转录组里能不能算出来」。本脚本回答的是另一个问题：

> **在一个真实的临床黑色素瘤 FFPE 样本上，用常规检测手段，这个特征拿不拿得到？**

经验代理（避免凭空判断）
------------------------
"常规 IO 面板"没有统一定义。本脚本用 **GSE215868 的 NanoString PanCancer
IO 360 面板（784 个靶标）** 作为代理——它是真实临床研究里实际在用的
免疫肿瘤面板，且本项目已经下载了它的基因清单（`scripts/17`）。
用它的覆盖率来分级，比"我觉得这个签名很成熟/很不成熟"要硬。

分级规则（定义性的，与结果无关；本项目主线结论为阴性，无动机调高等级）
------------------------------------------------------------------------
| 等级 | 判据 | 含义 |
|---|---|---|
| **A 常规可及** | 面板覆盖 ≥60%，且代表基因不依赖冷门组织 | FFPE + 商业面板/IHC 即可拿到 |
| **B 条件可及** | 面板覆盖 30–60% | 需覆盖该签名的定制面板，或 FFPE+RNA-seq 深度足够 |
| **C 难可及** | 面板覆盖 <30% | 需全转录组，且 FFPE 降解后仍可能失真 |
| **D 非常规可及** | 平台层面就不可及 | 需 WGS / scTCR / 单细胞等非 bulk 手段 |

输出
----
qc/feature_accessibility.tsv
qc/feature_accessibility.md
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCR = ROOT / "scripts"
GMT = SCR / "hallmark_2024.1.Hs.gmt"
PANEL_MATRIX = ROOT / "expr" / "GSE215868_symbol_matrix.tsv.gz"

# 非表达类特征：平台层面就不可及，直接 D
NON_EXPRESSION = {
    "GSE308433_cn_dosage": ("D 非常规可及", "低深度 WGS (lpWGS) 拷贝数", "常规 FFPE/IHC 无对应检测"),
    "GSE308435_tcr_clonality": ("D 非常规可及", "scTCR 文库测序", "需配对 TCR 捕获，常规组织检测拿不到"),
    "GSE308433_state_sd": ("C 难可及", "单核 RNA-seq 状态分散度", "需单核测序深度"),
}


def load_hallmarks() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(GMT, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            name = p[0].replace("HALLMARK_", "").lower()
            out[name] = p[2:]
    return out


def main() -> int:
    hall = load_hallmarks()
    with gzip.open(PANEL_MATRIX, "rt") as fh:
        panel = {ln.split("\t", 1)[0] for ln in fh if not ln.startswith("gene_symbol")}
    n_panel = len(panel)

    rows = []

    def grade(cov: float) -> str:
        if cov >= 0.60:
            return "A 常规可及"
        if cov >= 0.30:
            return "B 条件可及"
        return "C 难可及"

    for feat, genes in sorted(hall.items()):
        present = [g for g in genes if g in panel]
        cov = len(present) / len(genes)
        rows.append({
            "feature": feat,
            "kind": "Hallmark 签名",
            "n_genes": len(genes),
            "n_in_routine_panel": len(present),
            "panel_coverage": round(cov, 4),
            "grade": grade(cov),
            "measured_by": "bulk RNA-seq / 靶向转录组面板",
            "note": "",
        })

    for name, (g, how, note) in NON_EXPRESSION.items():
        rows.append({"feature": name, "kind": "非表达类", "n_genes": 0,
                     "n_in_routine_panel": 0, "panel_coverage": 0.0,
                     "grade": g, "measured_by": how, "note": note})

    # IMPRES：只有 15 个基因，是全部候选里最可及的一类
    sys.path.insert(0, str(SCR))
    from impres_pairs import FEATURE_IDX, idx_to_pair
    ig = sorted({g for i in FEATURE_IDX for g in idx_to_pair(i)})
    ip = [g for g in ig if g in panel]
    rows.append({
        "feature": "IMPRES_15pairs", "kind": "基因对打分",
        "n_genes": len(ig), "n_in_routine_panel": len(ip),
        "panel_coverage": round(len(ip) / len(ig), 4),
        "grade": grade(len(ip) / len(ig)),
        "measured_by": "靶向转录组面板 / 甚至 RT-PCR 即可",
        "note": "15 个基因全为已知免疫检查点相关基因，是本项目候选里"
                "**临床可及性最高**的一类",
    })

    df = pd.DataFrame(rows).sort_values(
        ["grade", "panel_coverage"], ascending=[True, False])
    out = ROOT / "qc" / "feature_accessibility.tsv"
    df.to_csv(out, sep="\t", index=False)

    log = ["# 候选特征的临床可及性分级（任务 11）\n",
           "生成脚本 `scripts/22_feature_accessibility.py`。",
           f"**经验代理**：NanoString PanCancer IO 360 面板——真实临床研究在用的 IO 面板。"
           f"面板共 784 个条目，其中 6 个阳性对照 + 8 个阴性对照不参与打分，"
           f"**实测基因 {n_panel} 个**，本项目已下载其清单（`scripts/17`）。\n",
           "分级规则见脚本头，**定义性判据、与结果无关**——"
           "本项目主线结论为阴性，没有把等级调高的动机。\n"]

    log.append("## 一、分布\n")
    vc = df["grade"].value_counts()
    log.append("| 等级 | 特征数 |")
    log.append("|---|---:|")
    for g, c in vc.items():
        log.append(f"| {g} | {c} |")
    log.append("")

    log.append("## 二、逐特征\n")
    log.append("| 特征 | 类型 | 基因数 | 面板覆盖 | 等级 | 测量手段 |")
    log.append("|---|---|---:|---:|---|---|")
    for _, r in df.iterrows():
        log.append(f"| `{r['feature']}` | {r['kind']} | {r['n_genes']} | "
                   f"{r['panel_coverage']:.1%} | **{r['grade']}** | "
                   f"{r['measured_by']} |")
    log.append("")

    log.append("## 三、对交付物的直接影响\n")
    a = df[df["grade"].str.startswith("A")]["feature"].tolist()
    b = df[df["grade"].str.startswith("B")]["feature"].tolist()
    c = df[df["grade"].str.startswith("C")]["feature"].tolist()
    d = df[df["grade"].str.startswith("D")]["feature"].tolist()
    log.append(f"- **A 级 {len(a)} 个**：常规 FFPE + 商业面板即可拿到，"
               f"落地无额外门槛")
    log.append(f"- **B 级 {len(b)} 个**：需定制面板或更好的 FFPE RNA-seq")
    log.append(f"- **C 级 {len(c)} 个**：需全转录组且对降解敏感，"
               f"**不应写进模型卡的「常规可测」清单**")
    log.append(f"- **D 级 {len(d)} 个**：平台层面不可及，"
               f"**必须从任何声称临床可部署的清单中剔除**")
    log.append("")
    log.append("### 一条必须写进结论的观察\n")
    imp = df[df["feature"] == "IMPRES_15pairs"].iloc[0]
    log.append(f"**IMPRES 的 15 个基因在本项目所有候选里临床可及性最高**"
               f"（面板覆盖 {imp['panel_coverage']:.1%}，"
               f"{imp['n_in_routine_panel']}/{imp['n_genes']} 个基因在常规 IO 面板上）。")
    log.append("")
    log.append("但 `docs/07` 的结论是：**它在 79 例预登记锁定的独立验证上"
               "AUC 0.505 / 0.483，与随机无法区分**。")
    log.append("")
    log.append("这两条事实并列，正好构成项目最有说服力的一处对比：")
    log.append("")
    log.append("> **一个最容易测、测起来最便宜的签名，恰恰是最没有判别力的那个。**")
    log.append("> 它的失败不能归咎于「测不准」「测不到」「样本量不够」——"
               "这三项在本项目里都已逐一排除（不变量校验 8/8、面板覆盖 86.7%、n=79）。")
    log.append("")
    log.append("**同时必须写明反向的边界**：可及性分级**不能**用来预测判别力。"
               "本表给的是测量成本，不是诊断价值；"
               "把它读成「C 级特征更值得做」是误读。\n")

    md = ROOT / "qc" / "feature_accessibility.md"
    md.write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    print(f"\n[OK] {out}\n[OK] {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())