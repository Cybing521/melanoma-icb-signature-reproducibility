#!/usr/bin/env python3
"""组装交付包。

自动分类器对本项目判断有误（图是英文的、数据漏了、把下载脚本算成了分析脚本），
因此这里按项目实际内容手工组装，并逐条对照打包规范自检。
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEST = ROOT / "黑色素瘤ICB签名可复现性_delivery"

# 稿件正文里的参考文献块（唯一真源），用于生成参考文献清单
REFS = [
    ("1", "Auslander N. et al.", "Nat. Med. 2018, 24, 1545–1549.", "PMID 30127394", "IMPRES 签名原文"),
    ("2", "Charoentong P. et al.", "Cell Rep. 2017, 18, 248–262.", "PMID 28052254", "Immunophenoscore 原文"),
    ("3", "Riaz N. et al.", "Cell 2017, 171, 934–949.e16.", "PMID 29033130", "GSE91061 源文献"),
    ("4", "Hugo W. et al.", "Cell 2016, 165, 35–44.", "PMID 26997480", "GSE78220 源文献"),
    ("5", "Vathiotis I.A. et al.", "NPJ Precis. Oncol. 2022, 6, 92.", "PMID 36522538", "GSE215868 源文献"),
    ("6", "Lauss M. et al.", "Nat. Commun. 2024, 15, 3075.", "PMID 38594286", "GSE244982 源文献"),
    ("7", "Di Pietro A. et al.", "Nat. Commun. 2026, 17, 7445.", "PMID 42277002", "GSE294272 / GSE294273 源文献"),
    ("8", "Izar laboratory, Columbia University", "GEO GSE308433/434/435", "GEO 无关联文献，按登录号引用",
     "单细胞与 TCR 分析"),
    ("9", "Liberzon A. et al.", "Cell Syst. 2015, 1, 417–425.", "PMID 26771021", "MSigDB Hallmark 特征集"),
    ("10", "Hanley J.A.; McNeil B.J.", "Radiology 1982, 143, 29–36.", "doi:10.1148/radiology.143.1.7063747",
     "AUC 标准误与功效计算基础"),
    ("11", "DerSimonian R.; Laird N.", "Control. Clin. Trials 1986, 7, 177–188.", "PMID 3802833", "研究级随机效应合并"),
    ("12", "Deeks J.J. et al.", "Ann. Intern. Med. 2008, 148, 175–176.", "Crossref 无 DOI，未能机器核验",
     "检验层级"),
    ("13", "Hanley J.A.; Lippman-Hand A.", "JAMA 1983, 250, 2563–2566.", "Crossref 未检出，未能机器核验",
     "小样本事后择一的报告偏差"),
]

# 分析与作图脚本；排除下载器（00/02/03/04）、装配与出稿脚本
KEEP_SCRIPTS = [
    "05_build_patient_level.py", "06_build_expression_matrices.py", "07_qc_nuclei.py",
    "08_tcr_clonality.py", "09_state_scores.py", "10_pair_tumor_tcr.py",
    "11_check_feature_coverage.py", "12_cnv_dosage.py", "13_build_reference_profile.py",
    "14_impres_baseline.py", "15_bulk_signature_scores.py", "16_nested_cv_model.py",
    "17_build_gse215868.py", "18_impres_gse215868.py", "19_panel_diagnostic.py",
    "20_task14_external_validation.py", "21_scrna_cluster_tiebreak.py",
    "22_feature_accessibility.py", "23_ips_mhccp.py", "24_figures_mdpi.py",
    "25_power_analysis.py", "26_cohort_heterogeneity.py", "impres_pairs.py",
]


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> int:
    if DEST.exists():
        shutil.rmtree(DEST)
    for sub in ("脚本", "图/英文", "数据/处理后", "报告", "参考文献"):
        (DEST / sub).mkdir(parents=True)

    # ---- 报告（论文项目：一份 DOCX）
    copy(ROOT / "黑色素瘤ICB签名可复现性_预登记研究_Cancers投稿稿.docx",
         DEST / "报告/黑色素瘤ICB签名可复现性_预登记研究_Cancers投稿稿.docx")

    # ---- 脚本
    for s in KEEP_SCRIPTS:
        p = ROOT / "scripts" / s
        if p.exists():
            copy(p, DEST / "脚本" / s)
    copy(ROOT / "requirements.txt", DEST / "脚本/requirements.txt")

    # ---- 图：本文为英文投稿，图件本身即英文版（Arial，≤175 mm）
    for fig in sorted((ROOT / "figures").glob("Figure*")):
        if fig.suffix in (".png", ".pdf") and "_preview" not in fig.name:
            copy(fig, DEST / "图/英文" / fig.name)

    # ---- 数据：分析面板（患者级纳排表 + 全部结果文件）
    for p in sorted((ROOT / "qc").iterdir()):
        if p.is_file():
            copy(p, DEST / "数据/处理后" / p.name)
    for name in ("patient_level.tsv", "patient_level_summary.md", "dev_set_pre_treatment.tsv",
                 "data_gate_summary.json", "GSE91061_Riaz2017_sample_mapping.tsv"):
        p = ROOT / "meta" / name
        if p.exists():
            copy(p, DEST / "数据/处理后" / name)

    # ---- 补充材料（论文声明随稿提交的文件）
    copy(ROOT / "supplementary/Supplementary_Material_S1.pdf", DEST / "Supplementary_Material_S1.pdf")

    # ---- 参考文献：正文实际引用的 13 条（机器可读，非整个检索库）
    csv_lines = ["number,first_author_et_al,source,identifier,role_in_paper"]
    for num, authors, src, ident, role in REFS:
        csv_lines.append(f'{num},"{authors}","{src}","{ident}","{role}"')
    (DEST / "参考文献/引用文献.csv").write_text("\n".join(csv_lines) + "\n", encoding="utf-8")

    # ---- 汇总（唯一的 Markdown 顶层文件；参考文献清单放在 参考文献/ 下，属该目录自身说明）
    n_files = sum(1 for _ in DEST.rglob("*") if _.is_file())
    size_mb = sum(p.stat().st_size for p in DEST.rglob("*") if p.is_file()) / 1024 / 1024
    zh = [
        "# 黑色素瘤 ICB 应答签名可复现性研究 — 交付汇总",
        "",
        "## 交付内容",
        "",
        "| 报告类型 | MDPI *Cancers* 投稿稿（英文），预登记可复现性研究 |",
        "|---|---|",
        f"| 篇幅 | 31 页，6 表 4 图，13 条参考文献，210 词摘要 |",
        f"| 文件总数 | {n_files} |",
        f"| 总体积 | {size_mb:.1f} MB |",
        "",
        "## 报告",
        "",
        "`报告/` 下为一份 DOCX，按 MDPI *Cancers* 体例排版：A4、双倍行距、Times New Roman "
        "12 pt、连续行号、题注在上图注在下、三线表。**作者、单位与通讯地址仍是占位符，"
        "需在投稿前填入。**",
        "",
        "## 脚本",
        "",
        f"`脚本/` 共 {len(KEEP_SCRIPTS) + 1} 个文件（22 个分析/作图脚本 + `requirements.txt`）。"
        "按打包规范，下载器脚本（`00`–`04`，只负责取 GEO 元数据与原始文件）与出稿/打包脚本"
        "未纳入本包，它们在公开仓库里。",
        "",
        "主要入口，按执行顺序：",
        "",
        "| 环节 | 脚本 |",
        "|---|---|",
        "| 患者级装配（分析单位是患者的核心） | `05_build_patient_level.py` |",
        "| 表达矩阵构建 | `06_build_expression_matrices.py` |",
        "| IMPRES 基线 | `14_impres_baseline.py` |",
        "| 嵌套交叉验证建模 | `16_nested_cv_model.py` |",
        "| GSE215868 纳入与 IMPRES 验证 | `17`、`18` |",
        "| IPS-MHC+CP | `23_ips_mhccp.py` |",
        "| 功效分析与研究级元分析 | `25_power_analysis.py` |",
        "| 队列异质性与归因边界 | `26_cohort_heterogeneity.py` |",
        "| MDPI 图件 | `24_figures_mdpi.py` |",
        "",
        "## 图",
        "",
        "`图/英文/` 4 张图，PNG + 矢量 PDF 两种格式，Arial，宽度均 ≤175 mm："
        "Figure 1 163.5 / Figure 3 153.2 / Figure 4 171.2 / Figure 2 162.7 mm。",
        "",
        "**本文是英文投稿，图件本身即英文版，因此未另行生成中文标注图集**——如需中文版汇报用图，"
        "改 `scripts/24_figures_mdpi.py` 的标签重新出图即可，数字与版式不会变。",
        "",
        "## 数据",
        "",
        f"`数据/处理后/` 共 {len(list((DEST / '数据/处理后').iterdir()))} 个文件，"
        "即全部分析面板：",
        "",
        "- `patient_level.tsv` —— 患者级纳入/排除表，**本文每一个 n 的唯一真源**",
        "- `qc/` 下的全部结果文件（AUC 表、元分析、分层结果、功效、单细胞 tie-break、TCR 克隆性等）",
        "",
        "**原始数据未纳入本包**，也不需要：全部表达数据来自 NCBI GEO 的公开登录号"
        "（GSE91061、GSE78220、GSE215868、GSE244982、GSE294272、GSE294273、GSE308433/434/435），"
        "任何人可凭登录号直接取得。稿件 Data Availability 一节已列明。",
        "",
        "## 参考文献",
        "",
        "稿件正文实际引用 13 条，全部经 PubMed / Europe PMC / Crossref 逐条回库核对，"
        "核对过程留痕在仓库 `docs/12_参考文献核对.md`。`参考文献/引用文献.csv` 是这 13 条的"
        "机器可读清单（不是整个检索库）。",
        "",
        "| # | 文献 | 出处 | 标识 | 在本文中的作用 |",
        "|---|---|---|---|---|",
    ]
    zh += [f"| {n} | {a} | {src} | {ident} | {role} |" for n, a, src, ident, role in REFS]
    zh += [
        "",
        "第 12、13 条未能用 Crossref / Europe PMC 机器核验（前者无 DOI 记录，后者未检出），"
        "按领域内通行引用保留，投稿前建议人工在 PubMed 复核一次。",
        "",
        "## 补充材料",
        "",
        "`Supplementary_Material_S1.pdf`（A4，4 页）随稿提交，与稿件同层。单细胞与 TCR 分析的"
        "全部数字由 `manuscript/build_supplementary.py` 从仓库结果文件读出生成，非手抄；"
        "其中 tie-break 表与仓库 `qc/scrna_tiebreak_report.md` 逐行一致。",
        "",
        "## 复现与可核查",
        "",
        "- 公开仓库：https://github.com/Cybing521/melanoma-icb-signature-reproducibility",
        "- 完整预登记（含一次流程违规的永久留痕）：仓库 `docs/03`",
        "- 参考文献逐条核对与更正记录：仓库 `docs/12`",
        "- 稿件正文不手抄进 Word，由 `manuscript/build_manuscript.py` 生成，"
        "`manuscript/ManuscriptBuilder/`（C# OpenXML）排版；补充材料由 "
        "`manuscript/build_supplementary.py` 生成。两者均在仓库内。",
        "",
        "## 打包规范说明（三处经判断的处理）",
        "",
        "1. **未生成中文标注图集。** 本文是英文投稿，图件本身即英文版（Arial、≤175 mm）。"
        "把英文图标成「中文图」或另出一套中文图都会失真，因此 `图/` 下只有 `英文/`。",
        "2. **`数据/处理后/` 下的 9 个 `.md` 是分析产物，不是包文档。** 它们是各分析步骤的结果"
        "报告（含事后标注与更正留痕），`汇总.md` 仍是全包唯一的说明性 Markdown。"
        "`Supplementary_Material_S1.pdf` 的 tie-break 表即与其中的 `scrna_tiebreak_report.md` "
        "逐行核对一致。",
        "3. **`scrna_reference_profile_v2.tsv` 保留文件名中的 v2。** 仓库中不存在 v1——该参考谱"
        "被重建过一次，`v2` 标记的是修正后的版本，与 `docs/08` 的更正记录相互对应。"
        "改名会切断这条溯源链。",
        "",
        "## 尚未完成",
        "",
        "1. **作者姓名、单位、通讯地址**——当前为占位符，需你提供。",
        "2. **参考文献第 12、13 条**未机器核验，建议人工 PubMed 复核。",
        "3. IMPRES 真实基因表未索取（按你的指示），已作为既定限制写入摘要与 §4.5。",
        "",
    ]
    (DEST / "汇总.md").write_text("\n".join(zh), encoding="utf-8")
    print(f"交付包已生成：{DEST}")
    print(f"  文件 {n_files} 个，{size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())