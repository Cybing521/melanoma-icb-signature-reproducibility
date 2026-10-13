#!/usr/bin/env python3
"""把参考文献从 13 条扩到 25 条，并按 MDPI 要求重排引用编号。

两件事一起做，顺序不能颠倒
----------------------------
1. **先插后编**。MDPI 要求参考文献按正文首次出现顺序编号。若先编好号再插新引用，
   新条目必然插不进正确位置。因此先用 `[[N1]]` 这类占位符把新引用写进正文，
   再做一次**同步**替换（同一趟里完成旧号→新号与占位符→新号），
   避免「1→2、2→3」这种级联错位。
2. **字段来自 PubMed**。12 条新条目的作者、刊名、年、卷、页由
   `refs/verify_refs.py` 回 efetch XML 取得，存在 `refs/pmid_records.json`；
   本脚本只负责排版，不重新解析题名。

用法：
    python refs/expand_refs.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONTENT = ROOT / "manuscript" / "content.json"
RECORDS = ROOT / "refs" / "pmid_records.json"

# ---------------------------------------------------------------- 旧 13 条
# 题名 / 刊卷页沿用稿件既有条目（已在 docs/12 逐条回库核对过），此处不重新解析。
OLD_REFS = {
    1: ("Auslander N.; Zhang G.; Lee J.S.; Frederick D.T.; Miao B.; Moll T.; Tian T.; Wei Z.; "
        "Madan S.; Sullivan R.J.; Boland G.; Flaherty K.; Herlyn M.; Ruppin E. "
        "Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. "
        "Nat. Med. 2018, 24, 1545–1549. "
        "(The IMPRES signature; note the publisher's correction of December 2018.)"),
    2: ("Charoentong P.; Finotello F.; Angelova M.; Mayer C.; Efremova M.; Rieder D.; Hackl H.; "
        "Trajanoski Z. Pan-cancer immunogenomic analyses reveal genotype–immunophenotype "
        "relationships and predictors of response to checkpoint blockade. "
        "Cell Rep. 2017, 18, 248–262. (The Immunophenoscore.)"),
    3: ("Riaz N.; Havel J.J.; Makarov V.; Desrichard A.; Urba W.J.; Sims J.S.; Hodi F.S.; "
        "Martín-Algarra S.; Mandal R.; Sharfman W.H.; Bhatia S.; Hwu W.J.; Gajewski T.F.; "
        "Slingluff C.L.; Chowell D.; Kendall S.M.; Chang H.; Shah R.; Kuo F.; Morris L.G.T.; "
        "Sidhom J.W.; Schneck J.P.; Horak C.E.; Weinhold N.; Chan T.A. "
        "Tumor and microenvironment evolution during immunotherapy with nivolumab. "
        "Cell 2017, 171, 934–949.e16. (Source of GSE91061.)"),
    4: ("Hugo W.; Zaretsky J.M.; Sun L.; Song C.; Moreno B.H.; Hu-Lieskovan S.; Berent-Maoz B.; "
        "Pang J.; Chmielowski B.; Cherry G.; Seja E.; Lomeli S.; Kong X.; Kelley M.C.; "
        "Sosman J.A.; Johnson D.B.; Ribas A.; Lo R.S. Genomic and transcriptomic features of "
        "response to anti-PD-1 therapy in metastatic melanoma. "
        "Cell 2016, 165, 35–44. (Source of GSE78220.)"),
    5: ("Vathiotis I.A.; Salichos L.; Martinez-Morilla S.; Gavrielatou N.; Aung T.N.; Shafi S.; "
        "Wong P.F.; Jessel S.; Kluger H.M.; Syrigos K.N.; Warren S.; Gerstein M.; Rimm D.L. "
        "Baseline gene expression profiling determines long-term benefit to programmed cell death "
        "protein 1 axis blockade. NPJ Precis. Oncol. 2022, 6, 92. "
        "(Source of GSE215868; native endpoint is 24-month long-term benefit.)"),
    6: ("Lauss M.; Phung B.; Borch T.H.; Harbst K.; Kaminska K.; Ebbesson A.; Hedenfalk I.; "
        "Yuan J.; Nielsen K.; Ingvar C.; Carneiro A.; Isaksson K.; Pietras K.; Svane I.M.; "
        "Donia M.; Jönsson G. Molecular patterns of resistance to immune checkpoint blockade in "
        "melanoma. Nat. Commun. 2024, 15, 3075. (Source of GSE244982.)"),
    7: ("Di Pietro A.; Au L.; Crock P.; Thio N.; Pizzolla A.; Nguyen T.N.; Macdonald S.; "
        "Chalmers H.; Zhu R.; Airaghi A.; Molden-Hauer T.; Bacac M.; Schwalie P.; Schlenker R.; "
        "Levesque M.P.; Mailer S.; Barnes-Cullen K.; Winch K.; Chan J.; Yeung G.A.; Spain L.; "
        "Rao A.D.; Sandhu S.; Gyorki D.E.; McArthur G.A.; Mackay L.K.; Neeson P.J. "
        "Tumor-resident T cells and dendritic cells form an in situ archetype during immunotherapy "
        "response in melanoma. Nat. Commun. 2026, 17, 7445. "
        "(Source of GSE294272 and GSE294273.)"),
    8: ("Tumor-intrinsic and extrinsic immune mechanisms define resistance to immune checkpoint "
        "blockade in metastatic melanoma. GEO Series GSE308433, GSE308434 and GSE308435 [Internet]. "
        "Columbia University, Izar laboratory; cited by accession because no peer-reviewed "
        "publication is linked to these accessions in GEO. Available at: "
        "https://www.ncbi.nlm.nih.gov/geo/ (accessed 7 October 2026). "
        "(Source of the single-cell and T-cell receptor analyses in Supplementary Material S1.)"),
    9: ("Liberzon A.; Birger C.; Thorvaldsdóttir H.; Ghandi M.; Mesirov J.P.; Tamayo P. "
        "The Molecular Signatures Database (MSigDB) hallmark gene set collection. "
        "Cell Syst. 2015, 1, 417–425. "
        "(Source of the feature annotations used in model development.)"),
    10: ("Hanley J.A.; McNeil B.J. The meaning and use of the area under a receiver operating "
         "characteristic (ROC) curve. Radiology 1982, 143, 29–36. "
         "(Standard error of the AUC, and the basis of the power calculation in Section 2.8.)"),
    11: ("DerSimonian R.; Laird N. Meta-analysis in clinical trials. "
         "Control. Clin. Trials. 1986, 7, 177–188. "
         "(Method for the study-level pooling in Section 2.8.)"),
    12: ("Deeks J.J.; Keating J.; Leeflang M.M.G. A hierarchy of diagnostic tests. "
         "Ann. Intern. Med. 2008, 148, 175–176."),
    13: ("Hanley J.A.; Lippman-Hand A. If a coin is tossed 3 times, is it fair? "
         "JAMA. 1983, 250, 2563–2566."),
}

# ---------------------------------------------------------------- 新 12 条
# (占位符, PMID, 用途, 正文注解)
NEW_REFS = [
    ("N1", "30127393", "Introduction — other published ICB signatures",
     "(Tumour Immune Dysfunction and Exclusion model.)"),
    ("N2", "28650338", "Introduction — interferon-γ-related gene-expression profile",
     "(The IFN-γ-related 28-gene profile.)"),
    ("N3", "31519915", "Introduction — B-cell-based signature from melanoma cohorts",
     "(An independent signature derived from melanoma directly.)"),
    ("N4", "30635339", "Discussion — the biological basis of immune exclusion",
     "(Mechanistic basis of the TIDE exclusion arm; cited in Section 4.3.)"),
    ("N5", "20010215", "Methods — framework for assessing discrimination and optimism",
     "(Framework for comparing AUC estimates across conventions.)"),
    ("N6", "18633328", "Introduction — inflation of published associations",
     "(Basis for treating a negative re-test as informative rather than anomalous.)"),
    ("N7", "16504092", "Methods — optimism of cross-validation used for model selection",
     "(Why the nested-CV conventions in Section 2.6 are all reported.)"),
    ("N8", "8870764", "Methods — sample size for accuracy studies",
     "(Sample-size framework underlying Table 4.)"),
    ("N9", "25569120", "Discussion — reporting standards for prediction models",
     "(TRIPOD; cited in Section 4.6.)"),
    ("N10", "16106022", "Discussion — reporting standards for tumour-marker studies",
     "(REMARK; cited in Section 4.6.)"),
    ("N11", "29531091", "Introduction — rationale for pre-registration",
     "(Rationale for fixing the protocol before computation.)"),
    ("N12", "25822800", "Supplementary — reference-profile deconvolution",
     "(CIBERSORT; the reference-profile method used in Supplementary Material S1.)"),
]

# -------------------------------------------------- 正文插入（占位符 → 片段）
INSERTS = [
    (4, " Both are used less as clinical tests than as comparators against which new models "
        "must demonstrate added value.",
        " Both are used less as clinical tests than as comparators against which new models "
        "must demonstrate added value. Several further signatures compete for the same role: "
        "the Tumour Immune Dysfunction and Exclusion model, which separates a tumour that is "
        "inaccessible to immune effectors from one whose effectors are dysfunctional [[N1]]; the "
        "interferon-\u03b3-related gene-expression profile [[N2]]; and a B-cell-based classifier "
        "derived from melanoma cohorts directly [[N3]]."),
    (5, " Retrospective replication can fail for reasons that have nothing to do with the "
        "signature, and a negative outcome becomes uninterpretable rather than informative.",
        " Retrospective replication can fail for reasons that have nothing to do with the "
        "signature, and a negative outcome becomes uninterpretable rather than informative. "
        "The concern is not hypothetical: associations between biomarkers and outcome are "
        "systematically inflated in the published literature [[N6]], which is precisely why "
        "re-testing under criteria fixed in advance is the only remedy that can produce an "
        "interpretable negative."),
    (6, " We then test two signatures selected to fail in distinguishable ways, and we report "
        "the boundary of what the available data can support.",
        " We then test two signatures selected to fail in distinguishable ways, and we report "
        "the boundary of what the available data can support. Fixing a protocol in advance is "
        "what converts a null result into evidence rather than an absence of evidence [[N11]]."),
    (28, " Model selection was nested: 100 outer folds, feature selection performed strictly "
         "within each training fold, a hard cap of 10 features, and 20 repeats of the entire "
         "procedure.",
         " Model selection was nested: 100 outer folds, feature selection performed strictly "
         "within each training fold, a hard cap of 10 features, and 20 repeats of the entire "
         "procedure. The nesting is what keeps the reported performance honest, since selecting "
         "features on the same data used to score them inflates apparent discrimination [[N7]]."),
    (33, " Standard errors of the AUC follow Hanley and McNeil [10].",
         " Standard errors of the AUC follow Hanley and McNeil [10], and the sample-size "
         "calculations follow the standard framework for accuracy studies [[N8]]."),
    (40, "Discrimination was quantified as the area under the receiver operating characteristic "
         "curve with 95% confidence intervals from 4000 bootstrap resamples of the patient-level "
         "score.",
         "Discrimination was quantified as the area under the receiver operating characteristic "
         "curve with 95% confidence intervals from 4000 bootstrap resamples of the patient-level "
         "score, following the standard framework for assessing prediction performance [[N5]]."),
    (42, " In the single-cell component, 88,715 cells passed initial loading",
         " Deconvolution used a reference-expression approach of the CIBERSORT type [[N12]]. "
         "In the single-cell component, 88,715 cells passed initial loading"),
    (96, " The natural reading is that the erratic signature is the one whose definition is "
         "least securely known, while the consistent signature is better specified and may be "
         "capturing a small amount of real information.",
         " The natural reading is that the erratic signature is the one whose definition is "
         "least securely known, while the consistent signature is better specified and may be "
         "capturing a small amount of real information. It is also consistent with what is known "
         "about what transcriptomic signatures in this field actually measure: the most widely "
         "used discriminators are dominated by a single axis of antigen presentation and "
         "interferon-\u03b3 response [[N2]], with immune exclusion driven by an upstream "
         "WNT/\u03b2-catenin programme rather than by the effector compartment itself [[N4]]. A "
         "signature built on that upstream programme would be expected to order patients "
         "inconsistently across cohorts whose exclusion burden differs, which is the failure mode "
         "reported here."),
    (105, " And they should be powered explicitly against the published comparator rather than "
          "against a conventional per-group sample size:",
          " They should also report their analyses to the standards already agreed for "
          "prediction models [[N9]] and tumour-marker studies [[N10]], so that a reader can tell "
          "a pre-specified comparison from one chosen after the fact. And they should be powered "
          "explicitly against the published comparator rather than against a conventional "
          "per-group sample size:"),
]


def build_new_ref_text(recs: dict, pmid: str, note: str) -> str:
    r = recs[pmid]
    authors = "; ".join(r["authors"]) + "."
    return f"{authors} {r['title']} {r['journal']} {r['year']}, {r['volume']}, {r['pages']}. {note}"


def main() -> int:
    recs = json.loads(RECORDS.read_text())
    doc = json.loads(CONTENT.read_text())
    blocks = doc["blocks"]

    # ---- 1. 插入新引用片段 -------------------------------------------------
    by_index = {i: b for i, b in enumerate(blocks)}
    for idx, old, new in INSERTS:
        b = by_index[idx]
        if old not in b["text"]:
            raise SystemExit(f"块 {idx} 未命中待替换片段:\n{old[:80]}")
        b["text"] = b["text"].replace(old, new, 1)
    print(f"已插入 {len(INSERTS)} 处新引用")

    # ---- 2. 按首次出现顺序定号 ---------------------------------------------
    # 扫描全部正文的占位符与旧方括号引用，得出首次出现次序
    order: list[str] = []
    for b in blocks:
        if b["type"] not in ("paragraph", "heading"):
            continue
        for m in re.finditer(r"\[\[N\d+\]\]|\[\d+\]", b.get("text", "")):
            tok = m.group(0)
            if tok not in order:
                order.append(tok)
    if len(order) != 25:
        raise SystemExit(f"期望 25 个引用标记，实得 {len(order)}: {order}")
    mapping = {tok: str(i + 1) for i, tok in enumerate(order)}
    print("编号映射:", ", ".join(f"{k}->{v}" for k, v in mapping.items()))

    # ---- 3. 同步替换（同一趟，避免级联）-----------------------------------
    pat = re.compile(r"\[\[N\d+\]\]|\[\d+\]")

    def repl(m: re.Match) -> str:
        return f"[{mapping[m.group(0)]}]"

    n = 0
    for b in blocks:
        if b["type"] in ("paragraph", "heading") and pat.search(b.get("text", "")):
            b["text"] = pat.sub(repl, b["text"])
            n += 1
    print(f"替换了 {n} 个块的引用标记")

    # ---- 4. 重写参考文献块 -------------------------------------------------
    head = next(i for i, b in enumerate(blocks)
                if b["type"] == "heading" and b["text"] == "References")
    tail = len(blocks)
    new_blocks: list[dict] = []
    note_by_pmid = {pmid: (tag, note) for tag, pmid, _u, note in NEW_REFS}
    for num, tok in enumerate(order, start=1):
        if tok.startswith("[["):
            pmid = next(p for t, p, _u, _n in NEW_REFS if f"[[{t}]]" == tok)
            note = next(n2 for t, p, _u, n2 in NEW_REFS if p == pmid)
            text = build_new_ref_text(recs, pmid, note)
        else:
            text = OLD_REFS[int(tok.strip("[]"))]
        new_blocks.append({"type": "paragraph", "text": f"{num}. {text}"})

    blocks[head + 1:tail] = new_blocks
    CONTENT.write_text(json.dumps(doc, ensure_ascii=False, indent=1))
    print(f"参考文献重写完毕：共 {len(new_blocks)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
