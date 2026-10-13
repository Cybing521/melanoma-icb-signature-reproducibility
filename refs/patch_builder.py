#!/usr/bin/env python3
"""把引用扩写与重编号固化进 `manuscript/build_manuscript.py`。

为什么不能只改 content.json
---------------------------
`build_manuscript.py` 是正文的**唯一**来源：它每次运行都会重新生成
`manuscript/content.json`。直接改 JSON 会在下一次 build 时被静默覆盖——
这正是本项目「数字必须可追回某次脚本运行」这条纪律要防的事。因此把新增
引用措辞、编号规则和 25 条参考文献固化进 builder。

做法
----
1. 新增 9 处引用措辞作为 `CITATION_INSERTS`（按正文出现顺序）；
2. 新增 `apply_citation_edits()`：插入后扫描首次出现顺序，
   **同步**替换旧号与占位符，避免「1→2、2→3」级联错位；
3. 把参考文献列表从 13 条换成 25 条，并按 MDPI 的首次出现顺序编号。

校验：运行本脚本前后分别 build 一次，比对 content.json 必须逐字节一致。

用法：
    python refs/patch_builder.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILDER = ROOT / "manuscript" / "build_manuscript.py"
CONTENT = ROOT / "manuscript" / "content.json"

INSERTS = [
    ("Both are used less as clinical tests than as comparators against which new models "
     "must demonstrate added value.",
     "Both are used less as clinical tests than as comparators against which new models "
     "must demonstrate added value. Several further signatures compete for the same role: "
     "the Tumour Immune Dysfunction and Exclusion model, which separates a tumour that is "
     "inaccessible to immune effectors from one whose effectors are dysfunctional [[N1]]; the "
     "interferon-\u03b3-related gene-expression profile [[N2]]; and a B-cell-based classifier "
     "derived from melanoma cohorts directly [[N3]]."),
    ("Retrospective replication can fail for reasons that have nothing to do with the "
     "signature, and a negative outcome becomes uninterpretable rather than informative.",
     "Retrospective replication can fail for reasons that have nothing to do with the "
     "signature, and a negative outcome becomes uninterpretable rather than informative. "
     "The concern is not hypothetical: associations between biomarkers and outcome are "
     "systematically inflated in the published literature [[N6]], which is precisely why "
     "re-testing under criteria fixed in advance is the only remedy that can produce an "
     "interpretable negative."),
    ("We then test two signatures selected to fail in distinguishable ways, and we report "
     "the boundary of what the available data can support.",
     "We then test two signatures selected to fail in distinguishable ways, and we report "
     "the boundary of what the available data can support. Fixing a protocol in advance is "
     "what converts a null result into evidence rather than an absence of evidence [[N11]]."),
    ("Model selection was nested: 100 outer folds, feature selection performed strictly "
     "within each training fold, a hard cap of 10 features, and 20 repeats of the entire "
     "procedure.",
     "Model selection was nested: 100 outer folds, feature selection performed strictly "
     "within each training fold, a hard cap of 10 features, and 20 repeats of the entire "
     "procedure. The nesting is what keeps the reported performance honest, since selecting "
     "features on the same data used to score them inflates apparent discrimination [[N7]]."),
    ("Standard errors of the AUC follow Hanley and McNeil [10].",
     "Standard errors of the AUC follow Hanley and McNeil [16], and the sample-size "
     "calculations follow the standard framework for accuracy studies [[N8]]."),
    ("Discrimination was quantified as the area under the receiver operating characteristic "
     "curve with 95% confidence intervals from 4000 bootstrap resamples of the patient-level "
     "score.",
     "Discrimination was quantified as the area under the receiver operating characteristic "
     "curve with 95% confidence intervals from 4000 bootstrap resamples of the patient-level "
     "score, following the standard framework for assessing prediction performance [[N5]]."),
    ("In the single-cell component, 88,715 cells passed initial loading",
     "Deconvolution used a reference-expression approach of the CIBERSORT type [[N12]]. "
     "In the single-cell component, 88,715 cells passed initial loading"),
    ("The natural reading is that the erratic signature is the one whose definition is "
     "least securely known, while the consistent signature is better specified and may be "
     "capturing a small amount of real information.",
     "The natural reading is that the erratic signature is the one whose definition is "
     "least securely known, while the consistent signature is better specified and may be "
     "capturing a small amount of real information. It is also consistent with what is known "
     "about what transcriptomic signatures in this field actually measure: the most widely "
     "used discriminators are dominated by a single axis of antigen presentation and "
     "interferon-\u03b3 response [4], with immune exclusion driven by an upstream "
     "WNT/\u03b2-catenin programme rather than by the effector compartment itself [[N4]]. A "
     "signature built on that upstream programme would be expected to order patients "
     "inconsistently across cohorts whose exclusion burden differs, which is the failure mode "
     "reported here."),
    ("And they should be powered explicitly against the published comparator rather than "
     "against a conventional per-group sample size:",
     "They should also report their analyses to the standards already agreed for "
     "prediction models [[N9]] and tumour-marker studies [[N10]], so that a reader can tell "
     "a pre-specified comparison from one chosen after the fact. And they should be powered "
     "explicitly against the published comparator rather than against a conventional "
     "per-group sample size:"),
]

HELPER = '''

# ─────────────────────── 引用插入与按出现顺序重编号 ───────────────────────
# 12 条新增参考文献的措辞改动集中在这里，便于逐条复核。
# 编号由脚本按**正文首次出现顺序**自动分配（MDPI 要求），因此新增引用只需
# 插入措辞，不必手工维护编号，也不会出现「1→2、2→3」的级联错位。
CITATION_INSERTS = [
'''

FIGPATH = {
    '"figures/Figure4_NestedCV.png"': '"figures/Figure2_NestedCV.png"',
    '"figures/Figure2_SignaturePerformance.png"': '"figures/Figure3_SignaturePerformance.png"',
    '"figures/Figure3_RuledOutExplanations.png"': '"figures/Figure4_RuledOutExplanations.png"',
}


def wrap(s: str, indent: str = "    ", width: int = 86) -> str:
    """折成 builder 里既有的多行隐式拼接写法。

    每行都自带一对引号（Python 的隐式字符串拼接），否则生成的源码不可解析：
    `        "1. Auslander N.; ... E. "` / `        "Robust prediction ... 1545–1549. (...)"`。
    """
    out, line = [], ""
    for tok in s.split():
        cand = f"{line} {tok}" if line else tok
        if len(cand) + 2 > width:
            out.append(f'{indent}"{line} "')
            line = tok
        else:
            line = cand
    if line:
        out.append(f'{indent}"{line}"')
    return "\n".join(out)


def main() -> int:
    src = BUILDER.read_text()
    doc = json.loads(CONTENT.read_text())
    refs_head = next(i for i, b in enumerate(doc["blocks"])
                     if b["type"] == "heading" and b["text"] == "References")
    refs = [b["text"] for b in doc["blocks"][refs_head + 1:]]
    assert len(refs) == 25, len(refs)

    # ---- 1. 插图文件路径（文件名已按正文顺序重排）-------------------------
    for old, new in FIGPATH.items():
        if old not in src:
            raise SystemExit(f"builder 中未找到图路径 {old}")
        src = src.replace(old, new)

    # ---- 2. 插入引用措辞表 + 重编号函数 ---------------------------------
    body = ["    (%r,\n     %r),\n" % (a, b) for a, b in INSERTS]
    helper = HELPER + "".join(body) + ''']


def apply_citation_edits(blocks: list[dict]) -> None:
    """就地插入新引用措辞，并按首次出现顺序重排全部编号。

    为什么重编号要自动做：MDPI 要求参考文献按正文首次出现顺序编号。
    手工维护编号时，新增一条引用会把后面全部条目挤动，极易出错；
    这里让脚本编号，规则只有一条——谁先出现谁在前。
    """
    for b in blocks:
        if b["type"] != "paragraph":
            continue
        for old, new in CITATION_INSERTS:
            if old in b["text"]:
                b["text"] = b["text"].replace(old, new, 1)
                break

    order: list[str] = []
    for b in blocks:
        if b["type"] not in ("paragraph", "heading"):
            continue
        for m in re.finditer(r"\\[\\[N\\d+\\]\\]|\\[\\d+\\]", b["text"]):
            if m.group(0) not in order:
                order.append(m.group(0))
    mapping = {tok: str(i + 1) for i, tok in enumerate(order)}

    pat = re.compile(r"\\[\\[N\\d+\\]\\]|\\[\\d+\\]")
    for b in blocks:
        if b["type"] in ("paragraph", "heading"):
            b["text"] = pat.sub(lambda m: f"[{mapping[m.group(0)]}]", b["text"])

'''
    src = src.replace("\ndef build() -> dict:", helper + "\ndef build() -> dict:", 1)
    if "apply_citation_edits" not in src:
        raise SystemExit("未能插入 apply_citation_edits")

    # ---- 3. 参考文献列表换成 25 条 ---------------------------------------
    start = src.index('    B.append(h(1, "References"))')
    end = src.index("        B.append(p(ref))", start)
    new_refs = ['    B.append(h(1, "References"))\n',
                "    for ref in [\n",
                "        # 25 条。原有 13 条的作者、刊名、卷页已于 2026-10-07 逐条回库核对，\n",
                "        # 记录见 docs/12_参考文献核对.md；新增 12 条由 refs/verify_refs.py 回 PubMed\n",
                "        # efetch XML 逐条取字段，并以 Europe PMC 交叉核对 DOI（见 docs/13）。\n",
                "        # 编号由 apply_citation_edits() 按正文首次出现顺序自动分配，此处顺序即定稿顺序。\n"]
    for i, r in enumerate(refs, start=1):
        new_refs.append(wrap(r, indent="        ") + ",\n")
    new_refs.append("    ]:\n")
    src = src[:start] + "".join(new_refs) + src[end:]

    # ---- 4. 在 build() 写盘前调用重编号 --------------------------------
    if "import re" not in src.split("def build")[0]:
        src = src.replace("import json\nimport subprocess", "import json\nimport re\nimport subprocess", 1)
    if "apply_citation_edits(B)" not in src:
        src = src.replace("    return {\n        \"meta\": {",
                          "    apply_citation_edits(B)\n\n    return {\n        \"meta\": {", 1)

    BUILDER.write_text(src)
    print("builder 已更新：25 条参考文献 + 9 处引用措辞 + 图路径重排")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
