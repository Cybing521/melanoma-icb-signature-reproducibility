#!/usr/bin/env python3
"""按 MDPI 制图规范出图。

MDPI 图件规范（https://www.mdpi.com/journal/cancers/authors/manuscript-preparation）
---------------------------------------------------------------------------
* 宽度须适配 **85 mm（单栏）或 175 mm（通栏）**
* 位图**最低 300 dpi**（线条图建议 600 dpi 以免放大发虚）
* 字体 **Arial**，最终尺寸下 **8–10 pt**
* 面板标号用**小写加粗括号** (a) (b) (c)，置于各面板左上角
* 色彩用 RGB；不使用 3D 效果、不使用会误导的坐标轴截断

本脚本只读已落盘的结果文件，不做任何再计算。
输出到 `figures/`，同时写一份 `figures/FIGURE_LEGENDS.md`（MDPI 要求图注单独成文）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy import stats as _st


def power_at(auc: float, n1: int, n2: int) -> float:
    """Hanley–McNeil 口径下检出「真实 AUC 高于 0.50」的单侧功效。

    与 `scripts/25_power_analysis.py` 同一公式，图形与数字必须同源。
    """
    import math
    a = min(max(auc, 1e-6), 1 - 1e-6)
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2)
    se = math.sqrt(max(var, 1e-12))
    return float(_st.norm.cdf((a - 0.50) / se - 1.959963985))
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
QC = ROOT / "qc"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------- MDPI 样式
MM = 1 / 25.4
FULL_W = 175 * MM          # 通栏宽
COL_W = 85 * MM            # 单栏宽
DPI = 600                  # 线条图 600，保证放大清晰
FONT_SIZE = 8

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "axes.titlesize": FONT_SIZE + 0.5,
    "xtick.labelsize": FONT_SIZE - 0.5,
    "ytick.labelsize": FONT_SIZE - 0.5,
    "legend.fontsize": FONT_SIZE - 0.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "pdf.fonttype": 42,      # 字体嵌入，供排版编辑
    "ps.fonttype": 42,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# 黑白友好 + 色彩友好并存的调色板（色盲安全）
C_SIGN1 = "#0072B2"   # IMPRES 蓝
C_SIGN2 = "#D55E00"   # IPS-MHCCP 橙
C_GREY = "#4D4D4D"
C_ACC = "#009E73"
C_NULL = "#BDBDBD"


def panel(ax, label: str) -> None:
    """MDPI 面板标号：小写加粗括号，左上角。"""
    ax.text(-0.16, 1.06, f"({label})", transform=ax.transAxes,
            fontsize=FONT_SIZE + 1.5, fontweight="bold", va="bottom", ha="left")


MAX_W_MM = 175.5          # MDPI 通栏上限 175 mm，留 0.5 mm 容差


def save(fig, name: str, _tries: int = 0) -> None:
    """按 MDPI 上限自校正宽度后落盘。

    `bbox_inches="tight"` 会因刻度标签、面板标号、表格 bbox 把成图撑得比
    `figsize` 更宽，因此**必须实测**，不能只看 figsize。这里反复收窄画布
    直到成图宽度 ≤ 175.5 mm，最多 6 轮。
    """
    png = OUT / f"{name}.png"
    fig.savefig(png, dpi=DPI)
    from PIL import Image
    w_mm = Image.open(png).size[0] / DPI * 25.4
    if w_mm > MAX_W_MM and _tries < 12:
        factor = MAX_W_MM / w_mm * 0.97   # 固定边距不随画布线性缩放，需过冲
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * factor, h * factor, forward=True)
        return save(fig, name, _tries + 1)
    for ext in ("pdf", "tiff"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=DPI)
    plt.close(fig)
    flag = "OK " if w_mm <= MAX_W_MM else "WIDE"
    print(f"  [{flag}] figures/{name}.{{pdf,tiff,png}}  成图宽 {w_mm:.1f} mm"
          + (f"（收窄 {_tries} 轮）" if _tries else ""))


# ================================================================ Figure 1
def fig1_design() -> None:
    """研究设计与队列流程（CONSORT 式）。

    布局用 2×2 而非 3 栏：三栏并排在 175 mm 通栏宽下每个面板仅约 55 mm，
    容纳不下这些表格，实测会互相压盖。2×2 每面板约 85 mm，实测无重叠。

    表格一律用 `ax.table()`：matplotlib 的 table 不自动换行，所以所有
    cellText 都**必须**手动 `\n` 断行，否则文字会撑破单元格。
    """
    fig = plt.figure(figsize=(FULL_W, 156 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.26, wspace=0.22)

    def hide(ax):
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def mk_table(ax, rows, headers, widths, fontsize, bbox,
                 highlight_rows=(), edge="#C8C8C8"):
        """headers 传 None 表示无表头（不生成表头行）。"""
        t = ax.table(cellText=rows,
                     colLabels=[""] * len(widths) if headers is None else headers,
                     colWidths=widths, loc="upper left", bbox=bbox, cellLoc="left")
        t.auto_set_font_size(False)
        t.set_fontsize(fontsize)
        for (r, c), cell in t.get_celld().items():
            cell.set_edgecolor(edge)
            cell.set_linewidth(0.4)
            cell.get_text().set_linespacing(1.35)
            if headers is None and r == 0:
                cell.set_visible(False)
                continue
            if r == 0:
                cell.set_facecolor("#E8E8E8")
                cell.get_text().set_fontweight("bold")
            else:
                cell.set_facecolor("#FDF0E6" if r in highlight_rows else "#FFFFFF")
        return t

    # ---- (a) 队列总览 ----
    ax = fig.add_subplot(gs[0, 0])
    hide(ax)
    panel(ax, "a")
    ax.set_title("Cohorts and their roles", loc="left", fontweight="bold", pad=4)
    mk_table(
        ax,
        [["GSE91061", "109 / 65", "development"],
         ["GSE78220", "28 / 26", "external validation"],
         ["GSE215868", "105 / 105", "large-scale\nvalidation"],
         ["GSE244982", "41 / 41", "direction\nconsistency"],
         ["GSE308433/4/5", "42 / 34", "paired\ntumour-TCR"]],
        ["cohort", "smp / pat", "role in this study"],
        [0.34, 0.24, 0.42], 6.2, [0, 0.02, 1.0, 0.90], highlight_rows={3})

    # ---- (b) 主分析纳排 ----
    ax = fig.add_subplot(gs[0, 1])
    hide(ax)
    panel(ax, "b")
    ax.set_title("Analysis unit: patient", loc="left", fontweight="bold", pad=4)
    steps = [
        ("65 patients with Pre/On paired biopsies", "#F2F2F2", "#4D4D4D"),
        ("43 in primary analysis\n(PRCR 14 / PD 29 / SD 19 / NE 3)", "#F2F2F2", "#4D4D4D"),
        ("33 with a pre-treatment sample\n= modelling set (PRCR 10 / PD 23)", "#E6F0F7", "#0072B2"),
        ("10 excluded: Pt13/32/35/68/69/\n80/81/93/105/109", "#FFFFFF", "#BDBDBD"),
    ]
    y, h = 0.96, 0.195
    for i, (txt, fc, ec) in enumerate(steps):
        ax.add_patch(Rectangle((0.04, y - h), 0.92, h, transform=ax.transAxes,
                               facecolor=fc, edgecolor=ec,
                               lw=1.4 if i == 2 else 0.7))
        ax.text(0.50, y - h / 2, txt, fontsize=6.8, va="center", ha="center",
                transform=ax.transAxes, linespacing=1.4,
                fontweight="bold" if i == 2 else "normal")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((0.5, y - h - 0.006), (0.5, y - h - 0.055),
                                         transform=ax.transAxes, arrowstyle="-|>",
                                         mutation_scale=7, lw=0.8, color="#4D4D4D"))
        y -= h + 0.070

    # ---- (c) 分析单位纪律 ----
    ax = fig.add_subplot(gs[1, 0])
    hide(ax)
    panel(ax, "c")
    ax.set_title("Analysis-unit checks performed", loc="left", fontweight="bold", pad=4)
    mk_table(
        ax,
        [["Repeated biopsies", "paired / mixed structures;\nnot independent samples"],
         ["Pt27A / Pt27B", "one patient, two sites;\ncollapsed by median"],
         ["Pt16 (GSE78220)", "on-treatment only;\ndropped by flag"],
         ["Mel-TIL-026 / HM026", "same patient in two\ncohorts; de-duplicated"],
         ["GSE215868 plate ID", "plate position,\nnot a patient ID"]],
        ["situation", "handling"],
        [0.42, 0.58], 6.2, [0, 0.02, 1.0, 0.90])

    # ---- (d) 队列间为何不合并 ----
    ax = fig.add_subplot(gs[1, 1])
    hide(ax)
    panel(ax, "d")
    ax.set_title("What may and may not cross cohorts", loc="left",
                 fontweight="bold", pad=4)
    mk_table(
        ax,
        [["1. Sampling time", "GSE91061 pairs pre-dose with\nC1D29; the others are pre-tx"],
         ["2. Label definitions", "verbatim RECIST strings vs.\nbest overall response vs. BOR"],
         ["3. Platform, units", "raw counts / FPKM / featureCounts /\nNanoString counts"]],
        None, [0.38, 0.62], 6.6, [0.0, 0.30, 1.0, 0.62], edge="#FFFFFF")
    ax.text(0.0, 0.20,
            "Forbidden: pooling patients or expression matrices across cohorts.\n"
            "Permitted: study-level meta-analysis of published per-cohort AUCs,\n"
            "which touches no patient-level record and cannot create one.",
            fontsize=6.6, style="italic", color="#4D4D4D", va="top", ha="left")

    save(fig, "Figure1_StudyDesign")


# ================================================================ Figure 2
def fig2_headline() -> None:
    """头条图：两个签名 × 三个队列的 AUC 与 95% CI。"""
    im = pd.read_csv(QC / "IMPRES_scores_GSE215868.tsv", sep="\t")
    im_auc = pd.read_csv(QC / "GSE215868_impres_auc.tsv", sep="\t")
    ip = pd.read_csv(QC / "IPS_MHCCP_auc.tsv", sep="\t")

    fig = plt.figure(figsize=(FULL_W * 0.86, 88 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 1.02], wspace=0.52)

    # ---- (a) IMPRES 三队列森林图 ----
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a")
    ax.set_title("1 — IMPRES", loc="left", fontweight="bold", fontsize=FONT_SIZE)
    rows = [
        ("GSE91061  n=33", 0.359, 0.172, 0.552),
        ("GSE78220  n=26", 0.298, 0.122, 0.521),
        ("GSE215868  n=79", 0.505, 0.380, 0.631),
    ]
    for i, (lab, a, lo, hi) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=C_SIGN1, lw=1.6, solid_capstyle="round")
        ax.plot(a, i, "o", ms=5.2, color=C_SIGN1, mec="white", mew=0.7)
    ax.axvline(0.5, ls="--", lw=0.9, color=C_NULL, zorder=0)
    ax.axvspan(0.70, 1.0, color="#FDECEA", zorder=0)
    ax.text(0.85, -0.42, "published\nrange", fontsize=6, ha="center",
            va="top", color="#C0392B")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, 2.95)
    ax.set_xlabel("AUC (95% bootstrap CI)")
    ax.invert_yaxis()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ---- (b) IPS-MHCCP 森林图 ----
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b")
    ax.set_title("2 — IPS-MHC+CP", loc="left", fontweight="bold", fontsize=FONT_SIZE)
    sub = ip[ip["orientation"] == "high_ips_responder"]
    names = ["GSE91061  n=33", "GSE78220  n=26", "GSE215868  n=79"]
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=C_SIGN2, lw=1.6,
                solid_capstyle="round")
        ax.plot(r["auc"], i, "o", ms=5.2, color=C_SIGN2, mec="white", mew=0.7)
    ax.axvline(0.5, ls="--", lw=0.9, color=C_NULL, zorder=0)
    ax.axvspan(0.70, 1.0, color="#FDECEA", zorder=0)
    ax.set_yticks(range(3))
    ax.set_yticklabels(names, fontsize=6.2)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlabel("AUC (95% bootstrap CI)")
    ax.invert_yaxis()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ---- (c) 两个签名的方向一致性对比 ----
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c")
    ax.set_title("3 — Direction consistency", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    imp_vals = [0.359, 0.298, 0.505]
    ips_vals = [0.604, 0.583, 0.591]
    x = np.arange(3)
    ax.plot(x, imp_vals, "o-", color=C_SIGN1, lw=1.5, ms=5,
            label="IMPRES (erratic)")
    ax.plot(x, ips_vals, "o-", color=C_SIGN2, lw=1.5, ms=5,
            label="IPS-MHC+CP (consistent)")
    ax.axhline(0.5, ls="--", lw=0.9, color=C_NULL)
    ax.axhline(0.70, ls=":", lw=0.9, color="#C0392B")
    ax.text(0.99, 0.705, "reproduction bar", fontsize=6, color="#C0392B",
            ha="right", va="bottom", transform=ax.get_yaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels(["91061", "78220", "215868"], fontsize=6.4)
    ax.set_xlabel("cohort", fontsize=6.6, labelpad=1)
    ax.set_ylabel("AUC", labelpad=1)
    ax.set_ylim(0.2, 0.78)
    ax.legend(frameon=False, loc="lower left", fontsize=6.0,
              bbox_to_anchor=(0.02, 0.02), ncol=1, handletextpad=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    save(fig, "Figure2_SignaturePerformance")


# ================================================================ Figure 3
def fig3_exclusions() -> None:
    """七项技术性辩解逐条排除 + 两项不可排除。

    图注与图必须一一对应：初版只画了 3 项，图注却写「七项已排除 + 两项不可排除」，
    属图文不符。本版按 2×3 面板补齐全部条目，底部横条写明不可排除的两项。
    """
    inv = pd.read_csv(QC / "GSE215868_impres_invariant_check.tsv", sep="\t")
    diag = pd.read_csv(QC / "GSE215868_panel_diagnostic.tsv", sep="\t")
    var = pd.read_csv(QC / "GSE215868_score_variance.tsv", sep="\t")
    strat = pd.read_csv(QC / "stratified_auc.tsv", sep="\t")

    fig = plt.figure(figsize=(FULL_W, 146 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.70, wspace=0.46,
                          left=0.055, right=0.985, top=0.945, bottom=0.235)

    # ================= (a) 1–2 实现与归一化 =================
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a")
    n_ok = int(inv["identical_to_raw"].sum())
    n_all = len(inv)
    ax.barh([0], [n_all], color="#E8E8E8", height=0.34)
    ax.barh([0], [n_ok], color=C_ACC, height=0.34)
    ax.text(n_all + 0.4, 0, f"{n_ok}/{n_all}", va="center", fontsize=15,
            fontweight="bold", color=C_ACC)
    ax.set_xlim(0, n_all * 1.7)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xlabel("normalisation × direction combinations", fontsize=6.4)
    ax.set_title("1–2. Implementation\nand normalisation", loc="left",
                 fontweight="bold", fontsize=FONT_SIZE, linespacing=1.5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    # ================= (b) 3 平台适用性 =================
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b")
    order = ["GSE91061", "GSE294272", "GSE215868", "GSE78220", "GSE244982"]
    ok = diag[diag["status"] == "可算"]
    good = [float(((ok[ok["cohort"] == c]["p1"] >= 0.15) &
                   (ok[ok["cohort"] == c]["p1"] <= 0.85)).sum() /
                  max(len(ok[ok["cohort"] == c]), 1)) for c in order]
    cols = [C_GREY, C_GREY, C_SIGN2, C_GREY, C_GREY]
    ax.bar(range(5), good, color=cols, width=0.64)
    ax.axhline(0.5, ls="--", lw=0.9, color=C_NULL)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["91061", "294272", "215868", "78220", "244982"],
                       fontsize=6.0, rotation=45, ha="right")
    ax.set_ylabel("informative features", fontsize=6.4)
    ax.set_ylim(0, 1.28)
    ax.set_title("3. Platform suitability", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    ax.text(2, 1.06, "targeted panel\nis not more degenerate",
            ha="center", fontsize=6.0, color=C_SIGN2, fontweight="bold",
            linespacing=1.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ================= (c) 3 支撑：分数未退化 =================
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c")
    ax.bar(range(len(var)), var["ratio"], color=C_ACC, width=0.58)
    ax.axhline(1.0, ls="--", lw=0.9, color=C_GREY)
    ax.set_xticks(range(len(var)))
    ax.set_xticklabels([v.replace("GSE", "") for v in var["cohort"]],
                       fontsize=6.0, rotation=45, ha="right")
    ax.set_ylabel("observed / expected SD", fontsize=6.4)
    ax.set_ylim(0, 2.35)
    ax.set_title("3. Score not\ndegenerate", loc="left",
                 fontweight="bold", fontsize=FONT_SIZE, linespacing=1.5)
    ax.text(0.5, -0.34, "1.0 = variance under independence",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
            color=C_GREY)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ================= (d) 5 活检时点 =================
    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "d")
    vis = strat[(strat["stratum_type"] == "visit") & (strat["stratum"] != "Pre+On")]
    vis = vis.set_index("stratum").loc[["Pre", "On"]]
    ys = [1, 0]
    for y, (idx, r) in zip(ys, vis.iterrows()):
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=C_SIGN1, lw=1.8,
                solid_capstyle="round")
        ax.plot(r["auc"], y, "o", ms=5.4, color=C_SIGN1, mec="white", mew=0.7)
    ax.axvline(0.5, ls="--", lw=0.9, color=C_NULL, zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{i}  n={int(vis.loc[i, 'n_PRCR'])}/{int(vis.loc[i, 'n_PD'])}"
                        for i in vis.index], fontsize=6.0)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("AUC (95% bootstrap CI)", fontsize=6.4)
    ax.set_title("5. Biopsy\ntimepoint", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE, linespacing=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ================= (e) 6 治疗方案 =================
    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "e")
    rx = strat[(strat["stratum_type"] == "regimen") &
               (strat["signature"] == "IMPRES g1_low") &
               (strat["n_PD"] > 0)]
    rx = rx.set_index("stratum").loc[["IPI+NIVO", "NIVO", "PEMBRO"]]
    ys = [2, 1, 0]
    for y, (idx, r) in zip(ys, rx.iterrows()):
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=C_SIGN1, lw=1.8,
                solid_capstyle="round")
        ax.plot(r["auc"], y, "o", ms=5.4, color=C_SIGN1, mec="white", mew=0.7)
    ax.axvline(0.5, ls="--", lw=0.9, color=C_NULL, zorder=0)
    ax.axvspan(0.70, 1.0, color="#FDECEA", zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{i}  n={int(rx.loc[i, 'n_PRCR'])}/{int(rx.loc[i, 'n_PD'])}"
                        for i in rx.index], fontsize=6.0)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlabel("AUC (95% bootstrap CI)", fontsize=6.4)
    ax.set_title("6. Treatment\nregimen", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE, linespacing=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ================= (f) 4 样本量（功效曲线） =================
    ax = fig.add_subplot(gs[1, 2])
    panel(ax, "f")
    grid = np.linspace(0.50, 0.85, 220)
    pw = np.array([power_at(g, 45, 34) for g in grid])
    ax.plot(grid, pw, color=C_SIGN2, lw=1.8)
    ax.axhline(0.8, ls=":", lw=0.9, color=C_GREY)
    ax.axvline(0.70, ls="--", lw=0.9, color=C_NULL)
    ax.axvline(0.77, ls="--", lw=0.9, color="#C0392B")
    ax.fill_between(grid, pw, 1.0, where=(pw >= 0.8), color=C_SIGN2, alpha=0.10)
    for xv, lab, col in ((0.70, "0.70", C_NULL), (0.77, "0.77", "#C0392B")):
        ax.text(xv + 0.004, 0.06, lab, fontsize=6.0, color=col, rotation=90,
                va="bottom")
    ax.text(0.50, 0.84, "0.8 power", fontsize=6.0, color=C_GREY, va="bottom")
    ax.set_xlabel("true AUC", fontsize=6.4)
    ax.set_ylabel("power (n = 45 / 34)", fontsize=6.4)
    ax.set_ylim(0, 1.04)
    ax.set_xlim(0.50, 0.85)
    ax.set_title("4. Sample size", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ================= 底部横条 =================
    band = (
        "Excluded outright (6 of 9); item 4 addressed but not fully resolved:\n"
        "1  implementation error      4  sample size         7  endpoint definition\n"
        "2  normalisation             5  biopsy timepoint     population composition   ✗\n"
        "3  platform unsuitability    6  treatment regimen   weakness in the\n"
        "                                                derivation cohorts   ✗\n"
        "✗ = the available data cannot exclude this. Item 4 is excluded against the published\n"
        "performance (power 0.93 at AUC 0.70) but not against a weak signal. Item 7 needed no\n"
        "data panel because the native endpoint of GSE215868 is 24-month long-term benefit\n"
        "(PFS-derived), not RECIST; response was re-derived from best overall response to match\n"
        "the endpoint on which IMPRES was reported."
    )
    fig.text(0.5, 0.118, band, ha="center", va="top", fontsize=6.0,
             color="#333333", linespacing=1.75,
             bbox=dict(boxstyle="round,pad=0.55", facecolor="#F6F6F6",
                       edgecolor="#CCCCCC"))

    save(fig, "Figure3_RuledOutExplanations")


# ================================================================ Figure 4
def fig4_cv() -> None:
    """嵌套交叉验证：两套特征集的折间分布。"""
    fig = plt.figure(figsize=(FULL_W, 72 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.0, 1.1], wspace=0.42)

    def load(tag):
        d = pd.read_csv(QC / f"{tag}_results.tsv", sep="\t")
        d["auc"] = pd.to_numeric(d["auc"], errors="coerce")
        return d.dropna(subset=["auc"])

    a = load("dev_cv")
    b = load("dev_cv_mainpanel")

    # ---- (a) 折 AUC 分布 ----
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a")
    ax.set_title("Per-fold AUC (100 folds)", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    bins = np.linspace(0, 1, 26)
    ax.hist([a["auc"], b["auc"]], bins=bins, color=[C_GREY, C_SIGN1],
            label=["17 features (as run)", "6 features (pre-registered)"],
            alpha=0.85, edgecolor="white", linewidth=0.3, rwidth=0.92)
    ax.axvline(0.5, color="#C0392B", lw=1.1, ls="--")
    ax.set_xlabel("AUC in one outer fold  (dashed line = chance)")
    ax.set_ylabel("number of folds")
    ax.legend(frameon=False, fontsize=6.0, loc="upper right",
              bbox_to_anchor=(1.0, 1.0))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ---- (b) 重复间均值 ----
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b")
    ax.set_title("Mean AUC per repeat", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    ra = a.groupby("repeat")["auc"].mean()
    rb = b.groupby("repeat")["auc"].mean()
    ax.scatter(range(1, len(ra) + 1), ra, s=14, color=C_GREY, label="17 features")
    ax.scatter(range(1, len(rb) + 1), rb, s=14, color=C_SIGN1, label="6 features")
    ax.axhline(0.5, color="#C0392B", lw=1.0, ls="--")
    ax.axhline(ra.mean(), color=C_GREY, lw=0.8, ls=":")
    ax.axhline(rb.mean(), color=C_SIGN1, lw=0.8, ls=":")
    ax.set_xlabel("CV repeat")
    ax.set_ylabel("mean AUC")
    ax.set_ylim(0.1, 0.95)
    ax.legend(frameon=False, fontsize=6.0, loc="lower left", ncol=2,
              columnspacing=0.9, handletextpad=0.35)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ---- (c) 三个口径 ----
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c")
    ax.set_title("Three conventions", loc="left", fontweight="bold",
                 fontsize=FONT_SIZE)
    labels = ["pooled\nOOF", "repeat\nmean", "fold\nmean"]
    v17 = [0.500, 0.516, 0.549]
    v6 = [0.417, 0.408, 0.347]
    x = np.arange(3)
    ax.bar(x - 0.19, v17, 0.36, color=C_GREY, label="17 features")
    ax.bar(x + 0.19, v6, 0.36, color=C_SIGN1, label="6 features")
    ax.axhline(0.5, color="#C0392B", lw=1.0, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.9)
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 0.75)
    ax.legend(frameon=False, fontsize=6.2, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # 说明文字放到坐标轴外，避免压在柱子上
    ax.text(0.5, -0.30,
            "pooled OOF merges all 20 repeats, so each patient is\n"
            "counted 20 times. Fold SD is 0.257 for the 17-feature set\n"
            "and 0.213 for the pre-registered panel.",
            transform=ax.transAxes, fontsize=6.0, style="italic", color=C_GREY,
            ha="center", va="top", linespacing=1.5)

    save(fig, "Figure4_NestedCV")


def main() -> int:
    print("按 MDPI 规范出图：")
    fig1_design()
    fig2_headline()
    fig3_exclusions()
    fig4_cv()
    print(f"\n全部输出至 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())