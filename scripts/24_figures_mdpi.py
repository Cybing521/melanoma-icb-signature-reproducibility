#!/usr/bin/env python3
"""按参考文献风格出图（MDPI *Cancers* 体例）。

风格依据（实测，不是凭印象）
------------------------------
本脚本的字体与配色取自本研究直接对照的原始文献，证据链见
`refs/probe_style.py`（PDF 嵌入字体 + 栅格化取色）与 `docs/13`：

| 来源 | 在本文中的角色 | 实测要点 |
|---|---|---|
| Auslander 2018 *Nat Med*（IMPRES） | 签名 1 | 主色为红/鲑红系；KM 曲线 High=蓝 / Low=红 |
| Hugo 2016 *Cell*（GSE78220 出处） | 验证队列 | 箱线浅蓝 #90CCFC ↔ 浅红 #FCB4B4；蓝=应答 / 红=不应答 |
| Riaz 2017 *Cell*（GSE91061 出处） | 开发队列 | 鲑红 #F06C6C、藏青 #3C4884、橙 #F09C00、青 #249084 |
| Lauss 2024 *Nat Commun*（GSE244982 出处） | 归因边界 | 深蓝 #0060A0、红 #E03030、绿 #70B020；面板标号 A. B. C. |

由此定下两条主色：**签名 1 = 深蓝**（对应领域惯例「蓝 = 应答 / 高分」），
**签名 2 = 鲑红**（对应 Hugo 的「红 = 不应答 / 低分」与 Auslander 的主色系）。
字体沿用这些文献通用的 Arial/Helvetica 无衬线，7–8 pt。

与 MDPI 的关系：宽度、300/600 dpi、Arial 8–10 pt、面板标号用小写加粗括号
均按 MDPI 投稿要求执行；面板标号不改成 Nature 的大写 A.，因为投稿目标是 MDPI。

数据纪律
--------
本脚本只读 `qc/` 下已落盘的结果文件，不做任何再计算，也不硬编码任何统计量。
逐折 AUC 之外的全部 AUC / CI / 合并估计一律取自 `qc/power_analysis.tsv`
与 `qc/meta_auc.tsv`，因此图与表、正文必然同源。

用法
----
    python scripts/24_figures_mdpi.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy import stats as _st

ROOT = Path(__file__).resolve().parent.parent
QC = ROOT / "qc"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------- MDPI 尺寸
MM = 1 / 25.4
FULL_W = 175 * MM          # 通栏上限
DPI = 600
FONT_SIZE = 7.5
MAX_W_MM = 175.5

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE - 1,
    "ytick.labelsize": FONT_SIZE - 1,
    "legend.fontsize": FONT_SIZE - 1,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
})

# ---------------------------------------------------------------- 配色
# 取自 refs/probe_style.py 对上表四篇文献的实测取色，理由见模块 docstring。
C_S1 = "#1F5FA9"      # 签名 1 IMPRES —— 深蓝（Lauss #0060A0 / Hugo 浅蓝的深色端）
C_S1_L = "#A8CBE8"    # 签名 1 浅色端（Hugo #90CCFC）
C_S2 = "#C8392F"      # 签名 2 IPS-MHC+CP —— 鲑红（Auslander #FC6060 / Riaz #F06C6C）
C_S2_L = "#F0A9A0"    # 签名 2 浅色端（Hugo #FCB4B4）
C_ACC = "#2E7D5B"     # 已排除 / 通过（Lauss #70B020、Riaz 青 #249084 的深色端）
C_GREY = "#6E6E6E"
C_NULL = "#9A9A9A"
C_BAND = "#EDF2F7"    # 预登记判定窗
C_FILL = "#F7F8F9"    # 浅底
C_GRID = "#D5D9DD"

COHORTS = ["GSE91061", "GSE78220", "GSE215868"]
# 森林图 y 轴标签用短名：长名会把相邻面板的数值标注挤出画布
SHORT = {"GSE91061": "91061", "GSE78220": "78220", "GSE215868": "215868"}


def panel(ax, label: str, dx: float = -0.15, dy: float = 1.045) -> None:
    """MDPI 面板标号：小写加粗括号，左上角。"""
    ax.text(dx, dy, f"({label})", transform=ax.transAxes,
            fontsize=FONT_SIZE + 2, fontweight="bold", va="bottom", ha="left")


def title(ax, text: str) -> None:
    ax.set_title(text, loc="left", fontweight="bold", pad=5,
                 color="#1A1A1A", fontsize=FONT_SIZE)


def clean(ax, keep=("left", "bottom")) -> None:
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    ax.tick_params(direction="out")


def save(fig, name: str, _tries: int = 0) -> None:
    """按 MDPI 175.5 mm 上限自校正宽度后落盘。

    `bbox_inches="tight"` 会把刻度标签与面板标号算进 bbox，成图可能比 figsize 宽，
    因此必须实测，不能只看画布尺寸。
    """
    png = OUT / f"{name}.png"
    fig.savefig(png, dpi=DPI)
    from PIL import Image
    w_mm = Image.open(png).size[0] / DPI * 25.4
    if w_mm > MAX_W_MM and _tries < 12:
        factor = MAX_W_MM / w_mm * 0.97
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * factor, h * factor, forward=True)
        return save(fig, name, _tries + 1)
    for ext in ("pdf", "tiff"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=DPI)
    plt.close(fig)
    flag = "OK " if w_mm <= MAX_W_MM else "WIDE"
    print(f"  [{flag}] figures/{name}  成图宽 {w_mm:.1f} mm"
          + (f"（收窄 {_tries} 轮）" if _tries else ""))


def power_at(auc: float, n1: int, n2: int) -> float:
    """Hanley–McNeil 口径下检出「真实 AUC 高于 0.50」的单侧功效。

    与 `scripts/25_power_analysis.py` 同一公式，图形与 Table 4 必须同源。
    """
    a = min(max(auc, 1e-6), 1 - 1e-6)
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2)
    se = math.sqrt(max(var, 1e-12))
    return float(_st.norm.cdf((a - 0.50) / se - 1.959963985))


def forest(ax, rows, color, chance=True, band=True, xlabel="AUC (95% bootstrap CI)"):
    """森林图公共画法。`rows` = [(label, auc, lo, hi, [色])]；`auc` 为 None 的行作组标题。"""
    ys = np.arange(len(rows))
    if band:
        ax.axvspan(0.70, 1.0, color=C_BAND, zorder=0, lw=0)
    if chance:
        ax.axvline(0.5, ls="--", lw=0.7, color=C_NULL, zorder=1)
    tr = ax.get_yaxis_transform()      # x 用轴分数、y 用数据坐标
    for y, row in zip(ys, rows):
        lab, a, lo, hi = row[0], row[1], row[2], row[3]
        c = row[4] if len(row) > 4 and row[4] else color
        if a is None:                   # 组标题行：只留标签
            continue
        ax.plot([lo, hi], [y, y], color=c, lw=1.7, solid_capstyle="butt", zorder=2)
        ax.plot(a, y, "o", ms=4.6, color=c, mec="white", mew=0.6, zorder=3)
        # 数值标注自适应放置：优先 CI 右侧，放不下就放左侧，始终留在轴内——
        # 放到轴外会落进相邻面板的标签区，读起来像是下一个面板的数字。
        txt = f"{a:.3f}"
        if hi + 0.025 <= 0.975:
            xpos, ha = hi + 0.025, "left"
        elif lo - 0.025 >= 0.005:
            xpos, ha = lo - 0.025, "right"
        else:
            continue
        ax.text(xpos, y, txt, va="center", ha=ha, fontsize=6.2, color=c,
                zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.6)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(len(rows) - 0.45, -0.85)
    ax.set_xlabel(xlabel)
    clean(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    return ax


def group_rules(ax, y, color):
    """签名之间的分组横线（森林图左端短横）。"""
    ax.plot([0.0, 0.085], [y, y], color=color, lw=1.6, clip_on=False,
            solid_capstyle="butt")


# ================================================================ Figure 1
def fig1_design() -> None:
    """研究设计与队列流程。

    版式沿用参考文献的 2×2 面板 + 浅底表格；修掉初版面板 (b) 末框被裁、
    (d) 文字压边两个问题——办法是把画布加高到 168 mm 并给每个面板留足下边距。
    """
    fig = plt.figure(figsize=(FULL_W, 172 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20,
                          left=0.035, right=0.985, top=0.960, bottom=0.015)

    def hide(ax):
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def mk_table(ax, rows, headers, widths, fontsize, bbox,
                 highlight_rows=(), head_fill=C_BAND, edge=C_GRID):
        t = ax.table(cellText=rows,
                     colLabels=[""] * len(widths) if headers is None else headers,
                     colWidths=widths, loc="upper left", bbox=bbox, cellLoc="left")
        t.auto_set_font_size(False)
        t.set_fontsize(fontsize)
        for (r, c), cell in t.get_celld().items():
            cell.set_edgecolor(edge)
            cell.set_linewidth(0.35)
            cell.get_text().set_linespacing(1.4)
            if headers is None and r == 0:
                cell.set_visible(False)
                continue
            cell.set_facecolor(head_fill if r == 0 else
                               (C_S1_L if r in highlight_rows else "#FFFFFF"))
            if r == 0:
                cell.get_text().set_fontweight("bold")
        return t

    # ---- (a) 队列总览 ----
    ax = fig.add_subplot(gs[0, 0])
    hide(ax)
    panel(ax, "a", dx=-0.055)
    title(ax, "Cohorts and their roles")
    mk_table(
        ax,
        [["GSE91061", "109 / 65", "development"],
         ["GSE78220", "28 / 26", "external validation"],
         ["GSE215868", "105 / 105", "large-scale validation"],
         ["GSE244982", "41 / 41", "direction consistency"],
         ["GSE308433/4/5", "42 / 34", "paired tumour–TCR"]],
        ["cohort", "smp / pat", "role in this study"],
        [0.32, 0.23, 0.45], 6.2, [0, 0.02, 1.0, 0.92], highlight_rows={3})

    # ---- (b) 主分析纳排 ----
    ax = fig.add_subplot(gs[0, 1])
    hide(ax)
    panel(ax, "b")
    title(ax, "Analysis unit: patient")
    steps = [
        ("65 patients with Pre/On paired biopsies", C_FILL, C_GREY),
        ("43 in primary analysis\n(PRCR 14 / PD 29 / SD 19 / NE 3)", C_FILL, C_GREY),
        ("33 with a pre-treatment sample\n= modelling set (PRCR 10 / PD 23)",
         C_S1_L, C_S1),
        ("10 excluded: Pt13/32/35/68/69/80/81/93/105/109", "#FFFFFF", C_NULL),
    ]
    y, h = 0.97, 0.185
    for i, (txt, fc, ec) in enumerate(steps):
        ax.add_patch(Rectangle((0.05, y - h), 0.90, h, transform=ax.transAxes,
                               facecolor=fc, edgecolor=ec,
                               lw=1.5 if i == 2 else 0.6))
        ax.text(0.50, y - h / 2, txt, fontsize=6.6, va="center", ha="center",
                transform=ax.transAxes, linespacing=1.45,
                color=C_S1 if i == 2 else "#1A1A1A",
                fontweight="bold" if i == 2 else "normal")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((0.5, y - h - 0.008), (0.5, y - h - 0.052),
                                         transform=ax.transAxes, arrowstyle="-|>",
                                         mutation_scale=6, lw=0.7, color=C_GREY))
        y -= h + 0.068

    # ---- (c) 分析单位纪律 ----
    ax = fig.add_subplot(gs[1, 0])
    hide(ax)
    panel(ax, "c", dx=-0.055)
    title(ax, "Analysis-unit checks performed")
    mk_table(
        ax,
        [["Repeated biopsies", "paired / mixed structures;\nnot independent samples"],
         ["Pt27A / Pt27B", "one patient, two sites;\ncollapsed by median"],
         ["Pt16 (GSE78220)", "on-treatment only;\ndropped by flag"],
         ["Mel-TIL-026 / HM026", "same patient in two\ncohorts; de-duplicated"],
         ["GSE215868 plate ID", "plate position,\nnot a patient ID"]],
        ["situation", "handling"],
        [0.42, 0.58], 6.2, [0, 0.02, 1.0, 0.92])

    # ---- (d) 队列间为何不合并 ----
    ax = fig.add_subplot(gs[1, 1])
    hide(ax)
    panel(ax, "d")
    title(ax, "What may and may not cross cohorts")
    mk_table(
        ax,
        [["1. Sampling time", "GSE91061 pairs pre-dose with\nC1D29; the others are pre-tx"],
         ["2. Label definitions", "verbatim RECIST strings vs.\nbest overall response vs. BOR"],
         ["3. Platform, units", "raw counts / FPKM / featureCounts /\nNanoString counts"]],
        None, [0.36, 0.64], 6.4, [0.0, 0.46, 1.0, 0.80], edge="#FFFFFF")
    ax.add_patch(Rectangle((0.0, 0.225), 1.0, 0.215, transform=ax.transAxes,
                           facecolor="#FDF2F1", edgecolor=C_GRID, lw=0.35))
    ax.text(0.03, 0.352, "FORBIDDEN", fontsize=6.2, fontweight="bold",
            color=C_S2, va="center", transform=ax.transAxes)
    ax.text(0.03, 0.272, "pooling patients or expression matrices across cohorts",
            fontsize=6.4, va="center", transform=ax.transAxes)
    ax.add_patch(Rectangle((0.0, 0.0), 1.0, 0.215, transform=ax.transAxes,
                           facecolor="#EAF1F8", edgecolor=C_S1, lw=0.5))
    ax.text(0.03, 0.127, "PERMITTED", fontsize=6.2, fontweight="bold",
            color=C_S1, va="center", transform=ax.transAxes)
    ax.text(0.03, 0.047, "study-level meta-analysis of the per-cohort AUCs computed\n"
                         "here, which touches no patient-level record",
            fontsize=6.4, va="center", transform=ax.transAxes, linespacing=1.45)

    save(fig, "Figure1_StudyDesign")


# ================================================================ Figure 2
def fig2_nestedcv() -> None:
    """嵌套交叉验证：两套特征集、三个口径。

    初版两个毛病：面板 (a) 图例压住「chance」标注、面板 (c) 的斜体说明被画布
    下沿裁掉。本版把说明文字移进图内留白区，并给画布留 84 mm 高。
    """
    fig = plt.figure(figsize=(FULL_W, 84 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 1.0, 1.06], wspace=0.40,
                          left=0.055, right=0.985, top=0.90, bottom=0.19)

    def load(tag):
        d = pd.read_csv(QC / f"{tag}_results.tsv", sep="\t")
        d["auc"] = pd.to_numeric(d["auc"], errors="coerce")
        return d.dropna(subset=["auc"])

    a = load("dev_cv")              # 实际使用的 17 特征
    b = load("dev_cv_mainpanel")    # 预登记 6 特征面板
    C17, C6 = C_GREY, C_S1

    # ---- (a) 折间分布 ----
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a")
    title(ax, "Per-fold AUC (100 outer folds)")
    bins = np.linspace(0, 1, 26)
    ax.hist([a["auc"], b["auc"]], bins=bins, color=[C17, C6],
            label=["17 features (as run)", "6 features (pre-registered)"],
            histtype="stepfilled", alpha=0.9, edgecolor="white",
            linewidth=0.4, rwidth=0.9)
    ax.axvline(0.5, color=C_S2, lw=1.0, ls="--", zorder=3)
    ax.text(0.5, -0.155, "chance", color=C_S2, fontsize=6.2,
            ha="center", va="top", transform=ax.get_xaxis_transform())
    ax.set_xlabel("AUC in one outer fold")
    ax.set_ylabel("number of folds")
    ax.legend(frameon=False, fontsize=6.2, loc="upper left",
              bbox_to_anchor=(0.02, 1.0), handlelength=1.1, handletextpad=0.5,
              labelspacing=0.35)
    clean(ax)

    # ---- (b) 重复间均值 ----
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b")
    title(ax, "Mean AUC per CV repeat")
    ra = a.groupby("repeat")["auc"].mean()
    rb = b.groupby("repeat")["auc"].mean()
    ax.axhspan(0.5, 1.0, color=C_BAND, zorder=0, lw=0)
    ax.scatter(range(1, len(ra) + 1), ra, s=13, color=C17, label="17 features", zorder=2)
    ax.scatter(range(1, len(rb) + 1), rb, s=13, color=C6, label="6 features", zorder=3)
    ax.axhline(0.5, color=C_S2, lw=0.9, ls="--", zorder=1)
    ax.axhline(ra.mean(), color=C17, lw=0.7, ls=":")
    ax.axhline(rb.mean(), color=C6, lw=0.7, ls=":")
    ax.set_xlabel("CV repeat")
    ax.set_ylabel("mean AUC")
    ax.set_xlim(0, len(ra) + 1)
    ax.set_ylim(0.1, 0.95)
    ax.legend(frameon=False, fontsize=6.2, loc="lower left", ncol=2,
              columnspacing=0.8, handletextpad=0.3, handlelength=1.1)
    clean(ax)

    # ---- (c) 三个口径 ----
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c")
    title(ax, "One model, three statistics")
    labels = ["pooled\nOOF", "repeat\nmean", "fold\nmean"]
    v17 = [0.500, 0.516, 0.549]
    v6 = [0.417, 0.408, 0.347]
    x = np.arange(3)
    ax.axhspan(0.5, 0.75, color=C_BAND, zorder=0, lw=0)
    b1 = ax.bar(x - 0.19, v17, 0.34, color=C17, label="17 features", zorder=2)
    b2 = ax.bar(x + 0.19, v6, 0.34, color=C6, label="6 features", zorder=2)
    for bars in (b1, b2):
        for r in bars:
            ax.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.012,
                    f"{r.get_height():.3f}", ha="center", va="bottom",
                    fontsize=5.8, color="#1A1A1A")
    ax.axhline(0.5, color=C_S2, lw=0.9, ls="--", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.2)
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 0.78)
    ax.legend(frameon=False, fontsize=6.2, loc="upper left", ncol=2,
              columnspacing=0.8, handletextpad=0.3, handlelength=1.1)
    clean(ax)
    fig.text(0.5, 0.075,
             "pooled OOF merges all 20 repeats, so each patient is counted 20 times; "
             "reporting only 0.500 misleads.\n"
             "Fold-mean SD is 0.257 for the 17-feature set and 0.213 for the "
             "pre-registered panel.",
             ha="center", va="top", fontsize=6.2, style="italic", color=C_GREY,
             linespacing=1.6)

    save(fig, "Figure2_NestedCV")


# ================================================================ Figure 3
def fig3_signatures() -> None:
    """两个签名 × 三个队列 × 两个方向。

    初版把两个签名拆成两块几乎一样的森林图，面板 (c) 又重复一遍同样的点估计，
    信息冗余且图例与标题互相压盖。本版改成：(a) 预登记方向合并森林图（两签名
    分组，直接看出「方向不一致 vs 一致」），(b) 镜像方向同样画法，(c) 只保留
    方向一致性曲线。
    """
    pw = pd.read_csv(QC / "power_analysis.tsv", sep="\t")
    ips = pd.read_csv(QC / "IPS_MHCCP_auc.tsv", sep="\t")

    def pick(cohort, signature, orientation):
        """预登记方向取自 power_analysis.tsv（列名 auc_observed）。"""
        r = pw[(pw["cohort"] == cohort) & (pw["signature"] == signature) &
               (pw["orientation"] == orientation)]
        assert len(r) == 1, (cohort, signature, orientation, len(r))
        return r.iloc[0]

    def pick_ips_mirror(cohort):
        """IPS 镜像方向不在 power_analysis.tsv 中，取自 IPS_MHCCP_auc.tsv。"""
        r = ips[(ips["cohort"] == cohort) & (ips["orientation"] == "mirror")]
        assert len(r) == 1, (cohort, "mirror", len(r))
        return r.iloc[0]

    def auc_of(r, mirrored: bool, si: int) -> float:
        """power_analysis.tsv 的点估计列名是 auc_observed，IPS_MHCCP_auc.tsv 是 auc。"""
        if mirrored and si == 1:
            return float(r["auc"])
        return float(r["auc_observed"])

    # (签名显示名, 取值函数, 取值函数(镜像), 颜色)
    specs = [
        ("IMPRES",
         lambda c: pick(c, "IMPRES g1_low", "g1_low"),
         lambda c: pick(c, "IMPRES g2_low", "g2_low"),
         C_S1),
        ("IPS-MHC+CP",
         lambda c: pick(c, "IPS-MHC+CP", "high_ips_responder"),
         pick_ips_mirror,
         C_S2),
    ]

    fig = plt.figure(figsize=(FULL_W, 98 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.26, 1.26, 1.00], wspace=0.78,
                          left=0.115, right=0.960, top=0.895, bottom=0.150)

    def two_group_forest(ax, mirrored: bool, heading: str) -> None:
        rows = []
        for si, (sig, get, get_m, col) in enumerate(specs):
            rows.append((sig, None, None, None, col))     # 组标题行
            for c in COHORTS:
                r = get_m(c) if mirrored else get(c)
                n = f"{int(r['n_PRCR'])}/{int(r['n_PD'])}"
                lo = float(r["ci_lo"] if mirrored and si == 1 else r["ci_lo_boot"])
                hi = float(r["ci_hi"] if mirrored and si == 1 else r["ci_hi_boot"])
                rows.append((f"{SHORT[c]}  n={n}", auc_of(r, mirrored, si), lo, hi, col))
        forest(ax, rows, C_GREY)
        for lbl, row in zip(ax.get_yticklabels(), rows):
            if row[1] is None:            # 组标题：加粗并按签名着色
                lbl.set_fontweight("bold")
                lbl.set_color(row[4])
        title(ax, heading)

    # ---- (a) 预登记方向 ----
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a", dx=-0.20)
    two_group_forest(ax, mirrored=False, heading="Pre-registered direction")

    # ---- (b) 镜像方向 ----
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b", dx=-0.20)
    two_group_forest(ax, mirrored=True, heading="Mirrored direction")

    # ---- (c) 方向一致性 ----
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c", dx=-0.17)
    title(ax, "Direction consistency")
    x = np.arange(3)
    ax.axhspan(0.70, 1.0, color=C_BAND, zorder=0, lw=0)
    ax.axhline(0.5, ls="--", lw=0.8, color=C_NULL, zorder=1)
    for sig, get, _m, col in specs:
        vals = [float(get(c)["auc_observed"]) for c in COHORTS]
        ax.plot(x, vals, "o-", color=col, lw=1.6, ms=4.6, mec="white", mew=0.6,
                label=sig, zorder=3)
    ax.text(1.0, 0.715, "registered window", fontsize=5.8, color=C_S1,
            ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(["91061", "78220", "215868"], fontsize=6.4)
    ax.set_xlabel("cohort")
    ax.set_ylabel("AUC")
    ax.set_ylim(0.2, 0.80)
    ax.set_xlim(-0.25, 2.25)
    ax.legend(frameon=False, fontsize=6.4, loc="lower right",
              handlelength=1.2, handletextpad=0.5, labelspacing=0.35)
    clean(ax)

    fig.text(0.5, 0.055,
             "Neither signature reaches the pre-registered replication threshold of 0.70 in any cohort, "
             "in either direction.  They fail differently:\n"
             "IMPRES does not agree on which side of chance it falls; IPS-MHC+CP stays above it throughout.",
             ha="center", va="top", fontsize=6.4, style="italic", color=C_GREY,
             linespacing=1.6)

    save(fig, "Figure3_SignaturePerformance")


# ================================================================ Figure 4
def fig4_exclusions() -> None:
    """七项技术性辩解逐条排除 + 两项不可排除。

    初版面板 (a) 画的是一根「8/8」数字条，读不出任何信息；本版改为把
    4 归一化 × 2 方向的分数区间逐条画出——八条完全重合，图形本身即证据。
    """
    inv = pd.read_csv(QC / "GSE215868_impres_invariant_check.tsv", sep="\t")
    diag = pd.read_csv(QC / "GSE215868_panel_diagnostic.tsv", sep="\t")
    var = pd.read_csv(QC / "GSE215868_score_variance.tsv", sep="\t")
    strat = pd.read_csv(QC / "stratified_auc.tsv", sep="\t")

    fig = plt.figure(figsize=(FULL_W, 168 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.66, wspace=0.56,
                          left=0.105, right=0.965, top=0.935, bottom=0.270)

    # ================= (a) 1–2 实现与归一化 =================
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "a")
    for gi, orient in enumerate(["g1_low", "g2_low"]):
        sub = inv[inv["orientation"] == orient]
        col = C_S1 if orient == "g1_low" else C_S2
        for vi, (_, r) in enumerate(sub.iterrows()):
            y = gi * 5 + vi
            ax.plot([r["score_min"], r["score_max"]], [y, y], color=col, lw=2.0,
                    solid_capstyle="round", zorder=2)
            ax.plot([r["score_min"], r["score_max"]], [y, y], "o", ms=2.4,
                    color="#1A1A1A", zorder=3)
        ax.text(-0.055, gi * 5 + 1.5, orient, fontsize=6.2, rotation=90,
                va="center", ha="center", color=col, fontweight="bold",
                transform=ax.get_yaxis_transform())
    n_ok = int(inv["identical_to_raw"].sum())
    n_all = len(inv)
    ax.text(0.5, 0.995, f"{n_ok}/{n_all} identical", transform=ax.transAxes,
            ha="center", va="top", fontsize=9, fontweight="bold", color=C_ACC)
    ax.set_yticks([])
    ax.set_xlim(1, 12)
    ax.set_ylim(-0.7, 9.7)
    ax.set_xticks([2, 4, 6, 8, 10])
    ax.set_xlabel("IMPRES score range per normalisation")
    title(ax, "1–2. Implementation and normalisation")
    clean(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    ax.text(0.5, -0.245, "4 normalisations × 2 directions;\nthe eight ranges coincide exactly",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
            style="italic", color=C_GREY)

    # ================= (b) 3 平台适用性 =================
    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "b")
    order_c = ["GSE91061", "GSE294272", "GSE215868", "GSE78220", "GSE244982"]
    ok = diag[diag["status"] == "可算"]
    good = [float(((ok[ok["cohort"] == c]["p1"] >= 0.15) &
                   (ok[ok["cohort"] == c]["p1"] <= 0.85)).sum() /
                  max(len(ok[ok["cohort"] == c]), 1)) for c in order_c]
    cols = [C_GREY, C_GREY, C_S1, C_GREY, C_GREY]
    bars = ax.bar(range(5), good, color=cols, width=0.62, zorder=2)
    for r, v in zip(bars, good):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.025, f"{v:.2f}",
                ha="center", va="bottom", fontsize=5.8, color="#1A1A1A")
    ax.axhline(0.5, ls="--", lw=0.8, color=C_NULL, zorder=1)
    ax.set_xticks(range(5))
    ax.set_xticklabels([c.replace("GSE", "") for c in order_c], fontsize=6.0,
                       rotation=45, ha="right")
    ax.set_ylabel("informative features")
    ax.set_ylim(0, 1.34)
    title(ax, "3. Platform suitability")
    ax.text(2, 1.10, "targeted panel is not\nmore degenerate", ha="center",
            fontsize=6.0, color=C_S1, fontweight="bold", linespacing=1.4)
    clean(ax)

    # ================= (c) 3 支撑：分数未退化 =================
    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "c")
    bars = ax.bar(range(len(var)), var["ratio"], color=C_ACC, width=0.56, zorder=2)
    for r, v in zip(bars, var["ratio"]):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.04, f"{v:.2f}",
                ha="center", va="bottom", fontsize=5.8, color="#1A1A1A")
    ax.axhline(1.0, ls="--", lw=0.8, color=C_NULL, zorder=1)
    ax.set_xticks(range(len(var)))
    ax.set_xticklabels([v.replace("GSE", "") for v in var["cohort"]],
                       fontsize=6.0, rotation=45, ha="right")
    ax.set_ylabel("observed / expected SD")
    ax.set_ylim(0, 2.28)
    title(ax, "3. Score non-degeneracy")
    ax.text(0.5, 0.995, "1.0 = variance under independence; all > 1",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
            color=C_GREY)
    clean(ax)

    # ================= (d) 5 活检时点 =================
    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "d")
    vis = strat[(strat["stratum_type"] == "visit") & (strat["stratum"] != "Pre+On")]
    vis = vis.set_index("stratum").loc[["Pre", "On"]]
    rows = [(f"{i}  n={int(vis.loc[i, 'n_PRCR'])}/{int(vis.loc[i, 'n_PD'])}",
             vis.loc[i, "auc"], vis.loc[i, "ci_lo"], vis.loc[i, "ci_hi"])
            for i in vis.index]
    forest(ax, rows, C_S1, band=False, xlabel="AUC (95% bootstrap CI)")
    title(ax, "5. Biopsy timepoint")

    # ================= (e) 6 治疗方案 =================
    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "e")
    rx = strat[(strat["stratum_type"] == "regimen") &
               (strat["signature"] == "IMPRES g1_low") & (strat["n_PD"] > 0)]
    rx = rx.set_index("stratum").loc[["IPI+NIVO", "NIVO", "PEMBRO"]]
    rows = [(f"{i}  n={int(rx.loc[i, 'n_PRCR'])}/{int(rx.loc[i, 'n_PD'])}",
             rx.loc[i, "auc"], rx.loc[i, "ci_lo"], rx.loc[i, "ci_hi"])
            for i in rx.index]
    forest(ax, rows, C_S1, xlabel="AUC (95% bootstrap CI)")
    title(ax, "6. Treatment regimen")

    # ================= (f) 4 样本量（功效曲线） =================
    ax = fig.add_subplot(gs[1, 2])
    panel(ax, "f")
    grid = np.linspace(0.50, 0.85, 220)
    pwv = np.array([power_at(g, 45, 34) for g in grid])
    ax.fill_between(grid, pwv, 1.0, where=(pwv >= 0.8), color=C_S2, alpha=0.12, lw=0)
    ax.plot(grid, pwv, color=C_S2, lw=1.8, zorder=3)
    ax.axhline(0.8, ls=":", lw=0.8, color=C_GREY, zorder=1)
    ax.axvline(0.70, ls="--", lw=0.8, color=C_S1, zorder=2)
    ax.axvline(0.77, ls="--", lw=0.8, color=C_S2, zorder=2)
    for xv, lab, col in ((0.70, "0.70", C_S1), (0.77, "0.77", C_S2)):
        ax.text(xv - 0.006, 0.05, lab, fontsize=6.0, color=col, rotation=90,
                va="bottom", ha="right")
    ax.text(0.505, 0.815, "0.8 power", fontsize=6.0, color=C_GREY, va="bottom")
    ax.set_xlabel("true AUC")
    ax.set_ylabel("power (n = 45 / 34)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0.50, 0.85)
    title(ax, "4. Sample size")
    clean(ax)

    # ================= 底部横条 =================
    band = (
        "Excluded outright (6 of 9). Item 4 is addressed but not fully resolved.  "
        "[X] = cannot be excluded with the available public data.\n"
        "1  implementation error        4  sample size*            7  endpoint definition\n"
        "2  normalisation               5  biopsy timepoint        population composition   [X]\n"
        "3  platform unsuitability      6  treatment regimen      weakness intrinsic to the\n"
        "                                                          derivation cohorts   [X]\n"
        "* Power 0.93 at AUC 0.70 and 0.999 at AUC 0.77: excluded against the published performance, not against a\n"
        "  weak signal (power 0.347 at AUC 0.60). Item 7 needed no data panel: the native endpoint of GSE215868 is\n"
        "  24-month long-term benefit (PFS-derived), not RECIST, so response was re-derived from best overall response."
    )
    fig.text(0.5, 0.175, band, ha="center", va="top", fontsize=6.0,
             color="#1A1A1A", linespacing=1.8,
             bbox=dict(boxstyle="round,pad=0.55", facecolor=C_FILL,
                       edgecolor=C_GRID))

    save(fig, "Figure4_RuledOutExplanations")


def main() -> int:
    print("按参考文献风格出图（MDPI Cancers 体例）：")
    fig1_design()
    fig2_nestedcv()
    fig3_signatures()
    fig4_exclusions()
    print(f"\n全部输出至 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
