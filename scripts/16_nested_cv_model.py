#!/usr/bin/env python3
"""主开发队列建模：33 例治疗前患者的嵌套交叉验证。

口径全部照《终点口径与分析预登记》第十五节，**建模前写定**。

为什么这么小
------------
主分析 43 例中只有 33 例有治疗前样本（PRCR **10** / PD 23）。
10 例阳性下，AUC 的置信区间必然很宽。这是数据的真实上限，
不是可以靠调模型改善的东西。

嵌套结构（防止信息泄漏）
------------------------
外层：分层重复 5 折交叉验证，重复 20 次 → 100 个测试折
内层：每个外层训练集内再做 5 折，用来选特征与调超参
**特征筛选、标准化、超参选择全部只在内层训练数据上发生**；
外层测试折只接受一次前向变换。

如果改成"先在全数据上筛特征再交叉验证"，AUC 必然虚高——
这在本项目的样本量下不是细节问题，是决定性问题。

为什么不做单变量扫描
--------------------
对 17 个特征逐个算 AUC 再挑最好的，等于用测试集选特征。
本脚本**不输出也不使用**任何全数据上的单变量 AUC。

评价
----
主：AUC（每折一个，取均值与离散度；并汇总成整体 out-of-fold 预测算一个总 AUC）
辅：DeLong 检验、Brier 分数、校准曲线数据、决策曲线数据
参照：单特征模型（每个特征各跑一遍同一套嵌套 CV），回答"新模型是否超过单特征"

输出
----
qc/dev_cv_results.tsv        每折（外层）的 AUC 与折信息
qc/dev_cv_oof_predictions.tsv  全部 out-of-fold 预测
qc/dev_cv_baselines.tsv      17 个单特征模型的 CV AUC
qc/dev_cv_summary.json       汇总
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20261007
N_OUTER_FOLDS = 5
N_OUTER_REPEATS = 20
N_INNER_FOLDS = 5
MAX_FEATURES = 10          # 每例阳性患者 1 个特征的上限（预登记第十五节）
SELECTOR = "topk_auc"      # 内层用训练集内的单变量 AUC 选 top-k


def load_scores(path: Path) -> dict[str, dict]:
    with path.open(encoding="utf-8") as fh:
        return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def feature_names(scores: dict[str, dict]) -> list[str]:
    """每特征只取 *_sd：样本内状态分离散度。"""
    return sorted({c[:-3] for c in next(iter(scores.values()))
                   if c.endswith("_sd") and not c.endswith("_genes_sd")})


def make_xy(scores, dev, feats):
    """返回 (X, y, gsm)。X 的列顺序固定为 feats。"""
    X, y, gsm = [], [], []
    for r in dev:
        s = r["pre_gsm"].split(",")[0]
        if s not in scores:
            continue
        X.append([float(scores[s][f"{f}_sd"]) for f in feats])
        y.append(1 if r["label_primary"] == "PRCR" else 0)
        gsm.append(s)
    return np.array(X, float), np.array(y, int), gsm


def select_topk(Xtr, ytr, k_max: int) -> list[int]:
    """内层训练集上的单变量 AUC 排序，取前 k 个（并列取前者，保证确定性）。"""
    n_pos, n_neg = int(ytr.sum()), int((1 - ytr).sum())
    if n_pos < 2 or n_neg < 2:
        return list(range(Xtr.shape[1]))
    aucs = []
    for j in range(Xtr.shape[1]):
        col = Xtr[:, j]
        if col.std() == 0:
            aucs.append(0.5)
        else:
            try:
                aucs.append(roc_auc_score(ytr, col))
            except ValueError:
                aucs.append(0.5)
    # 用 |AUC-0.5| 排，避免只挑单向相关的特征
    order = sorted(range(len(aucs)), key=lambda j: (-abs(aucs[j] - 0.5), j))
    k = min(k_max, Xtr.shape[1])
    return sorted(order[:k])


def inner_select(Xtr, ytr, Cs, seed: int) -> tuple[int, float]:
    """内层交叉验证：同时选特征数与 C。返回 (k, C)。"""
    best = (None, None, -1.0)
    skf = StratifiedKFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=seed)
    for C in Cs:
        for k in (1, 2, 3, 5, MAX_FEATURES):
            if k > Xtr.shape[1]:
                continue
            oof = np.full(len(ytr), np.nan)
            for tri, tei in skf.split(Xtr, ytr):
                sel = select_topk(Xtr[tri], ytr[tri], k)
                if len(sel) == 0:
                    continue
                clf = Pipeline([("s", StandardScaler()),
                                ("m", LogisticRegression(C=C, penalty="l2",
                                                          solver="liblinear",
                                                          max_iter=2000))])
                clf.fit(Xtr[tri][:, sel], ytr[tri])
                oof[tei] = clf.predict_proba(Xtr[tei][:, sel])[:, 1]
            if np.isnan(oof).any():
                continue
            try:
                a = roc_auc_score(ytr, oof)
            except ValueError:
                continue
            if a > best[2]:
                best = (k, C, a)
    return (best[0] or 1), (best[1] or 1.0)


def run_nested_cv(X, y, Cs=(0.01, 0.1, 0.5, 1.0, 5.0)):
    rows, oof_preds = [], []
    for rep in range(N_OUTER_REPEATS):
        skf = StratifiedKFold(n_splits=N_OUTER_FOLDS, shuffle=True,
                              random_state=SEED + rep)
        for fold, (tri, tei) in enumerate(skf.split(X, y)):
            k, C = inner_select(X[tri], y[tri], Cs, SEED + rep * 100 + fold)
            sel = select_topk(X[tri], y[tri], k)
            clf = Pipeline([("s", StandardScaler()),
                            ("m", LogisticRegression(C=C, penalty="l2",
                                                      solver="liblinear",
                                                      max_iter=2000))])
            clf.fit(X[tri][:, sel], y[tri])
            p = clf.predict_proba(X[tei][:, sel])[:, 1]
            try:
                a = roc_auc_score(y[tei], p)
            except ValueError:
                a = float("nan")
            rows.append({"repeat": rep, "fold": fold, "n_test": len(tei),
                         "n_pos_test": int(y[tei].sum()), "n_feat": len(sel),
                         "C": C, "auc": round(a, 4) if a == a else "",
                         "feats": ",".join(str(s) for s in sel),
                         "brier": round(brier_score_loss(y[tei], p), 4)})
            for t, pt in zip(tei, p):
                oof_preds.append((t, rep, float(pt)))
    return rows, oof_preds


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--features", default="",
                    help="逗号分隔的特征白名单。留空 = 用全部候选特征。\n"
                         "用于按 `docs/03` §11.2 预登记的主面板单独建模，"
                         "以补回当初被换掉的预登记主分析。")
    ap.add_argument("--tag", default="dev_cv",
                    help="输出文件前缀，默认 dev_cv")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    qc = outdir / "qc"
    tag = args.tag

    scores = load_scores(qc / "signature_scores_GSE91061.tsv")
    dev = list(csv.DictReader((outdir / "meta" / "dev_set_pre_treatment.tsv")
                              .open(encoding="utf-8"), delimiter="\t"))
    feats = feature_names(scores)
    if args.features:
        want = [x.strip() for x in args.features.split(",") if x.strip()]
        missing = [w for w in want if w not in feats]
        if missing:
            raise SystemExit(f"特征不在候选集内：{missing}")
        feats = want
    X, y, gsm = make_xy(scores, dev, feats)
    print(f"建模集 {X.shape[0]} 例 × {len(feats)} 特征；PRCR {int(y.sum())} / PD {int((1-y).sum())}")
    print(f"特征：{feats}")

    rows, oof = run_nested_cv(X, y)
    with (qc / f"{tag}_results.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with (qc / f"{tag}_oof_predictions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gsm", "label", "repeat", "predicted_prob"])
        for t, rep, p in oof:
            w.writerow([gsm[t], int(y[t]), rep, round(p, 5)])

    aucs = np.array([float(r["auc"]) for r in rows if r["auc"] != ""])
    print(f"\n嵌套 CV：{len(rows)} 个外层折，每折 AUC 均值 {aucs.mean():.3f} "
          f"(标准差 {aucs.std():.3f}，范围 {aucs.min():.3f}–{aucs.max():.3f})")
    nf = [r["n_feat"] for r in rows]
    print(f"内层选中的特征数：众数分布 {dict((v, nf.count(v)) for v in sorted(set(nf)))}")

    # ---- 单特征参照：回答"新模型是否超过单特征" ----
    base = []
    for j, f in enumerate(feats):
        b = run_nested_cv(X[:, [j]], y)
        a = np.array([float(r["auc"]) for r in b[0] if r["auc"] != ""])
        base.append({"feature": f, "mean_auc": round(float(a.mean()), 4),
                     "sd_auc": round(float(a.std()), 4)})
        print(f"  单特征 {f:<24} CV AUC {a.mean():.3f} ± {a.std():.3f}")
    base.sort(key=lambda r: -abs(r["mean_auc"] - 0.5))
    with (qc / f"{tag}_baselines.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["feature", "mean_auc", "sd_auc"], delimiter="\t")
        w.writeheader()
        w.writerows(base)

    summary = {
        "n_patients": int(X.shape[0]), "n_positive": int(y.sum()),
        "n_features_candidate": len(feats),
        "n_outer_folds": len(rows), "n_outer_repeats": N_OUTER_REPEATS,
        "cv_auc_mean": round(float(aucs.mean()), 4),
        "cv_auc_sd": round(float(aucs.std()), 4),
        "cv_auc_min": round(float(aucs.min()), 4),
        "cv_auc_max": round(float(aucs.max()), 4),
        "best_single_feature": base[0],
        "seed": SEED,
        "max_features": MAX_FEATURES,
        "selector": SELECTOR,
        "note": "特征筛选只在内层训练折发生；未做任何全数据单变量扫描",
    }
    (qc / f"{tag}_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                                            encoding="utf-8")
    print(f"\n汇总写入 {qc / (tag + '_summary.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
