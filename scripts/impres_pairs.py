"""从作者仓库的 IMPRES-codes 还原 15 个特征对。

索引空间（由 getIMPRESRAT.m / classifyImmuneCOMP.m 的枚举顺序决定）：
  28 个 IC 基因两两组合，C(28,2) = 378 对，分两个区块：
    区块 1  idx 1..378    : 按 i<j 枚举，g1=gene[j], g2=gene[i]
    区块 2  idx 379..756  : 同一枚举顺序，g1=gene[i], g2=gene[j]
  两个区块的 ratio 都是 exp_{g1}/exp_{g2}，故统一规则 F=1 当 exp_{g1} < exp_{g2}。

来源：https://github.com/noamaus/IMPRES-codes
  ADDITIONAL_CODES/Additional feature sets/FEATS.mat  → 15 个特征索引
  ADDITIONAL_CODES/Additional feature sets/CPall.mat  → 28 个 IC 基因（字母序）
"""

IC_GENES = [
    "BTLA", "C10orf54", "CD200", "CD200R1", "CD27", "CD274", "CD276", "CD28",
    "CD40", "CD80", "CD86", "CEACAM1", "CTLA4", "HAVCR2", "IDO1", "IL2RB",
    "LAG3", "PDCD1", "PDCD1LG2", "PVR", "PVRL2", "TIGIT", "TNFRSF14",
    "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9",
]
# ADDITIONAL_CODES/Additional feature sets/FEATS.mat（Matlab v5，uint16，1-based）
FEATURE_IDX = [31, 61, 128, 148, 169, 237, 493, 549, 567, 575, 603, 606, 619, 650, 710]


def _enumerate_pairs(n: int):
    """按 getIMPRESRAT.m 的枚举顺序产出 (i, j)：外层 i 升序、内层 j>i 升序，0-based。"""
    for i in range(n):
        for j in range(i + 1, n):
            yield i, j


def idx_to_pair(idx: int) -> tuple[str, str]:
    """Matlab 1-based 索引 → (g1, g2)，满足 F=1 当 exp_{g1} < exp_{g2}。"""
    n = len(IC_GENES)
    half = n * (n - 1) // 2
    if not 1 <= idx <= 2 * half:
        raise IndexError(f"索引 {idx} 超出 1..{2 * half} 的有序对空间")
    block1 = idx <= half
    k = idx if block1 else idx - half
    for pos, (i, j) in enumerate(_enumerate_pairs(n), start=1):
        if pos == k:
            # 区块 1：g1=gene[j], g2=gene[i]；区块 2：g1=gene[i], g2=gene[j]
            return ((IC_GENES[j], IC_GENES[i]) if block1
                    else (IC_GENES[i], IC_GENES[j]))
    raise IndexError(idx)


if __name__ == "__main__":
    print(f"{'#':>3} {'idx':>4} {'区块':>5}  {'g1 (低)':<12} {'g2 (高)':<12}")
    pairs = []
    for k, idx in enumerate(FEATURE_IDX, 1):
        g1, g2 = idx_to_pair(idx)
        block = "1" if idx <= 378 else "2"
        pairs.append((g1, g2))
        print(f"{k:>3} {idx:>4} {'blk' + block:>5}  {g1:<12} {g2:<12}")
    assert len(pairs) == 15
    assert len(set(map(frozenset, pairs))) == 15, "出现重复的无序基因对"
    genes = sorted({g for p in pairs for g in p})
    print(f"\n涉及 {len(genes)} 个基因：{genes}")
