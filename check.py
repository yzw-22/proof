# -*- coding: utf-8 -*-
# 判定 2^n - n (n < C) 是否为完全平方：偶数已证无解；奇支路 = 模筛(严格) + 幸存者精确复核
import math, sys, time
import numpy as np

def is_square(x):                       # 任意位长精确判定（标准库）
    if x < 0: return False
    s = math.isqrt(x); return s * s == x

def good_table(m):
    """布尔表 G over Z_L, L = lcm(m, ord_m(2))：G[t] = (2^t - t) mod m 是否为剩余"""
    o, cur = 1, 2 % m
    while cur != 1: cur = cur * 2 % m; o += 1
    L = m * o // math.gcd(m, o)
    small, c = np.empty(o, np.int64), 1 % m
    for i in range(o): small[i] = c; c = c * 2 % m
    t = np.arange(L, dtype=np.int64); pw = small[t % o]
    Q = np.zeros(m, bool)
    for x in range(m): Q[x * x % m] = True
    return L, Q[(pw - t) % m]

def build_wheel():                      # 一级轮子: n≡7(mod 8) + 模 9,5,7,11,17,13,23,19
    M, S = 8, np.array([7], dtype=np.int64)
    for m in (9, 5, 7, 11, 17, 13, 23, 19):
        L, G = good_table(m)
        steps = L // math.gcd(M, L)
        tvec = np.arange(steps, dtype=np.int64) * M
        out, CH = [], max(1, 3_000_000 // steps)
        for i in range(0, S.size, CH):
            cand = S[i:i+CH][:, None] + tvec[None, :]
            out.append(cand[G[cand % L]])
        S = np.concatenate(out); M *= steps
    return M, S                         # M1=2,677,114,440, |S1|=2,177,280, 密度 8e-4

def primes_upto(n):
    s = np.ones(n, bool); s[:2] = False
    for d in range(2, int(n**.5)+1):
        if s[d]: s[d*d::d] = False
    return [int(p) for p in np.nonzero(s)[0]]

M1, S1 = build_wheel()
FILT = [good_table(p) for p in primes_upto(430) if p > 23]   # 73 个二级素数表

def survivors_below(C):
    nq = (C + M1 - 1) // M1
    out = []
    for q in range(nq):                 # n = r + q*M1, r∈S1 —— 分片向量化
        N = S1 + q * M1
        if q == nq - 1: N = N[N < C]
        for L, G in FILT:
            N = N[G[N % L]]             # 一次取模+一次查表 = 一个素数的完整 QR 检验
            if N.size == 0: break
        if N.size: out.append(N)
    return np.concatenate(out) if out else np.zeros(0, np.int64)

def verify(C):
    sv = survivors_below(C)             # 被排除的每个 n 都已有素数证书证明非平方
    for n in map(int, sv):              # 幸存者精确复核（期望极少）
        assert is_square((1 << n) - n), f"n={n} 是疑似新解, 需人工精确复核!"
    return sorted([1] + list(map(int, sv)))    # n=1: 2^1-1=1=1^2; 偶数支路已证无解

if __name__ == "__main__":
    C = int(sys.argv[1]) if len(sys.argv) > 1 else 10**12
    t = time.perf_counter()
    print(f"n < {C:,}: 2^n-n 为完全平方的 n = {verify(C)}   [{time.perf_counter()-t:.1f}s]")
