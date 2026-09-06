# -*- coding: utf-8 -*-
# psieve.py —— 验证 n < C 内 2^n-n 是否为完全平方: 模筛 + 多进程并行
# 用法:  python psieve.py  [C]  [worker数, 默认=CPU核数]
# 例:    python psieve.py 1000000000000            # 自检, 应得 n ∈ {1, 7}
#        python psieve.py 1000000000000000 8        # 正式 10^15
import math, os, sys, time
import multiprocessing as mp
import numpy as np

# ---------- 数学内核 (与已实测的串行版一致) ----------
def good_table(m):
    o, cur = 1, 2 % m
    while cur != 1: cur = cur * 2 % m; o += 1
    L = m * o // math.gcd(m, o)
    small, c = np.empty(o, np.int64), 1 % m
    for i in range(o): small[i] = c; c = c * 2 % m
    t = np.arange(L, dtype=np.int64); pw = small[t % o]
    Q = np.zeros(m, bool)
    for x in range(m): Q[x * x % m] = True          # 0 必须计入 QR!
    return L, Q[(pw - t) % m]

def build_wheel():
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
    return M, S

def primes_upto(n):
    s = np.ones(n, bool); s[:2] = False
    for d in range(2, int(n**.5)+1):
        if s[d]: s[d*d::d] = False
    return [int(p) for p in np.nonzero(s)[0]]

def is_square(x):
    s = math.isqrt(x); return s * s == x

# 全局表: Linux fork 子进程写时复制继承; spawn 时由 initializer 重建
W_M1 = W_S1 = W_FILT = None

def init_tables():
    global W_M1, W_S1, W_FILT
    if W_M1 is not None: return
    W_M1, W_S1 = build_wheel()
    W_FILT = [good_table(p) for p in primes_upto(430) if p > 23]

# ---------- worker ----------
# ★ 修复点: 接收单个元组, 函数内解包 (imap_unordered 不自动解包) ★
def sweep_range(task):
    q0, q1, C = task
    out = []
    for q in range(q0, q1):
        N = W_S1 + q * W_M1                       # 本块候选 (2.18M 个)
        for L, G in W_FILT:                       # 73 筛
            N = N[G[N % L]]
            if N.size == 0: break                 # 杀空早退 (常态)
        if N.size:
            N = N[N < C]
            if N.size: out.append(N)
    return np.concatenate(out) if out else np.zeros(0, np.int64)

# ---------- 主控 ----------
def main():
    C = int(sys.argv[1]) if len(sys.argv) > 1 else 5*10**15
    workers = (int(sys.argv[2]) if len(sys.argv) > 2
               else len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
               else os.cpu_count() or 4)
    t0 = time.perf_counter()
    init_tables()
    qN = (C + W_M1 - 1) // W_M1                   # 总块数 (1e15 → 373,537)
    CHUNK = 32                                    # 每任务 32 块 ≈ 2.3s CPU
    tasks = [(q0, min(q0 + CHUNK, qN), C) for q0 in range(0, qN, CHUNK)]

    # macOS/Windows 用 spawn (安全); Linux 用 fork (表零成本共享)
    if sys.platform.startswith("linux") and "fork" in mp.get_all_start_methods():
        ctx = mp.get_context("fork")
    else:
        ctx = mp.get_context("spawn")
    print(f"C={C:,} | 块 {qN:,} | 任务 {len(tasks):,} | workers={workers} "
          f"({ctx.get_start_method()}) | 建表 {time.perf_counter()-t0:.2f}s", flush=True)

    out, done, t1 = [], 0, time.perf_counter()
    with ctx.Pool(workers, initializer=init_tables) as pool:
        for arr in pool.imap_unordered(sweep_range, tasks, chunksize=1):
            out.append(arr); done += 1
            step = max(1, len(tasks) // 20)
            if done % step == 0 or done == len(tasks):
                el = time.perf_counter() - t1
                print(f"  {done*100//len(tasks):>3}%  {el:7.0f}s  "
                      f"{qN/el:,.0f} 块/s  ETA {el/done*(len(tasks)-done):,.0f}s", flush=True)

    sv = np.concatenate(out) if out else np.zeros(0, np.int64)
    print(f"扫描完成: {time.perf_counter()-t1:,.0f}s  "
          f"(等效 {C/(time.perf_counter()-t1)/1e9:,.1f} G个n/s)", flush=True)

    sols = [1]                                    # n=1: 2-1=1=1²; 偶数支路已证无解
    for n in map(int, sv):                        # 幸存者精确复核
        v = (1 << n) - n
        ok = is_square(v)
        print(f"  幸存 n={n} -> "
              f"{'真解: ' + str(math.isqrt(v)) + '²' if ok else '假阳性(非平方)'}", flush=True)
        if ok: sols.append(n)
    print(f"结论: n < {C:,}: 2^n-n 为完全平方 ⟺ n ∈ {sorted(sols)}  "
          f"[总 {time.perf_counter()-t0:,.0f}s]", flush=True)

if __name__ == "__main__":
    main()

