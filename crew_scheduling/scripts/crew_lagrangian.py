# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt"
_lines = open(DATA).read().splitlines()
N, T = map(int, _lines[0].split())
tasks = [None] + [tuple(map(int, l.split())) for l in _lines[1:1+N]]
arc_cost = {}
for l in _lines[1+N:]:
    i, j, c = map(int, l.split()); arc_cost[(i, j)] = c
K_MIN = 27
EPS = 1e-7
M_DUMMY = 10**6

def enumerate_pool():
    adj = {}
    for (i, j), c in arc_cost.items():
        adj.setdefault(i, []).append((j, c))
    paths = []
    for start in range(1, N+1):
        s0 = tasks[start][0]
        stack = [(start, [start], 0)]
        while stack:
            u, seq, c = stack.pop()
            paths.append((tuple(seq), c, tasks[u][1]-s0))
            for v, cv in adj.get(u, []):
                span = tasks[v][1] - s0
                if span <= T:
                    stack.append((v, seq+[v], c+cv))
    return paths

paths = enumerate_pool()
P = len(paths)
pseq = [p[0] for p in paths]
pcost = [p[1] for p in paths]
pmask = []
for s2 in pseq:
    m2 = 0
    for i in s2:
        m2 |= (1 << i)
    pmask.append(m2)

def col_arcs(seq):
    a = []
    prev = 0
    for j in seq:
        a.append((prev, j))
        prev = j
    a.append((prev, N+1))
    return a

def run_lagrangian(verbose=True, max_iter=600):
    def rc_of(p, lam):
        s2 = 0.0
        mm = pmask[p]
        while mm:
            lb = mm & -mm
            s2 += lam[lb.bit_length()-1]
            mm -= lb
        return pcost[p] - s2
    lam = [0.0]*(N+1)
    rho = 2.0
    UB_anchor = 3139.0        # 已知原问题上界（池 IP 解），步长锚点须 >= 对偶最优值
    best_LB = -1e18
    best_UB = None
    best_routes = None
    no_improve = 0
    wall0 = time.time()
    for it in range(max_iter):
        neg = sorted((rc_of(p, lam), p) for p in range(P))
        S = [p for rc, p in neg[:K_MIN] if rc < -EPS]
        L = sum(lam[1:]) + sum(rc_of(p, lam) for p in S)
        if L > best_LB:
            best_LB = L
            no_improve = 0
        else:
            no_improve += 1
            if no_improve >= 30:
                rho /= 2.0
                no_improve = 0
        cnt = [0]*(N+1)
        for p in S:
            mm = pmask[p]
            while mm:
                lb = mm & -mm
                cnt[lb.bit_length()-1] += 1
                mm -= lb
        g = [0.0]*(N+1)
        for i in range(1, N+1):
            g[i] = 1.0 - cnt[i]
        # 贪心修复
        cov = [False]*(N+1)
        Sr = list(S)
        for p in Sr:
            mm = pmask[p]
            while mm:
                lb = mm & -mm
                cov[lb.bit_length()-1] = True
                mm -= lb
        for i in range(1, N+1):
            if cov[i]:
                continue
            best = None
            for p in range(P):
                if p in Sr or not ((pmask[p] >> i) & 1):
                    continue
                if len(Sr) >= K_MIN:
                    break
                if best is None or pcost[p] < pcost[best]:
                    best = p
            if best is None:
                break
            Sr.append(best)
            mm = pmask[best]
            while mm:
                lb = mm & -mm
                cov[lb.bit_length()-1] = True
                mm -= lb
        else:
            c = sum(pcost[p] for p in Sr)
            if best_UB is None or c < best_UB:
                best_UB = c
                best_routes = [pseq[p] for p in Sr]
        gnorm2 = sum(g[i]*g[i] for i in range(1, N+1))
        step = 0.0 if gnorm2 <= 1e-12 else min(rho*(UB_anchor - L)/gnorm2, 100.0)
        for i in range(1, N+1):
            lam[i] = min(1000.0, max(0.0, lam[i] + step*g[i]))
        if verbose and (it % 100 == 0 or it == max_iter-1):
            gap = (best_UB-best_LB)/best_LB*100 if (best_LB > 1e-9 and best_UB) else float('nan')
            print(f"iter {it+1:4d}: L={L:9.4f} best_LB={best_LB:9.4f} best_UB={best_UB if best_UB else 0:9.4f} gap={gap:.4f}% rho={rho:.4f}")
        if rho < 1e-5 or time.time()-wall0 > 110:
            if verbose:
                print("停机 at iter", it+1)
            break
    # 最终 MIP 修复
    m = cp_model.CpModel()
    yv = [m.NewBoolVar(f"y{p}") for p in range(P)]
    for i in range(1, N+1):
        m.Add(sum(yv[p] for p in range(P) if (pmask[p] >> i) & 1) >= 1)
    m.Add(sum(yv) <= K_MIN)
    m.Minimize(sum(pcost[p]*yv[p] for p in range(P)))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = 60
    st = s.Solve(m)
    routes = [pseq[p] for p in range(P) if s.Value(yv[p])]
    if verbose:
        print(f"MIP 修复: {s.StatusName(st)} obj={s.ObjectiveValue()} crew={len(routes)}")
        print(f"对偶下界 {round(best_LB,4)} vs 修复 {s.ObjectiveValue()} | gap {(s.ObjectiveValue()-best_LB)/best_LB*100 if best_LB>1e-9 else 0:.4f}%")
    return dict(lb=best_LB, greedy_ub=best_UB, ip=s.ObjectiveValue(), routes=routes,
                iterations=it+1, status=s.StatusName(st))

if __name__ == "__main__":
    run_lagrangian()
