import sys, math
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
import importlib.util
spec = importlib.util.spec_from_file_location("lbbd", "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts/lbbd.py")
# 只取数据部分：手动重建
DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
lines = open(DATA).read().splitlines()
custs = []
for l in lines:
    s = l.split()
    if len(s) == 7 and s[0].isdigit():
        custs.append((int(s[0]), int(s[1]), int(s[2]), int(s[3]), int(s[4]), int(s[5]), int(s[6])))
custs.sort()
depot = custs[0]; cs = custs[1:]
n = len(cs)
xc = [depot[1]] + [c[1] for c in cs]
yc = [depot[2]] + [c[2] for c in cs]
dem = [0] + [c[3] for c in cs]
ready = [0] + [c[4] for c in cs]
due = [0] + [c[5] for c in cs]
svc = [0] + [c[6] for c in cs]
depot_due = depot[5]
cap = 200
def dist(i, j):
    return math.hypot(xc[i]-xc[j], yc[i]-yc[j])

def tw_feasible(S):
    S = sorted(S, key=lambda i: due[i])
    if sum(dem[i] for i in S) > cap:
        return False
    memo = {}
    def dfs(last, mask, t):
        key = (last, mask)
        if key in memo:
            if t >= memo[key]:
                return False
        memo[key] = t
        if mask == (1 << len(S)) - 1:
            return t + svc[last] + dist(last, 0) <= depot_due + 1e-9
        for pos, j in enumerate(S):
            if mask & (1 << pos):
                continue
            arr = max(ready[j], t + svc[last] + dist(last, j))
            if arr <= due[j] and dfs(j, mask | (1 << pos), arr):
                return True
        return False
    for pos, j in enumerate(S):
        arr = max(ready[j], 0 + svc[0] + dist(0, j))
        if arr <= due[j] and dfs(j, 1 << pos, arr):
            return True
    return False

ang = sorted(range(1, n+1), key=lambda i: math.atan2(yc[i]-yc[0], xc[i]-xc[0]))
print("角度排序:", ang)
cum = [0]*(n+1)
for t_, i in enumerate(ang):
    cum[t_+1] = cum[t_] + dem[i]
print("前缀:", cum)
g1, g2, g3 = tuple(ang[:11]), tuple(ang[11:19]), tuple(ang[19:])
print("最优划分: g1=", g1, "dem=", sum(dem[i] for i in g1))
print("g2=", g2, "dem=", sum(dem[i] for i in g2))
print("g3=", g3, "dem=", sum(dem[i] for i in g3))
print("tw_feasible g1:", tw_feasible(g1))
print("tw_feasible g2:", tw_feasible(g2))
print("tw_feasible g3:", tw_feasible(g3))
# 单个最优路线顺序直接验证
r1 = (13,17,18,19,15,16,14,12)
t = 0.0; prev = 0
ok = True
for j in r1:
    t = max(ready[j], t + svc[prev] + dist(prev, j))
    print("  visit", j, "t=", round(t,3), "due=", due[j], "ok" if t <= due[j] else "FAIL")
    if t > due[j]: ok = False
    prev = j
print("返回:", round(t + svc[prev] + dist(prev, 0), 3), "<= 1236:", t + svc[prev] + dist(prev, 0) <= depot_due + 1e-9)
