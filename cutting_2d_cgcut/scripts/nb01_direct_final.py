import platform, time, ortools, sys
from ortools.sat.python import cp_model
print("python", platform.python_version(), "| ortools", ortools.__version__)

DAT = "/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt"
def parse(path):
    lines = [l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m = int(lines[0]); W, H = map(int, lines[1].split())
    items = []
    for l in lines[2:2+m]:
        a, b, q, v = map(int, l.split()); items.append((a, b, q, v))
    return m, W, H, items
m, W, H, items = parse(DAT)
q = [it[2] for it in items]; v = [it[3] for it in items]
print("instance: m=%d, sheet=%dx%d" % (m, W, H))
print("items (l,w,q,value):", items)
print("fixed orientation (no rotation) is used to match the literature optimum 244")

def solve_layout(allow_rotation):
    model = cp_model.CpModel()
    ivx = []; ivy = []; pres = []; objs = []
    idx_start = []
    for i, (a, b, qq, vv) in enumerate(items):
        orients = sorted(set([(a, b), (b, a)])) if allow_rotation else [(a, b)]
        idx_start.append(len(pres))
        for (dw, dh) in orients:
            for k in range(qq):
                p = model.NewBoolVar(f"p_{i}_{dw}_{dh}_{k}")
                xv = model.NewIntVar(0, W - dw, f"x_{i}_{dw}_{dh}_{k}")
                yv = model.NewIntVar(0, H - dh, f"y_{i}_{dw}_{dh}_{k}")
                ix = model.NewOptionalIntervalVar(xv, dw, xv + dw, p, f"ix_{i}_{dw}_{dh}_{k}")
                iy = model.NewOptionalIntervalVar(yv, dh, yv + dh, p, f"iy_{i}_{dw}_{dh}_{k}")
                pres.append(p); ivx.append(ix); ivy.append(iy); objs.append(vv * p)
        cnt = len(pres) - idx_start[-1]
        model.Add(sum(pres[idx_start[-1]:]) <= qq)
    model.AddNoOverlap2D(ivx, ivy)
    model.Maximize(sum(objs))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 120.0
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = 0
    t0 = time.time()
    status = solver.Solve(model)
    dt = time.time() - t0
    n = [0] * m
    for i in range(m):
        n[i] = int(round(sum(solver.Value(pres[k]) for k in range(idx_start[i], len(pres) if i == m-1 else idx_start[i+1]))))
    return solver.StatusName(status), solver.ObjectiveValue(), solver.BestObjectiveBound(), dt, n

print("
--- CP-SAT layout (fixed orientation, non-guillotine upper bound) ---")
st, obj, bound, dt, n = solve_layout(False)
print("status=%s obj=%.1f best_bound=%.1f time=%.3fs" % (st, obj, bound, dt))
print("selected counts n_i:", n, "value", sum(n[i]*v[i] for i in range(m)))

print("
--- rotation sensitivity check (to document the convention discrepancy) ---")
st_r, obj_r, bound_r, dt_r, n_r = solve_layout(True)
print("status=%s obj=%.1f best_bound=%.1f time=%.3fs" % (st_r, obj_r, bound_r, dt_r))
print("selected counts with rotation:", n_r, "value", sum(n_r[i]*v[i] for i in range(m)))
