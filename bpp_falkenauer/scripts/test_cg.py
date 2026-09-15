import math, time, datetime
from ortools.math_opt.python import mathopt

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes

name,C,n,best,sizes = parse()
print("instance", name, "C", C, "n", n, "file_best", best, "sum", sum(sizes))
print("ceil LB", math.ceil(sum(sizes)/C))

columns = [tuple([i]) for i in range(n)]
colset = set(columns)

def build_rmp(columns):
    model = mathopt.Model(name="rmp")
    lam = [model.add_variable(lb=0.0, is_integer=False, name=f"lam_{j}") for j in range(len(columns))]
    cons = []
    for i in range(n):
        expr = mathopt.fast_sum(lam[j] for j,col in enumerate(columns) if i in col)
        cons.append(model.add_linear_constraint(expr == 1.0, name=f"cov_{i}"))
    model.minimize(mathopt.fast_sum(lam))
    return model, lam, cons

def knapsack_exact(profits, weights, capacity):
    n = len(weights)
    dp = [0.0]*(capacity+1)
    keep = [[False]*(capacity+1) for _ in range(n)]
    for i in range(n):
        w = weights[i]; p = profits[i]
        for c in range(capacity, w-1, -1):
            cand = dp[c-w] + p
            if cand > dp[c] + 1e-12:
                dp[c] = cand
                keep[i][c] = True
    bestc = max(range(capacity+1), key=lambda c: dp[c])
    val = dp[bestc]
    items = []
    c = bestc
    for i in range(n-1, -1, -1):
        if keep[i][c]:
            items.append(i)
            c -= weights[i]
    return val, items

TOL = 1e-7
t0=time.time()
lp_solve_time=0.0
pricing_time=0.0
iters=0
obj=None
for it in range(2000):
    iters=it+1
    model, lam, cons = build_rmp(columns)
    ts=time.time()
    res = mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(enable_output=False))
    lp_solve_time += time.time()-ts
    if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
        print("LP not optimal", res.termination.reason); break
    obj = res.objective_value()
    pi = [res.dual_values()[cons[i]] for i in range(n)]
    ts=time.time()
    profit, items = knapsack_exact(pi, sizes, C)
    pricing_time += time.time()-ts
    rc = 1.0 - profit
    if rc < -TOL:
        newcol = tuple(sorted(items))
        if newcol in colset:
            print("duplicate column rc", rc, "items", len(items), "profit", profit); break
        columns.append(newcol); colset.add(newcol)
    else:
        break
    if it % 20 == 0:
        print("iter", it, "obj", round(obj,6), "ncol", len(columns), "rc", round(rc,9))

print("CG done: iters", iters, "ncol", len(columns), "LP obj", obj)
print("ceil(LP)", math.ceil(obj), "time", round(time.time()-t0,2), "lp_time", round(lp_solve_time,2), "pricing_time", round(pricing_time,2))

# integer MIP over pool
model = mathopt.Model(name="ip")
lam = [model.add_variable(lb=0.0, ub=float(n), is_integer=True, name=f"lam_{j}") for j in range(len(columns))]
for i in range(n):
    expr = mathopt.fast_sum(lam[j] for j,col in enumerate(columns) if i in col)
    model.add_linear_constraint(expr == 1.0, name=f"cov_{i}")
model.minimize(mathopt.fast_sum(lam))
ts=time.time()
res = mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
ip_time=time.time()-ts
print("IP status", res.termination.reason)
if res.has_primal_feasible_solution:
    print("IP obj", res.objective_value())
    try: print("IP bound", res.best_objective_bound())
    except Exception as e: print("no bound", e)
    vals = res.variable_values()
    used=[j for j in range(len(columns)) if vals[lam[j]]>0.5]
    total=sum(vals[lam[j]] for j in used)
    print("used patterns", len(used), "total bins", total)
    loads=[sum(sizes[i] for i in columns[j]) for j in used]
    print("min/max load", min(loads), max(loads), "all<=C", all(l<=C for l in loads))
print("IP time", round(ip_time,2))
