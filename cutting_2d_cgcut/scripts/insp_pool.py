# -*- coding: utf-8 -*-
import pickle, time, datetime
from ortools.math_opt.python import mathopt

D = "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts/"
data = pickle.load(open(D + "pool.pkl", "rb"))
items = data["items"]; m = data["m"]; W = data["W"]; H = data["H"]
patterns = data["patterns"]
v = [it[3] for it in items]
q = [it[2] for it in items]
print(f"m={m} 板材 {W}x{H} | 模式数 {len(patterns)} | 文献最优 244")

mm = mathopt.Model()
x = [mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{p}") for p in range(len(patterns))]
mm.maximize(mathopt.fast_sum([sum(v[i]*patterns[p][i] for i in range(m))*x[p] for p in range(len(patterns))]))
mm.add_linear_constraint(mathopt.fast_sum(x) == 1.0, name="one")
t0 = time.time()
res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                    params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
print(f"池 IP: {res.termination.reason.name} | 最优价值 {res.objective_value()} | {round(time.time()-t0,2)}s")
if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
    sel = [p for p in range(len(patterns)) if res.variable_values()[x[p]] > 0.5]
    if sel:
        pat = patterns[sel[0]]
        print("最优模式:", pat, "| 价值:", sum(v[i]*pat[i] for i in range(m)))
        print("件数 vs 上限:", [(pat[i], q[i]) for i in range(m)])
# LP 松弛
mm2 = mathopt.Model()
x2 = [mm2.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(len(patterns))]
mm2.maximize(mathopt.fast_sum([sum(v[i]*patterns[p][i] for i in range(m))*x2[p] for p in range(len(patterns))]))
mm2.add_linear_constraint(mathopt.fast_sum(x2) == 1.0, name="one")
res2 = mathopt.solve(mm2, mathopt.SolverType.GLOP)
print(f"池 LP(conv): {round(res2.objective_value(), 4)}")
