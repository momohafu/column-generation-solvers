import pickle, datetime, time
from ortools.math_opt.python import mathopt
data=pickle.load(open("/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts/pool.pkl","rb"))
m,W,H,items,P=data["m"],data["W"],data["H"],data["items"],data["patterns"]
q=[it[2] for it in items]; v=[it[3] for it in items]
print("pool",len(P))
model=mathopt.Model(name="pool_ip")
xv=[model.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{k}") for k in range(len(P))]
for i in range(m):
    model.add_linear_constraint(sum(P[k][i]*xv[k] for k in range(len(P))) <= q[i])
model.add_linear_constraint(sum(xv) <= 1)
model.maximize(sum(sum(P[k][j]*v[j] for j in range(m))*xv[k] for k in range(len(P))))
t=time.time()
res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120.0), enable_output=False))
dt=time.time()-t
print("status",res.termination.reason,"obj",res.objective_value(),"time",dt)
sol=[k for k in range(len(P)) if res.variable_values()[xv[k]]>0.5]
print("used",[P[k] for k in sol],"value",[sum(P[k][j]*v[j] for j in range(m)) for k in sol])
