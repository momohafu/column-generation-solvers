import math, time, sys
sys.path.insert(0,"/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg import build_data, enumerate_paths
from ortools.sat.python import cp_model
n,x,y,dem,ready,due,svc,cap,depot_due,dist = build_data()
t0=time.time()
paths,costs,masks,loads = enumerate_paths(n,dist,dem,ready,due,svc,cap,depot_due)
P=len(paths)
print("pool",P,"enumerate",time.time()-t0)
model=cp_model.CpModel()
yv=[model.NewBoolVar(f"y_{p}") for p in range(P)]
for i in range(1,n+1):
    model.Add(sum(yv[p] for p in range(P) if (masks[p]>>i)&1)>=1)
model.Add(sum(yv)<=25)
# scaled objective to cents
model.Minimize(sum(int(round(costs[p]*100))*yv[p] for p in range(P)))
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=60
solver.parameters.num_search_workers=8
t0=time.time()
status=solver.Solve(model)
print("status",solver.StatusName(status),"obj",solver.ObjectiveValue()/100,"bound",solver.BestObjectiveBound()/100,"time",time.time()-t0)
# extract routes
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    sel=[p for p in range(P) if solver.Value(yv[p])]
    print("selected",len(sel), sel[:10])
    for p in sel:
        print(paths[p], costs[p])
