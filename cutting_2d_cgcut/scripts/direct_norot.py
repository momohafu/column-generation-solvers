import time
from ortools.sat.python import cp_model
m=7; W,H=15,10
items=[(8,4,2,66),(3,7,1,35),(8,2,3,24),(3,4,5,17),(3,3,2,11),(3,2,2,8),(2,1,1,2)]
model=cp_model.CpModel()
ivx=[]; ivy=[]; pres=[]; objs=[]
for i,(l,w,q,v) in enumerate(items):
    dw,dh=l,w
    for k in range(q):
        p=model.NewBoolVar(f"p_{i}_{k}")
        xv=model.NewIntVar(0, W-dw, f"x_{i}_{k}")
        yv=model.NewIntVar(0, H-dh, f"y_{i}_{k}")
        iv=model.NewOptionalIntervalVar(xv,dw,xv+dw,p,f"ix_{i}_{k}")
        jv=model.NewOptionalIntervalVar(yv,dh,yv+dh,p,f"iy_{i}_{k}")
        pres.append(p); ivx.append(iv); ivy.append(jv); objs.append(v*p)
    # per type cap
    start=len(pres)-q
    model.Add(sum(pres[start:]) <= q)
model.AddNoOverlap2D(ivx,ivy)
model.Maximize(sum(objs))
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120
solver.parameters.num_search_workers=8
t=time.time(); status=solver.Solve(model); dt=time.time()-t
print(solver.StatusName(status), solver.ObjectiveValue(), dt, solver.BestObjectiveBound())
n=[0]*m; idx=0
for i in range(m):
    n[i]=sum(1 for p in pres[idx:idx+items[i][2]] if solver.Value(p)); idx+=items[i][2]
print("counts",n)
