import time
from ortools.sat.python import cp_model

def parse(path):
    lines=open(path).read().strip().splitlines()
    lines=[l.strip() for l in lines if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split())
        items.append((a,b,q,v))
    return m,W,H,items

path="/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt"
m,W,H,items=parse(path)
print("m,W,H",m,W,H,items)

def solve_direct():
    model=cp_model.CpModel()
    x=[]; y=[]; pres=[]; ivx=[]; ivy=[]; objs=[]
    for i,(l,w,q,v) in enumerate(items):
        orients=sorted(set([(l,w),(w,l)]))
        for (dw,dh) in orients:
            for k in range(q):
                p=model.NewBoolVar(f"p_{i}_{dw}_{dh}_{k}")
                xv=model.NewIntVar(0, W-dw, f"x_{i}_{dw}_{dh}_{k}")
                yv=model.NewIntVar(0, H-dh, f"y_{i}_{dw}_{dh}_{k}")
                iv=model.NewOptionalIntervalVar(xv, dw, xv+dw, p, f"ix_{i}_{dw}_{dh}_{k}")
                jv=model.NewOptionalIntervalVar(yv, dh, yv+dh, p, f"iy_{i}_{dw}_{dh}_{k}")
                pres.append(p); x.append(xv); y.append(yv); ivx.append(iv); ivy.append(jv); objs.append(v*p)
        model.Add(sum([p for (p,(_d,_h)) in []]) >=0)
    # constraints per type sum of pres
    idx=0
    for i,(l,w,q,v) in enumerate(items):
        orients=sorted(set([(l,w),(w,l)]))
        cnt=len(orients)*q
        model.Add(sum(pres[idx:idx+cnt]) <= q)
        idx+=cnt
    model.AddNoOverlap2D(ivx, ivy)
    model.Maximize(sum(objs))
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=120
    solver.parameters.num_search_workers=8
    solver.parameters.log_search_progress=False
    t=time.time()
    status=solver.Solve(model)
    dt=time.time()-t
    print("direct status",solver.StatusName(status),"obj",solver.ObjectiveValue(),"time",dt,"bound",solver.BestObjectiveBound())
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        n=[0]*m
        idx=0
        for i,(l,w,q,v) in enumerate(items):
            orients=sorted(set([(l,w),(w,l)]))
            cnt=len(orients)*q
            n[i]=sum(1 for p in pres[idx:idx+cnt] if solver.Value(p))
            idx+=cnt
        print("counts",n)

solve_direct()
