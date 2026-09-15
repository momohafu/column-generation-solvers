import math, time, sys
from ortools.sat.python import cp_model

DATA="/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
SCALE=10000
M=10**12

def parse(path=DATA):
    lines=open(path).read().splitlines()
    customers=[]
    for l in lines:
        s=l.split()
        if len(s)==7 and s[0].isdigit():
            no,x,y,dem,ready,due,svc=map(int,s)
            customers.append((no,x,y,dem,ready,due,svc))
    customers.sort(key=lambda c:c[0])
    depot=customers[0]
    return depot, customers[1:]

def dist(a,b):
    return math.hypot(a[1]-b[1], a[2]-b[2])

def build_solve(K, time_limit=120.0, print_info=True):
    depot, custs = parse()
    n=len(custs)
    xc=[depot[1]]+[c[1] for c in custs]
    yc=[depot[2]]+[c[2] for c in custs]
    demands=[0]+[c[3] for c in custs]
    ready=[0]+[c[4] for c in custs]
    due=[0]+[c[5] for c in custs]
    service=[0]+[c[6] for c in custs]
    cap=200
    due0=depot[5]
    d=[[0]*(n+1) for _ in range(n+1)]
    for i in range(n+1):
        for j in range(n+1):
            if i!=j:
                d[i][j]=round(math.hypot(xc[i]-xc[j], yc[i]-yc[j])*SCALE)
    model=cp_model.CpModel()
    x=[[[None]*(n+1) for _ in range(n+1)] for _ in range(K)]
    for k in range(K):
        for i in range(n+1):
            for j in range(n+1):
                x[k][i][j]=model.NewBoolVar(f"x_{k}_{i}_{j}")
    t=[[model.NewIntVar(0, due0*SCALE, f"t_{k}_{i}") for i in range(n+1)] for k in range(K)]
    q=[[model.NewIntVar(0, cap, f"q_{k}_{i}") for i in range(n+1)] for k in range(K)]
    for k in range(K):
        model.AddCircuit([(i, j, x[k][i][j]) for i in range(n+1) for j in range(n+1)])
    for i in range(1,n+1):
        model.Add(sum(x[k][i][i] for k in range(K))==K-1)
    for k in range(K-1):
        model.Add(x[k][0][0]<=x[k+1][0][0])
    for k in range(K):
        model.Add(t[k][0]==0)
        model.Add(q[k][0]==0)
        for i in range(1,n+1):
            model.Add(t[k][i]>=ready[i]*SCALE - M*x[k][i][i])
            model.Add(t[k][i]<=due[i]*SCALE + M*x[k][i][i])
            model.Add(t[k][i]<=M*(1-x[k][i][i]))
            model.Add(q[k][i]>=demands[i] - M*x[k][i][i])
            model.Add(q[k][i]<=M*(1-x[k][i][i]))
        for i in range(n+1):
            for j in range(1,n+1):
                if i==j: continue
                model.Add(t[k][j]>=t[k][i]+service[i]*SCALE+d[i][j]-M*(1-x[k][i][j]))
                model.Add(q[k][j]>=q[k][i]+demands[j]-M*(1-x[k][i][j]))
        for i in range(1,n+1):
            model.Add(t[k][i]+service[i]*SCALE+d[i][0]<=due0*SCALE+M*(1-x[k][i][0]))
    obj=[]
    for k in range(K):
        for i in range(n+1):
            for j in range(n+1):
                if i!=j:
                    obj.append(d[i][j]*x[k][i][j])
    model.Minimize(sum(obj))
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=time_limit
    solver.parameters.num_search_workers=8
    solver.parameters.log_search_progress=False
    t0=time.time()
    status=solver.Solve(model)
    elapsed=time.time()-t0
    obj_scaled = solver.ObjectiveValue() if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    bound_scaled = solver.BestObjectiveBound()
    if print_info:
        print("K",K,"status",solver.StatusName(status),
              "obj_scaled", obj_scaled, "obj2", (round(obj_scaled/SCALE,2) if obj_scaled is not None else None),
              "time",round(elapsed,2),"bound2",round(bound_scaled/SCALE,2))
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        routes=[]
        for k in range(K):
            route=[]
            cur=0
            while True:
                nxt=None
                for j in range(n+1):
                    if solver.Value(x[k][cur][j]):
                        nxt=j; break
                if nxt is None or nxt==0:
                    break
                route.append(nxt); cur=nxt
                if len(route)>n:
                    break
            if route:
                routes.append(route)
        # recompute exact full-precision cost
        def exact_cost(route):
            seq=[0]+route+[0]
            return sum(math.hypot(xc[seq[i]]-xc[seq[i+1]], yc[seq[i]]-yc[seq[i+1]]) for i in range(len(seq)-1))
        exact=sum(exact_cost(r) for r in routes)
        return status, obj_scaled/SCALE, routes, elapsed, bound_scaled/SCALE, exact
    return status, None, None, elapsed, bound_scaled/SCALE, None

if __name__=="__main__":
    for K in range(1,4):
        status,obj,routes,elapsed,bound,exact=build_solve(K,120)
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            print("routes",routes)
            print("exact total", round(exact,6), "rounded2", round(exact,2))
            break
