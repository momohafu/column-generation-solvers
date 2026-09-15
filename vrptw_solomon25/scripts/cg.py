import math, time, datetime
from ortools.math_opt.python import mathopt

DATA="/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"

def parse(path=DATA):
    lines=open(path).read().splitlines()
    customers=[]
    for l in lines:
        s=l.split()
        if len(s)==7 and s[0].isdigit():
            no,x,y,dem,ready,due,svc=map(int,s)
            customers.append((no,x,y,dem,ready,due,svc))
    customers.sort(key=lambda c:c[0])
    return customers[0], customers[1:]

def build_data():
    depot, custs = parse()
    n=len(custs)
    x=[depot[1]]+[c[1] for c in custs]
    y=[depot[2]]+[c[2] for c in custs]
    dem=[0]+[c[3] for c in custs]
    ready=[0]+[c[4] for c in custs]
    due=[0]+[c[5] for c in custs]
    svc=[0]+[c[6] for c in custs]
    cap=200
    depot_due=depot[5]
    def dist(i,j):
        return math.hypot(x[i]-x[j], y[i]-y[j])
    return n,x,y,dem,ready,due,svc,cap,depot_due,dist

def enumerate_paths(n,dist,dem,ready,due,svc,cap,depot_due):
    paths=[]
    costs=[]
    masks=[]
    loads=[]
    seq=[]
    def dfs(cur, t, load, mask, cost):
        for j in range(1,n+1):
            if mask & (1<<j):
                continue
            arr = t + svc[cur] + dist(cur,j)
            if arr < ready[j]:
                arr = ready[j]
            if arr > due[j]:
                continue
            if load + dem[j] > cap:
                continue
            if arr + svc[j] + dist(j,0) > depot_due + 1e-9:
                continue
            new_mask = mask | (1<<j)
            new_cost = cost + dist(cur,j)
            seq.append(j)
            paths.append(tuple(seq))
            costs.append(new_cost + dist(j,0))
            masks.append(new_mask)
            loads.append(load+dem[j])
            dfs(j, arr, load+dem[j], new_mask, new_cost)
            seq.pop()
    dfs(0, 0.0, 0, 0, 0.0)
    return paths, costs, masks, loads

def column_generation(verbose=True):
    n,x,y,dem,ready,due,svc,cap,depot_due,dist = build_data()
    t0=time.time()
    paths,costs,masks,loads = enumerate_paths(n,dist,dem,ready,due,svc,cap,depot_due)
    P=len(paths)
    if verbose: print("enumerated paths", P, "time", round(time.time()-t0,2))
    K=25
    # initial columns: single-customer paths (should exist)
    init=[]
    for p in range(P):
        if len(paths[p])==1:
            init.append(p)
    selected=list(init)
    sel_set=set(selected)
    iter_count=0
    lp_obj=None
    pricing_times=[]
    t_start=time.time()
    while iter_count < 200 and time.time()-t_start < 100:
        iter_count+=1
        # build RMP
        m=mathopt.Model(name=f"RMP_{iter_count}")
        vs=[m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x_{p}") for p in selected]
        cover_con=[]
        for i in range(1,n+1):
            expr=mathopt.fast_sum([vs[k] for k,p in enumerate(selected) if (masks[p]>>i)&1])
            cover_con.append(m.add_linear_constraint(expr>=1.0, name=f"cover_{i}"))
        veh_con=m.add_linear_constraint(mathopt.fast_sum(vs)<=K, name="vehicles")
        m.minimize(mathopt.fast_sum([costs[p]*vs[k] for k,p in enumerate(selected)]))
        res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(enable_output=False))
        if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
            print("RMP not optimal", res.termination.reason)
            break
        lp_obj=res.objective_value()
        duals=res.dual_values()
        pi=[0.0]*(n+1)
        for i in range(1,n+1):
            pi[i]=duals[cover_con[i-1]]
        mu=duals[veh_con]
        # pricing: min reduced cost over all paths
        tp=time.time()
        best_rc=float("inf")
        candidates=[]
        for p in range(P):
            if p in sel_set:
                continue
            # rc = cost - sum pi over customers - mu
            msk=masks[p]
            s=0.0
            # sum pi via path tuple is slower; use bitmask loop
            mm=msk
            while mm:
                lb=mm & -mm
                i=lb.bit_length()-1
                s+=pi[i]
                mm-=lb
            rc=costs[p]-s-mu
            if rc < -1e-8:
                candidates.append((rc,p))
        pricing_times.append(time.time()-tp)
        if not candidates:
            if verbose: print("CG converged at iter",iter_count,"lp_obj",round(lp_obj,6),"sel",len(selected))
            break
        # add up to 2000 most negative
        candidates.sort()
        add=[p for rc,p in candidates[:2000]]
        for p in add:
            sel_set.add(p)
            selected.append(p)
        if verbose and (iter_count<=5 or iter_count%10==0):
            print("iter",iter_count,"lp_obj",round(lp_obj,6),"cols",len(selected),"best_rc",round(candidates[0][0],6),"add",len(add))
    total_time=time.time()-t0
    return dict(paths=paths,costs=costs,masks=masks,loads=loads,n=n,K=K,selected=selected,sel_set=sel_set,lp_obj=lp_obj,iterations=iter_count,pricing_times=pricing_times,total_time=total_time,dist=dist,x=x,y=y)

def integer_recovery(info, time_limit=120.0, full_pool=False, verbose=True):
    n=info["n"]; K=info["K"]
    if full_pool:
        pool=list(range(len(info["paths"])))
    else:
        pool=info["selected"]
    m=mathopt.Model(name="IP_recovery")
    vs=[m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y_{p}") for p in pool]
    for i in range(1,n+1):
        m.add_linear_constraint(mathopt.fast_sum([vs[k] for k,p in enumerate(pool) if (info["masks"][p]>>i)&1])>=1.0, name=f"cover_{i}")
    m.add_linear_constraint(mathopt.fast_sum(vs)<=K, name="vehicles")
    m.minimize(mathopt.fast_sum([info["costs"][p]*vs[k] for k,p in enumerate(pool)]))
    t0=time.time()
    res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
    elapsed=time.time()-t0
    obj=res.objective_value() if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE) else None
    bound=res.best_dual_bound() if hasattr(res,"best_dual_bound") else None
    if verbose:
        print("IP reason",res.termination.reason,"obj",round(obj,6) if obj is not None else None,"time",round(elapsed,2),"bound",round(bound,6) if bound is not None else None)
    routes=[]
    if obj is not None:
        vals=res.variable_values()
        for k,p in enumerate(pool):
            if vals[vs[k]]>0.5:
                routes.append(info["paths"][p])
    return dict(obj=obj, bound=bound, time=elapsed, routes=routes, termination=str(res.termination.reason), ncols=len(pool))

if __name__=="__main__":
    info=column_generation(verbose=True)
    print("FINAL lp_obj", info["lp_obj"], "iterations", info["iterations"], "selected", len(info["selected"]), "total_time", round(info["total_time"],2))
    ip=integer_recovery(info, time_limit=120, full_pool=False)
    print("IP routes", ip["routes"])
    # exact cost of routes
    tot=0
    for r in ip["routes"]:
        seq=[0]+list(r)+[0]
        c=sum(info["dist"](seq[i],seq[i+1]) for i in range(len(seq)-1))
        print(r, round(c,6))
        tot+=c
    print("exact total", round(tot,6), "rounded2", round(tot,2))
