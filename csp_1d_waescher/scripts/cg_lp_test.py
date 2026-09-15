import time, math
from ortools.math_opt.python import mathopt
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
length=[l for l,d in items]; dem=[d for l,d in items]
patterns=[]
for i,(l,d) in enumerate(items):
    a=[0]*m; a[i]=min(d, L//l); patterns.append(a)
items01=[]
for i,(l,d) in enumerate(items):
    for c in range(d):
        items01.append((l,i))
n=len(items01)
def price(pi):
    dp=[-1e300]*(L+1); dp[0]=0.0
    keep=[bytearray(L+1) for _ in range(n)]
    for idx,(w,typ) in enumerate(items01):
        val=pi[typ]; k=keep[idx]
        for cap in range(L, w-1, -1):
            if dp[cap-w] > -1e299:
                nv=dp[cap-w]+val
                if nv > dp[cap]:
                    dp[cap]=nv; k[cap]=1
    best_w=max(range(L+1), key=lambda c: dp[c])
    best_val=dp[best_w]
    a=[0]*m; cap=best_w
    for idx in range(n-1,-1,-1):
        if keep[idx][cap]:
            w,typ=items01[idx]; a[typ]+=1; cap-=w
    return best_val,a,best_w
def solve_rmp(patterns):
    model=mathopt.Model(name='rmp')
    x=[model.add_variable(lb=0.0, is_integer=False, name=f'x_{j}') for j in range(len(patterns))]
    cons=[]
    for i in range(m):
        cons.append(model.add_linear_constraint(sum(patterns[j][i]*x[j] for j in range(len(patterns))) >= dem[i]))
    model.minimize(sum(x))
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(enable_output=False))
    if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
        return None,None
    duals=res.dual_values()
    pi=[duals[cons[i]] for i in range(m)]
    return res.objective_value(), pi
tol=1e-7
it=0; t0=time.time(); obj=None; pi=None
while True:
    obj,pi=solve_rmp(patterns)
    if obj is None:
        print('RMP fail'); break
    best_val,a,best_w=price(pi)
    rc=1.0-best_val
    it+=1
    if it<=10 or rc < -tol or it%100==0:
        print('it',it,'obj',obj,'rc',rc,'cols',len(patterns))
    if rc >= -tol:
        print('LP CONVERGED obj %.10f iters %d time %.2f'%(obj,it,time.time()-t0))
        break
    patterns.append(a)
    if it>500:
        print('max iter'); break
