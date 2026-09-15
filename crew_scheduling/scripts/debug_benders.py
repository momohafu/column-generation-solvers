import datetime
from ortools.math_opt.python import mathopt
# build same paths minimal
DATA='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(DATA).read().splitlines()
N,T=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs[(i,j)]=c
adj={i:[] for i in range(1,N+1)}
for (i,j) in arcs: adj[i].append(j)
paths=[]
for start in range(1,N+1):
    s0=tasks[start][0]
    stack=[(start,[start],0)]
    while stack:
        u,seq,c=stack.pop()
        paths.append({'tasks':list(seq),'cost':c,'span':tasks[u][1]-s0})
        for v in adj[u]:
            c2=c+arcs[(u,v)]
            if tasks[v][1]-s0<=T:
                stack.append((v,seq+[v],c2))
P=len(paths); K=27
cover=[[] for _ in range(N+1)]
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
# choose first 27 cheapest singleton paths
singletons=[p for p,o in enumerate(paths) if len(o['tasks'])==1]
sel=singletons[:K]
yhat=[0]*P
for p in sel: yhat[p]=1
m=mathopt.Model(name='sub')
x=[]; xcons=[]
for p in range(P):
    xv=m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'x{p}')
    x.append(xv); xcons.append(m.add_linear_constraint(xv <= yhat[p]))
sp=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sp{i}') for i in range(1,N+1)]
sm=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sm{i}') for i in range(1,N+1)]
eq=[]
for i in range(1,N+1):
    eq.append(m.add_linear_constraint(sum(x[p] for p in cover[i]) + sp[i-1] - sm[i-1] == 1))
m.minimize(sum(sp)+sum(sm))
res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
print('V',res.objective_value(),'term',res.termination.reason)
dv=res.dual_values()
pi=[dv[eq[i-1]] for i in range(1,N+1)]
alpha=[dv[xcons[p]] for p in range(P)]
print('sum pi',sum(pi),'pi minmax',min(pi),max(pi),'first pi',pi[:10])
print('alpha nonzero',sum(1 for a in alpha if abs(a)>1e-8),'sum alpha',sum(alpha),'first',alpha[:10])
# verify dual objective for selected y
print('dual obj computed',sum(pi)-sum(alpha[p]*yhat[p] for p in range(P)))
# Check dual feasibility for x constraints: sum_i a_ip pi_i - alpha_p should be <=0 for all p
viol=[]
for p in range(P):
    v=sum(pi[i-1] for i in paths[p]['tasks']) - alpha[p]
    if v>1e-6: viol.append((p,v,paths[p]['tasks']))
print('dual x constraint violations',len(viol),viol[:5])
# slack constraints |pi|<=1?
print('pi beyond +-1',sum(1 for pii in pi if abs(pii)>1+1e-6))
