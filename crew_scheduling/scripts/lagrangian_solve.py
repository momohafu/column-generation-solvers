import datetime, time, math, random
from ortools.math_opt.python import mathopt

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
K=27; P=len(paths)
cover=[[] for _ in range(N+1)]
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
# upper bound via full pool MIP
m=mathopt.Model(name='ub')
z=[m.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f'z{p}') for p in range(P)]
for i in range(1,N+1):
    m.add_linear_constraint(sum(z[p] for p in cover[i]) == 1)
m.add_linear_constraint(sum(z) == K)
m.minimize(sum(paths[p]['cost']*z[p] for p in range(P)))
t0=time.time()
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
UB=res.objective_value(); ub_wall=time.time()-t0
print('UB',UB,'ub_wall',round(ub_wall,3),'status',res.termination.reason)
# lagrangian subgradient
lam=[0.0]*(N+1)
best_LB=-1e18; best_lam=None
alpha=2.0; last_improve=0
MAXIT=500
hist=[]
random.seed(42)
start=time.time()
for it in range(1,MAXIT+1):
    # reduced costs
    rc=[]
    for p in range(P):
        r=paths[p]['cost'] - sum(lam[i] for i in paths[p]['tasks'])
        rc.append(r)
    order=sorted(range(P), key=lambda p: rc[p])
    sel=set(order[:K])
    L=sum(lam[i] for i in range(1,N+1)) + sum(rc[p] for p in sel)
    if L>best_LB:
        best_LB=L; last_improve=it
    # subgradient
    g=[0.0]*(N+1)
    for i in range(1,N+1):
        cnt=sum(1 for p in cover[i] if p in sel)
        g[i]=1.0-cnt
    gnorm2=sum(g[i]*g[i] for i in range(1,N+1))
    if gnorm2<1e-12:
        # feasible (unlikely with K selection but possible)
        print('zero subgradient at',it,'L',L); break
    step=alpha*(UB-L)/gnorm2
    if step<0: step=0
    for i in range(1,N+1):
        lam[i]+=step*g[i]
    hist.append((it,L,math.sqrt(gnorm2),alpha,step))
    if it%50==0:
        print(f'it {it}: L {L:.4f} best {best_LB:.4f} gap {(UB-best_LB)/UB*100:.4f}% step {step:.4f} alpha {alpha}', flush=True)
    if it - last_improve >= 30:
        alpha*=0.5; last_improve=it
wall=time.time()-start
gap=(UB-best_LB)/UB*100
print('Lagrangian done wall',round(wall,3),'best_LB',best_LB,'UB',UB,'gap%',gap,'final_it',it)
# print some trajectory
print('first hist',hist[:3]); print('last hist',hist[-3:])
