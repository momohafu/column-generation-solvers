import time, datetime
from ortools.math_opt.python import mathopt

p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
M=float(sum(costs)+1)
lp_params=mathopt.SolveParameters(enable_output=False)
mip_params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=10), enable_output=False)

def solve_sp(y):
    sp=mathopt.Model(name='sp')
    s=[sp.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f's{i}') for i in range(m)]
    x=[sp.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f'x{j}') for j in range(n)]
    row_cons=[]
    for i,row in enumerate(rows):
        row_cons.append(sp.add_linear_constraint(sum(x[j] for j in row) + s[i] >= 1.0, name=f'cov{i}'))
    xub=[sp.add_linear_constraint(x[j] <= float(y[j]), name=f'xub{j}') for j in range(n)]
    sp.minimize_linear_objective(sum(costs[j]*x[j] for j in range(n)) + M*sum(s[i] for i in range(m)))
    res=mathopt.solve(sp, mathopt.SolverType.GLOP, params=lp_params)
    pi=res.dual_values(row_cons)
    mu=res.dual_values(xub)
    if not isinstance(pi,list): pi=[pi[c] for c in row_cons]
    if not isinstance(mu,list): mu=[mu[c] for c in xub]
    alpha=sum(pi)
    lam=[-v for v in mu]  # standard nonnegative lambda
    obj=res.objective_value()
    return alpha, lam, obj

def solve_master(cuts):
    mp=mathopt.Model(name='master')
    y=[mp.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'y{j}') for j in range(n)]
    theta=mp.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='theta')
    for kk,(alpha,lam) in enumerate(cuts):
        mp.add_linear_constraint(theta + sum(lam[j]*y[j] for j in range(n)) >= alpha, name=f'bcut{kk}')
    mp.minimize_linear_objective(theta)
    res=mathopt.solve(mp, mathopt.SolverType.HIGHS, params=mip_params)
    yv=res.variable_values(y)
    ystar=[1 if yv[j]>0.5 else 0 for j in range(n)]
    th=res.variable_values([theta])[0]
    return ystar, th, res.termination.reason, res.objective_value()

t0=time.perf_counter()
current=[1]*n
cuts=[]
max_iter=30
for it in range(max_iter):
    alpha,lam,sp_obj=solve_sp(current)
    cuts.append((alpha,lam))
    ystar,th,mp_term,mp_obj=solve_master(cuts)
    print('it',it,'sp_obj',round(sp_obj,4),'alpha',round(alpha,4),'sum_lam',round(sum(lam),4),'mp_obj',round(mp_obj,4),'sel',sum(ystar),'mp_term',mp_term,'cuts',len(cuts))
    if ystar==current:
        print('same y, break'); break
    current=ystar
    if time.perf_counter()-t0>110:
        print('time limit'); break
print('wall', time.perf_counter()-t0, 'cuts', len(cuts), 'last_obj', round(mp_obj,4))
