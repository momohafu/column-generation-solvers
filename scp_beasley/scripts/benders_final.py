import time, datetime
from ortools.math_opt.python import mathopt

DATA='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(DATA).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)

M=float(sum(costs)+1)
lp_params=mathopt.SolveParameters(enable_output=False)
mip_params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False)

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
    lam=[-v for v in mu]  # MathOpt <= dual is nonpositive; convert to standard nonnegative lambda
    return alpha, lam, res.objective_value()

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
    return ystar, th, res.objective_value(), res.termination.reason

t0=time.perf_counter()
current=[1]*n
cuts=[]
max_iter=20
trace=[]
for it in range(max_iter):
    alpha,lam,sp_obj=solve_sp(current)
    cuts.append((alpha,lam))
    ystar,theta,mp_obj,mp_term=solve_master(cuts)
    trace.append((it+1, sp_obj, mp_obj, theta, sum(ystar), len(cuts), str(mp_term)))
    print(f'iter {it+1}: SP_obj={sp_obj:.4f}, MP_obj={mp_obj:.4f}, theta={theta:.4f}, selected={sum(ystar)}, cuts={len(cuts)}, MP_term={mp_term}')
    if ystar==current:
        print('master solution unchanged -> convergence')
        break
    current=ystar
    if time.perf_counter()-t0>110:
        print('time limit reached')
        break
wall=time.perf_counter()-t0
print('Benders wall', round(wall,3), 'iterations', len(trace), 'cuts', len(cuts), 'Benders lower bound', mp_obj)

# integer repair on full pool
t1=time.perf_counter()
mip=mathopt.Model(name='scp41_benders_repair')
x=[mip.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'x{j}') for j in range(n)]
mip.minimize_linear_objective(sum(costs[j]*x[j] for j in range(n)))
for i,row in enumerate(rows):
    mip.add_linear_constraint(sum(x[j] for j in row) >= 1.0, name=f'cov{i}')
mres=mathopt.solve(mip, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
mv=mres.variable_values(x)
sel=[j for j in range(n) if mv[j]>0.5]
print('integer repair', mres.termination.reason, mres.objective_value(), 'best_bound', mres.best_objective_bound(), 'cols', len(sel), 'wall', round(time.perf_counter()-t1,3))
print('selected_mip', sorted(sel))
