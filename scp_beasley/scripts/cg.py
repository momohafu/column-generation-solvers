import time, datetime
from ortools.math_opt.python import mathopt

p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split()
m,n=map(int,toks[:2])
costs=list(map(int,toks[2:2+n]))
idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1
    rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)

M=sum(costs)+1
model=mathopt.Model(name='scp41_rmp')
a=[model.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f'a{i}') for i in range(m)]
row_cons=[model.add_linear_constraint(a[i] >= 1.0, name=f'cov{i}') for i in range(m)]
xvar=[None]*n
added=set()
obj_terms=[(a[i], float(M)) for i in range(m)]
model.minimize_linear_objective(sum(coef*var for var,coef in obj_terms))
params=mathopt.SolveParameters(enable_output=False)
tol=1e-7
max_iter=2000
t0=time.perf_counter()
iters=0
lp_obj=None
hist=[]
for it in range(max_iter):
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=params)
    lp_obj=res.objective_value()
    pi=res.dual_values(row_cons)
    if not isinstance(pi, list):
        pi=[pi[c] for c in row_cons]
    # precompute row dual prefix? Just scan columns
    min_rc=float('inf'); best=-1
    for j in range(n):
        if j in added: continue
        rc=costs[j]
        for i in colrows[j]:
            rc -= pi[i]
        if rc < min_rc:
            min_rc=rc; best=j
    hist.append((it, lp_obj, min_rc, best, len(added)))
    if min_rc >= -tol:
        iters=it+1
        break
    v=model.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f'x{best}')
    xvar[best]=v
    for i in colrows[best]:
        row_cons[i].set_coefficient(v, 1.0)
    added.add(best)
    obj_terms.append((v, float(costs[best])))
    model.minimize_linear_objective(sum(coef*var for var,coef in obj_terms))
    iters=it+1
    if time.perf_counter()-t0 > 110:
        print('time limit in CG loop')
        break

cg_wall=time.perf_counter()-t0
print('CG iters', iters, 'added cols', len(added), 'lp_obj', lp_obj, 'min_rc_last', min_rc, 'wall', round(cg_wall,3))
# final RMP solution
res=mathopt.solve(model, mathopt.SolverType.GLOP, params=params)
vals=res.variable_values()
a_vals=[vals[a[i]] for i in range(m)]
print('artificial positive', sum(v>1e-7 for v in a_vals), 'max', max(a_vals))
print('first hist', hist[:5])
print('last hist', hist[-5:])

# integer repair on full pool
t1=time.perf_counter()
mip=mathopt.Model(name='scp41_repair')
x=[mip.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'x{j}') for j in range(n)]
mip.minimize_linear_objective(sum(costs[j]*x[j] for j in range(n)))
for i,row in enumerate(rows):
    mip.add_linear_constraint(sum(x[j] for j in row) >= 1.0, name=f'cov{i}')
mres=mathopt.solve(mip, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
mv=mres.variable_values(x)
sel=[j for j in range(n) if mv[j]>0.5]
rep_wall=time.perf_counter()-t1
print('repair termination', mres.termination.reason, 'obj', mres.objective_value(), 'best_bound', mres.best_objective_bound(), 'wall', round(rep_wall,3), 'sel', len(sel))
print('selected', sorted(sel))
