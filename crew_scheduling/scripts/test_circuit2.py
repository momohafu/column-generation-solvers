from ortools.sat.python import cp_model
m=cp_model.CpModel()
arcs=[]
for i,j in [(0,1),(1,2),(2,3),(3,0)]:
    arcs.append((i,j,m.NewBoolVar(f'a{i}_{j}')))
m.AddCircuit(arcs)
s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=5
st=s.Solve(m)
print('single', st, s.StatusName(st))
if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
    for i,j,b in arcs: print(' ',i,j,s.Value(b))
