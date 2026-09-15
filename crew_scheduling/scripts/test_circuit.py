from ortools.sat.python import cp_model
for name in ['AddCircuit','AddMultipleCircuit']:
    m=cp_model.CpModel()
    arcs=[]
    for i,j in [(0,1),(1,0),(2,3),(3,2)]:
        arcs.append((i,j,m.NewBoolVar(f'a{i}_{j}')))
    getattr(m,name)(arcs)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=5
    st=s.Solve(m)
    print(name, st, s.StatusName(st))
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        for i,j,b in arcs:
            print(' ',i,j,s.Value(b))
