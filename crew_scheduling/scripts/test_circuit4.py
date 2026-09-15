from ortools.sat.python import cp_model
for arcs in [[(0,1),(1,2),(2,0),(3,3)],[(0,1),(1,0),(2,2),(3,3)]]:
    m=cp_model.CpModel()
    a=[]
    for i,j in arcs:
        a.append((i,j,m.NewBoolVar(f'a{i}_{j}')))
    m.AddCircuit(a)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=5
    st=s.Solve(m)
    print('AddCircuit arcs',arcs,'->',s.StatusName(st))
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        for i,j,b in a: print(' ',i,j,s.Value(b))
