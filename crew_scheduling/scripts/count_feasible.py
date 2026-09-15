p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,L=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs=set()
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs.add((i,j))
for cond, ok in [('j>i', lambda i,j: i<j and tasks[i][1]<=tasks[j][0]), ('all', lambda i,j: i!=j and tasks[i][1]<=tasks[j][0])]:
    feas=[(i,j) for i in range(1,N+1) for j in range(1,N+1) if ok(i,j)]
    listed=[(i,j) for (i,j) in feas if (i,j) in arcs]
    print(cond, 'feasible', len(feas), 'listed', len(listed), 'missing', len(feas)-len(listed))
print('listed arcs with j>i', sum(1 for i,j in arcs if i<j), 'with j<i', sum(1 for i,j in arcs if i>j))
