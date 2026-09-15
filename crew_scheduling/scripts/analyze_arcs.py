p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,L=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs[(i,j)]=c
print('N',N)
for i in range(1,N+1):
    out=[j for j in range(1,N+1) if j!=i and tasks[i][1]<=tasks[j][0]]
    out_listed=[j for j in out if (i,j) in arcs]
    if len(out)!=len(out_listed):
        print(f'task {i} out feasible {len(out)} listed {len(out_listed)} missing {[j for j in out if j not in out_listed][:10]}')
# incoming missing
for j in range(1,N+1):
    inc=[i for i in range(1,N+1) if i!=j and tasks[i][1]<=tasks[j][0]]
    inc_listed=[i for i in inc if (i,j) in arcs]
    if len(inc)!=len(inc_listed):
        print(f'task {j} IN feasible {len(inc)} listed {len(inc_listed)}')
