p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,T=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
adj=[[] for _ in range(N+1)]; costs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); adj[i].append(j); costs[(i,j)]=c
paths=[]
for start in range(1,N+1):
    s0=tasks[start][0]
    stack=[(start,[start],0)]
    while stack:
        u,seq,c=stack.pop()
        paths.append((tuple(seq), c, tasks[u][1]-s0))
        for v in adj[u]:
            c2=c+costs[(u,v)]
            span=tasks[v][1]-s0
            if span<=T:
                stack.append((v,seq+[v],c2))
print('N',N,'T',T,'arcs',sum(len(a) for a in adj))
print('num paths',len(paths))
print('max len',max(len(s) for s,_,_ in paths),'min len',min(len(s) for s,_,_ in paths))
print('min cost',min(c for _,c,_ in paths),'max cost',max(c for _,c,_ in paths))
from collections import Counter
cnt=Counter(len(s) for s,_,_ in paths)
print('length dist',sorted(cnt.items()))
import json
with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/paths.json','w') as f:
    for s,c,span in paths:
        f.write(json.dumps({'tasks':list(s),'cost':c,'span':span})+chr(10))
print('saved paths.json')
