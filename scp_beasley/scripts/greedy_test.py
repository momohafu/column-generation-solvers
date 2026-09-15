import time
p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)
def greedy(selected, score='ratio'):
    sel=set(selected); covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    steps=0
    while not all(covered) and steps<1000:
        best=None
        for j in range(n):
            if j in sel: continue
            new=0
            for i in colrows[j]:
                if not covered[i]: new+=1
            if new==0: continue
            if score=='ratio': key=costs[j]/new
            elif score=='coverage': key=-new/costs[j]
            elif score=='cost': key=(costs[j], -new)
            if best is None or key<best[0]: best=(key,j,new)
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
        steps+=1
    return list(sel)
for score in ['ratio','coverage','cost']:
    sel=greedy([],score)
    print(score, sum(costs[j] for j in sel), len(sel))
