import time, math
p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)

def ratio_greedy(selected):
    sel=set(selected); covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    while not all(covered):
        best=None
        for j in range(n):
            if j in sel: continue
            new=sum(1 for i in colrows[j] if not covered[i])
            if new==0: continue
            key=costs[j]/new
            if best is None or key<best[0]: best=(key,j,new)
        if best is None: return None
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
    return list(sel)

def rc_greedy(pi, selected):
    sel=set(selected); covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    while not all(covered):
        best=None
        for j in range(n):
            if j in sel: continue
            has=False; s=0.0
            for i in colrows[j]:
                if not covered[i]: has=True
                s+=pi[i]
            if not has: continue
            rc=costs[j]-s
            key=(rc, costs[j], -sum(1 for i in colrows[j] if not covered[i]))
            if best is None or key<best[0]: best=(key,j)
        if best is None: return None
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
    return list(sel)

UB=sum(costs[j] for j in ratio_greedy([]))
print('initial UB', UB)
pi=[0.0]*m; lam=2.0; best_L=-1e9; best_pi=None; best_x=None; no_impr=0; max_iter=600
t0=time.perf_counter()
for it in range(max_iter):
    rc=[0.0]*n
    for j in range(n):
        s=0.0
        for i in colrows[j]: s+=pi[i]
        rc[j]=costs[j]-s
    x=[1 if rc[j]<0 else 0 for j in range(n)]
    L=sum(pi)+sum(min(0.0, rc[j]) for j in range(n))
    if L>best_L+1e-9:
        best_L=L; best_pi=pi[:]; best_x=x[:]; no_impr=0
    else: no_impr+=1
    g=[1.0]*m
    for j in range(n):
        if x[j]:
            for i in colrows[j]: g[i]-=1.0
    norm2=sum(gi*gi for gi in g)
    if it%50==0 or it==max_iter-1:
        sel=rc_greedy(pi, [j for j in range(n) if x[j]])
        if sel is not None:
            c=sum(costs[j] for j in sel)
            if c<UB: UB=c
    step=0.0 if norm2<1e-12 else lam*(UB-L)/norm2
    for i in range(m): pi[i]=max(0.0, pi[i]+step*g[i])
    if no_impr>=60: lam/=2.0; no_impr=0
selb=rc_greedy(best_pi, [j for j in range(n) if best_x[j]])
cb=sum(costs[j] for j in selb) if selb else None
if cb and cb<UB: UB=cb
print('best_L', best_L, 'UB', UB, 'greedy_best', cb, 'wall', time.perf_counter()-t0)
