# min path cover ignoring span
p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,T=map(int,lines[0].split())
arcs=[]
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs.append((i,j))
from collections import deque
# Hopcroft-Karp max matching on bipartite left i -> right j for arcs
adj=[[] for _ in range(N+1)]
for i,j in arcs: adj[i].append(j)
pairU=[0]*(N+1); pairV=[0]*(N+1); dist=[0]*(N+1)
INF=10**9
def bfs():
    q=deque()
    for u in range(1,N+1):
        if pairU[u]==0:
            dist[u]=0; q.append(u)
        else: dist[u]=INF
    dist[0]=INF
    while q:
        u=q.popleft()
        if dist[u]<dist[0]:
            for v in adj[u]:
                if dist[pairV[v]]==INF:
                    dist[pairV[v]]=dist[u]+1; q.append(pairV[v])
    return dist[0]!=INF
def dfs(u):
    if u:
        for v in adj[u]:
            if dist[pairV[v]]==dist[u]+1 and dfs(pairV[v]):
                pairV[v]=u; pairU[u]=v; return True
        dist[u]=INF; return False
    return True
matching=0
while bfs():
    for u in range(1,N+1):
        if pairU[u]==0 and dfs(u): matching+=1
print('N',N,'arcs',len(arcs),'max_matching',matching,'min_path_cover(no span)',N-matching)
