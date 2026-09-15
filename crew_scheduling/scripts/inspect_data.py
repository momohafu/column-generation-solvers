
import sys
p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,L=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs[(i,j)]=c
print('N',N,'L',L,'narcs',len(arcs))
print('task time min',min(s for s,f in tasks[1:]),'max',max(f for s,f in tasks[1:]))
# overlap max using sweep events half-open [start, finish)
ev=[]
for i,(s,f) in enumerate(tasks[1:],1):
    ev.append((s,1)); ev.append((f,-1))
ev.sort(key=lambda x:(x[0], -x[1]))
cur=0; mx=0; tmx=None
for t,d in ev:
    cur+=d
    if cur>mx: mx=cur; tmx=t
print('max overlap',mx,'at',tmx)
# total duration lower bound
total=sum((f-s) for s,f in tasks[1:])
print('total duration',total,'ceil/L', -(-total//L))
# check arcs i->j: all should satisfy finish_i <= start_j? and count
ok=0; bad=[]
for (i,j),c in arcs.items():
    if tasks[i][1] <= tasks[j][0]: ok+=1
    else: bad.append((i,j,tasks[i],tasks[j],c))
print('arcs feasible',ok,'infeasible',len(bad),'sample bad',bad[:5])
# how many ordered feasible pairs total?
feas=sum(1 for i in range(1,N+1) for j in range(1,N+1) if i!=j and tasks[i][1]<=tasks[j][0])
print('total ordered feasible pairs',feas)
# min/max arcs cost
costs=list(arcs.values())
print('cost min/max/sum',min(costs),max(costs),sum(costs))
