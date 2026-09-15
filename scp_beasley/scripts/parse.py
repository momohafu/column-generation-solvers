import sys
p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split()
print('total tokens', len(toks))
m,n=map(int,toks[:2]); print('m,n',m,n)
costs=list(map(int,toks[2:2+n]))
print('costs len',len(costs),'min max',min(costs),max(costs),'sum',sum(costs))
idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1
    cols=[int(t)-1 for t in toks[idx:idx+k]]
    idx+=k
    rows.append(cols)
print('parsed idx',idx,'remaining',len(toks)-idx)
print('rows len',len(rows))
print('row sizes min max', min(map(len,rows)), max(map(len,rows)))
print('first row', rows[0][:20], 'len', len(rows[0]))
colcounts=[0]*n
for r in rows:
    for j in r: colcounts[j]+=1
print('cols with coverage0', sum(c==0 for c in colcounts),'max cov', max(colcounts))
print('first costs', costs[:20])
