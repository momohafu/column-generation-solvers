import random, time, math
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
exp=[]
for i,(l,d) in enumerate(items): exp += [l]*d
total=sum(exp)
best=10**9; bestbins=None; seed=0
t0=time.time()
def ff(order):
    bins=[]; assigns=[]
    for it in order:
        b=-1
        for i,load in enumerate(bins):
            if load+it<=L:
                b=i; break
        if b<0:
            bins.append(it); assigns.append([it])
        else:
            bins[b]+=it; assigns[b].append(it)
    return bins,assigns
for trial in range(20000):
    random.shuffle(exp)
    bins,assigns=ff(exp)
    if len(bins)<best:
        best=len(bins); bestbins=(bins[:], [a[:] for a in assigns]); seed=trial
        print('trial',trial,'bins',best,'time',time.time()-t0)
        if best==28: break
print('best',best)
if bestbins:
    for i,b in enumerate(bestbins[0]): print(i,b,bestbins[1][i])
