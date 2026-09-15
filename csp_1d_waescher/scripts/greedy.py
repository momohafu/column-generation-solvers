p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
# expand items sorted descending for FFD
exp=[]
for l,d in items:
    exp += [l]*d
exp.sort(reverse=True)
bins=[]
for it in exp:
    best=None
    for i,b in enumerate(bins):
        if b+it<=L and (best is None or b>bins[best]):
            best=i
    if best is None:
        bins.append(it)
    else:
        bins[best]+=it
print('FFD bins',len(bins),'loads',bins)
# best fit
bins=[]
for it in exp:
    best=None
    for i,b in enumerate(bins):
        if b+it<=L and (best is None or b+it>bins[best]+0):
            best=i
    if best is None: bins.append(it)
    else: bins[best]+=it
print('BFD bins',len(bins),'loads',bins)
