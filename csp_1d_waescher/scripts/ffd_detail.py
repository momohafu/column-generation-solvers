p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
length=[l for l,d in items]
exp=[]
for i,(l,d) in enumerate(items): exp += [l]*d
exp.sort(reverse=True)
bins=[]; assign=[]
for it in exp:
    best=None
    for i,b in enumerate(bins):
        if b+it<=L and (best is None or b>bins[best]): best=i
    if best is None: bins.append(it); assign.append([it])
    else: bins[best]+=it; assign[best].append(it)
print('FFD bins',len(bins))
for k,b in enumerate(bins):
    print(k,'load',b,'slack',L-b,'contents',assign[k])
