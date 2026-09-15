p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
from collections import Counter
for mod in [2,5,10]:
    c=Counter()
    for l,d in items: c[l%mod]+=d
    print('mod',mod, dict(sorted(c.items())))
    total_mod=sum((l%mod)*d for l,d in items)%mod
    print('total sum mod',total_mod,'slack sum mod', (65)%mod)
