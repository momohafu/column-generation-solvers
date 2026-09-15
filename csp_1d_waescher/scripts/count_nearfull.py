p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
# DP count patterns with load in [L-65, L]
poly=[0]*(L+1); poly[0]=1
for l,d in items:
    new=poly[:]
    for w in range(1,d+1):
        wl=w*l
        if wl>L: break
        for cap in range(L, wl-1, -1):
            new[cap]+=poly[cap-wl]
    poly=new
cnt=sum(poly[L-65:])
print('near-full patterns count',cnt)
# count patterns with load <=65 too? no
# maybe number of patterns with load >= 9000 etc
for th in [9000,9500,9900,9935]:
    print(th, sum(poly[th:]))
