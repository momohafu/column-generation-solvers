p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
poly=[0]*(L+1); poly[0]=1
for l,d in items:
    new=poly[:]
    for w in range(1,d+1):
        wl=w*l
        if wl>L: break
        for cap in range(L, wl-1, -1):
            new[cap]+=poly[cap-wl]
    poly=new
count=sum(poly)
print('pattern count <=L',count)
print('top coeffs', max(poly), 'at', poly.index(max(poly)))
# count patterns with weight >= something? no
