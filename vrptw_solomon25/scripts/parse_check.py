import math
p="/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
lines=open(p).read().splitlines()
cust=[]
start=False
for l in lines:
    s=l.split()
    if len(s)==7 and s[0].isdigit():
        no,x,y,dem,ready,due,svc=map(int,s)
        cust.append((no,x,y,dem,ready,due,svc))
print("ncustomers",len(cust)-1)
depot=cust[0]
def dist(a,b):
    return round(math.hypot(a[1]-b[1], a[2]-b[2]),2)
D=sorted((dist(depot,c),c[0]) for c in cust[1:])
print(D)
# distance matrix
n=len(cust)
import itertools
for i in range(n):
    for j in range(n):
        if i!=j:
            pass
print("max depot",D[-1])
