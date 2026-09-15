import math
cust=[]
for l in open("/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"):
    s=l.split()
    if len(s)==7 and s[0].isdigit():
        no,x,y,dem,ready,due,svc=map(int,s); cust.append((no,x,y,dem,ready,due,svc))
D={c[0]:(c[1],c[2]) for c in cust}
def dd(a,b):
    x1,y1=D[a]; x2,y2=D[b]; return math.hypot(x1-x2,y1-y2)
routes=[[20,24,25,23,22,21],[13,17,18,19,15,16,14,12],[5,3,7,8,10,11,9,6,4,2,1]]
tot=0
for r in routes:
    seq=[0]+r+[0]
    s=0
    for i in range(len(seq)-1):
        v=dd(seq[i],seq[i+1])
        r2=round(v,2)
        ri=round(v*100)
        r_int=round(r2*100)
        print(seq[i],seq[i+1], repr(v), "round2",r2, "round2*100 int",int(r2*100), "round(v*100)",ri)
        s+=ri
    tot+=s
    print("route int sum",s, "cents")
print("total cents",tot, "value",tot/100)
