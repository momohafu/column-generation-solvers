import math
cust=[]
for l in open("/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"):
    s=l.split()
    if len(s)==7 and s[0].isdigit():
        no,x,y,dem,ready,due,svc=map(int,s); cust.append((no,x,y,dem,ready,due,svc))
D={c[0]:(c[1],c[2]) for c in cust}
def d(a,b):
    x1,y1=D[a]; x2,y2=D[b]; return math.hypot(x1-x2,y1-y2)
routes=[[20,24,25,23,22,21],[13,17,18,19,15,16,14,12],[5,3,7,8,10,11,9,6,4,2,1]]
total_r=0; total_f=0; total_t=0
for r in routes:
    seq=[0]+r+[0]
    sr=sum(round(d(seq[i],seq[i+1]),2) for i in range(len(seq)-1))
    sf=sum(math.floor(d(seq[i],seq[i+1])*100+1e-9)/100 for i in range(len(seq)-1))
    st=sum(d(seq[i],seq[i+1]) for i in range(len(seq)-1))
    total_r+=sr; total_f+=sf; total_t+=st
    print(r, "round", round(sr,4), "floor", round(sf,4), "full", round(st,4))
    print([(round(d(seq[i],seq[i+1]),2), math.floor(d(seq[i],seq[i+1])*100+1e-9)/100) for i in range(len(seq)-1)])
print("TOTALS round", round(total_r,4), "floor", round(total_f,4), "full", round(total_t,4))
