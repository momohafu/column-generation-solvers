import math, time, numpy as np, sys, pickle
sys.path.insert(0,"/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg import build_data, enumerate_paths

n,x,y,dem,ready,due,svc,cap,depot_due,dist = build_data()
t0=time.time()
paths,costs,masks,loads = enumerate_paths(n,dist,dem,ready,due,svc,cap,depot_due)
print("pool", len(paths), "enumerate time", round(time.time()-t0,2))
# min cost per mask (mask bits are 1..n; shift to 0..n-1 for SOS)
mincost={}
bestpath={}
for p,m,c in zip(paths,masks,costs):
    mm=m>>1
    if mm not in mincost or c < mincost[mm]:
        mincost[mm]=c
        bestpath[mm]=p
print("distinct masks", len(mincost))
N=1<<n
f=np.full(N, np.inf, dtype=np.float32)
t1=time.time()
for mm,c in mincost.items():
    f[mm]=np.float32(c)
print("fill time", round(time.time()-t1,2), "f MB", round(f.nbytes/1e6,1))
t1=time.time()
for b in range(n):
    step=1<<(b+1)
    half=1<<b
    f2=f.reshape(-1, step)
    f2[:, :half] = np.minimum(f2[:, :half], f2[:, half:])
print("sos time", round(time.time()-t1,2))
kept=[]
t1=time.time()
for mm,c in mincost.items():
    if f[mm] < c - 1e-6:
        continue
    kept.append((mm<<1, c, bestpath[mm]))
print("kept masks", len(kept), "filter check time", round(time.time()-t1,2))
print("total time", round(time.time()-t0,2))
opt_routes=[(13,17,18,19,15,16,14,12),(20,24,25,23,22,21),(5,3,7,8,10,11,9,6,4,2,1)]
for r in opt_routes:
    m=0
    for i in r:
        m |= (1<<i)
    print("opt mask", m, "in kept?", any(k[0]==m for k in kept), "mincost", mincost.get(m>>1))
with open("/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts/kept_pool.pkl","wb") as fh:
    pickle.dump(kept, fh)
