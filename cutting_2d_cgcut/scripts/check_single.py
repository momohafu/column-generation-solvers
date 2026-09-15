import sys, time
sys.setrecursionlimit(1000000)
from functools import lru_cache
m=7; W,H=15,10
items=[(8,4,2,66),(3,7,1,35),(8,2,3,24),(3,4,5,17),(3,3,2,11),(3,2,2,8),(2,1,1,2)]
calls=0
@lru_cache(maxsize=None)
def feas(w,h,cc):
    global calls; calls+=1
    if all(x==0 for x in cc): return True
    # single item base
    for i in range(m):
        if cc[i]==1 and all(cc[j]==0 for j in range(m) if j!=i):
            a,b,_,_=items[i]
            if a<=w and b<=h: return True
    # vertical splits
    for x in range(1,w):
        # enumerate c1 subcounts via recursion
        def gen(idx, cur):
            if idx==m:
                yield tuple(cur); return
            for val in range(cc[idx]+1):
                cur.append(val); yield from gen(idx+1,cur); cur.pop()
        for c1 in gen(0,[]):
            if all(v==0 for v in c1) or all(c1[j]==cc[j] for j in range(m)): continue
            c2=tuple(cc[j]-c1[j] for j in range(m))
            if feas(x,h,c1) and feas(w-x,h,c2): return True
    for y in range(1,h):
        def gen2(idx, cur):
            if idx==m:
                yield tuple(cur); return
            for val in range(cc[idx]+1):
                cur.append(val); yield from gen2(idx+1,cur); cur.pop()
        for c1 in gen2(0,[]):
            if all(v==0 for v in c1) or all(c1[j]==cc[j] for j in range(m)): continue
            c2=tuple(cc[j]-c1[j] for j in range(m))
            if feas(w,y,c1) and feas(w,h-y,c2): return True
    return False
cc=(2,1,1,2,1,1,0)
t=time.time(); f=feas(W,H,cc); print(f,time.time()-t,"calls",calls)
