import sys
sys.setrecursionlimit(1000000)
m=7; W,H=15,10
items=[(8,4,2,66),(3,7,1,35),(8,2,3,24),(3,4,5,17),(3,3,2,11),(3,2,2,8),(2,1,1,2)]
q=[it[2] for it in items]
orients=[]
for (a,b,qq,vv) in items:
    orients.append(sorted(set([(a,b),(b,a)])))
c=(2,1,3,0,1,1,1)
print("pattern area",sum(c[i]*(items[i][0]*items[i][1]) for i in range(m)),"sheet",W*H)
from functools import lru_cache
# return True and a description, with memo
calls=0
@lru_cache(maxsize=None)
def feas(w,h,cc):
    global calls
    calls+=1
    if all(x==0 for x in cc):
        return True, "waste"
    # base: one item and rest waste
    for i in range(m):
        if cc[i]==1 and all(cc[j]==0 for j in range(m) if j!=i):
            for (dw,dh) in orients[i]:
                if dw<=w and dh<=h:
                    return True, f"item{i}({dw}x{dh}) in {w}x{h}"
    # splits
    for x in range(1,w):
        # enumerate sub-counts c1
        for i in range(m):
            pass
        # generate all c1 componentwise <= cc
        ranges=[range(cc[j]+1) for j in range(m)]
        def gen(idx, cur):
            if idx==m:
                yield tuple(cur); return
            for val in ranges[idx]:
                cur.append(val)
                yield from gen(idx+1,cur)
                cur.pop()
        for c1 in gen(0,[]):
            if all(c1[j]==0 for j in range(m)) or all(c1[j]==cc[j] for j in range(m)):
                continue
            c2=tuple(cc[j]-c1[j] for j in range(m))
            f1,d1=feas(x,h,c1)
            if not f1: continue
            f2,d2=feas(w-x,h,c2)
            if f2:
                return True, f"V{x}: ({x}x{h})[{d1}] + ({w-x}x{h})[{d2}]"
    for y in range(1,h):
        ranges=[range(cc[j]+1) for j in range(m)]
        def gen2(idx, cur):
            if idx==m:
                yield tuple(cur); return
            for val in ranges[idx]:
                cur.append(val)
                yield from gen2(idx+1,cur)
                cur.pop()
        for c1 in gen2(0,[]):
            if all(c1[j]==0 for j in range(m)) or all(c1[j]==cc[j] for j in range(m)):
                continue
            c2=tuple(cc[j]-c1[j] for j in range(m))
            f1,d1=feas(w,y,c1)
            if not f1: continue
            f2,d2=feas(w,h-y,c2)
            if f2:
                return True, f"H{y}: ({w}x{y})[{d1}] + ({w}x{h-y})[{d2}]"
    return False, None
f,d=feas(W,H,c)
print("feasible",f,"calls",calls)
print(d)
