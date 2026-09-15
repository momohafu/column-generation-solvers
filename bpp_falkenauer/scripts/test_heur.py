import math, random
from collections import Counter

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes

name,C,n,best,sizes = parse()
items=sorted(range(n), key=lambda i:-sizes[i])

def first_fit(order):
    bins=[]; assign={}
    for i in order:
        s=sizes[i]; placed=False
        for b in range(len(bins)):
            if bins[b]>=s:
                bins[b]-=s; assign[i]=b; placed=True; break
        if not placed:
            assign[i]=len(bins); bins.append(C-s)
    return len(bins)

def best_fit(order):
    bins=[]; assign={}
    for i in order:
        s=sizes[i]
        # pick bin with min residual >= s
        best=-1; bestres=None
        for b in range(len(bins)):
            if bins[b]>=s and (bestres is None or bins[b]<bestres):
                best=b; bestres=bins[b]
        if best>=0:
            bins[best]-=s; assign[i]=best
        else:
            assign[i]=len(bins); bins.append(C-s)
    return len(bins)

print("FFD", first_fit(items), "BFD", best_fit(items))
print("FF increasing", first_fit(sorted(range(n), key=lambda i:sizes[i])))
print("BF increasing", best_fit(sorted(range(n), key=lambda i:sizes[i])))

# local search: random restarts with FFD-style + swap improvements
best_ub = 10**9
best_assign=None
random.seed(42)
for r in range(5000):
    order = list(range(n))
    random.shuffle(order)
    # try FF
    k = first_fit(order)
    if k < best_ub:
        best_ub = k
        print("random FF improved", k, "at r", r)
        if k <= 48: break
print("best random FF", best_ub)

