import math, time

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes
name,C,n,best,sizes=parse()

def knap(items, C, objective="load"):
    # items: list of (idx, size). maximize load <= C. tie-break: fewer items (MBS-like)
    n=len(items)
    dp=[0]*(C+1)
    keep=[[False]*(C+1) for _ in range(n)]
    cnt=[[0]*(C+1) for _ in range(n)]
    for i,(idx,s) in enumerate(items):
        for c in range(C,s-1,-1):
            cand=dp[c-s]+s
            if cand>dp[c]:
                dp[c]=cand; keep[i][c]=True
            # for tie: prefer fewer items handled separately
    bestc=max(range(C+1), key=lambda c: dp[c])
    chosen=[]; c=bestc
    for i in range(n-1,-1,-1):
        if keep[i][c]:
            chosen.append(items[i][0]); c-=items[i][1]
    return chosen, dp[bestc]

def knap_few(items, C):
    # DP that tracks (load, -count) lexicographic max: prefer max load, then fewer items
    n=len(items)
    # state: (load, count) maximize load then minimize count -> use tuple compare with (-count)
    dp=[(0,0)]*(C+1)
    keep=[[False]*(C+1) for _ in range(n)]
    for i,(idx,s) in enumerate(items):
        for c in range(C,s-1,-1):
            l2, cnt2 = dp[c-s]; l2+=s; cnt2+=1
            l1, cnt1 = dp[c]
            if (l2, -cnt2) > (l1, -cnt1):
                dp[c]=(l2,cnt2); keep[i][c]=True
    bestc=max(range(C+1), key=lambda c: (dp[c][0], -dp[c][1]))
    chosen=[]; c=bestc
    for i in range(n-1,-1,-1):
        if keep[i][c]:
            chosen.append(items[i][0]); c-=items[i][1]
    return chosen, dp[bestc][0], dp[bestc][1]

for variant in ["load_only","fewer_items"]:
    remaining=[(i,sizes[i]) for i in range(n)]
    nbins=0; loads=[]
    while remaining:
        if variant=="load_only":
            chosen, load = knap(remaining, C)
        else:
            chosen, load, cnt = knap_few(remaining, C)
        nbins+=1; loads.append(load)
        chosen=set(chosen)
        remaining=[it for it in remaining if it[0] not in chosen]
    print(variant, "bins", nbins, "min/max load", min(loads), max(loads), "avg", sum(loads)/len(loads))
