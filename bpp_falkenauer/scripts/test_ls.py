import random, math, time
from collections import defaultdict

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
C=150; n=120

def ffd_pack():
    order=sorted(range(n), key=lambda i:-sizes[i])
    bins=[]; assign={}
    for i in order:
        s=sizes[i]; placed=False
        for b in range(len(bins)):
            if bins[b]>=s:
                bins[b]-=s; assign[i]=b; placed=True; break
        if not placed:
            assign[i]=len(bins); bins.append(C-s)
    return assign, len(bins)

def assign_to_bins(assign, K):
    bins=[[] for _ in range(K)]
    for i,b in assign.items():
        bins[b].append(i)
    return bins

def load(b):
    return sum(sizes[i] for i in b)

def total_overflow(bins, K):
    return sum(max(0, load(b)-C) for b in bins)

# steepest-descent local search: move or swap to reduce overflow
def local_search(bins, K, max_iters=200000):
    # bins: list of lists of item indices
    # work with loads
    loads=[load(b) for b in bins]
    ov=sum(max(0,l-C) for l in loads)
    for it in range(max_iters):
        improved=False
        # try move: item from overloaded bin to another bin
        for bi in range(K):
            if loads[bi]<=C: continue
            for idx in list(bins[bi]):
                s=sizes[idx]
                # best target bin with min resulting overflow gain
                for bj in range(K):
                    if bj==bi: continue
                    delta = (max(0, loads[bj]+s-C)-max(0,loads[bj]-C)) + (max(0, loads[bi]-s-C)-max(0,loads[bi]-C))
                    if delta < -1e-9:
                        bins[bi].remove(idx); bins[bj].append(idx)
                        loads[bi]-=s; loads[bj]+=s
                        ov+=delta; improved=True
                        break
                if improved: break
            if improved: break
        if improved: continue
        # try swap: item in overloaded bin with item in other bin
        for bi in range(K):
            if loads[bi]<=C: continue
            for bj in range(K):
                if bj==bi: continue
                for ii in list(bins[bi]):
                    si=sizes[ii]
                    for jj in list(bins[bj]):
                        sj=sizes[jj]
                        if si<=sj: continue
                        new_li=loads[bi]-si+sj; new_lj=loads[bj]-sj+si
                        delta=(max(0,new_li-C)-max(0,loads[bi]-C))+(max(0,new_lj-C)-max(0,loads[bj]-C))
                        if delta < -1e-9:
                            bins[bi].remove(ii); bins[bj].remove(jj)
                            bins[bi].append(jj); bins[bj].append(ii)
                            loads[bi]=new_li; loads[bj]=new_lj
                            ov+=delta; improved=True
                            break
                    if improved: break
                if improved: break
            if improved: break
        if not improved:
            break
    return bins, ov

# try to eliminate one bin from FFD (49) and repair to 48
for seed in range(30):
    random.seed(seed)
    assign,K = ffd_pack()
    # pick a bin to empty: try smallest load bins
    bins0=assign_to_bins(assign,K)
    order_bins=sorted(range(K), key=lambda b: load(bins0[b]))
    for target in order_bins:
        bins=[list(b) for b in bins0]
        bins.pop(target)
        # move items of target into others greedily (best fit)
        K2=K-1
        bins=[list(b) for b in bins]
        loads=[load(b) for b in bins]
        items=list(bins0[target])
        ok=True
        for idx in sorted(items, key=lambda i:-sizes[i]):
            s=sizes[idx]
            # best fit among bins with room
            bestb=-1; bestres=None
            for b in range(K2):
                if loads[b]+s<=C and (bestres is None or loads[b]>bestres):
                    bestb=b; bestres=loads[b]
            if bestb>=0:
                bins[bestb].append(idx); loads[bestb]+=s
            else:
                # allow overflow: put into bin with least load
                bestb=min(range(K2), key=lambda b: loads[b])
                bins[bestb].append(idx); loads[bestb]+=s
        ov=sum(max(0,l-C) for l in loads)
        bins2, ov2 = local_search(bins, K2, max_iters=100000)
        if ov2<=1e-9:
            print("FOUND 48 packing seed", seed, "target", target)
            # verify
            for b in bins2:
                assert sum(sizes[i] for i in b)<=C
            print("bins", len(bins2))
            # save assignment
            asg={}
            for b in range(len(bins2)):
                for i in bins2[b]:
                    asg[i]=b
            with open("/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/scripts/packing48.txt","w") as f:
                f.write("\n".join(str(asg[i]) for i in range(n)))
            print("saved")
            raise SystemExit
    print("seed", seed, "no 48 found")
print("failed to find 48")
