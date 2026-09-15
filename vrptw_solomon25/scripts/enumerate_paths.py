import math, time
DATA="/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"

def parse(path=DATA):
    lines=open(path).read().splitlines()
    customers=[]
    for l in lines:
        s=l.split()
        if len(s)==7 and s[0].isdigit():
            no,x,y,dem,ready,due,svc=map(int,s)
            customers.append((no,x,y,dem,ready,due,svc))
    customers.sort(key=lambda c:c[0])
    return customers[0], customers[1:]

def enumerate_paths(verbose=True):
    depot, custs = parse()
    n=len(custs)
    # node 0 depot, 1..n customers
    x=[depot[1]]+[c[1] for c in custs]
    y=[depot[2]]+[c[2] for c in custs]
    dem=[0]+[c[3] for c in custs]
    ready=[0]+[c[4] for c in custs]
    due=[0]+[c[5] for c in custs]
    svc=[0]+[c[6] for c in custs]
    cap=200
    depot_due=depot[5]
    def dist(i,j):
        return math.hypot(x[i]-x[j], y[i]-y[j])
    paths=[]
    # DFS state: current node, current time (start of service at current), load, visited bitmask
    # We'll use recursion; n=25 small
    import sys
    sys.setrecursionlimit(10000)
    def dfs(cur, t, load, visited_mask):
        # try extend to each unvisited customer
        for j in range(1, n+1):
            if visited_mask & (1<<j):
                continue
            arr = t + svc[cur] + dist(cur, j)
            if arr < ready[j]:
                arr = ready[j]
            if arr > due[j]:
                continue
            if load + dem[j] > cap:
                continue
            # return feasibility
            if arr + svc[j] + dist(j,0) > depot_due + 1e-9:
                continue
            # record path ending at j (customers in path excluding depot)
            # path represented as tuple of visited customer indices
            new_mask = visited_mask | (1<<j)
            # We record only the sequence; but to avoid huge memory, append sequence list
            paths.append(list(seq) + [j])
            seq.append(j)
            dfs(j, arr, load+dem[j], new_mask)
            seq.pop()
    seq=[]
    dfs(0, 0.0, 0, 0)
    return paths, (x,y,dem,ready,due,svc,cap,depot_due)

if __name__=="__main__":
    t0=time.time()
    paths, info = enumerate_paths()
    print("num paths", len(paths), "time", round(time.time()-t0,2))
    from collections import Counter
    print("max len", max(len(p) for p in paths), "min", min(len(p) for p in paths))
    print("lens", Counter(len(p) for p in paths).most_common(10))
    # sample longest
    lp=max(paths, key=len)
    print("longest", len(lp), lp)
