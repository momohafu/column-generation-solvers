import math, time
from ortools.sat.python import cp_model

def parse_first_instance(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    P = int(lines[0].strip())
    idx = 1
    name = lines[idx].strip(); idx += 1
    parts = lines[idx].split(); idx += 1
    C = int(parts[0]); n = int(parts[1]); best = int(parts[2])
    sizes = [int(lines[idx+i]) for i in range(n)]
    idx += n
    return name, C, n, best, sizes

def ffd(C, sizes):
    items = sorted(range(len(sizes)), key=lambda i: -sizes[i])
    bins = []  # list of residual capacity
    assign = {}
    for i in items:
        s = sizes[i]
        placed = False
        for b in range(len(bins)):
            if bins[b] >= s:
                bins[b] -= s; assign[i] = b; placed = True; break
        if not placed:
            assign[i] = len(bins); bins.append(C - s)
    return len(bins), assign

name, C, n, best, sizes = parse_first_instance("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt")
print("instance", name, "C", C, "n", n, "file_best", best)
ub_ffd, _ = ffd(C, sizes)
lb = max((sum(sizes) + C - 1)//C, max(1, math.ceil(sum(sizes)/C)))
print("LB", lb, "FFD UB", ub_ffd)

K = ub_ffd
model = cp_model.CpModel()
x = [[model.NewBoolVar(f"x_{i}_{b}") for b in range(K)] for i in range(n)]
y = [model.NewBoolVar(f"y_{b}") for b in range(K)]

for i in range(n):
    model.Add(sum(x[i][b] for b in range(K)) == 1)
for b in range(K):
    model.Add(sum(sizes[i]*x[i][b] for i in range(n)) <= C*y[b])
for b in range(K-1):
    model.Add(y[b] >= y[b+1])
# symmetry: item i can go to bin b only if b <= i
for i in range(n):
    for b in range(i+1, K):
        model.Add(x[i][b] == 0)
model.Add(sum(y) >= lb)
model.Minimize(sum(y))

solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 120.0
solver.parameters.num_search_workers = 8
solver.parameters.log_search_progress = True
t0 = time.time()
status = solver.Solve(model)
elapsed = time.time() - t0
print("status", solver.StatusName(status), "elapsed", round(elapsed,2))
print("objective", solver.ObjectiveValue() if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None)
print("best_bound", solver.BestObjectiveBound())
print("wall", solver.WallTime())

