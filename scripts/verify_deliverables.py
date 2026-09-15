#!/usr/bin/env python3
"""交付验收：检查每个家族 5 个 notebook（已执行、无 error 输出）、README.md、results.json。"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # 套件根目录（跨平台）
FAMILIES = ["csp_1d_waescher", "bpp_falkenauer", "scp_beasley", "crew_scheduling", "vrptw_solomon25", "cutting_2d_cgcut"]
NBS = ["01_direct.ipynb", "02_column_generation.ipynb", "03_benders.ipynb", "04_lagrangian.ipynb", "05_lbbd.ipynb"]
METHOD_KEYS = ["direct", "column_generation", "benders", "lagrangian", "lbbd"]

problems = []
for fam in FAMILIES:
    d = os.path.join(ROOT, fam)
    for nb in NBS:
        p = os.path.join(d, nb)
        if not os.path.exists(p):
            problems.append(f"{fam}/{nb}: MISSING")
            continue
        try:
            book = json.load(open(p, encoding="utf-8"))
        except Exception as e:
            problems.append(f"{fam}/{nb}: bad JSON {e}")
            continue
        ncells = book.get("nbformat") and len(book.get("cells", []))
        code_cells = [c for c in book.get("cells", []) if c.get("cell_type") == "code"]
        executed = sum(1 for c in code_cells if c.get("execution_count"))
        errs = 0
        for c in code_cells:
            for o in c.get("outputs", []):
                if o.get("output_type") == "error":
                    errs += 1
                    problems.append(f"{fam}/{nb}: ERROR output {o.get('ename')}: {str(o.get('evalue'))[:80]}")
        if ncells is None:
            problems.append(f"{fam}/{nb}: not nbformat json")
        elif executed < max(1, len(code_cells) // 2):
            problems.append(f"{fam}/{nb}: only {executed}/{len(code_cells)} code cells executed")
        elif errs == 0:
            print(f"OK   {fam}/{nb}: {len(code_cells)} code cells, {executed} executed, no errors")
    # results.json
    rp = os.path.join(d, "results.json")
    if not os.path.exists(rp):
        problems.append(f"{fam}/results.json: MISSING")
    else:
        try:
            r = json.load(open(rp, encoding="utf-8"))
            missing = [k for k in METHOD_KEYS if k not in r.get("methods", {})]
            if missing:
                problems.append(f"{fam}/results.json: missing methods {missing}")
            else:
                print(f"OK   {fam}/results.json: instance={r.get('instance')} best_known={r.get('best_known')}")
        except Exception as e:
            problems.append(f"{fam}/results.json: bad JSON {e}")
    if not os.path.exists(os.path.join(d, "README.md")):
        problems.append(f"{fam}/README.md: MISSING")

print()
if problems:
    print("PROBLEMS:")
    for p in problems:
        print(" -", p)
    sys.exit(1)
print("ALL DELIVERABLES OK")
