#!/usr/bin/env python3
"""汇总 6 个家族的 results.json -> results_summary.csv + markdown 表格片段。"""
import json, csv, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILIES = ["csp_1d_waescher", "bpp_falkenauer", "scp_beasley", "crew_scheduling", "vrptw_solomon25", "cutting_2d_cgcut"]
METHODS = ["direct", "column_generation", "benders", "lagrangian", "lbbd"]

rows = []
for fam in FAMILIES:
    p = os.path.join(ROOT, fam, "results.json")
    if not os.path.exists(p):
        rows.append({"family": fam, "missing": True})
        continue
    d = json.load(open(p, encoding="utf-8"))
    row = {
        "family": fam,
        "instance": d.get("instance", ""),
        "best_known": d.get("best_known", ""),
        "proved_optimal": d.get("proved_optimal", ""),
    }
    for m in METHODS:
        mm = d.get("methods", {}).get(m, {})
        row[m + "_obj"] = mm.get("objective", "")
        row[m + "_time"] = mm.get("time_s", "")
        row[m + "_gap"] = mm.get("gap", "")
        row[m + "_opt"] = mm.get("optimal", "")
        row[m + "_note"] = mm.get("note", "")
    rows.append(row)

csv_path = os.path.join(ROOT, "results_summary.csv")
with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["family", "instance", "best_known", "proved_optimal"]
               + [f"{m}_{k}" for m in METHODS for k in ("obj", "time", "gap", "opt", "note")])
    for r in rows:
        w.writerow([r.get(k, "") for k in ["family", "instance", "best_known", "proved_optimal"]]
                   + [r.get(f"{m}_{k}", "") for m in METHODS for k in ("obj", "time", "gap", "opt", "note")])
print("written", csv_path)

# markdown snippet
print("\n| 家族 | 直接建模 | 列生成 | Benders | 拉格朗日 | LBBD | 基准 |")
print("|---|---|---|---|---|---|---|")
for r in rows:
    def cell(m):
        if m + "_obj" not in r or r.get(m + "_obj") == "":
            return "-"
        s = f"{r[m+'_obj']} / {r[m+'_time']}s / {r[m+'_gap']}%"
        if not r.get(m + "_opt"):
            s += " †"
        return s
    print(f"| {r['family']} | {cell('direct')} | {cell('column_generation')} | {cell('benders')} | {cell('lagrangian')} | {cell('lbbd')} | {r.get('best_known','')} |")
