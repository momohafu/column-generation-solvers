# -*- coding: utf-8 -*-
"""三方法最终结果一致性校验（完整版，正确容差）：
重跑三方法 + 首原理复算 + 交叉对比 + 与 BKS 对比。
容差说明：精确双精度值 = 191.8136197787，其 6 位小数表示 = 191.813620（文档口径），
2 位小数表示 = 191.81（BKS 口径）。跨方法比较用 1e-9，与 6 位小数字面量比较用 1e-6。
"""
import sys, math
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import parse, column_generation, ip_recovery
import benders, lbbd

depot, cs = parse()
xc = [depot[1]] + [c[1] for c in cs]
yc = [depot[2]] + [c[2] for c in cs]
def d(i, j):
    return math.hypot(xc[i]-xc[j], yc[i]-yc[j])
canon = [(13,17,18,19,15,16,14,12), (20,24,25,23,22,21), (5,3,7,8,10,11,9,6,4,2,1)]
total = 0.0
for r in canon:
    seq = [0] + list(r) + [0]
    total += sum(d(seq[i], seq[i+1]) for i in range(len(seq)-1))
print(f"首原理复算精确总距离 = {total:.10f}（6 位小数 {total:.6f}，2 位小数 {total:.2f}）")

info = column_generation(verbose=False)
ip = ip_recovery(info)
cg_obj, cg_routes = ip["exact"], set(frozenset(r) for r in ip["routes"])

rb = benders.run_benders("cg", verbose=False)
b_obj, b_routes = rb["exact"], set(frozenset(r) for r in rb["routes"])

rl = lbbd.run_lbbd(verbose=False)
l_obj, l_routes = rl["ub"], set(frozenset(r) for r in rl["routes"].values())

print()
print("========== 最终一致性校验 ==========")
c1 = abs(cg_obj-b_obj) < 1e-9 and abs(b_obj-l_obj) < 1e-9
c2 = cg_routes == b_routes == l_routes
c3 = abs(cg_obj-total) < 1e-9
c4 = abs(cg_obj-191.813620) < 1e-6            # 6 位小数文档口径
c5 = abs(round(cg_obj, 2)-191.81) < 1e-9      # BKS 两位小数口径
c6 = abs(info["lp_obj"]-cg_obj) < 1e-9        # CG: LP 下界 = 整数目标
c7 = abs(rb["lb"]-b_obj) < 1e-9               # Benders: 下界 = 修复
c8 = rl["proven"]                             # LBBD: 主问题不可行证明
print(f"02 列生成 : 目标 {cg_obj:.10f} | LP 下界 {info['lp_obj']:.6f}")
print(f"03 Benders: 目标 {b_obj:.10f} | 下界 {rb['lb']:.6f}")
print(f"05 LBBD   : 目标 {l_obj:.10f} | LB_K={rl['LB_K']} | proven={rl['proven']}")
print(f"[1] 三方法目标一致(1e-9): {c1}")
print(f"[2] 三方法最优路线集合一致: {c2}")
print(f"[3] 与首原理复算一致(1e-9): {c3}")
print(f"[4] 与文档口径 191.813620 一致(1e-6): {c4}")
print(f"[5] 两位小数 = BKS 191.81: {c5}")
print(f"[6] CG 证明 LP下界=整数目标: {c6}")
print(f"[7] Benders 证明 下界=整数修复: {c7}")
print(f"[8] LBBD 证明 主问题不可行: {c8}")
ok = all([c1, c2, c3, c4, c5, c6, c7, c8])
print(f"综合判定: {'全部一致 PASS' if ok else '存在不一致 FAIL'}")
