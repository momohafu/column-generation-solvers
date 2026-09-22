# -*- coding: utf-8 -*-
"""ROTATE=False 检验：固定朝向池 IP 是否 = 244。"""
import sys
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts")
import cg2_core as cc
cc.ROTATE = False
cc.ORIENTS = [tuple(sorted(set([(a, b)]))) for (a, b, q, v) in cc.ITEMS]
cc.feas.cache_clear()
pool = cc.build_pool()
P = len(pool)
cc.POOL = pool
cc.P = P
cc.val = [sum(cc.V[i]*pool[p][i] for i in range(cc.M)) for p in range(P)]
term, obj, sel, wt = cc.pool_ip()
print(f"固定朝向: 可行模式 {P} | 池 IP {term} = {obj} | 模式 {pool[sel[0]] if sel else None} | 文献 244 一致: {abs(obj-244)<1e-6}")
