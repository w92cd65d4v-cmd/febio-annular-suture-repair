#!/usr/bin/env python3
"""Metrics of tighten-knot-release runs (reference state t = 2 after release; F = 2000*(t-2)/9).
Same definitions as the main analysis: 8 NP-defect interface nodes (anchors of the 3-suture model excluded),
outward posterolateral direction, positive-only AUC, restoration vs the original defect/intact models,
rim von Mises over the 57 rim elements, and suture (knot) tension = sum of bite-spring forces.
Usage: python3 revision_tighten_metrics.py RUN_STEM [RUN_STEM ...]"""
import sys, re, json, numpy as np
import revision_meshsens_analyze as A
from revision_mesh_refine import parse as pfeb
BASE = {"pure": dict(defect=515.6, intact=154.7, dpk=0.358, ipk=0.107), "flex": dict(defect=556.5, intact=211.9, dpk=0.474, ipk=0.183)}
K, KNOT = 50.0, 0.1
def records(log):
    rec = {}
    for blk in re.split(r'Data Record #\d+', open(log, encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        k = {"x;y;z": "x", "sx;sy;sz;sxy;syz;sxz": "s"}.get(md.group(1))
        if not k: continue
        d = rec.setdefault(float(mt.group(1)), {"x": {}, "s": {}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            d[k][int(lm.group(1))] = np.array([float(v) for v in lm.group(2).split()])
    return rec
def tension(info, N, X):
    out = []
    for s in info:
        sp = [(a, b) for a in s["biteA"] for b in s["biteB"]]; kk = K / len(sp); u = None; T = 0.0
        a0, b0 = s["anchor"]; u = (N[b0] - N[a0]) / np.linalg.norm(N[b0] - N[a0])
        for a, b in sp:
            L0 = np.linalg.norm(N[b] - N[a]); Lt = np.linalg.norm((N[b] - s["closure"]/2*u) - (N[a] + s["closure"]/2*u))
            dl = max(L0 - Lt, 0.0) + KNOT
            if a in X and b in X: T += kk * max(np.linalg.norm(X[b] - X[a]) - L0 + dl, 0.0)
            else: return None
        out.append(T)
    return out
for stem in sys.argv[1:]:
    scen = "flex" if "_flex" in stem else "pure"
    rec = records(stem + ".log"); ts = sorted(t for t in rec if rec[t]["x"])
    ref = rec[min(ts, key=lambda t: abs(t - 2.0))]["x"]
    F, p = [], []
    for t in ts:
        if t < 2.0 - 1e-6: continue
        X = rec[t]["x"]; F.append(2000 * (t - 2) / 9)
        p.append(np.mean([(X[n] - ref[n])[:2] @ A.NOUT for n in A.MOUTH8]))
    F, p = np.array(F), np.array(p); auc = np.trapz(np.maximum(p, 0), F); b = BASE[scen]
    s_end = rec[ts[-1]]["s"]; vm = np.array([A.vm(s_end[e]) for e in A.RIM if e in s_end])
    info = json.load(open(stem + "_tighten.json"))
    N, _ = pfeb(open(stem + ".feb", encoding="ISO-8859-1").read())
    t_rel = tension(info, N, ref); t_end = tension(info, N, rec[ts[-1]]["x"])
    r = dict(run=stem.split("/")[-1], F_end=round(F[-1]), peak=round(p.max(), 3), F_peak=round(F[p.argmax()]), posAUC=round(auc, 1),
             rest_AUC=round(100 * (b["defect"] - auc) / (b["defect"] - b["intact"]), 1),
             rest_peak=round(100 * (b["dpk"] - p.max()) / (b["dpk"] - b["ipk"]), 1), endpoint=round(p[-1], 3),
             rim_max=round(vm.max(), 2), rim_p95=round(float(np.percentile(vm, 95)), 2),
             knot_after_release=None if t_rel is None else [round(x, 1) for x in t_rel],
             knot_at_end=None if t_end is None else [round(x, 1) for x in t_end])
    print(json.dumps(r))
