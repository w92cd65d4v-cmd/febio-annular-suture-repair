#!/usr/bin/env python3
"""Diagnose element distortion in a refined mesh-sensitivity run.
For every child element in the refinement zone, computes the corner Jacobian ratio J/J0 at the
last converged logged state and reports the worst elements, their parent, material and location
(relative to the defect centre / NP-defect interface).
Usage: python3 revision_meshsens_distortion.py meshsens/app12_pure_h2
"""
import sys, json, numpy as np, re
from revision_mesh_refine import parse as pfeb
from revision_meshsens_analyze import parse
NB = {0:(1,3,4),1:(2,0,5),2:(3,1,6),3:(0,2,7),4:(7,5,0),5:(4,6,1),6:(5,7,2),7:(6,4,3)}
def cornerJ(X):
    return np.array([np.dot(X[a]-X[i], np.cross(X[b]-X[i], X[c]-X[i])) for i,(a,b,c) in NB.items()])
stem = sys.argv[1]
N0, B = pfeb(open(stem + ".feb", encoding="ISO-8859-1").read())
E = {e: c for b in B for e, c in b["els"]}; NAME = {e: b["name"] for b in B for e, _ in b["els"]}
meta = json.load(open(stem + "_refine.json"))
kid2par = {c: int(p) for p, ch in meta["children"].items() for c in ch}
data = parse(stem + ".log"); t, d = data[-1]; X = d["x"]
C = np.array(json.load(open("region_sets.json"))["C"]); S = json.load(open("metric_sets.json"))
rimset = set(S["rim"])
rows = []
for c, p in kid2par.items():
    if not all(n in X for n in E[c]): continue
    J0 = cornerJ(np.array([N0[n] for n in E[c]])); J = cornerJ(np.array([X[n] for n in E[c]]))
    r = (J / J0).min()
    rows.append((r, c, p, NAME[c], np.linalg.norm(np.mean([N0[n] for n in E[c]], 0) - C), p in rimset))
rows.sort()
print(f"{stem}: last logged t={t:.4f} (F={2000*(t-1)/9:.0f} N); children checked {len(rows)}")
for r, c, p, nm, dist, rim in rows[:12]:
    print(f"  J/J0 {r:+.3f}  child {c} parent {p} {nm:14s} {dist:5.2f} mm from defect centre  rim-parent={rim}")
