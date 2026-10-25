#!/usr/bin/env python3
"""Anchor-pair gaps of a tightening run at given times (default t=1 end of tightening, t=2 after release).
Usage: python3 revision_tighten_gaps.py tighten/RUN [t1 t2 ...]"""
import sys, json, numpy as np
stem = sys.argv[1]; ts = [float(x) for x in sys.argv[2:]] or [1.0, 2.0]
pairs = json.load(open("suture_config.json"))["suture3"]
import re
def all_records(path):
    rec = {}
    for blk in re.split(r'Data Record #\d+', open(path, encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md) or md.group(1) != "x;y;z": continue
        d = rec.setdefault(float(mt.group(1)), {})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            d[int(lm.group(1))] = np.array([float(v) for v in lm.group(2).split()])
    return rec
R = all_records(stem + ".log")
out = {}
for t in ts:
    tt = min(R, key=lambda x: abs(x - t)); X = R[tt]
    g = [float(np.linalg.norm(X[a] - X[b])) for a, b in pairs]
    out[t] = g; print(f"t={tt:.3f}: gaps " + " ".join(f"{x:.3f}" for x in g) + f"  mean {np.mean(g):.3f}")
json.dump({str(k): v for k, v in out.items()}, open(stem + "_gaps.json", "w"))
