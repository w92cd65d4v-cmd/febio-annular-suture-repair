#!/usr/bin/env python3
"""Extract bite-node reaction forces at the end of the tightening step (t=1) from a pass-A log -> JSON."""
import sys, re, json
log, out = sys.argv[1], sys.argv[2]
best = None
for blk in re.split(r'Data Record #\d+', open(log, encoding="ISO-8859-1", errors="ignore").read()):
    mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
    if mt and md and md.group(1) in ("Rx;Ry;Rz", "react"):
        t = float(mt.group(1))
        R = {int(m.group(1)): [float(v) for v in m.group(2).split()] for m in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M)}
        if best is None or t >= best[0]: best = (t, R)
assert best and abs(best[0] - 1.0) < 1e-6, f"pass A did not reach t=1 (last {best and best[0]})"
json.dump(best[1], open(out, "w")); print(f"t={best[0]}: {len(best[1])} bite nodes, |R| max {max(sum(x*x for x in v)**.5 for v in best[1].values()):.2f} N")
