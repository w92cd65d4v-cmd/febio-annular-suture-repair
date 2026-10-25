#!/usr/bin/env python3
"""Consolidated re-analysis with the UNIFIED protrusion metric.

Protrusion reference state = start of compression (t~1.0, F=0), matching
compare_repair.py / make_sweep_plot.py (the manuscript's "load-induced" metric).
This script reports BOTH the 2000 N endpoint and path-based metrics
(peak protrusion, load at peak, area under protrusion-load curve), plus
restoration computed under each definition, for:
    intact / defect / suture1 / suture3   x   pure / flex   (suture E = 50 N/mm).
Pure forensic recomputation from the FEBio .log files; writes nothing but stdout
plus a machine-readable path_metrics.json.
"""
import re, math, json

S = json.load(open("metric_sets.json"))
SUT = set(n for cfg in json.load(open("suture_config.json")).values() for pr in cfg for n in pr)
MOUTH, RIM, CORE = set(S["mouth"]) - SUT, set(S["rim"]), set(S["core"])

base = open("comp_intact_pure.feb", encoding="ISO-8859-1").read()
nsec = base[base.index("<Nodes"):base.index("</Nodes>")]
init = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(','))
        for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec) if int(m.group(1)) in MOUTH}
cx, cy = -9.85, -17.12
mx = sum(p[0] for p in init.values())/len(init); my = sum(p[1] for p in init.values())/len(init)
nx, ny = mx-cx, my-cy; Ln = math.hypot(nx, ny); nout = (nx/Ln, ny/Ln)

def vm(s):
    sx,sy,sz,sxy,syz,sxz = s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def F_of_t(t): return max(0.0, 2000.0*(t-1)/9.0)

def series(path):
    """Return list of (F, protrusion, rimVM) over the compression phase, t~1 reference."""
    try: txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    except FileNotFoundError: return []
    rec = {}
    for blk in re.split(r'Data Record #\d+', txt):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        t = float(mt.group(1)); d = rec.setdefault(t, {"pos":{}, "rim":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            i = int(lm.group(1)); v = [float(x) for x in lm.group(2).split()]
            if md.group(1) == "x;y;z" and i in MOUTH and len(v)>=3: d["pos"][i]=v[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and i in RIM and len(v)>=6: d["rim"][i]=v[:6]
    ts = sorted(rec)
    if not ts: return []
    t0 = min(ts, key=lambda t: abs(t-1.0)); ref = dict(rec[t0]["pos"])
    out = []
    for t in ts:
        if t < t0 - 1e-9: continue
        d = rec[t]
        if not d["pos"]: continue
        pr = [(p[0]-ref[i][0])*nout[0]+(p[1]-ref[i][1])*nout[1] for i,p in d["pos"].items() if i in ref]
        prot = sum(pr)/len(pr) if pr else 0.0
        rs = max((vm(s) for s in d["rim"].values()), default=0.0)
        out.append((F_of_t(t), prot, rs))
    return out

def auc(rows):
    a = 0.0
    for (f0,p0,_),(f1,p1,_) in zip(rows, rows[1:]):
        a += 0.5*(p0+p1)*(f1-f0)
    return a

STATES = ["intact","defect","suture1","suture3"]
allm = {}
for scen in ("pure","flex"):
    print("="*88); print(f"SCENARIO: {scen} compression")
    M = {}
    for st in STATES:
        r = series(f"comp_{st}_{scen}.log")
        if not r: print(f"  {st}: no data"); continue
        F_end, p_end, rim_end = r[-1]
        p_peak, F_peak = max((p,f) for f,p,_ in r)
        rim_peak = max(rr for _,_,rr in r)
        a = auc(r)
        M[st] = dict(p_end=p_end, F_end=F_end, p_peak=p_peak, F_peak=F_peak,
                     rim_end=rim_end, rim_peak=rim_peak, auc=a)
    allm[scen] = M
    hdr = f"{'state':9s} {'end prot':>9s} {'peak prot':>10s} {'F@peak':>8s} {'AUC(mm·N)':>11s} {'rim_end':>8s} {'rim_peak':>9s}"
    print(hdr); print("-"*len(hdr))
    for st in STATES:
        if st in M:
            m = M[st]
            print(f"{st:9s} {m['p_end']:9.3f} {m['p_peak']:10.3f} {m['F_peak']:8.0f} {m['auc']:11.1f} {m['rim_end']:8.3f} {m['rim_peak']:9.3f}")
    # restoration under each definition
    if all(st in M for st in STATES):
        print("  -- restoration of the defect's deviation from intact (higher = better) --")
        for key,label in (("p_end","endpoint protrusion"),("p_peak","PEAK protrusion"),("auc","AUC protrusion")):
            di, dd = M["intact"][key], M["defect"][key]
            if abs(dd-di) < 1e-12: continue
            s1 = 100*(dd-M["suture1"][key])/(dd-di)
            s3 = 100*(dd-M["suture3"][key])/(dd-di)
            print(f"     {label:22s}:  1-suture {s1:6.1f}%   3-suture {s3:6.1f}%   (intact {di:.3f}, defect {dd:.3f})")
    print()

json.dump(allm, open("path_metrics.json","w"), indent=2)
print("wrote path_metrics.json")
