#!/usr/bin/env python3
"""Analyze the compression runs: nucleus protrusion through the defect, annulus
rim stress concentration, and intradiscal pressure vs applied axial load.

Applied axial force F(t) = 2000*(t-1)/9 N for t in [1,10] (linear load curve 2).
Protrusion = displacement of nucleus 'mouth' nodes projected on the outward
posterolateral (horizontal) direction, referenced to the START of compression
(t~1, F=0) so it is the load-induced movement only (unified with compare_repair.py).
NOTE: superseded by compare_repair.py (4-way) and path_metrics.py (path metrics);
retained only as an intact-vs-defect diagnostic.
Rim stress = max von Mises over the annulus rim elements.
Intradiscal pressure = -(sx+sy+sz)/3 averaged over nucleus core elements.
"""
import re, math, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

S = json.load(open("../validation_L45_defect/metric_sets.json"))
MOUTH, RIM, CORE = set(S["mouth"]), set(S["rim"]), set(S["core"])

# initial mouth-node coordinates + outward posterolateral direction
base = open("comp_intact_pure.feb", encoding="ISO-8859-1").read()
nsec = base[base.index("<Nodes"):base.index("</Nodes>")]
init = {}
for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec):
    nid = int(m.group(1))
    if nid in MOUTH: init[nid] = tuple(float(v) for v in m.group(2).split(','))
cx, cy = -9.85, -17.12
mx = sum(p[0] for p in init.values())/len(init); my = sum(p[1] for p in init.values())/len(init)
nout = (mx-cx, my-cy); L = math.hypot(*nout); nout = (nout[0]/L, nout[1]/L)  # horizontal outward

def F_of_t(t): return max(0.0, 2000.0*(t-1)/9.0)
def vm(s):  # s = sx,sy,sz,sxy,syz,sxz
    sx,sy,sz,sxy,syz,sxz = s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))

def parse(path):
    txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    # per time: gather node x;y;z (mouth), element stress (rim 6-comp, core 3-comp)
    rec = {}  # t -> dict
    for blk in re.split(r'Data Record #\d+', txt):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        t = float(mt.group(1)); d = rec.setdefault(t, {"pos":{}, "rim":{}, "core":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            iid = int(lm.group(1)); vals = [float(x) for x in lm.group(2).split()]
            if md.group(1) == "x;y;z" and iid in MOUTH and len(vals)>=3: d["pos"][iid]=vals[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and iid in RIM and len(vals)>=6: d["rim"][iid]=vals[:6]
            elif md.group(1)=="sx;sy;sz" and iid in CORE and len(vals)>=3: d["core"][iid]=vals[:3]
    # UNIFIED metric: reference = each node's position at the START of compression
    # (t~1.0, F=0), matching compare_repair.py / make_sweep_plot.py / path_metrics.py.
    # (Previously this script referenced the as-built initial geometry, which mixed
    #  the preload bulge into "load-induced" protrusion — superseded.)
    ts = sorted(rec)
    if not ts: return []
    t0 = min(ts, key=lambda t: abs(t-1.0))
    ref = dict(rec[t0]["pos"])
    rows = []
    for t in ts:
        if t < t0 - 1e-9: continue              # compression phase only
        d = rec[t]
        if not d["pos"]: continue
        prot = [ (p[0]-ref[nid][0])*nout[0] + (p[1]-ref[nid][1])*nout[1]
                 for nid, p in d["pos"].items() if nid in ref ]
        protrusion = sum(prot)/len(prot) if prot else 0.0
        rimstress = max((vm(s) for s in d["rim"].values()), default=0.0)
        idp = -sum((s[0]+s[1]+s[2])/3 for s in d["core"].values())/len(d["core"]) if d["core"] else 0.0
        rows.append((F_of_t(t), protrusion, rimstress, idp))
    return rows

cases = {}
for model in ("intact","defect"):
    for scen in ("pure","flex"):
        p = f"comp_{model}_{scen}.log"
        try: cases[(model,scen)] = parse(p)
        except FileNotFoundError: cases[(model,scen)] = []

print(f"outward posterolateral direction (horizontal) = ({nout[0]:.2f},{nout[1]:.2f})")
print("="*72)
print(f"{'case':18s} {'maxF(N)':>8s} {'protrusion(mm)':>15s} {'peakRimVM(MPa)':>15s} {'IDP(MPa)':>9s}")
print("-"*72)
for k,rows in cases.items():
    if rows:
        F,pr,rs,idp = rows[-1]
        print(f"{k[0]+'_'+k[1]:18s} {F:8.0f} {pr:15.3f} {rs:15.3f} {idp:9.3f}")
    else:
        print(f"{k[0]+'_'+k[1]:18s}   (no data)")

# figures: protrusion vs F, rim stress vs F  (pure & flex)
for scen,fname in (("pure","herniation_pure.png"),("flex","herniation_flex.png")):
    fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(13,4))
    for model,col in (("intact","#888"),("defect","#A32D2D")):
        r=cases.get((model,scen),[])
        if r:
            F=[x[0] for x in r]
            a1.plot(F,[x[1] for x in r],'-o',ms=2,color=col,label=model)
            a2.plot(F,[x[2] for x in r],'-o',ms=2,color=col,label=model)
            a3.plot(F,[x[3] for x in r],'-o',ms=2,color=col,label=model)
    a1.set_title("Nucleus protrusion through defect"); a1.set_ylabel("outward protrusion (mm)")
    a2.set_title("Peak annulus rim stress"); a2.set_ylabel("von Mises (MPa)")
    a3.set_title("Intradiscal pressure"); a3.set_ylabel("IDP (MPa)")
    for a in (a1,a2,a3): a.set_xlabel("axial compressive force (N)"); a.grid(alpha=.25); a.legend(fontsize=8)
    fig.suptitle(f"L4-L5 {scen} compression: intact vs defect",y=1.02); fig.tight_layout()
    fig.savefig(fname,dpi=140,bbox_inches="tight"); print("wrote",fname)
