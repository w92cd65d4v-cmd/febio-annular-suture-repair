#!/usr/bin/env python3
"""4-way comparison: intact vs defect vs 1-suture vs 3-suture repair.
Metrics vs axial load (0-2000 N): nucleus protrusion through the defect,
peak annulus rim von Mises, intradiscal pressure. Reports the fraction of the
defect's deficit that each repair restores toward intact."""
import re, math, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

S = json.load(open("metric_sets.json"))
# exclude suture anchor nodes from the protrusion set so the metric measures the
# nucleus bulge, not the nodes the sutures directly pull
SUT = set(n for cfg in json.load(open("suture_config.json")).values() for pr in cfg for n in pr)
MOUTH, RIM, CORE = set(S["mouth"]) - SUT, set(S["rim"]), set(S["core"])
print(f"protrusion measured on {len(MOUTH)} non-suture nucleus nodes (excluded anchors: {sorted(set(S['mouth'])&SUT)})")
STATES = ["intact", "defect", "suture1", "suture3"]
COL = {"intact":"#444", "defect":"#A32D2D", "suture1":"#1D9E75", "suture3":"#185FA5"}

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

def parse(path):
    try: txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    except FileNotFoundError: return []
    rec = {}
    for blk in re.split(r'Data Record #\d+', txt):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        t = float(mt.group(1)); d = rec.setdefault(t, {"pos":{}, "rim":{}, "core":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            i = int(lm.group(1)); v = [float(x) for x in lm.group(2).split()]
            if md.group(1) == "x;y;z" and i in MOUTH and len(v)>=3: d["pos"][i]=v[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and i in RIM and len(v)>=6: d["rim"][i]=v[:6]
            elif md.group(1)=="sx;sy;sz" and i in CORE and len(v)>=3: d["core"][i]=v[:3]
    # reference = each node's position at the START of compression (t~1.0, F=0),
    # so protrusion is the load-induced movement only ("arrow" metric)
    ts = sorted(rec)
    t0 = min(ts, key=lambda t: abs(t-1.0))
    ref = dict(rec[t0]["pos"])
    rows = []
    for t in sorted(rec):
        if t < t0 - 1e-9: continue   # compression phase only (drop preload points)
        d = rec[t]
        if not d["pos"]: continue
        pr = [ (p[0]-ref[i][0])*nout[0]+(p[1]-ref[i][1])*nout[1] for i,p in d["pos"].items() if i in ref ]
        prot = sum(pr)/len(pr) if pr else 0.0
        rs = max((vm(s) for s in d["rim"].values()), default=0.0)
        idp = -sum((s[0]+s[1]+s[2])/3 for s in d["core"].values())/len(d["core"]) if d["core"] else 0.0
        rows.append((F_of_t(t), prot, rs, idp))
    return rows

for scen in ("pure", "flex"):
    data = {st: parse(f"comp_{st}_{scen}.log") for st in STATES}
    print("="*70); print(f"SCENARIO: {scen} compression @ 2000 N")
    print(f"{'state':10s} {'protrusion(mm)':>15s} {'peakRimVM(MPa)':>15s} {'IDP(MPa)':>9s}")
    print("-"*70)
    last = {}
    for st in STATES:
        if data[st]:
            F,pr,rs,idp = data[st][-1]; last[st]=(pr,rs,idp)
            print(f"{st:10s} {pr:15.3f} {rs:15.3f} {idp:9.3f}")
    # % restoration of the defect deficit toward intact (protrusion & rim stress)
    if all(k in last for k in STATES):
        for metric,idx,name in ((0,0,"protrusion"),(1,1,"rim stress")):
            di = last["intact"][idx]; dd = last["defect"][idx]
            if abs(dd-di) > 1e-9:
                for rep in ("suture1","suture3"):
                    rr = last[rep][idx]
                    print(f"   {name}: {rep} restores {100*(dd-rr)/(dd-di):5.1f}% of the defect's deviation from intact")
    # figure
    fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(14,4.2))
    for st in STATES:
        r=data[st]
        if r:
            F=[x[0] for x in r]
            a1.plot(F,[x[1] for x in r],'-',color=COL[st],lw=1.8,label=st)
            a2.plot(F,[x[2] for x in r],'-',color=COL[st],lw=1.8,label=st)
            a3.plot(F,[x[3] for x in r],'-',color=COL[st],lw=1.8,label=st)
    a1.set_title("Nucleus protrusion (load-induced)"); a1.set_ylabel("outward protrusion (mm)")
    a2.set_title("Peak annulus rim stress"); a2.set_ylabel("von Mises (MPa)")
    a3.set_title("Intradiscal pressure"); a3.set_ylabel("IDP (MPa)")
    for a in (a1,a2,a3): a.set_xlabel("axial force (N)"); a.grid(alpha=.25); a.legend(fontsize=8)
    fig.suptitle(f"L4-L5 {scen} compression: intact vs defect vs suture repair",y=1.02)
    fig.tight_layout(); fig.savefig(f"repair_{scen}.png",dpi=140,bbox_inches="tight")
    print("wrote",f"repair_{scen}.png")
