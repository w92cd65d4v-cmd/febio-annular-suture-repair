#!/usr/bin/env python3
"""Issue #9 sensitivity analysis, using the UNIFIED t~1 protrusion metric.

(A) Defect-size:  intact / defectS(10 elem) / defect(18) / defectL(26)
(B) Suture number & placement:  defect / suture1 / suture2 / suture3 / suture3b
Metrics: endpoint (2000 N) & peak load-induced protrusion, AUC, endpoint peak rim VM.
Protrusion measured on the same mouth-node set for all states (suture anchors of every
config excluded) so values are directly comparable. Reads only the FEBio .log files.
"""
import re, math, json

S = json.load(open("metric_sets.json"))
RIM = set(S["rim"])

def anchors(feb):
    try: t = open(feb, encoding="ISO-8859-1").read()
    except FileNotFoundError: return set()
    m = re.search(r'<DiscreteSet name="SUTURE">(.*?)</DiscreteSet>', t, re.S)
    if not m: return set()
    return set(int(x) for d in re.finditer(r'<delem>(\d+),(\d+)</delem>', m.group(1)) for x in d.groups())

EXCL = set().union(*(anchors(f"comp_{s}_pure.feb") for s in ("suture1","suture2","suture3","suture3b")))
MOUTH = set(S["mouth"]) - EXCL
print(f"protrusion on {len(MOUTH)} mouth nodes (excluded {len(EXCL)} suture-anchor nodes across all configs)")

base = open("comp_intact_pure.feb", encoding="ISO-8859-1").read()
nsec = base[base.index("<Nodes"):base.index("</Nodes>")]
init = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(',')) for m in
        re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec) if int(m.group(1)) in MOUTH}
cx, cy = -9.85, -17.12
mx = sum(p[0] for p in init.values())/len(init); my = sum(p[1] for p in init.values())/len(init)
nx, ny = mx-cx, my-cy; Ln = math.hypot(nx, ny); nout = (nx/Ln, ny/Ln)

def vm(s):
    sx,sy,sz,sxy,syz,sxz = s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def F_of_t(t): return max(0.0, 2000.0*(t-1)/9.0)

def series(path):
    try: txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    except FileNotFoundError: return None
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
    if not ts: return None
    t0 = min(ts, key=lambda t: abs(t-1.0)); ref = dict(rec[t0]["pos"])
    rows = []
    for t in ts:
        if t < t0-1e-9 or not rec[t]["pos"]: continue
        pr = [(p[0]-ref[i][0])*nout[0]+(p[1]-ref[i][1])*nout[1] for i,p in rec[t]["pos"].items() if i in ref]
        prot = sum(pr)/len(pr) if pr else 0.0
        rs = max((vm(s) for s in rec[t]["rim"].values()), default=0.0)
        rows.append((F_of_t(t), prot, rs))
    return rows

def metrics(path):
    r = series(path)
    if not r: return None
    p_end = r[-1][1]; rim_end = r[-1][2]
    p_peak, F_peak = max((p, f) for f, p, _ in r)
    auc = sum(0.5*(a[1]+b[1])*(b[0]-a[0]) for a, b in zip(r, r[1:]))
    return dict(p_end=p_end, p_peak=p_peak, F_peak=F_peak, rim_end=rim_end, auc=auc)

def rest(dd, di, rr):  # % restoration of defect deviation toward intact
    return 100*(dd-rr)/(dd-di) if abs(dd-di) > 1e-9 else float('nan')

for scen in ("pure", "flex"):
    M = {st: metrics(f"comp_{st}_{scen}.log") for st in
         ("intact","defectS","defect","defectL","suture1","suture2","suture3","suture3b")}
    print("\n" + "="*94); print(f"SCENARIO: {scen} compression")
    print("(A) DEFECT-SIZE sensitivity")
    print(f"  {'model':10s} {'#elem':>6s} {'end prot':>9s} {'peak prot':>10s} {'F@peak':>7s} {'AUC':>8s} {'rim_end':>8s}")
    for st,ne in (("intact",0),("defectS",10),("defect",18),("defectL",26)):
        m = M[st]
        if m: print(f"  {st:10s} {ne:6d} {m['p_end']:9.3f} {m['p_peak']:10.3f} {m['F_peak']:7.0f} {m['auc']:8.1f} {m['rim_end']:8.3f}")
    print("(B) SUTURE number & placement (restoration of defect->intact at the 2000 N endpoint)")
    print(f"  {'model':10s} {'#sut':>4s} {'end prot':>9s} {'peak prot':>10s} {'rim_end':>8s} {'end%rest':>9s} {'peak%rest':>9s}")
    di, dd = M["intact"], M["defect"]
    for st,ns in (("defect",0),("suture1",1),("suture2",2),("suture3",3),("suture3b",3)):
        m = M[st]
        if not m: print(f"  {st:10s}  (no data)"); continue
        re_end = rest(dd['p_end'], di['p_end'], m['p_end']) if st!="defect" else 0.0
        re_pk  = rest(dd['p_peak'], di['p_peak'], m['p_peak']) if st!="defect" else 0.0
        tag = " (alt placement)" if st=="suture3b" else ""
        print(f"  {st:10s} {ns:4d} {m['p_end']:9.3f} {m['p_peak']:10.3f} {m['rim_end']:8.3f} {re_end:8.1f}% {re_pk:8.1f}%{tag}")
