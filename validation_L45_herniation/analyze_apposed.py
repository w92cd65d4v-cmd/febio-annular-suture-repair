#!/usr/bin/env python3
"""Compare intact / unrepaired defect / non-apposed bridge repair / apposed-anchor repair.
Metrics (load-induced 'arrow' protrusion, referenced to start of compression F=0, on the
nucleus mouth nodes with suture anchors excluded; consistent with compare_repair.py):
  endpoint protrusion @2000 N, peak protrusion, load at peak, AUC(protrusion vs F),
  peak rim von Mises @2000 N, intradiscal pressure @2000 N.
Restoration % computed vs the defect deficit toward intact (endpoint and AUC)."""
import re, math, json, os
S=json.load(open("metric_sets.json"))
SUT=set(n for cfg in json.load(open("suture_config.json")).values() for pr in cfg for n in pr)
MOUTH,RIM,CORE=set(S["mouth"])-SUT,set(S["rim"]),set(S["core"])
base=open("comp_intact_pure.feb",encoding="ISO-8859-1").read()
nsec=base[base.index("<Nodes"):base.index("</Nodes>")]
init={int(m.group(1)):tuple(float(v) for v in m.group(2).split(',')) for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',nsec) if int(m.group(1)) in MOUTH}
cx,cy=-9.85,-17.12
mx=sum(p[0] for p in init.values())/len(init); my=sum(p[1] for p in init.values())/len(init)
nx,ny=mx-cx,my-cy; Ln=math.hypot(nx,ny); nout=(nx/Ln,ny/Ln)
def vm(s):
    sx,sy,sz,sxy,syz,sxz=s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def F_of_t(t): return max(0.0,2000.0*(t-1)/9.0)
def series(path):
    if not os.path.exists(path): return None
    txt=open(path,encoding="ISO-8859-1",errors="ignore").read(); rec={}
    for blk in re.split(r'Data Record #\d+',txt):
        mt=re.search(r'Time\s*=\s*([\d.eE+-]+)',blk); md=re.search(r'Data\s*=\s*([a-zA-Z;]+)',blk)
        if not(mt and md): continue
        t=float(mt.group(1)); d=rec.setdefault(t,{"pos":{},"rim":{},"core":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$',blk,re.M):
            i=int(lm.group(1)); v=[float(x) for x in lm.group(2).split()]
            if md.group(1)=="x;y;z" and i in MOUTH and len(v)>=3: d["pos"][i]=v[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and i in RIM and len(v)>=6: d["rim"][i]=v[:6]
            elif md.group(1)=="sx;sy;sz" and i in CORE and len(v)>=3: d["core"][i]=v[:3]
    ts=sorted(rec); t0=min(ts,key=lambda t:abs(t-1.0)); ref=rec[t0]["pos"]
    out=[]
    for t in ts:
        if t<t0-1e-9: continue
        d=rec[t]
        if not d["pos"]: continue
        pr=[(d["pos"][i][0]-ref[i][0])*nout[0]+(d["pos"][i][1]-ref[i][1])*nout[1] for i in d["pos"] if i in ref]
        prot=sum(pr)/len(pr) if pr else 0.0
        rs=max((vm(s) for s in d["rim"].values()),default=0.0)
        idp=-sum((s[0]+s[1]+s[2])/3 for s in d["core"].values())/len(d["core"]) if d["core"] else 0.0
        out.append((F_of_t(t),prot,rs,idp))
    return out
def metrics(s):
    F=[x[0] for x in s]; P=[x[1] for x in s]
    end=s[-1]; ipk=max(range(len(P)),key=lambda i:P[i])
    auc=sum((P[i]+P[i-1])/2*(F[i]-F[i-1]) for i in range(1,len(s)))
    return dict(endpoint=end[1],peak=P[ipk],load_at_peak=F[ipk],auc=auc,rim=end[2],idp=end[3])

STATES=[("intact","comp_intact_%s.log"),("defect","comp_defect_%s.log"),
        ("bridge 1-suture","comp_suture1_%s.log"),("bridge 3-suture","comp_suture3_%s.log"),
        ("apposed 1-suture","comp_suture1_apposed_%s.log"),("apposed 3-suture","comp_suture3_apposed_%s.log")]
allres={}
for scen in ("pure","flex"):
    print("="*100); print(f"SCENARIO: {scen} compression @ 2000 N")
    print(f"{'state':18s} {'endpt(mm)':>10s} {'peak(mm)':>9s} {'load@pk(N)':>10s} {'AUC(mm·N)':>10s} {'rimVM(MPa)':>10s} {'IDP(MPa)':>9s} {'%rest(end/AUC)':>14s}")
    print("-"*100)
    M={}
    for name,patt in STATES:
        s=series(patt%scen)
        if s is None: print(f"{name:18s}   (pending)"); continue
        M[name]=metrics(s); allres[(scen,name)]=M[name]
    di,dd=M.get("intact"),M.get("defect")
    for name,_ in STATES:
        if name not in M: continue
        m=M[name]
        if di and dd and dd["endpoint"]!=di["endpoint"]:
            re_end=100*(dd["endpoint"]-m["endpoint"])/(dd["endpoint"]-di["endpoint"])
            re_auc=100*(dd["auc"]-m["auc"])/(dd["auc"]-di["auc"]) if dd["auc"]!=di["auc"] else float('nan')
            rest=f"{re_end:.0f}/{re_auc:.0f}"
        else: rest="—"
        print(f"{name:18s} {m['endpoint']:10.3f} {m['peak']:9.3f} {m['load_at_peak']:10.0f} {m['auc']:10.1f} {m['rim']:10.2f} {m['idp']:9.3f} {rest:>14s}")
json.dump({f"{s}|{n}":v for (s,n),v in allres.items()}, open("apposed_metrics.json","w"), indent=2)
print("\n-> apposed_metrics.json")
