#!/usr/bin/env python3
"""Pre-tension (closure-like apposition) sensitivity: compare the baseline no-preload
3-suture bridge with a pre-tensioned (0.2 mm shortened rest length, ~10 N pre-load)
3-suture bridge. Reports, on the nucleus mouth nodes (suture anchors excluded):
 - apposition during preload: outward position at F=0 vs the undeformed mesh
   (negative = edges drawn inward = closure-like apposition);
 - load-induced contained protrusion at 2000 N (arrow metric, referenced to F=0);
 - total outward position at 2000 N vs undeformed.
"""
import re, math, json, os
S=json.load(open("../validation_L45_defect/metric_sets.json"))
SUT=set(n for cfg in json.load(open("suture_config.json")).values() for pr in cfg for n in pr)
MOUTH=set(S["mouth"])-SUT
base=open("comp_intact_pure.feb",encoding="ISO-8859-1").read()
nsec=base[base.index("<Nodes"):base.index("</Nodes>")]
init={int(m.group(1)):tuple(float(v) for v in m.group(2).split(',')) for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',nsec) if int(m.group(1)) in MOUTH}
cx,cy=-9.85,-17.12
mx=sum(p[0] for p in init.values())/len(init); my=sum(p[1] for p in init.values())/len(init)
nx,ny=mx-cx,my-cy; Ln=math.hypot(nx,ny); nout=(nx/Ln,ny/Ln)
def positions(path):
    if not os.path.exists(path): return None
    txt=open(path,encoding="ISO-8859-1",errors="ignore").read(); rec={}
    for blk in re.split(r'Data Record #\d+',txt):
        mt=re.search(r'Time\s*=\s*([\d.eE+-]+)',blk); md=re.search(r'Data\s*=\s*([a-zA-Z;]+)',blk)
        if not(mt and md) or md.group(1)!="x;y;z": continue
        t=float(mt.group(1)); d=rec.setdefault(t,{})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$',blk,re.M):
            i=int(lm.group(1))
            if i in MOUTH: d[i]=tuple(float(x) for x in lm.group(2).split()[:3])
    return rec
def outward(pos, ref):  # mean outward (posterolateral) offset of mouth vs ref
    v=[(pos[i][0]-ref[i][0])*nout[0]+(pos[i][1]-ref[i][1])*nout[1] for i in pos if i in ref]
    return sum(v)/len(v) if v else None
print(f"{'scen':5s} {'model':16s} {'apposition@F0':>13s} {'loadProt@2kN':>13s} {'totalPos@2kN':>13s} {'%restore(arrow)':>15s}")
print("-"*82)
out={}
for scen in ("pure","flex"):
    intact=positions(f"comp_intact_{scen}.log"); defect=positions(f"comp_defect_{scen}.log")
    def metrics(rec):
        ts=sorted(rec); t0=min(ts,key=lambda t:abs(t-1.0)); t1=max(ts)
        f0=rec[t0]; f1=rec[t1]
        appos=outward(f0,init)              # F=0 vs undeformed
        total=outward(f1,init)              # 2000N vs undeformed
        load=outward(f1,f0)                 # arrow metric (vs F=0)
        return appos,total,load
    di=metrics(intact); dd=metrics(defect)
    for model,tag in (("defect",None),("3 sutures","comp_suture3_%s.log"),("3 sutures + pre-tension","comp_s3pre_%s.log")):
        rec = defect if model=="defect" else positions(tag%scen)
        if rec is None: print(f"{scen:5s} {model:16s}   (pending)"); continue
        ap,tot,load = metrics(rec)
        rest = 100*(dd[2]-load)/(dd[2]-di[2]) if dd[2]!=di[2] else float('nan')
        out[(scen,model)]=(ap,tot,load,rest)
        print(f"{scen:5s} {model:16s} {ap:13.3f} {load:13.3f} {tot:13.3f} {rest:14.1f}%")
