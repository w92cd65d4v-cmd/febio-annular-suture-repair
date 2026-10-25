#!/usr/bin/env python3
"""Closure-degree dose-response for the 3-SUTURE apposed repair.
States: intact, defect, open bridge, app12/app10b/app08/app06 (mean residual ~1.2/1.0/0.8/0.6 mm)
and the near-contact reference. Extended metrics + full protrusion-load curves, and an
automatic pick of the best 'physiologic apposed repair' (non-negative endpoint, restoration
not far >100%, large positive-AUC reduction vs bridge, better mesh quality than near)."""
import re, math, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
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
        out.append((F_of_t(t), sum(pr)/len(pr) if pr else 0.0,
                    max((vm(s) for s in d["rim"].values()),default=0.0),
                    -sum((s[0]+s[1]+s[2])/3 for s in d["core"].values())/len(d["core"]) if d["core"] else 0.0))
    return out
Q={d["tag"]:d for d in json.load(open("apposed3_levels_quality.json"))}
STATES=[("intact","comp_intact_%s.log","#444",None),("defect","comp_defect_%s.log","#A32D2D",None),
        ("bridge","comp_suture3_%s.log","#1D9E75",None),
        ("app17","comp_suture3_app17_%s.log","#C0DD97","app17"),("app15","comp_suture3_app15_%s.log","#97C459","app15"),
        ("app12","comp_suture3_app12_%s.log","#9FE1CB","app12"),("app10b","comp_suture3_app10b_%s.log","#5DCAA5","app10b"),
        ("app08","comp_suture3_app08_%s.log","#378ADD","app08"),("app06","comp_suture3_app06_%s.log","#185FA5","app06"),
        ("near","comp_suture3_apposed_%s.log","#0C447C","near")]
metr={}
fig,axes=plt.subplots(1,2,figsize=(14,5.4))
for ax,scen in zip(axes,("pure","flex")):
    print("="*120); print(f"SCENARIO: {scen} compression @ 2000 N  (3-suture closure-degree dose-response)")
    print(f"{'state':8s} {'meanGap':>7s} {'endpt':>7s} {'peak':>6s} {'minProt':>8s} {'signAUC':>8s} {'posAUC':>7s} {'%rest(posAUC)':>13s} {'rimVM':>6s} {'IDP':>6s} {'minSJ':>6s} {'#sJ<.2':>6s}")
    print("-"*120)
    base_pos=None; di=dd=None
    rows={}
    for name,patt,col,qt in STATES:
        s=series(patt%scen)
        if s is None: print(f"{name:8s}  (pending)"); continue
        F=[x[0] for x in s]; P=[x[1] for x in s]
        sAUC=sum((P[i]+P[i-1])/2*(F[i]-F[i-1]) for i in range(1,len(s)))
        pAUC=sum((max(P[i],0)+max(P[i-1],0))/2*(F[i]-F[i-1]) for i in range(1,len(s)))
        rows[name]=dict(end=P[-1],peak=max(P),mn=min(P),sAUC=sAUC,pAUC=pAUC,rim=s[-1][2],idp=s[-1][3])
        metr[(scen,name)]=rows[name]
        ax.plot(F,P,'-',color=col,lw=2,label=name)
    di,dd=rows.get("intact"),rows.get("defect")
    for name,patt,col,qt in STATES:
        if name not in rows: continue
        r=rows[name]; q=Q.get(qt,{})
        mg=q.get("mean_residual"); sj=q.get("min_scaled_jac"); nb=q.get("below_0p2")
        rest=100*(dd["pAUC"]-r["pAUC"])/(dd["pAUC"]-di["pAUC"]) if di and dd and dd["pAUC"]!=di["pAUC"] else float('nan')
        print(f"{name:8s} {(f'{mg:.2f}' if mg else '—'):>7s} {r['end']:7.3f} {r['peak']:6.3f} {r['mn']:8.3f} {r['sAUC']:8.1f} {r['pAUC']:7.1f} {rest:12.0f}% {r['rim']:6.2f} {r['idp']:6.3f} {(f'{sj:.3f}' if sj is not None else '—'):>6s} {(str(nb) if nb is not None else '—'):>6s}")
    ax.axhline(0,color='#999',lw=.8,ls=':'); ax.set_title(f"{scen} compression"); ax.set_xlabel("axial force (N)")
    ax.set_ylabel("load-induced protrusion (mm)"); ax.grid(alpha=.25); ax.legend(fontsize=8,title="closure ↑")
fig.suptitle("3-suture closure-degree dose-response: protrusion–load curves\nopen bridge → apposed mean residual 1.2 → 1.0 → 0.8 → 0.6 mm / near-contact",y=1.04)
fig.tight_layout(); fig.savefig("../apposed_repair_check/closure_doseresponse3_curves.png",dpi=140,bbox_inches="tight")
print("\nwrote apposed_repair_check/closure_doseresponse3_curves.png")

# pick best physiologic apposed repair
print("\n"+"="*120); print("BEST PHYSIOLOGIC APPOSED REPAIR (3-suture) — criteria: endpoint>=0 (both scen), restoration<=~110%, large posAUC drop vs bridge, mesh better than near")
cand=["app17","app15","app12","app10b","app08","app06","near"]
br={sc:metr[(sc,"bridge")]["pAUC"] for sc in ("pure","flex")}
intp={sc:metr[(sc,"intact")]["pAUC"] for sc in ("pure","flex")}
defp={sc:metr[(sc,"defect")]["pAUC"] for sc in ("pure","flex")}
near_sj=Q.get("near",{}).get("min_scaled_jac",0)
best=None
for c in cand:
    if (("pure",c) not in metr) or (("flex",c) not in metr): continue
    end_ok=all(metr[(sc,c)]["end"]>=-0.005 for sc in ("pure","flex"))
    rest={sc:100*(defp[sc]-metr[(sc,c)]["pAUC"])/(defp[sc]-intp[sc]) for sc in ("pure","flex")}
    rest_ok=all(rest[sc]<=112 for sc in ("pure","flex"))
    drop={sc:100*(br[sc]-metr[(sc,c)]["pAUC"])/br[sc] for sc in ("pure","flex")}
    drop_ok=all(drop[sc]>=15 for sc in ("pure","flex"))
    q=Q.get(c if c!="near" else "near",{}); mesh_ok=(q.get("min_scaled_jac",0) > near_sj+0.02) and (q.get("below_0p2",99)<=1)
    ok=end_ok and rest_ok and drop_ok and mesh_ok
    print(f"  {c:8s} endpoint>=0:{end_ok}  rest pure/flex={rest['pure']:.0f}/{rest['flex']:.0f}% (<=112:{rest_ok})  posAUC drop vs bridge={drop['pure']:.0f}/{drop['flex']:.0f}% (>=15:{drop_ok})  mesh_ok:{mesh_ok}  => {'PASS' if ok else 'no'}")
    if ok and best is None: best=c
print(f"\n  RECOMMENDED physiologic apposed 3-suture model: {best}")
json.dump({"metrics":{f"{s}|{n}":v for (s,n),v in metr.items()},"recommended":best}, open("apposed3_doseresponse.json","w"),indent=2)
