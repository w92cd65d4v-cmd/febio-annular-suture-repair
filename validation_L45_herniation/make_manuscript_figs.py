#!/usr/bin/env python3
"""Final manuscript figures for the apposed-anchor repair story.
Reads the validated cc01 dataset (logs + .feb + quality json) and writes figures into the
manuscript's home (validation_L45_herniation/). Terminology: 'contained protrusion'.
Figures:
  repair_geometry.png            open defect | bridge | apposed app12 | apposed near (upper bound)
  contained_protrusion_curves.png  intact/defect/bridge/app15/app12/near, pure & flex
  closure_doseresponse_summary.png mean residual gap vs peak protrusion, positive-AUC, mesh quality
  rim_stress_persistence.png       rim von Mises across states (riser persists despite protrusion control)
"""
import re, json, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patheffects as pe

CC="/Volumes/storage01/forCC/cc01/validation_L45_herniation"

if not os.path.isdir(CC): CC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nas_copy")  # local copy 2026-09-27
OUT="."
S=json.load(open(f"{CC}/metric_sets.json"))
SUT=set(n for cfg in json.load(open(f"{CC}/suture_config.json")).values() for pr in cfg for n in pr)
MOUTH,RIM,CORE=set(S["mouth"])-SUT,set(S["rim"]),set(S["core"])
PAIRS=json.load(open(f"{CC}/suture_config.json"))["suture3"]
QJ={d["tag"]:d for d in json.load(open(f"{CC}/apposed3_levels_quality.json"))}
def nodes_of(p):
    t=open(p,encoding="ISO-8859-1").read(); s=t[t.index("<Nodes"):t.index("</Nodes>")]
    return {int(m.group(1)):[float(v) for v in m.group(2).split(',')] for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',s)}
init={k:tuple(v) for k,v in nodes_of(f"{CC}/comp_intact_pure.feb").items() if k in MOUTH}
cx,cy=-9.85,-17.12
mx=sum(p[0] for p in init.values())/len(init); my=sum(p[1] for p in init.values())/len(init)
nx,ny=mx-cx,my-cy; Ln=math.hypot(nx,ny); nout=(nx/Ln,ny/Ln)
def vm(s):
    sx,sy,sz,sxy,syz,sxz=s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def F_of_t(t): return max(0.0,2000.0*(t-1)/9.0)
def series(tag,scen):
    fn={"intact":"comp_intact_%s.log","defect":"comp_defect_%s.log","bridge":"comp_suture3_%s.log",
        "near":"comp_suture3_apposed_%s.log"}.get(tag, f"comp_suture3_{tag}_%s.log")
    path=f"{CC}/{fn%scen}"
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
def met(tag,scen):
    s=series(tag,scen)
    if not s: return None
    P=[x[1] for x in s]; F=[x[0] for x in s]
    pAUC=sum((max(P[i],0)+max(P[i-1],0))/2*(F[i]-F[i-1]) for i in range(1,len(s)))
    return dict(end=P[-1],peak=max(P),pAUC=pAUC,rim=s[-1][2],idp=s[-1][3])

# ---- Figure 1: geometry (4 panels) ----
base=open(f"{CC}/comp_suture3_pure.feb",encoding="ISO-8859-1").read()
def elems_of(name):
    mm=re.search(r'<Elements type="hex8" mat="\d+" name="'+name+r'">(.*?)</Elements>',base,re.S)
    return {int(e.group(1)):[int(v) for v in e.group(2).split(',')] for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>',mm.group(1))}
conn={}
for q in ["ANT","L","POST","R"]:
    for lam in range(1,6): conn.update(elems_of(f"L45_{q}_LAM{lam}"))
conn.update(elems_of("L45_NUCLEUS"))
RS=json.load(open(f"{CC}/region_sets.json")); REG_E,DEFECT,C=set(RS["elems"]),set(RS["defect"]),RS["C"]
zc=C[2]
def cz(e,N): return sum(N[n][2] for n in conn[e])/8
def face(e,N): return [(N[n][0],N[n][1]) for n in conn[e][:4]]
N_orig=nodes_of(f"{CC}/comp_suture3_pure.feb")
slab=[e for e in REG_E if e in conn and abs(cz(e,N_orig)-zc)<=1.4]; hole=[e for e in DEFECT if e in conn and abs(cz(e,N_orig)-zc)<=1.4]
panels=[("open defect (no repair)",N_orig,False,None),
        ("open-gap bridge repair",N_orig,True,None),
        ("apposed app12 (main)\nmean gap ≈ 1.19 mm",nodes_of(f"{CC}/comp_suture3_app12_pure.feb"),True,"app12"),
        ("apposed near (upper bound)\nmean gap ≈ 0.63 mm",nodes_of(f"{CC}/comp_suture3_apposed_pure.feb"),True,"near")]
fig,axes=plt.subplots(1,4,figsize=(18,5.2))
for ax,(title,N,springs,qt) in zip(axes,panels):
    ax.add_collection(PolyCollection([face(e,N) for e in slab if e not in DEFECT],facecolors='#cdd7e5',edgecolors='#7c8aa0',linewidths=0.4))
    for e in hole:
        f=face(e,N)+[face(e,N)[0]]; ax.plot([p[0] for p in f],[p[1] for p in f],color='magenta',lw=1.1)
    if springs:
        for a,b in PAIRS:
            pa,pb=N[a],N[b]
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],color='yellow',lw=2.4,solid_capstyle='round',path_effects=[pe.Stroke(linewidth=4,foreground='black'),pe.Normal()],zorder=6)
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],'o',color='yellow',mec='black',ms=5,zorder=7)
        if qt:
            for a,b in PAIRS:
                g=math.dist(N[a],N[b]); ax.annotate(f"{g:.2f}",((N[a][0]+N[b][0])/2,(N[a][1]+N[b][1])/2),fontsize=7,bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=.7))
    ax.set_xlim(-16,2); ax.set_ylim(-15,0); ax.set_aspect('equal'); ax.set_title(title,fontsize=10); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (+post)")
fig.suptitle("Repair-model geometry: open annulotomy → tension-only suture bridge → apposed-anchor repair (3-suture, axial section)\n"
             "magenta = retained annulotomy (elements still deleted, no tissue restored); yellow = tension-only suture springs; numbers = residual anchor-pair gap (mm)",y=1.05,fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/repair_geometry.png",dpi=120,bbox_inches="tight"); print("wrote repair_geometry.png")

# ---- Figure 2: contained protrusion-load curves ----
CURVE=[("intact","#444"),("defect","#A32D2D"),("bridge","#1D9E75"),("app15","#97C459"),("app12","#185FA5"),("near","#0C447C")]
LBL={"app15":"apposed app15 (1.49 mm)","app12":"apposed app12 (1.19 mm)","near":"apposed near (0.63 mm)","bridge":"open-gap bridge"}
fig,axes=plt.subplots(1,2,figsize=(14,5.2))
for ax,scen in zip(axes,("pure","flex")):
    for tag,col in CURVE:
        s=series(tag,scen)
        if s: ax.plot([x[0] for x in s],[x[1] for x in s],'-',color=col,lw=2,label=LBL.get(tag,tag))
    ax.axhline(0,color='#999',lw=.8,ls=':'); ax.set_title(f"{scen} compression"); ax.set_xlabel("axial force (N)")
    ax.set_ylabel("contained protrusion (mm)"); ax.grid(alpha=.25); ax.legend(fontsize=8)
fig.suptitle("Contained protrusion vs axial load: intact, defect, open-gap bridge, and apposed-anchor repair levels",y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/contained_protrusion_curves.png",dpi=140,bbox_inches="tight"); print("wrote contained_protrusion_curves.png")

# ---- Figure 3: dose-response summary vs mean residual gap ----
bridge_gap=sum(math.dist(N_orig[a],N_orig[b]) for a,b in PAIRS)/3
LEVELS=[("bridge",bridge_gap,False),("app17",QJ.get("app17",{}).get("mean_residual"),False),("app15",QJ.get("app15",{}).get("mean_residual"),False),
        ("app12",QJ.get("app12",{}).get("mean_residual"),False),("app10b",QJ.get("app10b",{}).get("mean_residual"),False),
        ("app08",QJ.get("app08",{}).get("mean_residual"),False),("app06",QJ.get("app06",{}).get("mean_residual"),True),("near",QJ.get("near",{}).get("mean_residual"),True)]
fig,ax=plt.subplots(1,3,figsize=(16,4.8))
for sc,mk in (("pure","o"),("flex","s")):
    gp=[];pk=[];au=[];ms=[]
    for tag,g,over in LEVELS:
        m=met(tag,sc);
        if not m or g is None: continue
        gp.append(g); pk.append(m["peak"]); au.append(m["pAUC"]); ms.append(QJ.get(tag,{}).get("min_scaled_jac",None))
    ax[0].plot(gp,pk,'-'+mk,label=sc,color="#185FA5" if sc=="pure" else "#A32D2D")
    ax[1].plot(gp,au,'-'+mk,label=sc,color="#185FA5" if sc=="pure" else "#A32D2D")
# mesh quality (geometry only; same for pure/flex)
gq=[(QJ.get(t,{}).get("mean_residual"),QJ.get(t,{}).get("min_scaled_jac")) for t,_,_ in LEVELS if t!="bridge" and QJ.get(t,{}).get("min_scaled_jac") is not None]
ax[2].plot([x[0] for x in gq],[x[1] for x in gq],'-o',color="#0F6E56")
ax[2].axhline(0.2,ls=':',c='#888'); ax[2].text(0.65,0.21,"quality floor (0.2)",fontsize=8)
for a in ax[:2]:
    for tag,g,over in LEVELS:
        if over and g is not None: a.axvspan(g-0.03,g+0.03,color='#A32D2D',alpha=0.06)
ax[0].set_ylabel("peak contained protrusion (mm)"); ax[1].set_ylabel("positive-only AUC (mm·N)"); ax[2].set_ylabel("min scaled Jacobian")
for a in ax: a.set_xlabel("mean residual anchor gap (mm)"); a.grid(alpha=.25); a.invert_xaxis()
ax[0].legend(title="↓ closure →",fontsize=8); ax[1].legend(fontsize=8)
ax[0].set_title("peak protrusion"); ax[1].set_title("positive-only AUC"); ax[2].set_title("mesh quality (over-tightened ≤ ~0.6 mm)")
fig.suptitle("Closure-degree dose-response (3-suture): more apposition (smaller residual gap) → less peak & cumulative protrusion;\nmesh quality collapses at near-contact (≤ ~0.6 mm) = over-tightened upper bound. app12 (~1.19 mm) is the balanced model.",y=1.06,fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/closure_doseresponse_summary.png",dpi=140,bbox_inches="tight"); print("wrote closure_doseresponse_summary.png")

# ---- Figure 4: rim stress persistence ----
RST=[("intact","#444"),("defect","#A32D2D"),("bridge","#1D9E75"),("app15","#97C459"),("app12","#185FA5"),("near","#0C447C")]
fig,axes=plt.subplots(1,2,figsize=(13,4.6))
for ax,scen in zip(axes,("pure","flex")):
    names=[t for t,_ in RST]; vals=[met(t,scen)["rim"] for t,_ in RST]; cols=[c for _,c in RST]
    ax.bar(range(len(names)),vals,color=cols)
    ax.axhline(vals[0],ls='--',color='#444',lw=1); ax.text(len(names)-1,vals[0],"intact",va='bottom',ha='right',fontsize=8,color='#444')
    ax.set_xticks(range(len(names))); ax.set_xticklabels([n.replace('app','app ') for n in names],fontsize=8,rotation=20)
    ax.set_title(f"{scen} compression"); ax.set_ylabel("peak rim von Mises (MPa)"); ax.grid(axis='y',alpha=.25)
    for i,v in enumerate(vals): ax.text(i,v,f"{v:.2f}",ha='center',va='bottom',fontsize=7)
fig.suptitle("Rim stress persists despite protrusion control: every repair (bridge and apposed) keeps peak rim von Mises at or above the defect,\nwell above intact — protrusion control and rim-stress relief are distinct endpoints.",y=1.04,fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/rim_stress_persistence.png",dpi=140,bbox_inches="tight"); print("wrote rim_stress_persistence.png")
print("done")
