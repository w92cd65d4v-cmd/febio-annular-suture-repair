#!/usr/bin/env python3
"""Suture-stiffness sweep: protrusion (load-induced, arrow metric) and peak rim
stress at 2000 N vs suture stiffness (10/50/200 N/mm), for 1- and 3-suture repairs,
pure and flex. Reports % of the defect's herniation each stiffness restores."""
import re, math, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

S=json.load(open("metric_sets.json"))
SUT=set(n for cfg in json.load(open("suture_config.json")).values() for pr in cfg for n in pr)
MOUTH,RIM=set(S["mouth"])-SUT,set(S["rim"])
base=open("comp_intact_pure.feb",encoding="ISO-8859-1").read()
nsec=base[base.index("<Nodes"):base.index("</Nodes>")]
init={int(m.group(1)):tuple(float(v) for v in m.group(2).split(',')) for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',nsec) if int(m.group(1)) in MOUTH}
cx,cy=-9.85,-17.12
mx=sum(p[0] for p in init.values())/len(init); my=sum(p[1] for p in init.values())/len(init)
nx,ny=mx-cx,my-cy; Ln=math.hypot(nx,ny); nout=(nx/Ln,ny/Ln)
def vm(s):
    sx,sy,sz,sxy,syz,sxz=s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def final(path):
    if not __import__("os").path.exists(path): return None,None
    txt=open(path,encoding="ISO-8859-1",errors="ignore").read(); rec={}
    for blk in re.split(r'Data Record #\d+',txt):
        mt=re.search(r'Time\s*=\s*([\d.eE+-]+)',blk); md=re.search(r'Data\s*=\s*([a-zA-Z;]+)',blk)
        if not(mt and md): continue
        t=float(mt.group(1)); d=rec.setdefault(t,{"pos":{},"rim":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$',blk,re.M):
            i=int(lm.group(1)); v=[float(x) for x in lm.group(2).split()]
            if md.group(1)=="x;y;z" and i in MOUTH and len(v)>=3: d["pos"][i]=v[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and i in RIM and len(v)>=6: d["rim"][i]=v[:6]
    ts=sorted(rec); t0=min(ts,key=lambda t:abs(t-1.0)); t1=max(ts)
    ref=rec[t0]["pos"]
    pr=[(rec[t1]["pos"][i][0]-ref[i][0])*nout[0]+(rec[t1]["pos"][i][1]-ref[i][1])*nout[1] for i in rec[t1]["pos"] if i in ref]
    prot=sum(pr)/len(pr) if pr else None
    rs=max((vm(s) for s in rec[t1]["rim"].values()),default=None)
    return prot,rs

print(f"{'scen':5s} {'cfg':3s} {'E(N/mm)':>7s} {'prot(mm)':>9s} {'rim(MPa)':>9s} {'%restore':>9s}")
print("-"*52)
results={}
for scen in ("pure","flex"):
    ip,irs=final(f"comp_intact_{scen}.log"); dp,drs=final(f"comp_defect_{scen}.log")
    for cfg,tag in (("suture1","s1"),("suture3","s3")):
        for E in (10,50,200):
            path = f"comp_{cfg}_{scen}.log" if E==50 else f"comp_{tag}_E{E}_{scen}.log"
            p,rs=final(path)
            if p is None: print(f"{scen:5s} {tag:3s} {E:7d}   (pending)"); continue
            rest=100*(dp-p)/(dp-ip) if dp!=ip else float('nan')
            results.setdefault((scen,tag),[]).append((E,p,rs,rest))
            print(f"{scen:5s} {tag:3s} {E:7d} {p:9.3f} {rs:9.3f} {rest:8.1f}%")
    results[(scen,'_base')]=(ip,dp,irs,drs)

# plot restoration vs stiffness
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
for ax,scen in zip(axes,("pure","flex")):
    ip,dp,_,_=results[(scen,'_base')]
    for tag,col in (("s1","#1D9E75"),("s3","#185FA5")):
        pts=sorted(results.get((scen,tag),[]))
        if pts:
            E=[p[0] for p in pts]; rest=[p[3] for p in pts]
            ax.plot(E,rest,'-o',color=col,label=f"{tag} ({'1' if tag=='s1' else '3'} suture{'s' if tag=='s3' else ''})")
    ax.axhline(100,ls=':',c='#444'); ax.text(11,101,"full restoration (=intact)",fontsize=8)
    ax.set_xscale('log'); ax.set_xticks([10,50,200]); ax.set_xticklabels(['10','50','200'])
    ax.set_xlabel("suture stiffness (N/mm)"); ax.set_ylabel("% contained protrusion restored toward intact")
    ax.set_title(f"{scen} compression"); ax.grid(alpha=.25); ax.legend(fontsize=9)
fig.suptitle("Suture-stiffness sweep: how much of the defect's load-induced contained protrusion is restored",y=1.02)
fig.tight_layout(); fig.savefig("sweep_restoration.png",dpi=140,bbox_inches="tight"); print("\nwrote sweep_restoration.png")
