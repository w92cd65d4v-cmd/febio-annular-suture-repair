#!/usr/bin/env python3
"""Primary quantitative Results figures R1, R2, R4, and R5.
Reads verified metrics from the FEBio .log files (local + cc01). No simulation reruns.
Exports vector (PDF, SVG) + raster (PNG) for each figure.
Metric definitions match analyze_levels3.py / path_metrics.py exactly."""
import re, math, json, os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"font.family":"DejaVu Sans","axes.linewidth":0.8,
                     "svg.fonttype":"none","pdf.fonttype":42})
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","results_figures")
os.makedirs(OUT,exist_ok=True)
CC="/Volumes/storage01/forCC/cc01/validation_L45_herniation"
if not os.path.isdir(CC): CC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nas_copy")  # local copy 2026-09-27
S=json.load(open(os.path.join(HERE,"metric_sets.json")))
SUT=set(n for cfg in json.load(open(os.path.join(HERE,"suture_config.json"))).values() for pr in cfg for n in pr)
MOUTH,RIM,CORE=set(S["mouth"])-SUT,set(S["rim"]),set(S["core"])
base=open(os.path.join(HERE,"comp_intact_pure.feb"),encoding="ISO-8859-1").read()
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
        t=float(mt.group(1)); d=rec.setdefault(t,{"pos":{},"rim":{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$',blk,re.M):
            i=int(lm.group(1)); v=[float(x) for x in lm.group(2).split()]
            if md.group(1)=="x;y;z" and i in MOUTH and len(v)>=3: d["pos"][i]=v[:3]
            elif md.group(1).startswith("sx;sy;sz;sxy") and i in RIM and len(v)>=6: d["rim"][i]=v[:6]
    ts=sorted(rec); t0=min(ts,key=lambda t:abs(t-1.0)); ref=rec[t0]["pos"]
    out=[]
    for t in ts:
        if t<t0-1e-9: continue
        d=rec[t]
        if not d["pos"]: continue
        pr=[(d["pos"][i][0]-ref[i][0])*nout[0]+(d["pos"][i][1]-ref[i][1])*nout[1] for i in d["pos"] if i in ref]
        out.append((F_of_t(t),sum(pr)/len(pr) if pr else 0.0,max((vm(s) for s in d["rim"].values()),default=0.0)))
    return out
def logp(state,scen):
    loc={"intact":"comp_intact_%s.log","defect":"comp_defect_%s.log",
         "suture1":"comp_suture1_%s.log","bridge":"comp_suture3_%s.log",
         "app12":"comp_suture3_app12_%s_map.log"}
    if state in loc: return os.path.join(HERE,loc[state]%scen)
    m={"near":"comp_suture3_apposed_%s.log"}.get(state,f"comp_suture3_{state}_%s.log")
    return f"{CC}/{m%scen}"
COL={'intact':'#555555','defect':'#B0281A','suture1':'#D98F39','bridge':'#1F9E77','app17':'#bcd98f',
     'app15':'#8FBF4C','app12':'#2E6FB0','app10b':'#57b6c4','app08':'#3A8FD0',
     'app06':'#12417A','near':'#0A2E55'}
GAP={"app17":1.690,"app15":1.485,"app12":1.186,"app10b":0.982,"app08":0.796,"app06":0.591,"near":0.628}
SJ={"app17":0.2196,"app15":0.2196,"app12":0.2196,"app10b":0.2196,"app08":0.2196,"app06":0.0359,"near":0.0720}
def save(fig,name):
    for ext in ("pdf","svg","png"):
        fig.savefig(os.path.join(OUT,f"{name}.{ext}"),dpi=600,bbox_inches="tight")
    fig.savefig(os.path.join(OUT,f"{name}.tiff"),dpi=600,bbox_inches="tight",
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig); print("wrote",name,"(pdf/svg/png/tiff)")

# ---------- R1: contained protrusion-load ----------
def R1():
    states=["intact","defect","suture1","bridge","app12"]
    lab={"intact":"intact","defect":"unrepaired defect",
         "suture1":"load-targeted open-gap 1-suture bridge",
         "bridge":"load-targeted open-gap 3-suture bridge",
         "app12":"3-suture partial apposition (app12; gap 1.19 mm)"}
    ls={"intact":"-","defect":"-","suture1":"-","bridge":"-","app12":"-"}
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
    for ax,scen,tt in zip(axes,("pure","flex"),("pure axial compression","flexion + axial compression")):
        for st in states:
            r=series(logp(st,scen))
            if not r: continue
            F=[x[0] for x in r]; P=[x[1] for x in r]
            ax.plot(F,P,ls[st],color=COL[st],lw=2.2,label=lab[st],zorder=3)
            pk=max(range(len(P)),key=lambda i:P[i])
            ax.plot(F[pk],P[pk],'o',color=COL[st],ms=5,mec='white',mew=0.8,zorder=4)
            ax.plot(F[-1],P[-1],'s',color=COL[st],ms=4.5,mec='white',mew=0.7,zorder=4)
        ax.axhline(0,color='#999',lw=0.8,ls=':')
        ax.set_title(tt,fontsize=11); ax.set_xlabel("axial compressive force (N)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("load-induced contained protrusion (mm)")
    axes[0].legend(fontsize=8,loc="upper left",framealpha=0.9)
    axes[1].text(0.98,0.03,"● = path maximum; ■ = 2000-N endpoint\npeak + positive-only AUC are primary because paths are non-monotonic",
                 transform=axes[1].transAxes,fontsize=7.5,ha="right",va="bottom",color="#444",style="italic")
    fig.suptitle("Contained protrusion–load response",fontsize=12,x=0.09,ha="left")
    fig.tight_layout(); save(fig,"figR1_protrusion_load")

# ---------- R2: residual-gap dose-response ----------
def R2():
    order=["app17","app15","app12","app10b","app08"]  # clean monotonic branch (>=0.8mm)
    ub=["app06","near"]  # upper-bound distorted
    MT={}
    for scen in ("pure","flex"):
        for st in order+ub+["intact","defect","bridge"]:
            r=series(logp(st,scen)); P=[x[1] for x in r]; F=[x[0] for x in r]
            pAUC=sum((max(P[i],0)+max(P[i-1],0))/2*(F[i]-F[i-1]) for i in range(1,len(r)))
            MT[(scen,st)]=dict(peak=max(P),pAUC=pAUC)
    fig,axes=plt.subplots(1,3,figsize=(15,4.6))
    scn_style={"pure":("-","o"),"flex":("--","s")}
    for scen in ("pure","flex"):
        ls,mk=scn_style[scen]
        gx=[GAP[s] for s in order];
        axes[0].plot(gx,[MT[(scen,s)]["peak"] for s in order],ls,marker=mk,color="#2E6FB0",label=f"{scen}")
        axes[1].plot(gx,[MT[(scen,s)]["pAUC"] for s in order],ls,marker=mk,color="#2E6FB0",label=f"{scen}")
        # upper-bound points (open markers)
        for s in ub:
            axes[0].plot(GAP[s],MT[(scen,s)]["peak"],mk,mfc='none',mec="#B0281A",ms=7)
            axes[1].plot(GAP[s],MT[(scen,s)]["pAUC"],mk,mfc='none',mec="#B0281A",ms=7)
        # intact + bridge references
        axes[0].axhline(MT[(scen,"intact")]["peak"],color=COL['intact'],lw=1,ls=":",alpha=0.7)
        axes[1].axhline(MT[(scen,"intact")]["pAUC"],color=COL['intact'],lw=1,ls=":",alpha=0.7)
    # mesh quality panel
    gx=[GAP[s] for s in order+ub]
    axes[2].plot(gx,[SJ[s] for s in order+ub],'o',color="#444")
    axes[2].axhline(0.2,color="#B0281A",ls="--",lw=1); axes[2].text(1.7,0.21,"0.2 screening threshold",fontsize=8,color="#B0281A")
    # annotate app12 + upper-bound on each
    for ax,key,yl in ((axes[0],"peak","peak contained protrusion (mm)"),(axes[1],"pAUC","positive-only AUC (mm·N)"),(axes[2],None,"min scaled Jacobian")):
        ax.set_xlabel("mean residual anchor-pair gap (mm)"); ax.set_ylabel(yl); ax.grid(alpha=0.25); ax.invert_xaxis()
        if key:
            ax.annotate("app12 (representative)",xy=(GAP["app12"],MT[("pure","app12")][key]),
                        xytext=(-6,-30),textcoords="offset points",fontsize=7.5,color="#2E6FB0",ha="right",
                        arrowprops=dict(arrowstyle="->",color="#2E6FB0"))
    axes[0].text(0.02,0.03,"open-gap bridge (no imposed gap),\noff-scale: pure "+f'{MT[("pure","bridge")]["peak"]:.3f}'+" / flex "+f'{MT[("flex","bridge")]["peak"]:.3f}'+" mm peak",
                 transform=axes[0].transAxes,fontsize=6.8,color="#1F9E77",va="bottom")
    axes[2].text(0.03,0.42,"app06 / near =\ndistorted upper-bound\nsensitivity models",transform=axes[2].transAxes,fontsize=7.5,color="#B0281A")
    axes[0].legend(title="loading",fontsize=8); axes[1].legend(title="loading",fontsize=8)
    fig.suptitle("Residual-gap dose–response for the 3-suture apposition series",fontsize=11,x=0.09,ha="left")
    fig.tight_layout(); save(fig,"figR2_dose_response")

# ---------- R4: persistent rim stress ----------
def R4():
    states=["intact","defect","suture1","bridge","app15","app12","app06","near"]
    lab={"intact":"intact","defect":"defect","suture1":"1 suture","bridge":"3 sutures","app15":"app15","app12":"app12 (repr.)",
         "app06":"app06 (UB)","near":"near (UB)"}
    with open(os.path.join(HERE,"..","RESULTS_MASTER_TABLE.csv"),newline="") as fh:
        master={(r["loading"],r["state"]):r for r in csv.DictReader(fh)}
    with open(os.path.join(HERE,"..","RESULTS_SECONDARY_TABLE.csv"),newline="") as fh:
        secondary=list(csv.DictReader(fh))

    def rim_value(scen,state):
        if state=="suture1":
            row=next(r for r in secondary if r["loading"]==scen
                     and r["sensitivity_group"]=="suture number"
                     and r["variant"].startswith("1 suture"))
            return float(row["rim_VM_MPa"])
        return float(master[(scen,state)]["rim_VM_peak_MPa"])

    fig,ax=plt.subplots(figsize=(9,4.8))
    x=range(len(states)); w=0.36
    for k,(scen,off,mk) in enumerate((("pure",-w/2,"o"),("flex",w/2,"s"))):
        vals=[rim_value(scen,st) for st in states]
        ax.bar([i+off for i in x],vals,w,label=f"{scen}",
               color=["#888" if st=="intact" else ("#B0281A" if st=="defect" else ("#c9a227" if "UB" in lab[st] else "#2E6FB0")) for st in states],
               edgecolor="black",lw=0.4,alpha=0.6 if scen=="flex" else 0.95)
        for i,v in zip(x,vals): ax.text(i+off,v+0.03,f"{v:.2f}",ha="center",fontsize=6.5)
    ax.axhline(rim_value("pure","intact"),color="#888",ls=":",lw=1)
    ax.set_xticks(list(x)); ax.set_xticklabels([lab[s] for s in states],fontsize=8,rotation=15)
    ax.set_ylabel("peak annular-rim von Mises (MPa)"); ax.grid(axis='y',alpha=0.25)
    ax.set_ylim(0,6.3); ax.legend(title="loading (solid=pure, faded=flex)",fontsize=8,loc="lower left")
    ax.set_title("Annular-rim stress was not restored to intact levels after repair\n(UB = distorted upper-bound sensitivity model)",fontsize=10.5,loc="left")
    fig.tight_layout(); save(fig,"figR4_rim_stress")

# ---------- R5: suture-number and stiffness sensitivities ----------
def R5():
    secondary_path=os.path.join(HERE,"..","RESULTS_SECONDARY_TABLE.csv")
    with open(secondary_path,newline="") as fh:
        rows=list(csv.DictReader(fh))

    def row_for(loading,group,needle=None,stiffness=None):
        for row in rows:
            if row["loading"]!=loading or row["sensitivity_group"]!=group:
                continue
            variant=row["variant"]
            if needle is not None and not variant.startswith(f"{needle} suture") and not variant.startswith(f"{needle} sutures"):
                continue
            if stiffness is not None and f"E{stiffness}" not in variant:
                continue
            return row
        raise KeyError((loading,group,needle,stiffness))

    def value(row,key): return float(row[key])

    number_data={}
    for loading in ("pure","flex"):
        for needle in (1,2,3):
            number_data[(loading,needle)]=row_for(loading,"suture number",needle=needle)

    stiff_data={}
    for loading in ("pure","flex"):
        for needle in (1,3):
            for stiffness in (10,50,200):
                if stiffness==50:
                    row=number_data[(loading,needle)]
                else:
                    row=row_for(loading,"suture stiffness",needle=needle,stiffness=stiffness)
                stiff_data[(loading,needle,stiffness)]=row

    fig,axes=plt.subplots(2,2,figsize=(10.8,7.2))
    load_style={"pure":("#2E6FB0","o","pure axial compression"),
                "flex":("#D55E00","s","flexion + axial compression")}
    metrics=(("peak_prot_mm","peak contained protrusion (mm)"),
             ("positive_AUC_mmN","positive-only AUC (mm·N)"))

    for ci,(key,ylabel) in enumerate(metrics):
        ax=axes[0,ci]
        for loading in ("pure","flex"):
            color,marker,label=load_style[loading]
            yy=[value(number_data[(loading,n)],key) for n in (1,2,3)]
            ax.plot((1,2,3),yy,color=color,marker=marker,lw=1.8,ms=5,label=label)
        ax.set_xticks((1,2,3)); ax.set_xlabel("number of load-targeted sutures")
        top=max(value(number_data[(loading,n)],key)
                for loading in ("pure","flex") for n in (1,2,3))
        ax.set_ylabel(ylabel); ax.set_ylim(0,top*1.08); ax.grid(alpha=0.22)
        max_change=max(
            100*abs(value(number_data[(loading,3)],key)-value(number_data[(loading,1)],key))
            / value(number_data[(loading,1)],key)
            for loading in ("pure","flex")
        )
        ax.text(0.98,0.06,f"largest 1-to-3 change: {max_change:.1f}%",
                transform=ax.transAxes,ha="right",fontsize=8,color="#444")
        if ci==0: ax.legend(fontsize=8,loc="lower left")

    line_style={("pure",1):("#2E6FB0","o","-"),("pure",3):("#2E6FB0","s","--"),
                ("flex",1):("#D55E00","o","-"),("flex",3):("#D55E00","s","--")}
    xcat=(0,1,2); stiffnesses=(10,50,200)
    for ci,(key,ylabel) in enumerate(metrics):
        ax=axes[1,ci]
        for loading in ("pure","flex"):
            for needle in (1,3):
                color,marker,ls=line_style[(loading,needle)]
                yy=[value(stiff_data[(loading,needle,E)],key) for E in stiffnesses]
                ax.plot(xcat,yy,color=color,marker=marker,ls=ls,lw=1.6,ms=5,
                        label=f"{load_style[loading][2]}, {needle} suture{'s' if needle>1 else ''}")
        ax.set_xticks(xcat,stiffnesses); ax.set_xlabel("axial suture stiffness (N/mm)")
        top=max(value(stiff_data[(loading,needle,E)],key)
                for loading in ("pure","flex") for needle in (1,3) for E in stiffnesses)
        ax.set_ylabel(ylabel); ax.set_ylim(0,top*1.08); ax.grid(alpha=0.22)
        max_change=max(
            100*(max(value(stiff_data[(loading,needle,E)],key) for E in stiffnesses)
                 -min(value(stiff_data[(loading,needle,E)],key) for E in stiffnesses))
            / value(stiff_data[(loading,needle,50)],key)
            for loading in ("pure","flex") for needle in (1,3)
        )
        ax.text(0.98,0.20 if ci==1 else 0.06,f"largest 10-to-200 range: {max_change:.1f}%",
                transform=ax.transAxes,ha="right",fontsize=8,color="#444")
        if ci==1: ax.legend(fontsize=7.2,loc="lower left",ncol=2)

    for label,ax in zip("abcd",axes.flat):
        ax.text(-0.12,1.04,label,transform=ax.transAxes,fontsize=11,fontweight="bold")
    axes[0,0].set_title("Suture-number sensitivity",fontsize=10,loc="left")
    axes[0,1].set_title("Suture-number sensitivity",fontsize=10,loc="left")
    axes[1,0].set_title("Suture-stiffness sensitivity",fontsize=10,loc="left")
    axes[1,1].set_title("Suture-stiffness sensitivity",fontsize=10,loc="left")
    fig.suptitle("Contained protrusion showed limited sensitivity to suture number and axial suture stiffness",
                 fontsize=12,x=0.08,ha="left")
    fig.text(0.5,0.015,
             "Deterministic outputs from one finite-element specimen; no error bars or inferential statistics. Stiffness is not failure strength or pull-out resistance.",
             ha="center",fontsize=8,color="#333")
    fig.tight_layout(rect=(0.03,0.045,1,0.95)); save(fig,"figR5_suture_number_stiffness")

if __name__=="__main__":
    R1(); R2(); R4(); R5()
    print("done R1,R2,R4,R5")
