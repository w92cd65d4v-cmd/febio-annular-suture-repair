"""Reviewer R1-2/R1-8/R1-9, R2-3/R2-4/R2-7: robustness of protrusion and rim-stress metrics.
Uses only existing *_map.log files (591 region-node positions + 374 element stresses per step).
No reruns. Outputs REVISION_METRIC_ROBUSTNESS.csv and REVISION_NODE_LAYOUT.csv."""
import re, json, math, csv, numpy as np
from revision_morph_diagnostics import load, dN, fibers
S=json.load(open("metric_sets.json")); R=json.load(open("region_sets.json"))
SUT=set(n for c in json.load(open("suture_config.json")).values() for p in c for n in p)
MOUTH10=list(S["mouth"]); MOUTH8=[n for n in MOUTH10 if n not in SUT]
RIM=S["rim"]; REGION_E=R["elems"]; DEF=R["defect"]
Ni,Ei,NAMEi,AX=load("comp_intact_pure.feb")
nuc={n for e,c in Ei.items() if 'NUCLEUS' in NAMEi[e] for n in c}
dn={n for e in DEF for n in Ei[e]}; WALL=sorted(dn-nuc)            # 46 defect-margin AF nodes
FACES=[(0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
FAC=[]                                                             # exposed NP facets (outward normal)
for e,c in Ei.items():
    if 'NUCLEUS' in NAMEi[e]:
        ctr=np.mean([Ni[n] for n in c],0)
        for f in FACES:
            q=[c[i] for i in f]
            if set(q)<=set(MOUTH10):
                X=np.array([Ni[n] for n in q]); nrm=np.cross(X[2]-X[0],X[3]-X[1]); A=np.linalg.norm(nrm)/2; nrm/=np.linalg.norm(nrm)
                if nrm@(X.mean(0)-ctr)<0: nrm=-nrm
                FAC.append((q,A,nrm))
cx,cy=-9.85,-17.12
m=np.mean([Ni[n][:2] for n in MOUTH8],0); nout=np.r_[(m-[cx,cy])/np.linalg.norm(m-[cx,cy]),0]
def vm(s):
    sx,sy,sz,sxy,syz,sxz=s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))
def parse(path):
    rec={}
    for blk in re.split(r'Data Record #\d+',open(path,encoding="ISO-8859-1",errors="ignore").read()):
        mt=re.search(r'Time\s*=\s*([\d.eE+-]+)',blk); md=re.search(r'Data\s*=\s*([a-zA-Z;]+)',blk)
        if not(mt and md): continue
        k={'x;y;z':'x','sx;sy;sz;sxy;syz;sxz':'s'}.get(md.group(1))
        if not k: continue
        d=rec.setdefault(float(mt.group(1)),{'x':{},'s':{}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$',blk,re.M):
            d[k][int(lm.group(1))]=np.array([float(v) for v in lm.group(2).split()])
    ts=sorted(t for t in rec if rec[t]['x']); t0=min(ts,key=lambda t:abs(t-1)); return [(t,rec[t]) for t in ts if t>=t0-1e-9]
def curve_metrics(F,p):
    F=np.array(F);p=np.array(p); i=int(np.argmax(p))
    return dict(peak=p[i],F_peak=F[i],posAUC=np.trapz(np.maximum(p,0),F),signedAUC=np.trapz(p,F),endpoint=p[-1])
# reference geometry per model for fiber directions
def fib_geo(feb):
    N,E,NAME,AX=load(feb); return N,E,NAME,AX
STATES={"intact":("comp_intact_%s_map.log","comp_intact_pure.feb"),"defect":("comp_defect_%s_map.log","comp_defect_pure.feb"),
        "suture1":("comp_suture1_%s_map.log","comp_suture1_pure.feb"),"suture3":("comp_suture3_%s_map.log","comp_suture3_pure.feb"),
        "app12":("comp_suture3_app12_%s_map.log","comp_suture3_app12_pure_map.feb")}
rows=[]
for st,(lg,feb) in STATES.items():
    N,E,NAME,AX=fib_geo(feb)
    rimE=[e for e in RIM if e in E]; regAF=[e for e in REGION_E if e in E and 'LAM' in NAME[e]]
    for scen in ("pure","flex"):
        data=parse(lg%scen); ref=data[0][1]['x']
        Fs=[];P={k:[] for k in ('base8','all10','defect_rel','surface','surface_out','wall_only')}; LOO={n:[] for n in MOUTH8}
        for t,d in data:
            u=lambda n:(d['x'][n]-ref[n])
            Fs.append(max(0.0,2000*(t-1)/9))
            pn={n:u(n)@nout for n in MOUTH10}
            P['base8'].append(np.mean([pn[n] for n in MOUTH8])); P['all10'].append(np.mean(list(pn.values())))
            w=np.mean([u(n)@nout for n in WALL]); P['wall_only'].append(w); P['defect_rel'].append(P['base8'][-1]-w)
            P['surface'].append(sum(A*np.mean([u(n) for n in q],0)@nrm for q,A,nrm in FAC)/sum(A for _,A,_ in FAC))
            P['surface_out'].append(sum(A*np.mean([u(n) for n in q],0)@nout for q,A,nrm in FAC)/sum(A for _,A,_ in FAC))
            for n in MOUTH8: LOO[n].append(np.mean([pn[k] for k in MOUTH8 if k!=n]))
        row=dict(state=st,scen=scen)
        for k,p in P.items():
            for mk,v in curve_metrics(Fs,p).items(): row[f"{k}_{mk}"]=v
        lo=[curve_metrics(Fs,p) for p in LOO.values()]
        row.update(loo_peak_min=min(x['peak'] for x in lo),loo_peak_max=max(x['peak'] for x in lo),
                   loo_posAUC_min=min(x['posAUC'] for x in lo),loo_posAUC_max=max(x['posAUC'] for x in lo))
        # rim/regional stress + fiber measures at the 2000-N endpoint
        dE=data[-1][1]; X=dE['x']
        vr=np.array([vm(dE['s'][e]) for e in rimE]); vg=np.array([vm(dE['s'][e]) for e in regAF])
        row.update(rim_max=vr.max(),rim_top5=np.sort(vr)[-5:].mean(),rim_p95=np.percentile(vr,95),rim_p90=np.percentile(vr,90),rim_mean=vr.mean(),
                   regAF_p95=np.percentile(vg,95),regAF_mean=vg.mean(),n_rim=len(vr),n_regAF=len(vg))
        lam=[];sff=[]
        for e in rimE:
            c=E[e]
            if not all(n in X for n in c): continue
            X0=np.array([N[n] for n in c]); X1=np.array([X[n] for n in c])
            Fg=(X1.T@dN)@np.linalg.inv(X0.T@dN)
            sx,sy,sz,sxy,syz,sxz=dE['s'][e]; sig=np.array([[sx,sxy,sxz],[sxy,sy,syz],[sxz,syz,sz]])
            for a0 in fibers(X0,AX[NAME[e]]):
                a=Fg@a0; l=np.linalg.norm(a); a/=l; lam.append(l); sff.append(a@sig@a)
        lam=np.array(lam);sff=np.array(sff)
        row.update(fib_lam_max=lam.max(),fib_lam_p95=np.percentile(lam,95),fib_sff_max=sff.max(),fib_sff_p95=np.percentile(sff,95),n_fib=len(lam))
        rows.append(row)
with open("REVISION_METRIC_ROBUSTNESS.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
    for r in rows: w.writerow({k:(round(float(v),4) if not isinstance(v,str) else v) for k,v in r.items()})
# node layout (defect-local coordinates)
sp=json.load(open("suture_pairs.json")); C=np.array(sp["C"]); r=np.r_[sp["rad"],0]; tg=np.r_[sp["tang"],0]
with open("REVISION_NODE_LAYOUT.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["node","tangential_mm","radial_mm","axial_z_mm","suture_anchor","used_in_metric"])
    for n in MOUTH10:
        d=Ni[n]-C; w.writerow([n,round(d@tg,3),round(d@r,3),round(d[2],3),n in SUT,n in MOUTH8])
print("facets",len(FAC),"area",round(sum(A for _,A,_ in FAC),3),"wall nodes",len(WALL),"nout",nout.round(3))
