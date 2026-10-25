"""Reviewer R1-1 / R2-5: quantify what the apposition morph does to the mesh.
Compares the open-gap 3-suture model with each morphed (app*) model.
Morph displacement below 1e-4 mm is treated as coordinate-rounding noise."""
import re, sys, json, numpy as np, csv
TOL=1e-4
def load(p):
    t=open(p,encoding="ISO-8859-1").read()
    ns=t[t.index("<Nodes"):t.index("</Nodes>")]
    N={int(m.group(1)):np.array([float(v) for v in m.group(2).split(',')]) for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',ns)}
    E={};NAME={}
    for b in re.finditer(r'<Elements type="hex8" mat="(\d+)" name="([^"]+)">(.*?)</Elements>',t,re.S):
        for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>',b.group(3)):
            E[int(e.group(1))]=[int(v) for v in e.group(2).split(',')];NAME[int(e.group(1))]=b.group(2)
    mat=t[t.index("<Material>"):t.index("</Material>")]
    AX={m.group(1):m.group(2) for m in re.finditer(r'<material id="\d+" name="([^"]+)"[^>]*>\s*<mat_axis type="local">([^<]+)</mat_axis>',mat)}
    return N,E,NAME,AX
# hex8 shape-function derivatives at centroid
XI=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)
dN=XI/8.0
TETS=[(0,1,3,4),(1,2,3,6),(1,4,5,6),(3,4,6,7),(1,3,4,6)]
def vol(X): return sum(abs(np.dot(X[b]-X[a],np.cross(X[c]-X[a],X[d]-X[a])))/6 for a,b,c,d in TETS)
def fibers(X,axis,ang=(45,135)):
    i,j,k=[int(v)-1 for v in axis.split(',')]
    a=X[j]-X[i]; d=X[k]-X[i]
    e1=a/np.linalg.norm(a); e3=np.cross(a,d); e3/=np.linalg.norm(e3)
    return [np.sin(np.radians(g))*e1+np.cos(np.radians(g))*e3 for g in ang]
def run(ref,mor):
    N0,E,NAME,AX=load(ref); N1,_,_,_=load(mor)
    moved={n for n in N0 if np.linalg.norm(N1[n]-N0[n])>TOL}
    aff=[e for e,c in E.items() if any(n in moved for n in c)]
    out=dict(model=mor,moved_nodes=len(moved),max_disp=max(np.linalg.norm(N1[n]-N0[n]) for n in moved))
    grp={'annulus':[e for e in aff if 'LAM' in NAME[e]],'nucleus':[e for e in aff if 'NUCLEUS' in NAME[e]]}
    grp['other']=[e for e in aff if e not in set(grp['annulus'])|set(grp['nucleus'])]
    for g,L in grp.items(): out[f'n_{g}']=len(L)
    dv=[];ang=[];lam=[];J=[]
    for e in aff:
        X0=np.array([N0[n] for n in E[e]]);X1=np.array([N1[n] for n in E[e]])
        dv.append(100*(vol(X1)-vol(X0))/vol(X0))
        if e in set(grp['annulus']):
            F=(X1.T@dN)@np.linalg.inv(X0.T@dN)
            for f0,f1 in zip(fibers(X0,AX[NAME[e]]),fibers(X1,AX[NAME[e]])):
                ang.append(np.degrees(np.arccos(min(1,abs(f0@f1)))))
                lam.append(np.linalg.norm(F@f0))   # geometric fiber stretch discarded by stress-free reset
    dv=np.array(dv);ang=np.array(ang);lam=np.array(lam)
    q=lambda a,p:float(np.percentile(a,p))
    out.update(vol_med=q(dv,50),vol_p5=q(dv,5),vol_p95=q(dv,95),vol_min=dv.min(),vol_max=dv.max(),
               vol_gt10=int(np.sum(abs(dv)>10)),vol_gt20=int(np.sum(abs(dv)>20)),n_elem=len(dv),
               ang_med=q(ang,50),ang_p95=q(ang,95),ang_max=ang.max(),ang_gt5=int(np.sum(ang>5)),n_fib=len(ang),
               lam_p5=q(lam,5),lam_med=q(lam,50),lam_p95=q(lam,95),lam_min=lam.min(),lam_max=lam.max(),
               lam_gt5=int(np.sum(abs(lam-1)>0.05)),lam_gt10=int(np.sum(abs(lam-1)>0.10)))
    # NP-AF interface
    nn={n for e,c in E.items() if 'NUCLEUS' in NAME[e] for n in c}; an={n for e,c in E.items() if 'LAM' in NAME[e] for n in c}
    ifc=[np.linalg.norm(N1[n]-N0[n]) for n in (nn&an)]
    ifm=[x for x in ifc if x>TOL]
    out.update(iface_total=len(ifc),iface_moved=len(ifm),iface_max=max(ifm) if ifm else 0,iface_gt01=sum(x>0.1 for x in ifm))
    # overall AF volume
    V0=sum(vol(np.array([N0[n] for n in E[e]])) for e in E if 'LAM' in NAME[e]); V1=sum(vol(np.array([N1[n] for n in E[e]])) for e in E if 'LAM' in NAME[e])
    out['AF_total_vol_change_pct']=100*(V1-V0)/V0
    V0=sum(vol(np.array([N0[n] for n in E[e]])) for e in E if 'NUCLEUS' in NAME[e]); V1=sum(vol(np.array([N1[n] for n in E[e]])) for e in E if 'NUCLEUS' in NAME[e])
    out['NP_total_vol_change_pct']=100*(V1-V0)/V0
    return out
if __name__=="__main__":
    ref="comp_suture3_pure.feb"
    import glob
    mods=sorted(glob.glob("comp_suture3_app*_pure_map.feb"))
    rows=[run(ref,m) for m in mods]
    with open("REVISION_MORPH_DIAGNOSTICS.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();[w.writerow(r) for r in rows]
    for r in rows: print({k:(round(float(v),3) if not isinstance(v,(str,int)) or isinstance(v,float) else v) for k,v in r.items()})
