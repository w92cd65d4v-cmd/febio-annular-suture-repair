#!/usr/bin/env python3
"""Locate non-smooth (sub-element) deformation modes in a refined run: for each refined parent, the max
deviation of child-node displacement from the trilinear interpolation of the parent's corner displacements.
Usage: python3 revision_meshsens_mode.py RUN_STEM REFINE_JSON BASE_FEB"""
import sys, json, numpy as np
from revision_mesh_refine import parse as pfeb, XI
from revision_meshsens_analyze import parse
stem, rj, basefeb = sys.argv[1:4]
N0,B=pfeb(open(stem+".feb",encoding="ISO-8859-1").read()); E={e:c for b in B for e,c in b["els"]}; NAME={e:b["name"] for b in B for e,_ in b["els"]}
ch={int(p):v for p,v in json.load(open(rj))["children"].items()}
Nb,Bb=pfeb(open(basefeb,encoding="ISO-8859-1").read()); Eb={e:c for b in Bb for e,c in b["els"]}
data=dict(parse(stem+".log")); ts=sorted(data); C=np.array(json.load(open("region_sets.json"))["C"])
D=set(json.load(open("region_sets.json"))["defect"]); RIM=set(json.load(open("metric_sets.json"))["rim"])
XIa=np.array(XI,float)
def dev_at(t):
    X=data[t]["x"]; X0=data[ts[0]]["x"]; out=[]
    for p,kids in ch.items():
        cn=Eb[p]
        if not all(n in X for n in cn): continue
        U=np.array([X[n]-X0[n] for n in cn]); P0=np.array([Nb[n] for n in cn]); m=0
        for n in {n for k in kids for n in E[k]}:
            if n not in X: continue
            x=N0[n]; xi=np.zeros(3)
            for _ in range(8):
                Nf=np.prod(1+XIa*xi,axis=1)/8
                dN=np.array([[a[0]*(1+a[1]*xi[1])*(1+a[2]*xi[2]),a[1]*(1+a[0]*xi[0])*(1+a[2]*xi[2]),a[2]*(1+a[0]*xi[0])*(1+a[1]*xi[1])] for a in XIa])/8
                xi-=np.linalg.solve(P0.T@dN, Nf@P0-x)
            m=max(m,np.linalg.norm((X[n]-X0[n])-(np.prod(1+XIa*xi,axis=1)/8)@U))
        out.append((m,p,NAME[kids[0]],np.linalg.norm(P0.mean(0)-C),p in RIM,p in D))
    out.sort(reverse=True); return out
for t in (ts[len(ts)//3], ts[-1]):
    o=dev_at(t); print(f"t={t:.3f} F={2000*(t-1)/9:.0f}N median dev {np.median([x[0] for x in o])*1000:.1f} um; top:")
    for m,p,nm,d,rim,de in o[:8]: print(f"   {m*1000:6.1f} um parent {p} {nm:14s} {d:4.1f} mm rim={rim} defect-site={de}")
