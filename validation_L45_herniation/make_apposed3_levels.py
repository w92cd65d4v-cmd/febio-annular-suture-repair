#!/usr/bin/env python3
"""Closure-degree sensitivity for the 3-SUTURE apposed repair, at looser residual gaps
to remove the endpoint overcorrection seen at near-contact. Same self-balancing
distance-decaying morph as make_apposed.py, fixed sigma=2.2 (the spread that was valid at
near-contact for this 3-pair config). Closure is calibrated to a target MEAN anchor-pair
residual gap; every pair gap is reported. Element quality = min corner Jacobian + min
SCALED Jacobian + counts below thresholds. Targets: mean residual 1.2/1.0/0.8/0.6 mm.
Outputs comp_suture3_app{12,10,08,06}_{pure,flex}.feb."""
import re, json, math
CFG=json.load(open("suture_config.json")); SIGMA=2.2; MAXR=6.0
TARGETS=[(1.7,"app17"),(1.5,"app15"),(1.2,"app12"),(1.0,"app10b"),(0.8,"app08"),(0.6,"app06")]   # suture3 mean residual
CORN=[(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
def jac_corners(xyz8):
    mnd=1e18; mns=1e18
    for (xi,et,ze) in CORN:
        g=[[0.0]*3 for _ in range(3)]
        for i,(xj,ej,zj) in enumerate(CORN):
            dN=(0.125*xj*(1+et*ej)*(1+ze*zj),0.125*ej*(1+xi*xj)*(1+ze*zj),0.125*zj*(1+xi*xj)*(1+et*ej))
            x=xyz8[i]
            for a in range(3):
                for k in range(3): g[a][k]+=dN[a]*x[k]
        det=(g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])-g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])+g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]))
        nrm=math.sqrt(sum(c*c for c in g[0]))*math.sqrt(sum(c*c for c in g[1]))*math.sqrt(sum(c*c for c in g[2]))
        mnd=min(mnd,det); mns=min(mns,det/nrm if nrm>0 else 0)
    return mnd,mns
def load_nodes(txt):
    ns=txt[txt.index("<Nodes"):txt.index("</Nodes>")]
    return {int(m.group(1)):[float(v) for v in m.group(2).split(',')] for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',ns)}
def all_elems(txt):
    c={}
    for blk in re.finditer(r'<Elements type="hex8"[^>]*>(.*?)</Elements>',txt,re.S):
        for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>',blk.group(1)): c[int(e.group(1))]=[int(v) for v in e.group(2).split(',')]
    return c

src=open("comp_suture3_pure.feb",encoding="ISO-8859-1").read()
N0=load_nodes(src); conn=all_elems(src)
units=[]
for a,b in CFG["suture3"]:
    d=[N0[b][k]-N0[a][k] for k in range(3)]; L=math.hypot(*d); units.append((a,b,[di/L for di in d],L))
affected=[n for n,x in N0.items() if any(math.dist(x,N0[a])<MAXR or math.dist(x,N0[b])<MAXR for a,b,_,_ in units)]
Kc={n:[(math.exp(-(math.dist(N0[n],N0[a])/SIGMA)**2)-math.exp(-(math.dist(N0[n],N0[b])/SIGMA)**2)) for a,b,_,_ in units] for n in affected}
aff_el=[e for e,ns in conn.items() if any(n in set(affected) for n in ns)]
def disp(scale):
    out={}
    for n in affected:
        u=[0,0,0]
        for j,(a,b,e,L) in enumerate(units):
            A=scale*L/2.0
            for k in range(3): u[k]+=e[k]*A*Kc[n][j]
        out[n]=u
    return out
def gaps(d): return [math.dist([N0[a][k]+d[a][k] for k in range(3)],[N0[b][k]+d[b][k] for k in range(3)]) for a,b,_,_ in units]
report=[]
for tgt,tag in TARGETS:
    s=0.0; prev=1e9; chosen=None
    while s<4.0:
        g=gaps(disp(s)); mg=sum(g)/len(g)
        if mg>prev+1e-6: break
        if mg<=tgt and chosen is None: chosen=s; break
        prev=mg; s+=0.005
    scale=chosen if chosen is not None else s
    d=disp(scale); Nn={n:[N0[n][k]+d.get(n,[0,0,0])[k] for k in range(3)] for n in N0}
    g=gaps(d); mg=sum(g)/len(g)
    js=[jac_corners([Nn[n] for n in conn[e]]) for e in aff_el]
    minj=min(j[0] for j in js); mins=min(j[1] for j in js); below=lambda th: sum(1 for j in js if j[1]<th)
    print(f"[suture3 {tag}] target mean={tgt}mm -> mean residual={mg:.3f}mm scale={scale:.3f} | pairs={[round(x,3) for x in g]} | minJac={minj:.4f} minScaledJac={mins:.4f} sJ<0.2:{below(0.2)} <0.1:{below(0.1)} <=0:{below(0.0)}")
    for scen in ("pure","flex"):
        txt=open(f"comp_suture3_{scen}.feb",encoding="ISO-8859-1").read()
        a=txt.index("<Nodes"); b=txt.index("</Nodes>")
        nsec=re.sub(r'<node id="(\d+)">\s*[^<]+</node>', lambda m:f'<node id="{m.group(1)}">{Nn[int(m.group(1))][0]:.7g},{Nn[int(m.group(1))][1]:.7g},{Nn[int(m.group(1))][2]:.7g}</node>', txt[a:b])
        open(f"comp_suture3_{tag}_{scen}.feb","w",encoding="ISO-8859-1").write(txt[:a]+nsec+txt[b:])
    report.append(dict(tag=tag,target_mean=tgt,mean_residual=round(mg,3),pairs=[round(x,3) for x in g],scale=round(scale,3),
                       minjac=round(minj,4),min_scaled_jac=round(mins,4),below_0p2=below(0.2),below_0p1=below(0.1),inverted=below(0.0)))
# score existing near-apposed suture3 mesh for the table
txt=open("comp_suture3_apposed_pure.feb",encoding="ISO-8859-1").read(); Nn=load_nodes(txt)
js=[jac_corners([Nn[n] for n in conn[e]]) for e in aff_el]
g=[math.dist(Nn[a],Nn[b]) for a,b,_,_ in units]
print(f"[suture3 near] mean residual={sum(g)/3:.3f}mm pairs={[round(x,3) for x in g]} | minJac={min(j[0] for j in js):.4f} minScaledJac={min(j[1] for j in js):.4f} sJ<0.2:{sum(1 for j in js if j[1]<0.2)} <=0:{sum(1 for j in js if j[1]<=0)}")
report.append(dict(tag="near",mean_residual=round(sum(g)/3,3),pairs=[round(x,3) for x in g],minjac=round(min(j[0] for j in js),4),min_scaled_jac=round(min(j[1] for j in js),4),below_0p2=sum(1 for j in js if j[1]<0.2),inverted=sum(1 for j in js if j[1]<=0)))
json.dump(report,open("apposed3_levels_quality.json","w"),indent=2); print("\n-> apposed3_levels_quality.json")
