#!/usr/bin/env python3
"""Build the knot-tightened APPOSED-anchor suture repair models.

Geometric rest-mesh morph (defect elements stay deleted; tension-only springs kept).
Opposing defect-edge anchor nodes are drawn toward each pair midpoint by a SMOOTH,
SELF-BALANCING, distance-decaying displacement field so closure is distributed:
   u(n) = sum_pairs e_p * A_p * ( K(|n-a|) - K(|n-b|) ),  K(d)=exp(-(d/sigma)^2)
e_p = unit a->b. A_p = scale*(gap0_p)/2 (scale calibrated by a forward sweep on the
CLOSING branch — before nodes cross over). The morphed mesh is the new stress-free
reference; the existing springs then have ~residual rest length.

Because the slot-adjacent lamella elements are thin (min corner Jacobian ~0.08 in the
as-built mesh), exact co-location inverts them. We therefore select the TIGHTEST residual
gap that keeps every affected element's min corner Jacobian positive (with margin), and
report the residual-vs-Jacobian tradeoff. Outputs comp_<cfg>_apposed_<scen>.feb.
"""
import re, json, math
CFG = json.load(open("suture_config.json"))
SIGMA = 1.6; MAXR = 5.0
RES_TARGETS = [0.05,0.1,0.2,0.3,0.45,0.6,0.8,1.0,1.3]   # try tightest first
CORN=[(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
def min_jac(xyz8):
    mn=1e18
    for (xi,et,ze) in CORN:
        J=[[0.0]*3 for _ in range(3)]
        for i,(xj,ej,zj) in enumerate(CORN):
            dN=(0.125*xj*(1+et*ej)*(1+ze*zj),0.125*ej*(1+xi*xj)*(1+ze*zj),0.125*zj*(1+xi*xj)*(1+et*ej))
            x=xyz8[i]
            for r in range(3):
                for c in range(3): J[r][c]+=dN[c]*x[r]
        det=(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]))
        mn=min(mn,det)
    return mn
def load_nodes(txt):
    nsec=txt[txt.index("<Nodes"):txt.index("</Nodes>")]
    return {int(m.group(1)):[float(v) for v in m.group(2).split(',')] for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec)}
def all_elems(txt):
    conn={}
    for blk in re.finditer(r'<Elements type="hex8"[^>]*>(.*?)</Elements>', txt, re.S):
        for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', blk.group(1)):
            conn[int(e.group(1))]=[int(v) for v in e.group(2).split(',')]
    return conn

def build(cfg):
    pairs=CFG[cfg]
    src=open(f"comp_{cfg}_pure.feb",encoding="ISO-8859-1").read()
    N0=load_nodes(src); conn=all_elems(src)
    units=[]
    for a,b in pairs:
        d=[N0[b][k]-N0[a][k] for k in range(3)]; L=math.hypot(*d); units.append((a,b,[di/L for di in d],L))
    # auto-widen the morph spread (sigma) until a valid (non-inverting) tight apposition is found
    for SIGMA in (1.6,2.2,2.8,3.6,4.6):
        r=_try_sigma(cfg,N0,conn,units,SIGMA)
        if r is not None and r.get("viable"): return r
        last=r
    print(f"[{cfg}] no valid geometric apposition at any tested sigma")
    return last

def _try_sigma(cfg,N0,conn,units,SIGMA):
    affected=[n for n,xyz in N0.items() if any(math.dist(xyz,N0[a])<MAXR or math.dist(xyz,N0[b])<MAXR for a,b,_,_ in units)]
    Kc={n:[(math.exp(-(math.dist(N0[n],N0[a])/SIGMA)**2)-math.exp(-(math.dist(N0[n],N0[b])/SIGMA)**2)) for a,b,_,_ in units] for n in affected}
    def disp_at(scale):
        out={}
        for n in affected:
            u=[0.0,0.0,0.0]
            for j,(a,b,e,L) in enumerate(units):
                A=scale*L/2.0; w=Kc[n][j]
                for k in range(3): u[k]+=e[k]*A*w
            out[n]=u
        return out
    def worst_res(disp):
        return max(math.dist([N0[a][k]+disp[a][k] for k in range(3)],[N0[b][k]+disp[b][k] for k in range(3)]) for a,b,_,_ in units)
    # forward sweep on the closing branch: scale 0 -> until worst residual reaches its minimum
    scales=[i*0.01 for i in range(0,400)]
    curve=[]  # (scale, worst_res)
    prev=1e9
    for s in scales:
        wr=worst_res(disp_at(s))
        curve.append((s,wr))
        if wr>prev+1e-6: break   # passed the minimum (overshoot) -> stop
        prev=wr
    min_achievable=curve[-1][1]
    # affected elements (those touching any moved node)
    moved=set(affected)
    aff_el=[e for e,ns in conn.items() if any(n in moved for n in ns)]
    jac0=min(min_jac([N0[n] for n in conn[e]]) for e in aff_el)
    # choose tightest residual target with positive Jacobian margin
    chosen=None
    table=[]
    for tgt in RES_TARGETS:
        if tgt < min_achievable - 1e-6:   # cannot close that far on closing branch
            continue
        # smallest scale on closing branch achieving worst_res<=tgt
        s=next((sc for sc,wr in curve if wr<=tgt), curve[-1][0])
        disp=disp_at(s)
        Nn={n:[N0[n][k]+disp.get(n,[0,0,0])[k] for k in range(3)] for n in N0}
        jmin=min(min_jac([Nn[n] for n in conn[e]]) for e in aff_el)
        ninv=sum(1 for e in aff_el if min_jac([Nn[n] for n in conn[e]])<=0)
        table.append((tgt,s,worst_res(disp),jmin,ninv))
        if ninv==0 and jmin>0.003 and chosen is None:
            chosen=(tgt,s,Nn,worst_res(disp),jmin)
    print(f"[{cfg}] sigma={SIGMA} affected_nodes={len(affected)} aff_elems={len(aff_el)} as-built minJac={jac0:.4f} min-achievable-gap={min_achievable:.3f}mm")
    print(f"    residual-vs-Jacobian tradeoff (target, scale, worstGap, minJac, #inverted):")
    for tgt,s,wr,jm,ni in table: print(f"      tgt={tgt:.2f}  scale={s:.2f}  gap={wr:.3f}  minJac={jm:+.4f}  inverted={ni}")
    if chosen is None:
        print("    !! no valid residual found on closing branch -> geometric morph not viable");
        return dict(cfg=cfg,viable=False,asbuilt_minjac=jac0,min_gap=min_achievable)
    tgt,s,Nn,gap,jmin=chosen
    print(f"    CHOSEN residual gap={gap:.3f} mm (scale={s:.2f}) keeping minJac={jmin:.4f} (no inverted elements)")
    for a,b,_,L0 in units:
        print(f"    pair {a}-{b}: {L0:.3f} mm -> {math.dist(Nn[a],Nn[b]):.3f} mm")
    def write(scen):
        txt=open(f"comp_{cfg}_{scen}.feb",encoding="ISO-8859-1").read()
        a=txt.index("<Nodes"); b=txt.index("</Nodes>"); head,nsec,tail=txt[:a],txt[a:b],txt[b:]
        nsec=re.sub(r'<node id="(\d+)">\s*[^<]+</node>', lambda m:f'<node id="{m.group(1)}">{Nn[int(m.group(1))][0]:.7g},{Nn[int(m.group(1))][1]:.7g},{Nn[int(m.group(1))][2]:.7g}</node>', nsec)
        open(f"comp_{cfg}_apposed_{scen}.feb","w",encoding="ISO-8859-1").write(head+nsec+tail); print(f"    wrote comp_{cfg}_apposed_{scen}.feb")
    write("pure"); write("flex")
    return dict(cfg=cfg,viable=True,residual_gap=gap,scale=s,minjac=jmin,asbuilt_minjac=jac0,
                pairs={f"{a}-{b}":[round(L0,3),round(math.dist(Nn[a],Nn[b]),3)] for a,b,_,L0 in units})

if __name__=="__main__":
    res=[build(c) for c in ("suture1","suture3")]
    json.dump(res, open("apposed_build_report.json","w"), indent=2); print("\n-> apposed_build_report.json")
