#!/usr/bin/env python3
"""Model-construction figure panels (Figure 1A-1D) for the manuscript.

Renders directly from the intact FEBio mesh (comp_intact_flex_map.feb) with a
painter's-algorithm opaque surface renderer. No mechanical experiment is rerun.

 1A  whole intact L4-L5 FSU, posterolateral view, opaque, restrained anatomy colors
 1B  same camera; posterior bony arch + facets faded so the posterolateral annulus
     shows; annulotomy footprint (18 elements = L45_POST_LAM1..5) drawn on the annulus
 1C  axial section through the defect: 5 affected lamellae + removed-element slot + dims
 1D  boundary-condition / loading schematic (coronal silhouette + arrows + phase timeline)

Outputs individual PNGs, a combined figure, and figure_render_settings.txt.
"""
from __future__ import annotations
import json, re, math
from collections import defaultdict
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.patheffects as pe

HERE = Path(__file__).resolve().parent
FEB = HERE / "comp_intact_flex_map.feb"
REGION = json.loads((HERE / "region_sets.json").read_text())
DEFECT = set(REGION["defect"])
SP = json.loads((HERE / "suture_pairs.json").read_text())
TANG = np.array(SP["tang"]); RAD = np.array(SP["rad"])

FACE_IDX = ((0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
LIGHT = np.array([-0.5,-0.8,0.9]); LIGHT = LIGHT/np.linalg.norm(LIGHT)

# ---- restrained anatomical palette ----
C_ANN="#3a6f9c"; C_NUC="#6bb3b0"; C_FC="#bcd0d8"; C_CEND="#cfd6d2"
C_VEND="#cdbf9f"; C_POST="#c3ae83"; C_CORT="#dccfb0"; C_TRAB="#e7ddc5"; C_OTHER="#e0d8c4"
LAM_COLORS=["#9ec4e0","#6fa3cf","#4a83b8","#33689c","#204f7c"]  # LAM1(inner)->LAM5(outer)

def parse():
    txt=FEB.read_text(encoding="ISO-8859-1",errors="ignore")
    nsec=txt[txt.index("<Nodes"):txt.index("</Nodes>")]
    nodes={int(m.group(1)):np.array([float(v) for v in m.group(2).split(",")])
           for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>',nsec)}
    elems={}; ename={}
    for b in re.finditer(r'<Elements\s+type="hex8"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',txt,re.S):
        nm=b.group(1)
        for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>',b.group(2)):
            eid=int(e.group(1)); elems[eid]=[int(v) for v in e.group(2).split(",")]; ename[eid]=nm
    return nodes,elems,ename

def dcolor(name):
    if "NUCLEUS" in name: return C_NUC
    if "LAM" in name: return C_ANN
    if "_FC" in name: return C_FC
    if "CEND" in name: return C_CEND
    if "VEND" in name: return C_VEND
    if name.endswith("POST"): return C_POST
    if "CORT" in name: return C_CORT
    if "TRAB" in name: return C_TRAB
    return C_OTHER

def is_posterior(name):  # occluding structures faded in 1B
    return name.endswith("POST") or "_FC" in name

def exterior_faces(nodes,elems,ename,keep):
    counts=defaultdict(int); owner={}
    for eid in keep:
        conn=elems[eid]
        for idx in FACE_IDX:
            face=tuple(conn[i] for i in idx); key=tuple(sorted(face))
            counts[key]+=1; owner[key]=(eid,face)
    out=[]
    for key,c in counts.items():
        if c==1:
            eid,face=owner[key]
            if all(n in nodes for n in face): out.append((eid,[nodes[n] for n in face],ename[eid]))
    return out

def camera(nodes):
    pts=np.array(list(nodes.values())); center=pts.mean(0)
    view=np.array([-0.82,-1.0,0.18]); view/=np.linalg.norm(view)
    right=np.cross([0,0,1.0],view); right/=np.linalg.norm(right)
    up=np.cross(view,right); up/=np.linalg.norm(up)
    return center,view,right,up

def project(face,center,view,right,up):
    p=np.array(face)-center
    return list(zip(p@right,p@up)), float((p@view).mean())

def fnormal(face):
    p0,p1,p2=(np.array(x) for x in face[:3]); n=np.cross(p1-p0,p2-p1); nn=np.linalg.norm(n)
    return n/nn if nn else n

def shade(color,face):
    base=np.array(plt.matplotlib.colors.to_rgb(color))
    s=0.72+0.28*max(0.0,float(fnormal(face)@LIGHT))
    return np.clip(base*s,0,1)

def draw_triad(ax,x0,y0,L=10,labels=("+X","+Z")):
    ax.annotate("",xy=(x0+L,y0),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",color="#333",lw=1.4))
    ax.annotate("",xy=(x0,y0+L),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",color="#333",lw=1.4))
    ax.text(x0+L+1.5,y0,labels[0],fontsize=8,va="center",color="#333")
    ax.text(x0,y0+L+1.5,labels[1],fontsize=8,ha="center",color="#333")

# ---------------- 3D panels (1A, 1B) ----------------
def render_3d(ax,nodes,elems,ename,fade_posterior=False,footprint=False):
    keep=[e for e,n in ename.items() if n.startswith(("L4","L5","L45_"))]
    faces=exterior_faces(nodes,elems,ename,keep)
    center,view,right,up=camera(nodes)
    proj=[]
    for eid,face,name in faces:
        poly,depth=project(face,center,view,right,up)
        proj.append((depth,eid,poly,name,face))
    proj.sort(key=lambda t:t[0],reverse=True)
    polys=[p for _,_,p,_,_ in proj]; fcs=[]
    for _,_,_,name,face in proj:
        rgb=shade(dcolor(name),face)
        a=0.14 if (fade_posterior and is_posterior(name)) else 1.0
        fcs.append((rgb[0],rgb[1],rgb[2],a))
    ax.add_collection(PolyCollection(polys,facecolors=fcs,edgecolors=(0.4,0.4,0.4,0.35),linewidths=0.04,closed=True))
    if footprint:
        dfaces=exterior_faces(nodes,elems,ename,[e for e in DEFECT if e in elems])
        for eid,face,name in dfaces:
            poly,depth=project(face,center,view,right,up)
            xs,ys=zip(*(poly+[poly[0]]))
            ax.fill(xs,ys,facecolor="#d81ba6",alpha=0.55,zorder=6)
            ax.plot(xs,ys,color="#8a006a",lw=0.8,zorder=7)
    allp=np.array([pt for poly in polys for pt in poly])
    pad=4
    ax.set_xlim(allp[:,0].min()-pad,allp[:,0].max()+pad)
    ax.set_ylim(allp[:,1].min()-pad,allp[:,1].max()+pad)
    ax.set_aspect("equal"); ax.axis("off")
    return allp

def orient_arrows(ax,allp):
    x0=allp[:,0].min()+6; y0=allp[:,1].min()+7
    ax.annotate("",xy=(x0,y0+15),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",color="#444",lw=1.4))
    ax.annotate("",xy=(x0+15,y0),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",color="#444",lw=1.4))
    ax.text(x0,y0+17.5,"Sup",ha="center",va="bottom",fontsize=9,color="#444")
    ax.text(x0+17.5,y0,"Post",ha="left",va="center",fontsize=9,color="#444")

# ---------------- 1C axial section ----------------
def render_axial(ax,nodes,elems,ename,bottom_caption=True):
    dz=[nodes[n][2] for e in DEFECT for n in elems[e]]
    zc=(min(dz)+max(dz))/2
    band=[e for e,n in ename.items() if ("LAM" in n or "NUCLEUS" in n)
          and (min(nodes[k][2] for k in elems[e])<=zc<=max(nodes[k][2] for k in elems[e]))]
    C=np.array(REGION["C"])
    nucC=np.mean([nodes[n] for e,nm in ename.items() if "NUCLEUS" in nm for n in elems[e]],axis=0)[:2]
    def topface(e): return [(nodes[n][0],nodes[n][1]) for n in elems[e][:4]]
    polys=[]; cols=[]
    for e in band:
        nm=ename[e]
        if "NUCLEUS" in nm: col=C_NUC
        else:
            lam=int(nm[-1]) if nm[-1].isdigit() else 3
            col=LAM_COLORS[min(max(5-lam,0),4)]   # LAM5=inner→light, LAM1=outer→dark (verified radial order)
        polys.append(topface(e)); cols.append(col)
    ax.add_collection(PolyCollection(polys,facecolors=cols,edgecolors=(0.35,0.35,0.35,0.5),linewidths=0.25))
    # removed elements: magenta empty (light hatch so the radial arrow reads over it)
    dpoly=[topface(e) for e in DEFECT if e in elems and e in band] or [topface(e) for e in DEFECT if e in elems]
    ax.add_collection(PolyCollection(dpoly,facecolors="none",edgecolors="#c0007a",linewidths=1.8,hatch="//"))
    dc=np.array([C[0],C[1]])
    ax.set_xlim(dc[0]-12,dc[0]+12); ax.set_ylim(dc[1]-12,dc[1]+12); ax.set_aspect("equal"); ax.axis("off")
    dn=np.array([nodes[n] for e in DEFECT for n in elems[e]])
    tv=(dn[:,:2]-C[:2])@TANG; rv=(dn[:,:2]-C[:2])@RAD; tmid=(tv.min()+tv.max())/2
    # RADIAL arrow — through the defect centreline, from the nucleus-facing (inner) edge outward through the 5 lamellae
    e_in =C[:2]+TANG*tmid+RAD*rv.min(); e_out=C[:2]+TANG*tmid+RAD*rv.max()
    if np.hypot(*(e_in-nucC))>np.hypot(*(e_out-nucC)): e_in,e_out=e_out,e_in   # ensure inner=closer to nucleus
    ax.annotate("",xy=e_out,xytext=e_in,zorder=8,
                arrowprops=dict(arrowstyle="-|>",color="#111",lw=2.2,shrinkA=0,shrinkB=0))
    rlab=(e_in+e_out)/2+TANG*3.4
    ax.text(rlab[0],rlab[1],"radial depth\n~6.5 mm\n(5 lamellae,\ninner→outer)",fontsize=8,ha="center",va="center",color="#111",
            bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.85))
    ax.text(*(e_in-RAD*1.2),"inner\n(nucleus)",fontsize=7,ha="center",va="center",color="#0a6",
            bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.8))
    # CIRCUMFERENTIAL arrow — spans the slot width, offset just outside the outer edge along the true tangent
    off=RAD*(rv.max()+2.2)
    a=C[:2]+TANG*tv.min()+off; b=C[:2]+TANG*tv.max()+off
    ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="<|-|>",color="#111",lw=1.6))
    clab=(a+b)/2+RAD*1.6
    ax.text(clab[0],clab[1],"circumferential ~4.8 mm",fontsize=8,ha="center",va="center",color="#111",
            bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.85))
    # orientation: in-plane X/Y + out-of-plane superior/inferior (Z)
    draw_triad(ax,dc[0]-10.5,dc[1]-10.5,6,labels=("+X","+Y"))
    ax.text(dc[0]-10.5,dc[1]-6.2,"⊙ +Z (Sup)",fontsize=7.5,ha="left",color="#444")
    if bottom_caption:
        ax.text(0.5,0.015,"out-of-plane height ≈ 4.0 mm (superior–inferior)",
                transform=ax.transAxes,fontsize=8.5,ha="center",color="#333")

# ---------------- 1D BC / loading schematic ----------------
def render_bc(ax,nodes,elems,ename,mode="pure"):
    # coronal silhouette: project exterior faces to (x,z)
    keep=[e for e,n in ename.items() if n.startswith(("L4","L5","L45_"))]
    faces=exterior_faces(nodes,elems,ename,keep)
    polys=[[ (p[0],p[2]) for p in face] for _,face,_ in faces]
    cols=[]
    for _,face,name in faces:
        cols.append(C_ANN if "LAM" in name else C_NUC if "NUCLEUS" in name else "#ddd3bf")
    ax.add_collection(PolyCollection(polys,facecolors=cols,edgecolors=(0.6,0.6,0.6,0.25),linewidths=0.03))
    allp=np.array([pt for poly in polys for pt in poly])
    xmin,xmax=allp[:,0].min(),allp[:,0].max(); zmin,zmax=allp[:,1].min(),allp[:,1].max()
    xc=(xmin+xmax)/2
    rpad=62 if mode=="flex" else 16
    ax.set_xlim(xmin-16,xmax+rpad); ax.set_ylim(zmin-30,zmax+33); ax.set_aspect("equal"); ax.axis("off")
    # inferior fixed: hatched ground under L5
    gy=zmin-3
    ax.add_patch(Rectangle((xmin-6,gy-3),(xmax-xmin)+12,3,facecolor="none",edgecolor="#333",hatch="///",lw=1.0))
    ax.plot([xmin-6,xmax+6],[gy,gy],color="#333",lw=1.5)
    ax.text(xc,gy-7,"L5 inferior body: fixed in all six DOF",fontsize=9.5,ha="center",color="#222",
            bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.9))
    # superior endplate bar
    ax.plot([xmin-4,xmax+4],[zmax+2,zmax+2],color="#7a5a2a",lw=3,solid_capstyle="round")
    ax.text(xmin-4,zmax+4,"L4 rigid endplate",fontsize=9.5,ha="left",color="#5a3d18")
    # SINGLE central resultant axial-force arrow
    ax.annotate("",xy=(xc,zmax+2.5),xytext=(xc,zmax+15),
                arrowprops=dict(arrowstyle="-|>",color="#b02020",lw=3.0))
    ax.text(xc+3,zmax+15,r"$F_z = 0\to2000\ \mathrm{N}$",fontsize=12,ha="left",va="center",color="#b02020",fontweight="bold")
    if mode=="flex":
        arc=FancyArrowPatch((xmax+6, zmax-3),(xmax+6, zmin+7),
                            connectionstyle="arc3,rad=-0.5",arrowstyle="-|>",
                            mutation_scale=15,lw=2.0,color="#1f6f3f")
        ax.add_patch(arc)
        ax.text(xmax+27,(zmin+zmax)/2,r"$M_x = 7.5\ \mathrm{N\cdot m}$"+"\n(flexion)",fontsize=11,ha="left",va="center",
                color="#1f6f3f",fontweight="bold")
        title="Flexion + axial compression"
        sub="L4 superior body — fixed: x, y, Ry, Rz\nfree/loaded: z, Rx"
    else:
        title="Pure axial compression"
        sub="L4 superior body — fixed: x, y, Rx, Ry, Rz\nfree/loaded: z"
    ax.text(xc,zmax+27,title,fontsize=13,ha="center",fontweight="bold")
    ax.text(xc,gy-12,sub,fontsize=9.5,ha="center",va="top",color="#444",style="italic")
    # coordinate triad + out-of-plane Y indicator (right-handed: x right, z up → +Y out of page)
    draw_triad(ax,xmin-13,zmin+5,7,labels=("+X","+Z"))
    ax.text(xmin-13,zmin+0.5,"⊙ +Y (out of page)",fontsize=7.5,ha="left",color="#444")

def phase_timeline(ax):
    ax.set_xlim(0,10); ax.set_ylim(0,1); ax.axis("off")
    ax.add_patch(Rectangle((0,0.3),1,0.4,facecolor="#cfe0ef",edgecolor="#333"))
    ax.add_patch(Rectangle((1,0.3),9,0.4,facecolor="#f3ddd0",edgecolor="#333"))
    ax.text(0.5,0.5,"Phase 1",ha="center",va="center",fontsize=10,fontweight="bold")
    ax.text(5.5,0.5,"Phase 2",ha="center",va="center",fontsize=10,fontweight="bold")
    ax.text(0.5,0.12,"expansive nucleus-boundary\npressure preload 0 → 1 MPa (held)",ha="center",va="top",fontsize=9.5)
    ax.text(5.5,0.12,"mechanical load ramp: $F_z$ (and $M_x$) 0 → max",ha="center",va="top",fontsize=9.5)
    ax.text(0.5,0.88,"t=0",ha="center",fontsize=7); ax.text(1,0.88,"t=1",ha="center",fontsize=7)
    ax.text(10,0.88,"t=10",ha="center",fontsize=7)
    ax.annotate("",xy=(10,0.72),xytext=(0,0.72),arrowprops=dict(arrowstyle="-|>",color="#333",lw=1))
    ax.text(5,0.8,"pseudo-time",ha="center",fontsize=7,color="#333")

# ---------------- build all ----------------
def main():
    nodes,elems,ename=parse()
    # 1A
    fig,ax=plt.subplots(figsize=(6.4,7.2)); allp=render_3d(ax,nodes,elems,ename); orient_arrows(ax,allp)
    ax.set_title("A  Intact L4–L5 functional spinal unit",fontsize=12,loc="left")
    fig.savefig(HERE/"fig1A_whole_model.png",dpi=300,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("wrote fig1A_whole_model.png")
    # 1B
    fig,ax=plt.subplots(figsize=(6.4,7.2)); allp=render_3d(ax,nodes,elems,ename,fade_posterior=True,footprint=True); orient_arrows(ax,allp)
    ax.set_title("B  Annulotomy footprint on posterolateral annulus",fontsize=12,loc="left")
    # (explanatory sentence moved to the figure legend; only the essential magenta label kept)
    ax.text(0.02,0.02,"magenta = annulotomy footprint (18 hex elements)",
            transform=ax.transAxes,fontsize=8,color="#8a006a")
    fig.savefig(HERE/"fig1B_annulotomy_location.png",dpi=300,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("wrote fig1B_annulotomy_location.png")
    # 1C
    fig,ax=plt.subplots(figsize=(6.8,6.4)); render_axial(ax,nodes,elems,ename)
    ax.set_title("C  Defect geometry (axial section)",fontsize=12,loc="left")
    fig.savefig(HERE/"fig1C_defect_geometry.png",dpi=300,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("wrote fig1C_defect_geometry.png")
    # 1D
    fig=plt.figure(figsize=(12.5,6.6))
    gs=fig.add_gridspec(2,2,height_ratios=[5,1],hspace=0.05,wspace=0.05)
    axp=fig.add_subplot(gs[0,0]); render_bc(axp,nodes,elems,ename,"pure")
    axf=fig.add_subplot(gs[0,1]); render_bc(axf,nodes,elems,ename,"flex")
    axt=fig.add_subplot(gs[1,:]); phase_timeline(axt)
    fig.suptitle("D  Boundary conditions & loading",fontsize=12,x=0.09,ha="left")
    fig.savefig(HERE/"fig1D_boundary_conditions.png",dpi=300,bbox_inches="tight",pad_inches=0.05); plt.close(fig)
    print("wrote fig1D_boundary_conditions.png")
    # combined — rows A|B, C|D, and E (revision R1: suture placement and apposition methods) spanning the width
    import subprocess, sys as _s
    subprocess.run([_s.executable, str(HERE/"make_fig1E.py")], check=True, cwd=str(HERE))
    fig=plt.figure(figsize=(17,21))
    gs=fig.add_gridspec(3,2,height_ratios=[1.15,1.0,0.62],width_ratios=[1.0,1.35],hspace=0.04,wspace=0.03)
    for spec,im in ((gs[0,0],"fig1A_whole_model.png"),(gs[0,1],"fig1B_annulotomy_location.png"),
                    (gs[1,0],"fig1C_defect_geometry.png"),(gs[1,1],"fig1D_boundary_conditions.png"),
                    (gs[2,:],"fig1E_suture_apposition.png")):
        a=fig.add_subplot(spec); a.imshow(plt.imread(HERE/im)); a.axis("off")
    fig.savefig(HERE/"fig1_model_construction_combined.png",dpi=220,bbox_inches="tight",pad_inches=0.05); plt.close(fig)
    print("wrote fig1_model_construction_combined.png")
    # settings doc
    (HERE/"figure_render_settings.txt").write_text(
"""Figure 1 (model construction) — render settings, provenance & suggested legend text
Data source (all panels): validation_L45_herniation/comp_intact_flex_map.feb — the intact model, an ALL-HEXAHEDRAL
structured mesh of 33,532 eight-node hexahedral (hex8) solid elements and 38,566 nodes (NOT a hex/tet hybrid;
"hybrid" is only the source file-variant name). Region/defect: region_sets.json (18 IDs), suture_pairs.json (tang/rad).

Camera (1A,1B): posterolateral, view dir (-0.82,-1.0,0.18) normalised, up=+Z; orthographic painter's-algorithm.
Lighting: directional L=(-0.5,-0.8,0.9); shade=0.72+0.28*max(0,n·L). Opaque except 1B faded structures.
Colors: annulus #3a6f9c, nucleus #6bb3b0, facet cart #bcd0d8, cartilage endplate #cfd6d2, bony endplate #cdbf9f,
  posterior arch #c3ae83, cortical #dccfb0, trabecular #e7ddc5. (1C) lamella blue ramp: LAM5 inner=light → LAM1 outer=dark.

1A whole model: all L4/L5/L45 domains, exterior faces, opaque.
1B annulotomy location: SAME camera; posterior arch (L4POST/L5POST) + facet cartilages (_FC) faded to alpha 0.14;
   disc/endplates/bodies opaque; magenta = exterior faces of the 18 defect elements (all in L45_POST_LAM1..5),
   verified on the posterolateral annulus (not bone/facet).
1C defect geometry: axial (X-Y) section of annulus+nucleus elements crossing the defect mid-height (z≈15.4).
   VERIFIED radial order LAM1=outer .. LAM5=inner. Radial arrow runs from the nucleus-facing (inner) edge OUTWARD
   through the 5 lamellae along the defect centreline; circumferential (tangential) ~4.8 mm; out-of-plane height ~4.0 mm (S-I).
1D boundary conditions: coronal (X-Z) silhouette. L5 inferior body fixed all 6 DOF (hatched). L4 superior rigid endplate:
   single central resultant Fz 0→2000 N (red); constraints labelled per panel; flex adds Mx=7.5 N·m (green). Out-of-plane
   axis ⊙ +Y. Phase 1 (t0→1): prescribed follower-pressure preload 0→1 MPa (held); phase 2 (t1→10): mechanical ramp.
Note: 1A/1B share the 3D posterolateral camera; 1C is an axial section, 1D a coronal schematic; each has its own triad.

SUGGESTED FIGURE LEGEND (technical sentences moved out of panels B & C):
 B — "Posterior bony arch and facet cartilages are rendered semi-transparent; the intervertebral disc is opaque. The
      magenta patch marks the 18 removed hex elements (posterior annulus lamellae L45_POST_LAM1–5) on the posterolateral
      annulus surface."
 C — "Axial section through the defect mid-height. The annular defect removes 18 elements spanning all five posterior
      lamellae (LAM1 outer → LAM5 inner). Measured extents: ~4.8 mm circumferential, ~6.5 mm radial, ~4.0 mm superior–inferior."
 Preload — describe as "a model-prescribed 1-MPa follower-pressure preload inherited from the source model" (no osmotic claim).
""")
    print("wrote figure_render_settings.txt")

if __name__=="__main__":
    main()
