#!/usr/bin/env python3
"""Schematic of the idealized tension-only suture-bridge repair model (axial section).
Shows: (a) the annulotomy — annulus elements deleted and NOT restored; (b) the suture
modelled as discrete tension-only spring elements between selected edge nodes; (c) the
retained nucleus; (d) a note that the baseline repair has no pre-tension and no
geometric closure. No solver run — purely the as-built mesh + suture definition."""
import re, json, math
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patheffects as pe

S = json.load(open("region_sets.json"))
REG_E, DEFECT, C = set(S["elems"]), set(S["defect"]), S["C"]
SUT = json.load(open("suture_config.json"))["suture3"]
base = open("comp_intact_pure_map.feb", encoding="ISO-8859-1").read()
nsec = base[base.index("<Nodes"):base.index("</Nodes>")]
init = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(',')) for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec)}
def elems_of(name):
    mm = re.search(r'<Elements type="hex8" mat="\d+" name="'+name+r'">(.*?)</Elements>', base, re.S)
    return {int(e.group(1)): [int(v) for v in e.group(2).split(',')] for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', mm.group(1))}
ann={}
for q in ["ANT","L","POST","R"]:
    for lam in range(1,6): ann.update(elems_of(f"L45_{q}_LAM{lam}"))
nuc=elems_of("L45_NUCLEUS"); conn={**ann,**nuc}
def cz(e): return sum(init[n][2] for n in conn[e])/8
def face(e): return [(init[n][0],init[n][1]) for n in conn[e][:4]]
zc=C[2]; sl=lambda d:[e for e in d if abs(cz(e)-zc)<=1.4]

fig,ax=plt.subplots(figsize=(9.5,7.2))
# annulus (light blue-grey) and nucleus (pale yellow), kept elements only
ax.add_collection(PolyCollection([face(e) for e in sl(ann) if e not in DEFECT],facecolors='#cdd7e5',edgecolors='#7c8aa0',linewidths=0.4))
ax.add_collection(PolyCollection([face(e) for e in sl(nuc) if e not in DEFECT],facecolors='#f3ecc6',edgecolors='#c9bd86',linewidths=0.4))
# annulotomy: removed elements drawn hatched/empty with magenta dashed outline
ax.add_collection(PolyCollection([face(e) for e in sl(DEFECT)],facecolors='white',edgecolors='magenta',linewidths=1.4,linestyle='--',hatch='xx'))
# suture springs (3) between edge nodes
for a,b in SUT:
    pa,pb=init[a],init[b]
    ax.plot([pa[0],pb[0]],[pa[1],pb[1]],color='yellow',lw=3,solid_capstyle='round',
            path_effects=[pe.Stroke(linewidth=5,foreground='black'),pe.Normal()],zorder=6)
    ax.plot([pa[0],pb[0]],[pa[1],pb[1]],'o',color='yellow',mec='black',ms=8,zorder=7)
# legend proxies
ax.add_patch(plt.Rectangle((0,0),0,0,fc='#cdd7e5',ec='#7c8aa0',label="annulus elements (retained)"))
ax.add_patch(plt.Rectangle((0,0),0,0,fc='#f3ecc6',ec='#c9bd86',label="nucleus (retained, solid)"))
ax.add_patch(plt.Rectangle((0,0),0,0,fc='white',ec='magenta',ls='--',hatch='xx',label="annulotomy: elements deleted (NOT restored)"))
ax.plot([],[],color='yellow',lw=3,path_effects=[pe.Stroke(linewidth=5,foreground='black'),pe.Normal()],marker='o',mec='black',
        label="tension-only spring 'suture bridge'\nbetween selected edge nodes")
ax.legend(loc='upper left',fontsize=9,framealpha=.92)
# callouts
mid=lambda p:((init[p[0]][0]+init[p[1]][0])/2,(init[p[0]][1]+init[p[1]][1])/2)
sx,sy=mid(SUT[0])
ax.annotate("discrete tension-only spring\n(carries load only when the\nedge nodes separate)",
            xy=(sx,sy),xytext=(-15.5,-1.5),fontsize=9,color='black',
            arrowprops=dict(arrowstyle="->",color='black'),
            bbox=dict(boxstyle="round",fc="#fffbe6",ec="#caa"))
hx=sum(init[n][0] for e in sl(DEFECT) for n in conn[e][:4])/(4*len(sl(DEFECT)))
hy=sum(init[n][1] for e in sl(DEFECT) for n in conn[e][:4])/(4*len(sl(DEFECT)))
ax.annotate("annular defect: 18 elements removed;\nno tissue regenerated; gap remains open",
            xy=(hx,hy),xytext=(-15.5,-13.3),fontsize=9,color='black',
            arrowprops=dict(arrowstyle="->",color='magenta'),
            bbox=dict(boxstyle="round",fc="#fff0f6",ec="#caa"))
ax.text(0.5,-0.5,"Baseline repair: NO suture pre-tension,\nNO geometric closure — springs only\nresist further opening under load.",
        transform=ax.transAxes if False else ax.transData, ha='left', va='top', fontsize=9,
        bbox=dict(boxstyle="round",fc="#eef7ee",ec="#9c9"))
ax.set_xlim(-16.5,4); ax.set_ylim(-15,1.5); ax.set_aspect('equal')
ax.set_xlabel("x (mm)"); ax.set_ylabel("y (+posterior)")
ax.set_title("Idealized tension-only suture-bridge repair model (schematic, axial section through the defect)\n"
             "3-suture configuration shown; the 1-suture model uses only the upper spring",fontsize=11)
fig.tight_layout(); fig.savefig("suture_schematic.png",dpi=140,bbox_inches="tight")
print("wrote suture_schematic.png")
