#!/usr/bin/env python3
"""Regenerate the defect-location figure from the CURRENT 18-element annulotomy.

Reads the intact L4-L5 mesh (../validation_L45/L45_flexion.feb), highlights the
18 deleted posterior-annulus elements (DEFECT_IDS, identical to defect_spec.json /
make_defect.py), and renders the axial (top) and unrolled-posterior views.
Slot dimensions (height, circumferential width, radial depth) are MEASURED from
the deleted elements' node coordinates and written into the title — no hard-coded
numbers. Output: defect_location.png (overwrites the stale 56-element figure).
"""
import re, math, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SRC = "../validation_L45/L45_flexion.feb"
DEFECT_IDS = set(json.load(open("../validation_L45/defect_spec.json"))["ids"])

txt = open(SRC, encoding="ISO-8859-1").read()

# --- nodes ---
nsec = txt[txt.index("<Nodes"):txt.index("</Nodes>")]
node = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
        for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec)}

# --- disc element parts (region label per Elements block) ---
def region_of(name):
    if "NUCLEUS" in name: return "nucleus"
    if "_ANT_" in name:   return "anterior"
    if "_POST_" in name:  return "posterior"
    if "_L_LAM" in name:  return "left"
    if "_R_LAM" in name:  return "right"
    return None

elems = {}          # eid -> (region, centroid xyz)
for m in re.finditer(r'<Elements type="hex8" mat="\d+" name="([^"]+)">(.*?)</Elements>', txt, re.S):
    reg = region_of(m.group(1))
    if reg is None: continue
    for em in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', m.group(2)):
        eid = int(em.group(1)); ns = [int(v) for v in em.group(2).split(',')]
        c = np.mean([node[n] for n in ns], axis=0)
        elems[eid] = (reg, c, ns)

# disc centre (nucleus centroid) for radial/circumferential geometry
nuc = np.mean([c for r,c,_ in elems.values() if r == "nucleus"], axis=0)
cx, cy = nuc[0], nuc[1]

# --- measure the slot in a LOCAL radial/circumferential frame at the defect ---
# (avoids the angular-span x radius inflation; gives true tangential width &
#  radial depth as the node extents along the local outward-radial and
#  circumferential-tangential unit vectors.)
dnodes = sorted({n for eid in DEFECT_IDS for n in elems[eid][2]})
P = np.array([node[n] for n in dnodes])
dc = np.mean([elems[e][1] for e in DEFECT_IDS], axis=0)   # defect centroid (global)
rhat = np.array([dc[0]-cx, dc[1]-cy]); rhat = rhat/np.hypot(*rhat)   # outward radial
that = np.array([-rhat[1], rhat[0]])                                  # circumferential
Pxy = P[:, :2] - np.array([cx, cy])
height = P[:,2].max() - P[:,2].min()                      # z extent (mm)
radial_depth = (Pxy @ rhat).max() - (Pxy @ rhat).min()    # inner->outer (mm)
width = (Pxy @ that).max() - (Pxy @ that).min()           # tangential width (mm)

COL = {"anterior":"#7FB069","left":"#5AA9C9","posterior":"#9B8BC4","right":"#E0A458","nucleus":"#D7D2CB"}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

# axial (top) view: x-y
for reg in COL:
    pts = np.array([c[:2] for e,(r,c,_) in elems.items() if r==reg and e not in DEFECT_IDS])
    if len(pts): ax1.scatter(pts[:,0], pts[:,1], s=10, c=COL[reg], label=reg, edgecolors='none')
dp = np.array([elems[e][1][:2] for e in DEFECT_IDS])
ax1.scatter(dp[:,0], dp[:,1], s=42, marker='s', facecolors='none', edgecolors='#C0392B', linewidths=1.6, label='defect (18 elem)')
ax1.scatter([cx],[cy], marker='+', c='k', s=60)
ax1.set_aspect('equal'); ax1.set_xlabel("x (mm)"); ax1.set_ylabel("y (mm)")
ax1.set_title("Axial (top) view — L4-L5 disc"); ax1.legend(fontsize=7, loc='upper right'); ax1.grid(alpha=.2)

# unrolled posterior view: circumferential angle vs z
for reg in COL:
    pts = [(math.degrees(math.atan2(c[1]-cy, c[0]-cx)), c[2]) for e,(r,c,_) in elems.items() if r==reg and e not in DEFECT_IDS]
    if pts:
        a = np.array(pts); ax2.scatter(a[:,0], a[:,1], s=10, c=COL[reg], edgecolors='none')
dpa = np.array([(math.degrees(math.atan2(elems[e][1][1]-cy, elems[e][1][0]-cx)), elems[e][1][2]) for e in DEFECT_IDS])
ax2.scatter(dpa[:,0], dpa[:,1], s=42, marker='s', facecolors='none', edgecolors='#C0392B', linewidths=1.6)
ax2.set_xlabel("circumferential angle (deg)"); ax2.set_ylabel("z height (mm)")
ax2.set_title("Unrolled view (angle vs height)"); ax2.grid(alpha=.2)

fig.suptitle(
    f"Posterolateral-left annulotomy: 18 elements, "
    f"~{height:.1f} mm high x ~{width:.1f} mm wide x ~{radial_depth:.1f} mm radial (full thickness)\n"
    f"defect centroid (global) x={dc[0]:.1f}, y={dc[1]:.1f}, z={dc[2]:.1f} mm",
    y=1.04, fontsize=10)
fig.tight_layout()
fig.savefig("defect_location.png", dpi=140, bbox_inches="tight")
print(f"measured slot: height={height:.2f} mm, width(arc@mid-r)={width:.2f} mm, radial_depth={radial_depth:.2f} mm")
print(f"disc centre (nucleus centroid) global = ({cx:.2f}, {cy:.2f}); defect centroid = ({dc[0]:.2f},{dc[1]:.2f},{dc[2]:.2f})")
print("wrote validation_L45_defect/defect_location.png")
