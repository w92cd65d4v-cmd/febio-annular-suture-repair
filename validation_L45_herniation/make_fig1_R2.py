#!/usr/bin/env python3
"""Figure 1 for revision round 2 (Reviewer 1, minor point 3: balanced panels, legible labels at publication size).

All six panels are drawn in ONE matplotlib figure at the final print size (Frontiers full width, 180 mm), so every
font size below is a true point size (minimum 7 pt). Geometry and data are unchanged from the round-1 figure:
the 3D renderer and mesh parsing are imported from make_model_construction_figs.py, and panel E/F data come from
the same files as make_fig1E.py. The old scripts are not modified.

  A  intact L4-L5 functional spinal unit (posterolateral view)
  B  defect location (posterior arch faded; inset zooms on the 18 removed elements)
  C  defect geometry (axial section)
  D  boundary conditions and loading (pure compression, flexion + compression, two-phase timeline)
  E  20 candidate suture anchor pairs on the defect wall
  F  repair representations (open-gap bridge, morphed apposition, tightened + knotted)

Outputs: fig1_R2.png (600 dpi), fig1_R2.pdf, fig1_R2.tiff
"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch, Rectangle, Polygon
import make_model_construction_figs as M
from revision_morph_diagnostics import load

HERE = os.path.dirname(os.path.abspath(__file__))
MM = 1 / 25.4
W, H = 180 * MM, 222 * MM
FS, FS_S, FS_T, FS_L = 7.5, 7.0, 8.0, 10.0          # body, small (minimum), panel title, panel letter (pt)
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.linewidth": 0.6, "mathtext.default": "regular"})
INK, MUTED, GREY = "#1f1f1e", "#555550", "#a3a39d"
C1, C2, NPC, WALL, MAG = "#2a78d6", "#eb6834", "#6bb3b0", "#3a6f9c", "#c0007a"

nodes, elems, ename = M.parse()
fig = plt.figure(figsize=(W, H))
def ax_in(x, y, w, h):                     # axes placed in millimetres from the bottom-left corner
    return fig.add_axes([x / 180, y / 222, w / 180, h / 222])
def panel(x, y, letter, title):            # panel letter + title at a position in millimetres
    fig.text(x / 180, y / 222, letter, fontsize=FS_L, fontweight="bold", va="top", ha="left")
    fig.text((x + 4.2) / 180, (y - 0.4) / 222, title, fontsize=FS_T, va="top", ha="left")

# ---------------------------------------------------------------- A, B: 3D renders
def arrows_sup_post(ax, allp, L=13):
    x0 = allp[:, 0].min() + 4; y0 = allp[:, 1].min() + 4
    for dx, dy in ((0, L), (L, 0)):
        ax.annotate("", xy=(x0 + dx, y0 + dy), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.9, mutation_scale=7))
    ax.text(x0, y0 + L + 1.5, "Sup", ha="center", va="bottom", fontsize=FS_S, color="#444")
    ax.text(x0 + L + 1.5, y0, "Post", ha="left", va="center", fontsize=FS_S, color="#444")

TOP1 = 220
axA = ax_in(0, 150, 60, 64); allp = M.render_3d(axA, nodes, elems, ename); arrows_sup_post(axA, allp)
panel(0, TOP1, "A", "Intact L4–L5 segment")
axB = ax_in(60, 150, 60, 64); allp = M.render_3d(axB, nodes, elems, ename, fade_posterior=True, footprint=True); arrows_sup_post(axB, allp)
panel(60, TOP1, "B", "Defect location")
# inset: zoom on the removed elements
c, v, r, u = M.camera(nodes)
dpts = np.array([pt for e in M.DEFECT for pt in M.project([nodes[n] for n in elems[e]], c, v, r, u)[0]])
dcx, dcy = dpts.mean(0)
axBi = ax_in(92, 152, 27, 27); M.render_3d(axBi, nodes, elems, ename, fade_posterior=True, footprint=True)
axBi.set_xlim(dcx - 9, dcx + 9); axBi.set_ylim(dcy - 9, dcy + 9); axBi.axis("on")
axBi.set_xticks([]); axBi.set_yticks([])
for s in axBi.spines.values(): s.set_edgecolor(MAG); s.set_linewidth(0.8)
axB.add_patch(Rectangle((dcx - 9, dcy - 9), 18, 18, fill=False, ec=MAG, lw=0.8, zorder=9))
axBi.text(0.5, 0.04, "18 removed elements", transform=axBi.transAxes, ha="center", va="bottom", fontsize=FS_S, color="#8a006a",
          bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))

# ---------------------------------------------------------------- C: axial section
axC = ax_in(122, 150, 58, 64); panel(121, TOP1, "C", "Defect geometry (axial section)")
dz = [nodes[n][2] for e in M.DEFECT for n in elems[e]]; zc = (min(dz) + max(dz)) / 2
band = [e for e, n in ename.items() if ("LAM" in n or "NUCLEUS" in n) and min(nodes[k][2] for k in elems[e]) <= zc <= max(nodes[k][2] for k in elems[e])]
Cc = np.array(M.REGION["C"]); nucC = np.mean([nodes[n] for e, nm in ename.items() if "NUCLEUS" in nm for n in elems[e]], axis=0)[:2]
top = lambda e: [(nodes[n][0], nodes[n][1]) for n in elems[e][:4]]
cols = [M.C_NUC if "NUCLEUS" in ename[e] else M.LAM_COLORS[min(max(5 - (int(ename[e][-1]) if ename[e][-1].isdigit() else 3), 0), 4)] for e in band]
axC.add_collection(PolyCollection([top(e) for e in band], facecolors=cols, edgecolors=(0.35, 0.35, 0.35, 0.5), linewidths=0.2))
axC.add_collection(PolyCollection([top(e) for e in M.DEFECT if e in elems and e in band], facecolors="none", edgecolors=MAG, linewidths=1.2, hatch="//"))
dc = Cc[:2]; R = 11
axC.set_xlim(dc[0] - R, dc[0] + R); axC.set_ylim(dc[1] - R * 64 / 58, dc[1] + R * 64 / 58); axC.set_aspect("equal"); axC.axis("off")
dn = np.array([nodes[n] for e in M.DEFECT for n in elems[e]])
tv = (dn[:, :2] - Cc[:2]) @ M.TANG; rv = (dn[:, :2] - Cc[:2]) @ M.RAD; tmid = (tv.min() + tv.max()) / 2
e_in = Cc[:2] + M.TANG * tmid + M.RAD * rv.min(); e_out = Cc[:2] + M.TANG * tmid + M.RAD * rv.max()
if np.hypot(*(e_in - nucC)) > np.hypot(*(e_out - nucC)): e_in, e_out = e_out, e_in
axC.annotate("", xy=e_out, xytext=e_in, zorder=8, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3, shrinkA=0, shrinkB=0, mutation_scale=8))
box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.88)
mid = (e_in + e_out) / 2
axC.annotate("radial depth ≈ 6.5 mm\n(5 lamellae)", xy=mid, xytext=(0.02, 0.86), textcoords="axes fraction", fontsize=FS_S, ha="left", va="center",
             color=INK, bbox=box, zorder=9, arrowprops=dict(arrowstyle="-", color=INK, lw=0.5, shrinkA=0, shrinkB=2))
axC.text(*(e_in - M.RAD * 1.6), "nucleus", fontsize=FS_S, ha="center", va="center", color="#0a6655", bbox=box, zorder=9)
off = M.RAD * (rv.max() + 1.6); a = Cc[:2] + M.TANG * tv.min() + off; b = Cc[:2] + M.TANG * tv.max() + off
axC.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.0, mutation_scale=7))
axC.text(*((a + b) / 2 + M.RAD * 1.9), "circumferential ≈ 4.8 mm", fontsize=FS_S, ha="center", va="center", color=INK, bbox=box)
tx, ty = dc[0] - R + 1.0, dc[1] - R * 64 / 58 + 3.2
for dx, dy, t, ha, va in ((4.5, 0, "+X", "left", "center"), (0, 4.5, "+Y", "center", "bottom")):
    axC.annotate("", xy=(tx + dx, ty + dy), xytext=(tx, ty), arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.8, mutation_scale=6))
    axC.text(tx + dx * 1.15, ty + dy * 1.12, t, fontsize=FS_S, ha=ha, va=va, color="#333")
axC.text(tx + 7.5, ty, "+Z (Sup) out of page", fontsize=FS_S, ha="left", va="center", color="#333")
axC.text(0.5, 0.0, "height ≈ 4.0 mm (superior–inferior)", transform=axC.transAxes, fontsize=FS_S, ha="center", va="bottom", color=INK, bbox=box)

# ---------------------------------------------------------------- D: boundary conditions and loading
panel(0, 147, "D", "Boundary conditions and loading")
keep = [e for e, n in ename.items() if n.startswith(("L4", "L5", "L45_"))]
faces = M.exterior_faces(nodes, elems, ename, keep)
spolys = [[(p[0], p[2]) for p in face] for _, face, _ in faces]
scols = [WALL if "LAM" in nm else NPC if "NUCLEUS" in nm else "#ddd3bf" for _, _, nm in faces]
sp_all = np.array([pt for poly in spolys for pt in poly]); xmin, xmax = sp_all[:, 0].min(), sp_all[:, 0].max(); zmin, zmax = sp_all[:, 1].min(), sp_all[:, 1].max()
xc = (xmin + xmax) / 2
def bc(ax, mode):
    ax.add_collection(PolyCollection(spolys, facecolors=scols, edgecolors=(0.6, 0.6, 0.6, 0.25), linewidths=0.02))
    rpad = 58 if mode == "flex" else 8
    ax.set_xlim(xmin - 30, xmax + rpad); ax.set_ylim(zmin - 13, zmax + 14); ax.set_aspect("equal"); ax.axis("off")
    gy = zmin - 2
    ax.add_patch(Rectangle((xmin - 5, gy - 3), (xmax - xmin) + 10, 3, facecolor="none", edgecolor="#333", hatch="////", lw=0.6))
    ax.plot([xmin - 5, xmax + 5], [gy, gy], color="#333", lw=0.9)
    ax.annotate("L5 inferior body fixed (all six DOF)", xy=(xc, gy - 3), xytext=(0, -2), textcoords="offset points", fontsize=FS_S, ha="center", va="top", color=INK)
    ax.plot([xmin - 3, xmax + 3], [zmax + 2, zmax + 2], color="#7a5a2a", lw=2.0, solid_capstyle="round")
    ax.annotate("L4 rigid\nendplate", xy=(xmin - 3, zmax + 2), xytext=(-3, 0), textcoords="offset points", fontsize=FS_S, ha="right", va="center", color="#5a3d18")
    ax.annotate("", xy=(xc, zmax + 2.5), xytext=(xc, zmax + 13), arrowprops=dict(arrowstyle="-|>", color="#b02020", lw=1.6, mutation_scale=9))
    ax.annotate(r"$F_z$ = 0 → 2,000 N", xy=(xc, zmax + 9), xytext=(5, 0), textcoords="offset points", fontsize=FS, ha="left", va="center", color="#b02020", fontweight="bold")
    if mode == "flex":
        ax.add_patch(FancyArrowPatch((xmax + 5, zmax - 3), (xmax + 5, zmin + 7), connectionstyle="arc3,rad=-0.5", arrowstyle="-|>",
                                     mutation_scale=9, lw=1.3, color="#1f6f3f"))
        ax.text(xmax + 31, (zmin + zmax) / 2, r"$M_x$ = 7.5 N m" + "\n(flexion)", fontsize=FS, ha="left", va="center", color="#1f6f3f", fontweight="bold")
        sub = "L4: fixed x, y, Ry, Rz; free z, Rx"
    else:
        sub = "L4: fixed x, y, Rx, Ry, Rz; free z"
    ax.annotate(sub, xy=(xc, gy - 3), xytext=(0, -11), textcoords="offset points", fontsize=FS_S, ha="center", va="top", color=MUTED, style="italic")
    tx, ty = xmin - 26, zmin + 8
    for dx, dy, t, ha, va in ((9, 0, "+X", "left", "center"), (0, 9, "+Z", "center", "bottom")):
        ax.annotate("", xy=(tx + dx, ty + dy), xytext=(tx, ty), arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.7, mutation_scale=6))
        ax.text(tx + dx * 1.15, ty + dy * 1.1, t, fontsize=FS_S, ha=ha, va=va, color="#333")
    ax.annotate("+Y out of page", xy=(tx, ty), xytext=(0, -3), textcoords="offset points", fontsize=FS_S, ha="left", va="top", color="#333")
axDp = ax_in(0, 93, 80, 44); bc(axDp, "pure"); fig.text(40 / 180, 141 / 222, "Pure axial compression", fontsize=FS_T, fontweight="bold", ha="center", va="top")
axDf = ax_in(80, 93, 100, 44); bc(axDf, "flex"); fig.text(118 / 180, 141 / 222, "Flexion + axial compression", fontsize=FS_T, fontweight="bold", ha="center", va="top")
axT = ax_in(6, 76, 168, 13); axT.set_xlim(0, 10); axT.set_ylim(0, 1); axT.axis("off")
axT.add_patch(Rectangle((0, 0.42), 1.6, 0.36, facecolor="#cfe0ef", edgecolor="#333", lw=0.6))
axT.add_patch(Rectangle((1.6, 0.42), 8.4, 0.36, facecolor="#f3ddd0", edgecolor="#333", lw=0.6))
axT.text(0.8, 0.6, "Phase 1", ha="center", va="center", fontsize=FS, fontweight="bold")
axT.text(5.8, 0.6, "Phase 2", ha="center", va="center", fontsize=FS, fontweight="bold")
axT.text(0.8, 0.36, "1-MPa expansive nucleus-\nboundary pressure preload", ha="center", va="top", fontsize=FS_S)
axT.text(5.8, 0.36, r"mechanical load ramp: $F_z$ (and $M_x$) from 0 to maximum", ha="center", va="top", fontsize=FS_S)
for x_, t_ in ((0, "t = 0"), (1.6, "t = 1"), (10, "t = 10")): axT.text(x_, 0.84, t_, ha="center", va="bottom", fontsize=FS_S, color=MUTED)
axT.text(5.8, 0.84, "pseudo-time", ha="center", va="bottom", fontsize=FS_S, color=MUTED)

# ---------------------------------------------------------------- E: candidate anchor pairs
panel(0, 72, "E", "Candidate suture anchor pairs (20)")
sp = json.load(open(os.path.join(HERE, "suture_pairs.json"))); P = sp["pairs"]; Cs = np.array(sp["C"]); rd = np.r_[sp["rad"], 0]
S3 = [tuple(p) for p in json.load(open(os.path.join(HERE, "suture_config.json")))["suture3"]]
OUT = sorted({0, 6, 11, 16})                         # union of the outer-half layouts: outer1 [11], outer2 [6, 16], outer3 [0, 11, 16]
Nd, *_ = load(os.path.join(HERE, "comp_defect_pure_map.feb"))
pos = [((Nd[a] + Nd[b]) / 2 - Cs) for a, b in P]; zz = np.array([p[2] for p in pos]); rr = np.array([p @ rd for p in pos])
axE = ax_in(14, 16, 64, 50)
axE.axvspan(-2.6, -1.6, color=NPC, alpha=0.18, lw=0); axE.text(-2.55, 1.42, "deep\n(NP side)", fontsize=FS_S, color="#2f6f6c", va="top")
axE.axvspan(0.8, 2.0, color=C2, alpha=0.08, lw=0); axE.text(1.95, 1.42, "outer half", fontsize=FS_S, color="#a5461c", va="top", ha="right")
oth = [i for i in range(len(P)) if tuple(P[i]) not in S3 and i not in OUT]
axE.scatter(rr[oth], zz[oth], s=14, facecolor="white", edgecolor=GREY, lw=0.8, label="other candidate pairs", zorder=3)
axE.scatter(rr[OUT], zz[OUT], s=20, marker="s", color=C2, edgecolor="white", lw=0.5, label="outer-half layouts (1–3 sutures)", zorder=4)
tg = [P.index(list(p)) if list(p) in P else P.index([p[1], p[0]]) for p in S3]
axE.scatter(rr[tg], zz[tg], s=42, color=C1, edgecolor="white", lw=0.5, label="load-targeted pairs (rank 1–3)", zorder=5)
for k, i in enumerate(tg): axE.text(rr[i], zz[i], str(k + 1), ha="center", va="center", color="white", fontsize=6.0, fontweight="bold", zorder=6)
axE.set_xlim(-2.6, 2.0); axE.set_ylim(-2.45, 1.5); axE.set_yticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
axE.set_xlabel("Radial position on defect wall (mm)\nNP side ←      → outer surface", fontsize=FS)
axE.set_ylabel("Cranio-caudal position (mm)", fontsize=FS)
for s_ in ("top", "right"): axE.spines[s_].set_visible(False)
axE.grid(color="#e6e6e2", lw=0.4); axE.tick_params(labelsize=FS_S, colors=MUTED, width=0.5, length=2.5)
axE.legend(loc="lower center", fontsize=FS_S, frameon=True, framealpha=0.9, edgecolor="none", borderpad=0.3, handletextpad=0.3,
           labelspacing=0.25, bbox_to_anchor=(0.5, 0.0), ncol=1)

# ---------------------------------------------------------------- F: repair representations
panel(84, 72, "F", "Repair representations")
bx = ax_in(84, 2, 96, 64); bx.set_xlim(0, 33); bx.set_ylim(-2.4, 11.4); bx.axis("off")
def wall(x0, w, side, bulge=0.0):
    ys = np.linspace(1, 8, 30); face = x0 - side * bulge * np.sin(np.pi * (ys - 1) / 7)
    poly = [(x0 - w, 1), *zip(face, ys), (x0 - w, 8)] if side < 0 else [(x0 + w, 1), *zip(face, ys), (x0 + w, 8)]
    bx.add_patch(Polygon(poly, closed=True, fc=WALL, ec="white", lw=0.6, alpha=0.9))
def cap(x, t): bx.text(x, 0.4, t, ha="center", va="top", fontsize=FS_S, color=INK)
def note(x, t): bx.text(x, 9.1, t, ha="center", va="bottom", fontsize=FS_S, color=MUTED)
wall(4.0, 2.6, -1); wall(7.0, 2.6, +1)
bx.plot([3.4, 7.6], [5.5, 5.5], color=C2, lw=1.4); bx.plot([3.4, 7.6], [5.5, 5.5], "o", color=C2, ms=3)
bx.annotate("", xy=(4.0, 8.5), xytext=(7.0, 8.5), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.7, mutation_scale=6))
note(5.5, "open gap"); cap(5.5, "Open-gap bridge:\ntension-only spring,\nno pretension")
wall(15.6, 2.6, -1, bulge=-0.9); wall(17.4, 2.6, +1, bulge=-0.9)
bx.plot([15.9, 17.1], [4.5, 4.5], color=C2, lw=1.4); bx.plot([15.9, 17.1], [4.5, 4.5], "o", color=C2, ms=3)
note(16.5, "rest geometry morphed,\nreset as stress-free"); cap(16.5, "Morphed apposition:\nprimary series\n(app17–app06)")
wall(27.4, 2.6, -1, bulge=-0.75); wall(29.2, 2.6, +1, bulge=-0.75)
ba = [(27.75, 4.0), (27.85, 4.6), (27.75, 5.2)]; bb = [(28.85, 4.2), (28.75, 4.8)]
for p in ba:
    for q in bb: bx.plot([p[0], q[0]], [p[1], q[1]], color=C2, lw=0.7)
bx.plot(*zip(*ba), "o", color=C2, ms=2.6); bx.plot(*zip(*bb), "o", color=C2, ms=2.6)
for y in (2.0, 7.0):
    bx.annotate("", xy=(27.95, y), xytext=(26.0, y), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=6))
    bx.annotate("", xy=(28.65, y), xytext=(30.6, y), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=6))
note(28.3, "bites pulled together,\nknotted, released;\nwall contact enabled"); cap(28.3, "Tightened + knotted:\nforce-based check\n(app17 and app12 gaps)")
for x in (11.0, 22.4): bx.plot([x, x], [0.6, 9.6], color="#e0dfd8", lw=0.7)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {}), ("tiff", dict(dpi=600, pil_kwargs={"compression": "tiff_lzw"}))):
    fig.savefig(os.path.join(HERE, f"fig1_R2.{ext}"), **kw)
print("wrote fig1_R2.png/.pdf/.tiff", f"({W/MM:.0f} x {H/MM:.0f} mm)")
