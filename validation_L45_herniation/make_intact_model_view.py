#!/usr/bin/env python3
"""Clean intact L4-L5 model view for manuscript figure panels.

This post-processing script renders the complete intact L4/L5 functional spinal
unit without stress values. It uses the original intact FEBio mesh and marks the
planned annulotomy site only in the annotated export.
"""
from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
import matplotlib.patheffects as pe


HERE = Path(__file__).resolve().parent
FEB = HERE / "comp_intact_flex_map.feb"
REGION = json.loads((HERE / "region_sets.json").read_text())
DEFECT_IDS = set(REGION["defect"])

OUT_CLEAN = HERE / "intact_L45_model_posterolateral_clean.png"
OUT_MARKED = HERE / "intact_L45_model_posterolateral_defect_site.png"
OUT_OUTLINE = HERE / "intact_L45_model_posterolateral_planned_site_outline.png"

FACE_IDX = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def parse_feb(path: Path):
    txt = path.read_text(encoding="ISO-8859-1", errors="ignore")
    nsec = txt[txt.index("<Nodes") : txt.index("</Nodes>")]
    nodes = {
        int(m.group(1)): tuple(float(v) for v in m.group(2).split(","))
        for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec)
    }
    elems = {}
    elem_name = {}
    for block in re.finditer(
        r'<Elements\s+type="hex8"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt,
        re.S,
    ):
        name, body = block.group(1), block.group(2)
        for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', body):
            eid = int(e.group(1))
            elems[eid] = [int(v) for v in e.group(2).split(",")]
            elem_name[eid] = name
    return nodes, elems, elem_name


def selected_domain(name: str) -> bool:
    return name.startswith("L4") or name.startswith("L5") or name.startswith("L45_")


def domain_color(name: str):
    if name.startswith("L45_NUCLEUS"):
        return "#65bfc1"
    if name.startswith("L45_"):
        return "#2f6f9f"
    if "CART" in name or "END" in name:
        return "#d6d1c6"
    return "#e7e3da"


def exterior_faces(nodes, elems, elem_name, keep_ids):
    counts = defaultdict(int)
    owner = {}
    for eid in keep_ids:
        conn = elems[eid]
        for idx in FACE_IDX:
            face = tuple(conn[i] for i in idx)
            key = tuple(sorted(face))
            counts[key] += 1
            owner[key] = (eid, face)
    out = []
    for key, count in counts.items():
        if count == 1:
            eid, face = owner[key]
            if all(n in nodes for n in face):
                out.append((eid, [nodes[n] for n in face], elem_name[eid]))
    return out


def make_camera(nodes):
    pts = np.array(list(nodes.values()))
    center = pts.mean(axis=0)
    # Camera from left-posterolateral side. In this mesh, the annulotomy is in
    # the negative-x / negative-y quadrant, so this view exposes the future defect.
    view = np.array([-0.82, -1.0, 0.18], dtype=float)
    view /= np.linalg.norm(view)
    up_hint = np.array([0.0, 0.0, 1.0])
    right = np.cross(up_hint, view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    up /= np.linalg.norm(up)
    return center, view, right, up


def project(face, center, view, right, up):
    pts = np.array(face) - center
    x = pts @ right
    y = pts @ up
    depth = pts @ view
    return list(zip(x, y)), float(depth.mean())


def face_normal(face):
    p0, p1, p2 = (np.array(p) for p in face[:3])
    n = np.cross(p1 - p0, p2 - p1)
    norm = np.linalg.norm(n)
    return n / norm if norm else n


def draw_model(mark_defect: bool, label: bool = False):
    nodes, elems, elem_name = parse_feb(FEB)
    keep = [eid for eid, name in elem_name.items() if selected_domain(name)]
    faces = exterior_faces(nodes, elems, elem_name, keep)
    center, view, right, up = make_camera(nodes)

    projected = []
    for eid, face, name in faces:
        poly, depth = project(face, center, view, right, up)
        # Simple view-dependent shading; still opaque, just not flat-looking.
        n = face_normal(face)
        light = np.array([-0.5, -0.8, 0.9], dtype=float)
        light /= np.linalg.norm(light)
        shade = 0.72 + 0.28 * max(0.0, float(n @ light))
        projected.append((depth, eid, poly, name, shade))
    projected.sort(key=lambda t: t[0], reverse=True)

    polys = [p for _, _, p, _, _ in projected]
    facecolors = []
    for _, _, _, name, shade in projected:
        base = np.array(plt.matplotlib.colors.to_rgb(domain_color(name)))
        facecolors.append(np.clip(base * shade, 0, 1))

    fig, ax = plt.subplots(figsize=(7.2, 8.0))
    pc = PolyCollection(
        polys,
        facecolors=facecolors,
        edgecolors=(0.42, 0.42, 0.42, 0.45),
        linewidths=0.045,
        closed=True,
    )
    ax.add_collection(pc)

    if mark_defect:
        defect_faces = [f for f in projected if f[1] in DEFECT_IDS]
        # Draw only the visible outer patch outlines. The intact elements remain present.
        for _, _, poly, _, _ in defect_faces:
            xs, ys = zip(*(poly + [poly[0]]))
            ax.plot(
                xs,
                ys,
                color="#d000a8",
                lw=1.1,
                path_effects=[pe.Stroke(linewidth=2.2, foreground="white"), pe.Normal()],
                zorder=5,
            )
        if defect_faces and label:
            allp = np.array([pt for _, _, poly, _, _ in defect_faces for pt in poly])
            anchor = allp.mean(axis=0)
            ax.annotate(
                "planned annulotomy site",
                xy=anchor,
                xytext=(anchor[0] - 14, anchor[1] + 16),
                fontsize=10,
                color="#9b007a",
                arrowprops=dict(arrowstyle="-|>", color="#9b007a", lw=1.2),
                path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            )

    pts2 = np.array([pt for poly in polys for pt in poly])
    pad = 4.0
    ax.set_xlim(pts2[:, 0].min() - pad, pts2[:, 0].max() + pad)
    ax.set_ylim(pts2[:, 1].min() - pad, pts2[:, 1].max() + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    out = OUT_MARKED if (mark_defect and label) else OUT_OUTLINE if mark_defect else OUT_CLEAN
    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    draw_model(mark_defect=False)
    draw_model(mark_defect=True, label=False)
    draw_model(mark_defect=True, label=True)


if __name__ == "__main__":
    main()
