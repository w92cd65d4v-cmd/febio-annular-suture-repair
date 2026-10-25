#!/usr/bin/env python3
"""Posterior whole-FSU context render for manuscript review.

This is a post-processing figure only. It uses the already completed FEBio
models/logs: L4/L5 geometry is rendered as a translucent context surface, while
available von Mises output is mapped on the annulus/defect-region elements.
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
import matplotlib.patheffects as pe
import numpy as np
from matplotlib import cm, colors
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


HERE = Path(__file__).resolve().parent
CC01 = Path("/Volumes/storage01/forCC/cc01/validation_L45_herniation")
if not CC01.is_dir(): CC01 = Path(__file__).resolve().parent / "nas_copy"  # local copy 2026-09-27
OUT = HERE / "whole_fsu_stress_context_flex.png"
OUT2D = HERE / "whole_fsu_stress_posterior_flex.png"

REGION = json.loads((HERE / "region_sets.json").read_text())
DEFECT_IDS = set(REGION["defect"])
REGION_IDS = set(REGION["elems"])
SUTURE_CONFIG = json.loads((HERE / "suture_config.json").read_text())

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


def vm(vals):
    sx, sy, sz, sxy, syz, sxz = vals[:6]
    return math.sqrt(
        0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + sxz**2)
    )


def final_vm(path: Path):
    txt = path.read_text(encoding="ISO-8859-1", errors="ignore")
    out = {}
    blocks = re.split(r"Data Record #\d+", txt)
    for block in blocks:
        md = re.search(r"Data\s*=\s*sx;sy;sz;sxy;syz;sxz", block)
        if not md:
            continue
        local = {}
        for line in block[md.end() :].splitlines():
            parts = line.split()
            if len(parts) >= 7 and parts[0].isdigit():
                local[int(parts[0])] = vm([float(x) for x in parts[1:7]])
        if local:
            out = local
    return out


def surface_faces(elems, nodes, keep_ids):
    owner = {}
    counts = defaultdict(int)
    for eid in keep_ids:
        conn = elems.get(eid)
        if not conn:
            continue
        for idx in FACE_IDX:
            face = tuple(conn[i] for i in idx)
            key = tuple(sorted(face))
            counts[key] += 1
            owner[key] = (eid, face)
    faces = []
    for key, count in counts.items():
        if count == 1:
            eid, face = owner[key]
            if all(n in nodes for n in face):
                faces.append((eid, [nodes[n] for n in face]))
    return faces


def face_normal_y(face):
    p0, p1, p2 = (np.array(p) for p in face[:3])
    n = np.cross(p1 - p0, p2 - p1)
    return float(n[1])


def project_face(face):
    # Show posterior anatomy as an x-z view with shallow y-depth, which avoids the
    # unreadable occlusion that a fully opaque 3D posterior mesh creates.
    return [(x, z, -0.16 * y) for x, y, z in face]


def project_point(point):
    x, y, z = point
    return x, z, -0.16 * y


def elem_center(elems, nodes, eid):
    pts = [nodes[n] for n in elems[eid]]
    return np.mean(np.array(pts), axis=0)


def name_is_bone(name):
    return (name.startswith("L4") or name.startswith("L5")) and not name.startswith("L45")


def name_is_disc(name):
    return name.startswith("L45_")


def add_collection(ax, faces, *, color=None, colors_list=None, alpha=1.0, lw=0.05, edge="#555"):
    poly = Poly3DCollection(
        [project_face(f) for f in faces],
        facecolors=colors_list if colors_list is not None else color,
        edgecolors=edge,
        linewidths=lw,
        alpha=alpha,
    )
    ax.add_collection3d(poly)
    return poly


def project2(face):
    return [(x, z) for x, _y, z in face]


def add_poly2d(ax, faces, *, color=None, colors_list=None, alpha=1.0, lw=0.05, edge="#555"):
    pc = PolyCollection(
        [project2(f) for f in faces],
        facecolors=colors_list if colors_list is not None else color,
        edgecolors=edge,
        linewidths=lw,
        alpha=alpha,
    )
    ax.add_collection(pc)
    return pc


def draw_panel(ax, title, feb, log, *, suture_pairs=None, apposed=False, cmap_norm=None):
    nodes, elems, elem_name = parse_feb(feb)
    stress = final_vm(log)

    bone_ids = [eid for eid, name in elem_name.items() if name_is_bone(name)]
    disc_ids = [eid for eid, name in elem_name.items() if name_is_disc(name)]
    logged_disc_ids = [eid for eid in disc_ids if eid in stress]
    region_disc_ids = [eid for eid in disc_ids if eid in REGION_IDS and eid not in DEFECT_IDS]

    bone_faces = surface_faces(elems, nodes, bone_ids)
    # Keep mainly posterior-facing bone surface for a clean posterior view.
    bone_faces = [f for _, f in bone_faces if face_normal_y(f) > -1e-6]
    if len(bone_faces) > 3400:
        step = math.ceil(len(bone_faces) / 3400)
        bone_faces = bone_faces[::step]

    disc_faces = surface_faces(elems, nodes, region_disc_ids or disc_ids)
    disc_faces = [(eid, f) for eid, f in disc_faces if face_normal_y(f) > -1e-6]
    stressed_faces = [(eid, f) for eid, f in disc_faces if eid in stress]
    quiet_faces = [f for eid, f in disc_faces if eid not in stress]

    add_collection(ax, bone_faces, color=(0.76, 0.76, 0.76, 0.13), alpha=0.13, lw=0.02, edge="#999")
    if quiet_faces:
        add_collection(ax, quiet_faces, color=(0.64, 0.68, 0.70, 0.22), alpha=0.22, lw=0.03, edge="#777")

    mapper = cm.ScalarMappable(norm=cmap_norm, cmap="turbo")
    if stressed_faces:
        face_cols = [mapper.to_rgba(stress[eid]) for eid, _ in stressed_faces]
        add_collection(ax, [f for _, f in stressed_faces], colors_list=face_cols, alpha=0.98, lw=0.08, edge="#222")

    # Draw the original defect footprint. In intact it is a planned site; in
    # operated models it marks the annulotomy boundary.
    intact_nodes, intact_elems, _ = parse_feb(HERE / "comp_intact_flex_map.feb")
    defect_faces = surface_faces(intact_elems, intact_nodes, DEFECT_IDS)
    defect_faces = [f for _, f in defect_faces if face_normal_y(f) > -1e-6]
    for f in defect_faces:
        xs, zs, ys = zip(*[project_point(p) for p in f + [f[0]]])
        if title.lower().startswith("intact"):
            ax.plot(xs, zs, ys, color="#333", lw=0.7, ls=(0, (4, 3)), zorder=8)
        else:
            ax.plot(xs, zs, ys, color="#d000a8", lw=0.9, zorder=8)

    if suture_pairs:
        for a, b in suture_pairs:
            if a in nodes and b in nodes:
                pa, pb = nodes[a], nodes[b]
            elif a in intact_nodes and b in intact_nodes:
                pa, pb = intact_nodes[a], intact_nodes[b]
            else:
                continue
            xs, zs, ys = zip(*[project_point(pa), project_point(pb)])
            ax.plot(
                xs,
                zs,
                ys,
                color="yellow",
                lw=2.2,
                solid_capstyle="round",
                path_effects=[pe.Stroke(linewidth=3.8, foreground="black"), pe.Normal()],
                zorder=10,
            )
            ax.scatter(xs, zs, ys, s=12, color="yellow", edgecolor="black", zorder=11)

    centers = np.array([elem_center(elems, nodes, eid) for eid in (bone_ids + disc_ids) if eid in elems])
    mins, maxs = centers.min(axis=0), centers.max(axis=0)
    ax.set_xlim(mins[0] - 3, maxs[0] + 3)
    ax.set_ylim(mins[2] - 3, maxs[2] + 3)
    ax.set_zlim(-0.16 * maxs[1] - 2, -0.16 * mins[1] + 2)
    ax.view_init(elev=8, azim=-90)
    ax.set_box_aspect((1.0, 1.55, 0.25))
    ax.set_title(title, fontsize=10, pad=6)
    ax.set_axis_off()
    ax.text2D(
        0.02,
        0.03,
        f"colored elements: {len(logged_disc_ids)}",
        transform=ax.transAxes,
        fontsize=7,
        color="#333",
    )
    if apposed:
        ax.text2D(0.02, 0.095, "apposed-edge geometry", transform=ax.transAxes, fontsize=7, color="#333")
    return stress


def draw_panel_2d(ax, title, feb, log, *, suture_pairs=None, apposed=False, cmap_norm=None):
    nodes, elems, elem_name = parse_feb(feb)
    stress = final_vm(log)

    bone_ids = [eid for eid, name in elem_name.items() if name_is_bone(name)]
    disc_ids = [eid for eid, name in elem_name.items() if name_is_disc(name)]
    region_disc_ids = [eid for eid in disc_ids if eid in REGION_IDS and eid not in DEFECT_IDS]
    logged_disc_ids = [eid for eid in disc_ids if eid in stress]

    bone_surface = surface_faces(elems, nodes, bone_ids)
    if bone_surface:
        cy = np.array([np.mean([p[1] for p in f]) for _, f in bone_surface])
        posterior_cut = np.percentile(cy, 45)
        bone_faces = [f for (_, f), y in zip(bone_surface, cy) if y <= posterior_cut]
    else:
        bone_faces = []
    if len(bone_faces) > 4800:
        step = math.ceil(len(bone_faces) / 4800)
        bone_faces = bone_faces[::step]

    disc_surface = surface_faces(elems, nodes, region_disc_ids or disc_ids)
    if disc_surface:
        cy = np.array([np.mean([p[1] for p in f]) for _, f in disc_surface])
        posterior_cut = np.percentile(cy, 60)
        disc_surface = [(ef, y) for ef, y in zip(disc_surface, cy) if y <= posterior_cut]
    disc_faces = [(eid, f) for (eid, f), _y in disc_surface]
    stressed_faces = [(eid, f) for eid, f in disc_faces if eid in stress]
    quiet_faces = [f for eid, f in disc_faces if eid not in stress]

    add_poly2d(ax, bone_faces, color=(0.72, 0.72, 0.72, 0.16), alpha=0.16, lw=0.02, edge="#a8a8a8")
    if quiet_faces:
        add_poly2d(ax, quiet_faces, color=(0.60, 0.66, 0.68, 0.28), alpha=0.28, lw=0.03, edge="#777")

    mapper = cm.ScalarMappable(norm=cmap_norm, cmap="turbo")
    if stressed_faces:
        face_cols = [mapper.to_rgba(stress[eid]) for eid, _ in stressed_faces]
        add_poly2d(ax, [f for _, f in stressed_faces], colors_list=face_cols, alpha=0.98, lw=0.08, edge="#222")

    intact_nodes, intact_elems, _ = parse_feb(HERE / "comp_intact_flex_map.feb")
    defect_faces = surface_faces(intact_elems, intact_nodes, DEFECT_IDS)
    if defect_faces:
        cy = np.array([np.mean([p[1] for p in f]) for _, f in defect_faces])
        posterior_cut = np.percentile(cy, 70)
        defect_faces = [f for (_, f), y in zip(defect_faces, cy) if y <= posterior_cut]
    for f in defect_faces:
        xs, zs = zip(*project2(f + [f[0]]))
        if title.lower().startswith("intact"):
            ax.plot(xs, zs, color="#333", lw=0.7, ls=(0, (4, 3)), zorder=8)
        else:
            ax.plot(xs, zs, color="#d000a8", lw=0.9, zorder=8)

    if suture_pairs:
        for a, b in suture_pairs:
            if a in nodes and b in nodes:
                pa, pb = nodes[a], nodes[b]
            elif a in intact_nodes and b in intact_nodes:
                pa, pb = intact_nodes[a], intact_nodes[b]
            else:
                continue
            ax.plot(
                [pa[0], pb[0]],
                [pa[2], pb[2]],
                color="yellow",
                lw=2.0,
                solid_capstyle="round",
                path_effects=[pe.Stroke(linewidth=3.5, foreground="black"), pe.Normal()],
                zorder=10,
            )
            ax.scatter([pa[0], pb[0]], [pa[2], pb[2]], s=18, color="yellow", edgecolor="black", zorder=11)

    pts = np.array(list(nodes.values()))
    ax.set_xlim(np.percentile(pts[:, 0], 0.5) - 2.5, np.percentile(pts[:, 0], 99.5) + 2.5)
    ax.set_ylim(np.percentile(pts[:, 2], 0.5) - 2.5, np.percentile(pts[:, 2], 99.5) + 2.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    ax.text(
        0.02,
        0.03,
        f"colored elements: {len(logged_disc_ids)}",
        transform=ax.transAxes,
        fontsize=7,
        color="#333",
    )
    if apposed:
        ax.text(0.02, 0.085, "apposed-edge geometry", transform=ax.transAxes, fontsize=7, color="#333")


def main():
    panels = [
        (
            "Intact",
            HERE / "comp_intact_flex_map.feb",
            HERE / "comp_intact_flex_map.log",
            None,
            False,
        ),
        (
            "Defect",
            HERE / "comp_defect_flex_map.feb",
            HERE / "comp_defect_flex_map.log",
            None,
            False,
        ),
        (
            "3-suture bridge",
            HERE / "comp_suture3_flex_map.feb",
            HERE / "comp_suture3_flex_map.log",
            SUTURE_CONFIG["suture3"],
            False,
        ),
        (
            "3-suture apposed app12",
            CC01 / "comp_suture3_app12_flex.feb",
            CC01 / "comp_suture3_app12_flex.log",
            SUTURE_CONFIG["suture3"],
            True,
        ),
    ]
    all_vals = []
    for _, _, log, _, _ in panels:
        all_vals.extend(final_vm(log).values())
    norm = colors.Normalize(vmin=0, vmax=np.percentile(all_vals, 98))

    fig = plt.figure(figsize=(16, 7.2))
    axes = [fig.add_subplot(1, 4, i + 1, projection="3d") for i in range(4)]
    for ax, (title, feb, log, sutures, apposed) in zip(axes, panels):
        draw_panel(ax, title, feb, log, suture_pairs=sutures, apposed=apposed, cmap_norm=norm)

    sm = cm.ScalarMappable(norm=norm, cmap="turbo")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.55, pad=0.015, fraction=0.025)
    cbar.set_label("von Mises stress in logged annulus/rim elements (MPa)")
    fig.suptitle(
        "Posterior whole-L4/L5 context view with defect-region von Mises stress at 2000 N flexion + compression\n"
        "grey = full vertebral/disc context; dashed = planned defect in intact; magenta = annulotomy footprint; yellow = suture springs",
        fontsize=11,
        y=0.98,
    )
    fig.savefig(OUT, dpi=220, bbox_inches="tight")
    print(f"wrote {OUT}")

    fig2, axes2 = plt.subplots(1, 4, figsize=(15.5, 5.8), constrained_layout=True)
    for ax, (title, feb, log, sutures, apposed) in zip(axes2, panels):
        draw_panel_2d(ax, title, feb, log, suture_pairs=sutures, apposed=apposed, cmap_norm=norm)
    sm2 = cm.ScalarMappable(norm=norm, cmap="turbo")
    sm2.set_array([])
    cbar2 = fig2.colorbar(sm2, ax=axes2, shrink=0.78, pad=0.012, fraction=0.025)
    cbar2.set_label("von Mises stress in logged annulus/rim elements (MPa)")
    fig2.suptitle(
        "Posterior full-L4/L5 projection with defect-region von Mises stress at 2000 N flexion + compression\n"
        "grey = vertebral context; dashed = planned defect in intact; magenta = annulotomy footprint; yellow = suture springs",
        fontsize=11,
    )
    fig2.savefig(OUT2D, dpi=240, bbox_inches="tight")
    print(f"wrote {OUT2D}")


if __name__ == "__main__":
    main()
