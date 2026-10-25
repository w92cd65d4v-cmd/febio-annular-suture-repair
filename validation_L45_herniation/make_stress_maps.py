#!/usr/bin/env python3
"""Axial-section annular von Mises maps for the Results figure.

Rows are pure compression and flexion plus compression. Columns are intact,
unrepaired defect, load-targeted open-gap one-suture and three-suture bridges,
and app12 partial apposition. Each panel uses its own FEBio mesh. The section
plane intersects the hex8 elements geometrically, so intact retains the 18
annulotomy elements while defect and repair panels show the actual empty
opening. No mechanical simulation is rerun.
"""
import json
import math
import os
import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.patheffects as pe
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_figures")
REGION = json.load(open(os.path.join(HERE, "region_sets.json")))
REG_E = set(REGION["elems"])
DEF = set(REGION["defect"])
ANN_E = REG_E | DEF
Z_SECTION = float(REGION["C"][2])
SUTURE_CONFIG = json.load(open(os.path.join(HERE, "suture_config.json")))
SUT1 = SUTURE_CONFIG["suture1"]
SUT3 = SUTURE_CONFIG["suture3"]

HEX_EDGES = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
)
LEVELS = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.42])
CMAP = ListedColormap([
    "#2445A2", "#356FD0", "#32A5E6", "#28C8C2", "#55D879",
    "#9ADF45", "#E1E33A", "#F7B53B", "#F0642E", "#A50026",
])
NORM = BoundaryNorm(LEVELS, CMAP.N, clip=True)


def von_mises(s):
    sx, sy, sz, sxy, syz, sxz = s
    return math.sqrt(
        0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
        + 3 * (sxy ** 2 + syz ** 2 + sxz ** 2)
    )


def parse_feb(path):
    text = open(path, encoding="ISO-8859-1", errors="ignore").read()
    node_section = text[text.index("<Nodes"):text.index("</Nodes>")]
    nodes = {
        int(m.group(1)): np.array([float(v) for v in m.group(2).split(",")])
        for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', node_section)
    }
    elems = {}
    names = {}
    for block in re.finditer(
        r'<Elements\s+type="hex8"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        text,
        re.S,
    ):
        name = block.group(1)
        for item in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', block.group(2)):
            eid = int(item.group(1))
            elems[eid] = [int(v) for v in item.group(2).split(",")]
            names[eid] = name
    return nodes, elems, names


def endpoint_vm(path):
    text = open(path, encoding="ISO-8859-1", errors="ignore").read()
    records = {}
    for block in re.split(r"Data Record #\d+", text):
        mt = re.search(r"Time\s*=\s*([\d.eE+-]+)", block)
        md = re.search(r"Data\s*=\s*([a-zA-Z;]+)", block)
        if not (mt and md) or not md.group(1).startswith("sx;sy;sz;sxy"):
            continue
        data = records.setdefault(float(mt.group(1)), {})
        for line in re.finditer(r"^\s*(\d+)\s+([-\d.eE+ ]+)$", block, re.M):
            values = [float(x) for x in line.group(2).split()]
            if len(values) >= 6:
                data[int(line.group(1))] = von_mises(values[:6])
    return records[max(records)]


def section_polygon(connectivity, nodes, z0, tol=1e-8):
    """Intersect one hex8 with z=z0 and return its ordered x-y polygon."""
    points = []
    xyz = [nodes[n] for n in connectivity]
    for i, j in HEX_EDGES:
        a, b = xyz[i], xyz[j]
        da, db = a[2] - z0, b[2] - z0
        if abs(da) <= tol and abs(db) <= tol:
            points.extend((a[:2], b[:2]))
        elif abs(da) <= tol:
            points.append(a[:2])
        elif abs(db) <= tol:
            points.append(b[:2])
        elif da * db < 0:
            t = (z0 - a[2]) / (b[2] - a[2])
            points.append((a + t * (b - a))[:2])
    unique = {}
    for p in points:
        unique[(round(float(p[0]), 7), round(float(p[1]), 7))] = np.asarray(p)
    if len(unique) < 3:
        return None
    poly = np.array(list(unique.values()))
    center = poly.mean(axis=0)
    order = np.argsort(np.arctan2(poly[:, 1] - center[1], poly[:, 0] - center[0]))
    return poly[order]


def opening_boundary(nodes, intact_connectivity):
    """Return only the external section edges of the deleted element chain."""
    edges = Counter()
    original = {}
    for eid in DEF:
        poly = section_polygon(intact_connectivity[eid], nodes, Z_SECTION)
        if poly is None:
            continue
        for a, b in zip(poly, np.roll(poly, -1, axis=0)):
            ka = (round(float(a[0]), 5), round(float(a[1]), 5))
            kb = (round(float(b[0]), 5), round(float(b[1]), 5))
            key = tuple(sorted((ka, kb)))
            edges[key] += 1
            original[key] = (a, b)
    return [original[key] for key, count in edges.items() if count == 1]


BASE_NODES, BASE_ELEMS, _ = parse_feb(os.path.join(HERE, "comp_intact_pure.feb"))


def draw_panel(ax, model, feb, log_path):
    nodes, elems, names = parse_feb(feb)
    stress = endpoint_vm(log_path)

    # Nucleus is neutral context; only annulus elements carry the stress scale.
    nucleus_polys = []
    for eid, conn in elems.items():
        if "NUCLEUS" not in names[eid]:
            continue
        poly = section_polygon(conn, nodes, Z_SECTION)
        if poly is not None:
            nucleus_polys.append(poly)
    if nucleus_polys:
        ax.add_collection(PolyCollection(
            nucleus_polys,
            facecolors="#eeeeee",
            edgecolors="#cccccc",
            linewidths=0.25,
            zorder=0,
        ))

    polys = []
    values = []
    missing = []
    for eid in sorted(ANN_E):
        if eid not in elems or "LAM" not in names[eid]:
            continue
        poly = section_polygon(elems[eid], nodes, Z_SECTION)
        if poly is None:
            continue
        polys.append(poly)
        if eid not in stress:
            missing.append(eid)
        else:
            values.append(stress[eid])
    if missing:
        raise RuntimeError(
            f"{model}: stress logger omitted section elements {sorted(missing)}; "
            "missing values must never be rendered as zero"
        )
    colors = [CMAP(NORM(value)) for value in values]
    ax.add_collection(PolyCollection(
        polys,
        facecolors=colors,
        edgecolors=(0.10, 0.10, 0.10, 0.48),
        linewidths=0.28,
        zorder=2,
    ))

    # Only non-intact panels receive an opening boundary. The interior stays
    # white because those elements are genuinely absent from each repair mesh.
    if model != "intact":
        for a, b in opening_boundary(nodes, BASE_ELEMS):
            ax.plot(
                [a[0], b[0]], [a[1], b[1]],
                color="#d000b8", lw=1.45, zorder=6,
            )

    suture_pairs = {
        "suture1": SUT1,
        "bridge": SUT3,
        "app12": SUT3,
    }.get(model, [])
    if suture_pairs:
        for a, b in suture_pairs:
            if a in nodes and b in nodes:
                pa, pb = nodes[a], nodes[b]
                ax.plot(
                    [pa[0], pb[0]], [pa[1], pb[1]],
                    color="#ffd400", lw=1.9, solid_capstyle="round", zorder=8,
                    path_effects=[pe.Stroke(linewidth=3.1, foreground="#222222"), pe.Normal()],
                )

    ax.set_xlim(-15.5, 1.5)
    ax.set_ylim(-14.5, -1.0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.55)


def main():
    labels = {
        "intact": "Intact",
        "defect": "Unrepaired defect",
        "suture1": "Load-targeted open-gap\n1-suture bridge",
        "bridge": "Load-targeted open-gap\n3-suture bridge",
        "app12": "Load-targeted partial apposition\n(app12; mean gap 1.186 mm)",
    }
    febs = {
        "intact": "comp_intact_%s.feb",
        "defect": "comp_defect_%s.feb",
        "suture1": "comp_suture1_%s_map.feb",
        "bridge": "comp_suture3_%s.feb",
        "app12": "comp_suture3_app12_%s_map.feb",
    }
    logs = {
        "intact": "comp_intact_%s_fullstress.log",
        "defect": "comp_defect_%s_map.log",
        "suture1": "comp_suture1_%s_map.log",
        "bridge": "comp_suture3_%s_map.log",
        "app12": "comp_suture3_app12_%s_map.log",
    }
    columns = ["intact", "defect", "suture1", "bridge", "app12"]
    rows = ["pure", "flex"]
    fig, axes = plt.subplots(2, 5, figsize=(11.5, 5.5))

    for ri, scenario in enumerate(rows):
        for ci, model in enumerate(columns):
            draw_panel(
                axes[ri, ci],
                model,
                os.path.join(HERE, febs[model] % scenario),
                os.path.join(HERE, logs[model] % scenario),
            )
            if ri == 0:
                axes[ri, ci].set_title(labels[model], fontsize=8.5)
            if ci == 0:
                axes[ri, ci].text(
                    -0.11, 0.5,
                    "Pure axial\ncompression" if scenario == "pure" else "Flexion +\naxial compression",
                    transform=axes[ri, ci].transAxes,
                    rotation=90,
                    va="center",
                    ha="center",
                    fontsize=8.5,
                    fontweight="bold",
                )

    fig.suptitle(
        "Annular von Mises stress in an axial section through the defect at the 2000-N endpoint",
        fontsize=10.5,
        y=0.985,
    )
    sm = ScalarMappable(norm=NORM, cmap=CMAP)
    cax = fig.add_axes([0.25, 0.120, 0.50, 0.024])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=LEVELS, boundaries=LEVELS)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label(
        "element-averaged von Mises stress (MPa) | common discrete scale; highest range shown in red",
        fontsize=8.5,
    )
    fig.text(
        0.5, 0.012,
        "Disc interior is at the bottom and the posterolateral exterior at the top. Magenta = actual opening boundary; yellow = actual 3D discrete-spring path(s) projected into the axial plane; intact contains tissue at this site.",
        ha="center", fontsize=6.8, color="#222222",
    )
    fig.subplots_adjust(left=0.062, right=0.992, top=0.89, bottom=0.225, wspace=0.08, hspace=0.08)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"figR3_stress_maps.{ext}"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figR3_stress_maps.tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)
    print("wrote figR3_stress_maps (pdf/svg/png/tiff)")


if __name__ == "__main__":
    main()
