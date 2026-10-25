#!/usr/bin/env python3
"""Compare intact vs posterolateral-annulotomy L4-L5 ROM (7.5 N.m, 6 load cases).
Reads material-11 rigid-body data from both the intact runs (../validation_L45)
and the defect runs (here), reports ROM and % change, and overlays the
moment-rotation curves (intact dashed vs defect solid)."""
import re, math, os

LOADS = ["flexion", "extension", "latbend_pos", "latbend_neg", "axrot_pos", "axrot_neg"]
AXISOF = {"flexion": 0, "extension": 0, "latbend_pos": 1, "latbend_neg": 1,
          "axrot_pos": 2, "axrot_neg": 2}
PLANE = {"flexion": "Flexion-extension", "extension": "Flexion-extension",
         "latbend_pos": "Lateral bending", "latbend_neg": "Lateral bending",
         "axrot_pos": "Axial rotation", "axrot_neg": "Axial rotation"}

def series(path, axis):
    if not os.path.exists(path): return []
    txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    q, m = {}, {}
    for blk in re.split(r'Data Record #\d+', txt):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk)
        md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        mv = re.search(r'^\s*11\s+([-\d.eE+ ]+)$', blk, re.M)
        if not (mt and md and mv): continue
        t = float(mt.group(1)); vals = [float(x) for x in mv.group(1).split()]
        if md.group(1).startswith("qx"): q[t] = vals[:4]
        elif md.group(1).startswith("Mx"): m[t] = vals[:3]
    out = []
    for t in sorted(set(q) & set(m)):
        qx, qy, qz, qw = q[t]; qw = max(-1, min(1, qw))
        rot = math.degrees(2*math.atan2((qx, qy, qz)[axis], qw))
        out.append((m[t][axis]/1000.0, rot))
    return out

def rom(s): return abs(max(s, key=lambda p: abs(p[0]))[1]) if s else None

intact, defect = {}, {}
for ld in LOADS:
    ax = AXISOF[ld]
    intact[ld] = series(f"../validation_L45/L45_{ld}.log", ax)
    defect[ld] = series(f"L45d_{ld}.log", ax)

print("="*70)
print(f"{'load case':16s} {'plane':18s} {'intact':>8s} {'defect':>8s} {'Δ%':>7s}")
print("-"*70)
for ld in LOADS:
    ri, rd = rom(intact[ld]), rom(defect[ld])
    if ri and rd:
        print(f"{ld:16s} {PLANE[ld]:18s} {ri:8.3f} {rd:8.3f} {100*(rd-ri)/ri:+7.1f}")
    else:
        print(f"{ld:16s} {PLANE[ld]:18s} {'?':>8s} {'?':>8s}  (pending)")

# per-plane combined
def comb(d):
    return {"Flexion-extension": (rom(d['flexion']) or 0)+(rom(d['extension']) or 0),
            "Lateral bending": ((rom(d['latbend_pos']) or 0)+(rom(d['latbend_neg']) or 0))/2,
            "Axial rotation": ((rom(d['axrot_pos']) or 0)+(rom(d['axrot_neg']) or 0))/2}
ci, cd = comb(intact), comb(defect)
print("\n" + "="*70)
print(f"{'plane (combined)':22s} {'intact':>8s} {'defect':>8s} {'Δ%':>7s}")
print("-"*70)
for k in ci:
    if ci[k] and cd[k]:
        print(f"{k:22s} {ci[k]:8.3f} {cd[k]:8.3f} {100*(cd[k]-ci[k])/ci[k]:+7.1f}")

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    planes = [("Flexion-extension", 'flexion', 'extension'),
              ("Lateral bending", 'latbend_pos', 'latbend_neg'),
              ("Axial rotation", 'axrot_pos', 'axrot_neg')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (title, vp, vn) in zip(axes, planes):
        for v in (vp, vn):
            si, sd = intact[v], defect[v]
            if si: ax.plot([r for _, r in si], [m for m, _ in si], '--', color='#888', lw=1.5)
            if sd: ax.plot([r for _, r in sd], [m for m, _ in sd], '-', color='#A32D2D', lw=1.8)
        ax.axhline(0, color='#bbb', lw=.6); ax.axvline(0, color='#bbb', lw=.6)
        ax.set_title(title); ax.set_xlabel("rotation (deg)"); ax.set_ylabel("moment (N·m)"); ax.grid(alpha=.25)
    axes[0].plot([], [], '--', color='#888', label='intact'); axes[0].plot([], [], '-', color='#A32D2D', label='defect'); axes[0].legend(fontsize=9)
    fig.suptitle("L4–L5 intact vs posterolateral annulotomy (7.5 N·m)", y=1.02)
    fig.tight_layout(); fig.savefig("intact_vs_defect_L45.png", dpi=150, bbox_inches="tight")
    print("\nwrote intact_vs_defect_L45.png")
except Exception as e:
    print("plot skipped:", e)
