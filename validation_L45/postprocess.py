#!/usr/bin/env python3
"""Extract moment-rotation / range-of-motion from the 6 FEBio runs and compare
to representative in-vitro L1-L2 data.

Reads each <variant>.log, parsing the rigid-body data records for material 11
(the loaded superior endplate): the rotation quaternion (qx;qy;qz;qw) and the
reaction moment (Mx;My;Mz) at every converged time step. Rotation angle about
the loaded axis is recovered from the quaternion; the applied moment is taken
from the recorded reaction moment about that axis.
"""
import re, math, os, csv, sys

P = sys.argv[1] if len(sys.argv) > 1 else "L12"   # level prefix (L12, L45, ...)

AXIS = {  # loaded axis per variant -> component index (0=x,1=y,2=z)
    f"{P}_flexion": 0, f"{P}_extension": 0,
    f"{P}_latbend_pos": 1, f"{P}_latbend_neg": 1,
    f"{P}_axrot_pos": 2, f"{P}_axrot_neg": 2,
}
PLANE = {f"{P}_flexion": "Flexion-extension", f"{P}_extension": "Flexion-extension",
         f"{P}_latbend_pos": "Lateral bending", f"{P}_latbend_neg": "Lateral bending",
         f"{P}_axrot_pos": "Axial rotation", f"{P}_axrot_neg": "Axial rotation"}

def parse_log(path, axis):
    """Return list of (moment_Nm_signed, rotation_deg_signed) per time step."""
    txt = open(path, encoding="ISO-8859-1", errors="ignore").read()
    # split into Data Record blocks; capture Time, Data type, value line
    blocks = re.split(r'Data Record #\d+', txt)
    q_by_t, m_by_t = {}, {}
    for blk in blocks:
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk)
        md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        mv = re.search(r'^\s*11\s+([-\d.eE+ ]+)$', blk, re.M)
        if not (mt and md and mv):
            continue
        t = float(mt.group(1)); vals = [float(x) for x in mv.group(1).split()]
        if md.group(1).startswith("qx") and len(vals) >= 4:
            q_by_t[t] = vals[:4]
        elif md.group(1).startswith("Mx") and len(vals) >= 3:
            m_by_t[t] = vals[:3]
    out = []
    for t in sorted(set(q_by_t) & set(m_by_t)):
        qx, qy, qz, qw = q_by_t[t]
        qw = max(-1.0, min(1.0, qw))
        qvec = (qx, qy, qz)
        # swing-twist: signed rotation ABOUT the loaded axis (not total magnitude),
        # so coupled off-axis motion does not inflate the reported primary ROM
        rot_deg = math.degrees(2.0 * math.atan2(qvec[axis], qw))
        mom_Nm = m_by_t[t][axis] / 1000.0   # N.mm -> N.m, about loaded axis
        out.append((mom_Nm, rot_deg))
    return out

def final_rom(series):
    """ROM magnitude at peak |moment| (end of loading)."""
    if not series: return None, None
    mom, rot = max(series, key=lambda p: abs(p[0]))
    return abs(rot), mom

results = {}
for v in AXIS:
    p = v + ".log"
    if os.path.exists(p):
        results[v] = parse_log(p, AXIS[v])

print("="*64)
print(f"{'variant':18s} {'plane':18s} {'M(N.m)':>8s} {'ROM(deg)':>9s} {'steps':>6s}")
print("-"*64)
rom = {}
for v in AXIS:
    s = results.get(v, [])
    r, m = final_rom(s)
    rom[v] = r
    if r is not None:
        print(f"{v:18s} {PLANE[v]:18s} {m:8.3f} {r:9.3f} {len(s):6d}")
    else:
        print(f"{v:18s} {PLANE[v]:18s} {'--':>8s} {'(no data)':>9s}")

# combined per-plane ROM
fe = (rom.get(f'{P}_flexion') or 0) + (rom.get(f'{P}_extension') or 0)
lb_l, lb_r = rom.get(f'{P}_latbend_pos'), rom.get(f'{P}_latbend_neg')
ar_l, ar_r = rom.get(f'{P}_axrot_pos'), rom.get(f'{P}_axrot_neg')

# representative in-vitro L1-L2 ROM at ~7.5 N.m (confirm vs your target studies)
LIT = {  # (low, high) degrees
    "Flexion (flex only)": (2.5, 6.5), "Extension (ext only)": (1.5, 5.0),
    "Flex-ext (total)": (4.0, 11.0),
    "Lateral bending (per side)": (3.0, 6.5),
    "Axial rotation (per side)": (1.0, 3.0),
}
print("\n" + "="*64)
print(f"MODEL ({P}) vs representative in-vitro lumbar (7.5 N.m)")
print("-"*64)
def chk(name, val, key):
    lo, hi = LIT[key]
    flag = "in range" if (val is not None and lo <= val <= hi) else "CHECK"
    vs = f"{val:.2f}" if val is not None else "n/a"
    print(f"  {name:30s} model={vs:>6s}  lit=[{lo:.1f},{hi:.1f}]  {flag}")
chk("Flexion", rom.get(f'{P}_flexion'), "Flexion (flex only)")
chk("Extension", rom.get(f'{P}_extension'), "Extension (ext only)")
chk("Flex-ext total", fe if fe else None, "Flex-ext (total)")
chk("Lateral bending (avg/side)", (lb_l+lb_r)/2 if lb_l and lb_r else None, "Lateral bending (per side)")
chk("Axial rotation (avg/side)", (ar_l+ar_r)/2 if ar_l and ar_r else None, "Axial rotation (per side)")

# write curves CSV
with open("moment_rotation_curves.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["variant", "plane", "moment_Nm", "rotation_deg"])
    for v in AXIS:
        for mom, rot in results.get(v, []):
            w.writerow([v, PLANE[v], f"{mom:.5f}", f"{rot:.5f}"])
print("\nwrote moment_rotation_curves.csv")

# plot if matplotlib available
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    planes = [("Flexion-extension", f'{P}_flexion', f'{P}_extension'),
              ("Lateral bending", f'{P}_latbend_pos', f'{P}_latbend_neg'),
              ("Axial rotation", f'{P}_axrot_pos', f'{P}_axrot_neg')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (title, vp, vn) in zip(axes, planes):
        for vv, col in ((vp, '#534AB7'), (vn, '#1D9E75')):
            s = results.get(vv, [])
            if s:
                xs = [r for _, r in s]; ys = [m for m, _ in s]
                ax.plot(xs, ys, '-o', ms=3, color=col, label=vv.replace(f'{P}_', ''))
        ax.axhline(0, color='#888', lw=.6); ax.axvline(0, color='#888', lw=.6)
        ax.set_title(title); ax.set_xlabel("rotation (deg)"); ax.set_ylabel("moment (N·m)")
        ax.grid(alpha=.25); ax.legend(fontsize=8)
    label = "–".join(["L"+c for c in P[1:]]) if P[0] == "L" and P[1:].isdigit() else P
    fig.suptitle(f"Intact {label} moment–rotation (49F, FEBio 4.12)", y=1.02)
    fig.tight_layout()
    fig.savefig(f"moment_rotation_{P}.png", dpi=150, bbox_inches="tight")
    print(f"wrote moment_rotation_{P}.png")
except Exception as e:
    print("plot skipped:", e)
