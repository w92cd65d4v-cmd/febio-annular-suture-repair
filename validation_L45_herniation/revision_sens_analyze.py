#!/usr/bin/env python3
"""Analyse revision sensitivity runs in sens/ (preload, independent placement) on the original mesh.
Metrics as in the manuscript (8 NP-defect interface nodes, positive-only AUC, 57-element rim VM) plus rim p95,
IDP (8 core elements) at F=0 and 2000 N, and L4 rigid-body axial displacement (FSU/disc deformation).
Output: REVISION_SENS_RESULTS.csv (+ printed summary)."""
import os, re, glob, json, csv, math, numpy as np
import revision_meshsens_analyze as A
HERE = os.path.dirname(os.path.abspath(__file__)); SD = os.path.join(HERE, "sens")

def parse(path):
    rec = {}
    for blk in re.split(r'Data Record #\d+', open(path, encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        rows = {int(m.group(1)): np.array([float(v) for v in m.group(2).split()]) for m in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M)}
        d = rec.setdefault(float(mt.group(1)), {"x": {}, "s": {}, "rb": None})
        if md.group(1) == "x;y;z":
            if list(rows) == [11]: d["rb"] = rows[11]
            else: d["x"].update(rows)
        elif md.group(1).startswith("sx;sy;sz"): d["s"].update(rows)
    ts = sorted(t for t in rec if rec[t]["x"])
    return [(t, rec[t]) for t in ts]

def metrics(path, t0=1.0):
    data = parse(path); data = [(t, d) for t, d in data if t >= t0 - 1e-9]
    ref = data[0][1]; F, p = [], []
    for t, d in data:
        F.append(max(0.0, 2000*(t - t0)/9))
        p.append(np.mean([(d["x"][n] - ref["x"][n])[:2] @ A.NOUT for n in A.MOUTH8]))
    F, p = np.array(F), np.array(p)
    end = data[-1][1]; s = end["s"]
    rim = np.array([A.vm(s[e]) for e in A.RIM if e in s])
    idp = lambda d: -np.mean([(d["s"][e][0]+d["s"][e][1]+d["s"][e][2])/3 for e in A.CORE if e in d["s"]])
    out = dict(F_end=F[-1], complete=F[-1] > 1999, peak=p.max(), F_peak=F[p.argmax()], posAUC=np.trapz(np.maximum(p, 0), F),
               signedAUC=np.trapz(p, F), endpoint=p[-1], rim_max=rim.max(), rim_p95=np.percentile(rim, 95),
               IDP_F0=idp(ref), IDP_end=idp(end))
    if ref["rb"] is not None and end["rb"] is not None:
        out["axial_disp_mm"] = float(end["rb"][2] - ref["rb"][2])
    return out

if __name__ == "__main__":
    rows = []
    for log in sorted(glob.glob(os.path.join(SD, "*.log"))):
        name = os.path.basename(log)[:-4]
        if "N O R M A L" not in open(log, errors="ignore").read()[-5000:]: continue
        r = dict(model=name); r.update(metrics(log)); rows.append(r)
    R = {r["model"]: r for r in rows}
    def rest(r, base, scen):
        d, i = R.get(f"{base}defect_{scen}"), R.get(f"{base}intact_{scen}")
        if not (d and i): return {}
        return {f"rest_{k}": 100*(d[k]-r[k])/(d[k]-i[k]) for k in ("peak", "posAUC")}
    for r in rows:
        m = re.match(r'(pl\d{3}_)(\w+?)_(pure|flex)$', r["model"])
        if m and m.group(2) not in ("intact", "defect"): r.update(rest(r, m.group(1), m.group(3)))
        m = re.match(r'(single_p\d\d|outer\d)_(pure|flex)$', r["model"])
        if m:  # placement runs: reference = baseline-preload (1.0 MPa) intact/defect from the original logs
            pass
    keys = ["model"] + [k for k in rows[0] if k != "model"] + ["rest_peak", "rest_posAUC"]
    with open(os.path.join(HERE, "REVISION_SENS_RESULTS.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({k: (round(float(r[k]), 4) if isinstance(r.get(k), (float, np.floating)) else r.get(k, "")) for k in keys})
    for r in rows:
        print(f"{r['model']:20s} peak {r['peak']:.3f} @{r['F_peak']:4.0f} posAUC {r['posAUC']:6.1f} end {r['endpoint']:+.3f} "
              f"rim {r['rim_max']:.2f}/p95 {r['rim_p95']:.2f} IDP {r['IDP_F0']:+.3f}->{r['IDP_end']:.3f} ax {r.get('axial_disp_mm', float('nan')):+.3f}"
              + (f" | rest pk {r['rest_peak']:.0f}% AUC {r['rest_posAUC']:.0f}%" if "rest_peak" in r else ""))
