#!/usr/bin/env python3
"""Reviewer R1-10 / R2-1 / R2-7: evaluate the local mesh-sensitivity runs (h, h/2, h/3).

h   = original *_map.log (validation_L45_herniation/)
h/n = meshsens/<state>_<scen>_h<n>.log (from revision_mesh_refine.py)
Protrusion uses the same 8 original NP-defect interface vertices at every level (same physical points),
plus an 'all interface nodes' variant that includes the new refined interface nodes (anchors excluded).
Rim stress: (i) raw max over child elements; (ii) parent-averaged (children's stress tensors averaged
per original rim element, then von Mises) max / top-5 mean / 95th percentile -> same spatial support
as the original 57-element metric. IDP: -(sx+sy+sz)/3 averaged over the 8 nucleus core elements.
Output: REVISION_MESHSENS.csv
"""
import re, json, math, csv, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(HERE, "metric_sets.json")))
SUT = set(n for c in json.load(open(os.path.join(HERE, "suture_config.json"))).values() for p in c for n in p)
MOUTH8 = [n for n in S["mouth"] if n not in SUT]; RIM = S["rim"]; CORE = S["core"]
cx, cy = -9.85, -17.12

def init_xy():
    t = open(os.path.join(HERE, "comp_intact_pure.feb"), encoding="ISO-8859-1").read()
    ns = t[t.index("<Nodes"):t.index("</Nodes>")]
    P = {int(m.group(1)): [float(v) for v in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', ns)}
    m = np.mean([P[n][:2] for n in MOUTH8], 0); v = m - [cx, cy]; return v / np.linalg.norm(v)
NOUT = init_xy()

def vm(s):
    sx, sy, sz, sxy, syz, sxz = s
    return math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(sxy**2+syz**2+sxz**2))

def parse(path):
    rec = {}
    for blk in re.split(r'Data Record #\d+', open(path, encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if not (mt and md): continue
        k = {"x;y;z": "x", "sx;sy;sz;sxy;syz;sxz": "s"}.get(md.group(1))
        if not k: continue
        d = rec.setdefault(float(mt.group(1)), {"x": {}, "s": {}})
        for lm in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
            d[k][int(lm.group(1))] = np.array([float(v) for v in lm.group(2).split()])
    ts = sorted(t for t in rec if rec[t]["x"]); t0 = min(ts, key=lambda t: abs(t - 1))
    return [(t, rec[t]) for t in ts if t >= t0 - 1e-9]

def metrics(log, refine_json=None, master_idp=None):
    data = parse(log)
    ch = {}; iface = []
    if refine_json:
        meta = json.load(open(refine_json)); ch = {int(k): v for k, v in meta["children"].items()}; iface = meta["iface_new_nodes"]
    allif = MOUTH8 + iface
    ref = data[0][1]["x"]; F = []; p8 = []; pall = []
    for t, d in data:
        F.append(max(0.0, 2000*(t-1)/9))
        pr = lambda ns: np.mean([(d["x"][n]-ref[n])[:2] @ NOUT for n in ns if n in d["x"]])
        p8.append(pr(MOUTH8)); pall.append(pr(allif))
    F = np.array(F); p8 = np.array(p8); pall = np.array(pall)
    out = dict(peak=p8.max(), F_peak=F[p8.argmax()], posAUC=np.trapz(np.maximum(p8, 0), F), endpoint=p8[-1],
               peak_allif=pall.max(), posAUC_allif=np.trapz(np.maximum(pall, 0), F), n_iface=len(allif))
    s = data[-1][1]["s"]
    raw, pav = [], []
    for e in RIM:
        kids = ch.get(e, [e])
        ks = [s[k] for k in kids if k in s]
        raw += [vm(x) for x in ks]; pav.append(vm(np.mean(ks, 0)))
    pav = np.array(pav)
    out.update(rim_raw_max=max(raw), rim_raw_p95=np.percentile(raw, 95), rim_pav_max=pav.max(),
               rim_pav_top5=np.sort(pav)[-5:].mean(), rim_pav_p95=np.percentile(pav, 95), rim_pav_mean=pav.mean())
    core = [s[e] for e in CORE if e in s]
    out["IDP"] = -np.mean([(x[0]+x[1]+x[2])/3 for x in core]) if core else master_idp
    return out

if __name__ == "__main__":
    master = {(r["loading"], r["state"]): r for r in csv.DictReader(open(os.path.join(HERE, "..", "RESULTS_MASTER_TABLE.csv")))}
    MKEY = {"intact": "intact", "defect": "defect", "suture3": "bridge", "app12": "app12"}
    BASE = {"intact": "comp_intact_%s_map.log", "defect": "comp_defect_%s_map.log",
            "suture3": "comp_suture3_%s_map.log", "app12": "comp_suture3_app12_%s_map.log"}
    rows = []
    for st in BASE:
        for scen in ("pure", "flex"):
            m = master.get((scen, MKEY[st]))
            idp = float(m["IDP_endpoint_MPa"]) if m else None
            for lev in ("h", "h2", "h3"):
                if lev == "h":
                    log, rj = os.path.join(HERE, BASE[st] % scen), None
                else:
                    log = os.path.join(HERE, "meshsens", f"{st}_{scen}_{lev}.log"); rj = log.replace(".log", "_refine.json")
                    if not (os.path.exists(log) and "N O R M A L" in open(log, errors="ignore").read()[-3000:]):
                        continue
                r = dict(state=st, scen=scen, level=lev); r.update(metrics(log, rj, idp)); rows.append(r)
    # restoration per level
    for r in rows:
        g = {x["state"]: x for x in rows if x["scen"] == r["scen"] and x["level"] == r["level"]}
        if "intact" in g and "defect" in g and r["state"] not in ("intact", "defect"):
            for k in ("peak", "posAUC", "posAUC_allif"):
                r[f"rest_{k}"] = 100*(g["defect"][k]-r[k])/(g["defect"][k]-g["intact"][k])
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(HERE, "REVISION_MESHSENS.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    for r in rows:
        print(f"{r['state']:8s} {r['scen']:4s} {r['level']:3s} peak {r['peak']:.3f} posAUC {r['posAUC']:6.1f} allif {r['posAUC_allif']:6.1f} "
              f"end {r['endpoint']:+.3f} | rim raw {r['rim_raw_max']:.2f} pav {r['rim_pav_max']:.2f} top5 {r['rim_pav_top5']:.2f} p95 {r['rim_pav_p95']:.2f} | IDP {r['IDP']}"
              + (f" | rest pk {r['rest_peak']:.0f}% AUC {r['rest_posAUC']:.0f}%" if 'rest_peak' in r else ""))
