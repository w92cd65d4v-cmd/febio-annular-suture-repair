#!/usr/bin/env python3
"""Mesh sensitivity over the common attainable load range (0-756 N, t<=4.4).
Refined models (h/2, h/3) develop a local instability near 750-900 N (NP bulging into the defect /
outer posterior lamella at the defect site) and cannot be driven to 2000 N, so all levels are compared
on the same truncated path. h = original *_map.log truncated at t=4.4.
Output: REVISION_MESHSENS_PARTIAL.csv
"""
import os, json, csv, numpy as np
import revision_meshsens_analyze as A
HERE = os.path.dirname(os.path.abspath(__file__)); P = os.path.join(HERE, "meshsens", "partial")
import sys
TEND = float(sys.argv[1]) if len(sys.argv) > 1 else 4.4
OUTCSV = f"REVISION_MESHSENS_PARTIAL_t{TEND:g}.csv"
BASE = {"intact": "comp_intact_%s_map.log", "defect": "comp_defect_%s_map.log",
        "suture3": "comp_suture3_%s_map.log", "app12": "comp_suture3_app12_%s_map.log"}

def run(log, rj=None):
    data = [(t, d) for t, d in A.parse(log) if t <= TEND + 1e-6]
    ch, iface = {}, []
    if rj:
        m = json.load(open(rj)); ch = {int(k): v for k, v in m["children"].items()}; iface = m["iface_new_nodes"]
    ref = data[0][1]["x"]; F, p8, pall = [], [], []
    for t, d in data:
        F.append(max(0.0, 2000*(t-1)/9))
        pr = lambda ns: np.mean([(d["x"][n]-ref[n])[:2] @ A.NOUT for n in ns if n in d["x"]])
        p8.append(pr(A.MOUTH8)); pall.append(pr(A.MOUTH8 + iface))
    F, p8, pall = map(np.array, (F, p8, pall))
    s = data[-1][1]["s"]; raw, pav = [], []
    for e in A.RIM:
        ks = [s[k] for k in ch.get(e, [e]) if k in s]
        raw += [A.vm(x) for x in ks]; pav.append(A.vm(np.mean(ks, 0)))
    pav = np.array(pav)
    core = [s[e] for e in A.CORE if e in s]
    return dict(t_end=data[-1][0], F_end=F[-1], prot_end=p8[-1], peak=p8.max(), F_peak=F[p8.argmax()],
                posAUC=np.trapz(np.maximum(p8, 0), F), prot_end_allif=pall[-1], posAUC_allif=np.trapz(np.maximum(pall, 0), F),
                rim_pav_max=pav.max(), rim_pav_top5=np.sort(pav)[-5:].mean(), rim_pav_p95=np.percentile(pav, 95),
                rim_raw_max=max(raw), IDP=-np.mean([(x[0]+x[1]+x[2])/3 for x in core]) if core else np.nan)

rows = []
for scen in ("pure", "flex"):
    for st in BASE:
        for lev in ("h", "h2", "h3"):
            if lev == "h":
                r = run(os.path.join(HERE, BASE[st] % scen))
            else:
                log = os.path.join(P, f"{st}_{scen}_{lev}.log")
                if not os.path.exists(log): continue
                r = run(log, log.replace(".log", "_refine.json"))
            r.update(state=st, scen=scen, level=lev); rows.append(r)
for r in rows:
    g = {x["state"]: x for x in rows if x["scen"] == r["scen"] and x["level"] == r["level"]}
    if {"intact", "defect"} <= set(g) and r["state"] in ("suture3", "app12"):
        for k in ("posAUC", "prot_end", "posAUC_allif"):
            d, i = g["defect"][k], g["intact"][k]
            r["rest_" + k] = 100*(d - r[k])/(d - i) if d != i else np.nan
keys = ["scen", "state", "level"] + [k for k in rows[0] if k not in ("scen", "state", "level")] + ["rest_posAUC", "rest_prot_end", "rest_posAUC_allif"]
with open(os.path.join(HERE, OUTCSV), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow({k: (round(float(r[k]), 4) if isinstance(r.get(k), (float, np.floating)) else r.get(k, "")) for k in keys})
for r in rows:
    inc = r["t_end"] < TEND - 1e-3
    print(f"{r['scen']:4s} {r['state']:7s} {r['level']:2s} F_end {r['F_end']:4.0f}{'*' if inc else ' '} prot_end {r['prot_end']:.4f} (allif {r['prot_end_allif']:.4f}) "
          f"posAUC {r['posAUC']:6.1f} | rim pav {r['rim_pav_max']:.3f} p95 {r['rim_pav_p95']:.3f} raw {r['rim_raw_max']:.3f} | IDP {r['IDP']:.3f}"
          + (f" | rest AUC {r['rest_posAUC']:.0f}% end {r['rest_prot_end']:.0f}% allifAUC {r['rest_posAUC_allif']:.0f}%" if "rest_posAUC" in r else ""))
