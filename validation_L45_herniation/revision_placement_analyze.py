#!/usr/bin/env python3
"""Independent suture placement (R1-3, R2-6): one 50 N/mm suture at each of the 20 candidate anchor pairs, plus
geometry-only outer-surface layouts (1, 2, 3 sutures).  Restoration relative to the baseline (1.0 MPa) intact and
defect.  Pair descriptors: position on the slot (axial z, radial depth; + = outer) and load-induced opening in the
unrepaired defect (0 -> 2000 N).  Output: REVISION_PLACEMENT.csv"""
import os, json, csv, numpy as np
from revision_sens_analyze import metrics, parse
from revision_morph_diagnostics import load
HERE = os.path.dirname(os.path.abspath(__file__)); SD = os.path.join(HERE, "sens")
sp = json.load(open(os.path.join(HERE, "suture_pairs.json"))); P = sp["pairs"]; C = np.array(sp["C"]); rd = np.r_[sp["rad"], 0]
S3 = {tuple(sorted(p)) for p in json.load(open(os.path.join(HERE, "suture_config.json")))["suture3"]}
REF = {"pure": dict(intact=(0.107, 154.7), defect=(0.358, 515.6), suture1=(0.281, 398.4), suture3=(0.282, 405.7)),
       "flex": dict(intact=(0.183, 211.9), defect=(0.474, 556.5), suture1=(0.301, 450.0), suture3=(0.304, 456.6))}
N, *_ = load(os.path.join(HERE, "comp_defect_pure_map.feb"))
opening = {}
for scen in ("pure", "flex"):
    d = parse(os.path.join(HERE, f"comp_defect_{scen}_map.log")); d = [x for x in d if x[0] >= 1 - 1e-9]
    X0, Xe = d[0][1]["x"], d[-1][1]["x"]
    opening[scen] = [float(np.linalg.norm(Xe[a]-Xe[b]) - np.linalg.norm(X0[a]-X0[b])) for a, b in P]
rows = []
for scen in ("pure", "flex"):
    (pi, ai), (pd, ad) = REF[scen]["intact"], REF[scen]["defect"]
    rest = lambda pk, auc: (100*(pd-pk)/(pd-pi), 100*(ad-auc)/(ad-ai))
    for i, (a, b) in enumerate(P):
        log = os.path.join(SD, f"single_p{i:02d}_{scen}.log")
        if not os.path.exists(log) or "N O R M A L" not in open(log, errors="ignore").read()[-5000:]: continue
        m = metrics(log); rp, ra = rest(m["peak"], m["posAUC"]); mid = (N[a]+N[b])/2 - C
        rows.append(dict(config=f"single_p{i:02d}", scen=scen, n_sutures=1, pairs=f"{a}-{b}", z_mm=round(mid[2], 2),
                         radial_mm=round(float(mid @ rd), 2), opening_mm=round(opening[scen][i], 3), in_S3=tuple(sorted((a, b))) in S3,
                         peak=round(m["peak"], 4), posAUC=round(m["posAUC"], 1), rest_peak=round(rp, 1), rest_AUC=round(ra, 1),
                         rim_max=round(m["rim_max"], 2), rim_p95=round(m["rim_p95"], 2)))
    lay = json.load(open(os.path.join(SD, "sens_manifest.json")))["outer_layouts"]
    for tag, prs in lay.items():
        log = os.path.join(SD, f"{tag}_{scen}.log")
        if not os.path.exists(log) or "N O R M A L" not in open(log, errors="ignore").read()[-5000:]: continue
        m = metrics(log); rp, ra = rest(m["peak"], m["posAUC"])
        idx = [P.index(p) for p in prs]
        rows.append(dict(config=tag, scen=scen, n_sutures=len(prs), pairs=";".join(f"{a}-{b}" for a, b in prs), z_mm="", radial_mm="",
                         opening_mm=round(sum(opening[scen][i] for i in idx), 3), in_S3=False, peak=round(m["peak"], 4),
                         posAUC=round(m["posAUC"], 1), rest_peak=round(rp, 1), rest_AUC=round(ra, 1),
                         rim_max=round(m["rim_max"], 2), rim_p95=round(m["rim_p95"], 2)))
with open(os.path.join(HERE, "REVISION_PLACEMENT.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); [w.writerow(r) for r in rows]
for scen in ("pure", "flex"):
    R = [r for r in rows if r["scen"] == scen]
    print(f"== {scen}  (reference: suture1 targeted rest AUC {100*(REF[scen]['defect'][1]-REF[scen]['suture1'][1])/(REF[scen]['defect'][1]-REF[scen]['intact'][1]):.0f}%, suture3 {100*(REF[scen]['defect'][1]-REF[scen]['suture3'][1])/(REF[scen]['defect'][1]-REF[scen]['intact'][1]):.0f}%)")
    for r in sorted([r for r in R if r["n_sutures"] == 1 and r["config"].startswith("single")], key=lambda r: -r["opening_mm"]):
        print(f"  {r['config']} z {r['z_mm']:+.2f} rad {r['radial_mm']:+.2f} opening {r['opening_mm']:+.3f} {'S3' if r['in_S3'] else '  '} -> rest AUC {r['rest_AUC']:5.1f}% peak {r['rest_peak']:5.1f}% rim {r['rim_max']:.2f}")
    sing = [r for r in R if r["config"].startswith("single")]
    if len(sing) > 3:
        o = np.array([r["opening_mm"] for r in sing]); ra = np.array([r["rest_AUC"] for r in sing])
        print(f"  Spearman-like rank corr(opening, rest AUC) = {np.corrcoef(np.argsort(np.argsort(o)), np.argsort(np.argsort(ra)))[0,1]:.2f}; n={len(sing)}")
    for r in [r for r in R if r["config"].startswith("outer")]:
        print(f"  {r['config']} ({r['n_sutures']} outer sutures) -> rest AUC {r['rest_AUC']:.1f}% peak {r['rest_peak']:.1f}% rim {r['rim_max']:.2f}")
