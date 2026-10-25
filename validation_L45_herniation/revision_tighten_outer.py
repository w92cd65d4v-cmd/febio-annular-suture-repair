#!/usr/bin/env python3
"""Tightened OUTER-surface sutures (follow-up to R1-3 / R2-6, 2026-10-08).
Same tighten-knot-release + slot-wall-contact model as revision_tighten_driver.py, but with the
geometry-only outer-half layouts (revision_gen_sens.OUTER): outer1 = pair 11; outer3 = pairs 0, 11, 16.
Each suture is tightened so that its anchor-pair gap after release matches the app12 MEAN residual gap
(1.186 mm), i.e. the same degree of apposition as the load-targeted app12 tightening.

  python3 revision_tighten_outer.py calibrate LAYOUT          (pure compression; iterates closures)
  python3 revision_tighten_outer.py full LAYOUT SCEN c1 [c2 c3]
Runs go to tighten_outer/.
"""
import sys, os, json, subprocess, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, "tighten_outer"); os.makedirs(T, exist_ok=True)
FEBIO = "/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4"
ALL = json.load(open(os.path.join(HERE, "suture_pairs.json")))["pairs"]
LAYOUT = {"outer1": [11], "outer3": [0, 11, 16]}
TARGET = 1.186
ENV = dict(os.environ, OMP_NUM_THREADS="2")

def febio(stem):
    for ext in (".log", ".xplt"):
        if os.path.exists(stem + ext): os.remove(stem + ext)
    subprocess.run([FEBIO, "-i", stem + ".feb", "-silent"], cwd=T, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV)
    return "N O R M A L" in open(stem + ".log", errors="ignore").read()[-5000:]

def gaps(stem, t, pairs):
    best = {}
    for blk in re.split(r'Data Record #\d+', open(stem + ".log", encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if mt and md and md.group(1) == "x;y;z" and abs(float(mt.group(1)) - t) < 1e-6:
            for m in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
                best[int(m.group(1))] = np.array([float(v) for v in m.group(2).split()])
    return [float(np.linalg.norm(best[a] - best[b])) for a, b in pairs]

def run(layout, scen, tag, cl, full):
    idx = LAYOUT[layout]; pairs = [ALL[i] for i in idx]
    gen = os.path.join(HERE, "revision_gen_tighten.py"); cs = [f"{c:.4f}" for c in cl]; pf = ["--pairs", ",".join(map(str, idx))]
    a = os.path.join(T, f"{tag}_A"); b = os.path.join(T, tag)
    subprocess.run(["python3", gen, scen, a + ".feb", *cs, "--passA", *pf], cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    if not febio(a): return None, "passA failed"
    subprocess.run(["python3", os.path.join(HERE, "revision_tighten_react.py"), a + ".log", a + "_react.json"], cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["python3", gen, scen, b + ".feb", *cs, "--react", a + "_react.json", *pf] + ([] if full else ["--short"]), cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    ok = febio(b)
    if not ok and not full: return None, "passB failed"
    try: g = gaps(b, 2.0, pairs)
    except KeyError: g = None
    return g, ("ok" if ok else "passB incomplete")

if __name__ == "__main__":
    mode, layout = sys.argv[1], sys.argv[2]
    pairs = [ALL[i] for i in LAYOUT[layout]]
    if mode == "calibrate":
        from revision_mesh_refine import parse
        N, _ = parse(open(os.path.join(HERE, "comp_suture3_pure_map.feb"), encoding="ISO-8859-1").read())
        g_open = np.array([np.linalg.norm(N[a] - N[b]) for a, b in pairs]); tgt = np.full(len(pairs), TARGET)
        cl = np.clip(g_open - tgt, 0, None); hist = []
        for it in range(6):
            g, st = run(layout, "pure", f"cal_{layout}_pure_it{it}", cl, False)
            print(f"it{it} closures {np.round(cl,3).tolist()} -> gaps(t=2) {None if g is None else np.round(g,3).tolist()} [{st}]", flush=True)
            if g is None: break
            g = np.array(g); hist.append((cl.copy(), g))
            if np.max(np.abs(g - tgt)) < 0.03: break
            cl = np.clip(cl + (g - tgt), 0, None)
        json.dump(dict(open_gaps=g_open.tolist(), target=TARGET, history=[(c.tolist(), g.tolist()) for c, g in hist]),
                  open(os.path.join(T, f"cal_{layout}_pure.json"), "w"))
        if hist: print("FINAL", layout, "closures", np.round(hist[-1][0], 4).tolist(), "gaps", np.round(hist[-1][1], 3).tolist())
    elif mode == "full":
        scen = sys.argv[3]; cl = [float(x) for x in sys.argv[4:]]
        g, st = run(layout, scen, f"tight_{layout}_{scen}", cl, True); print(f"tight_{layout}_{scen}", "gaps(t=2)", g, st, flush=True)
