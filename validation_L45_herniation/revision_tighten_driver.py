#!/usr/bin/env python3
"""Drive the tighten-knot-release model: pass A (tightening only, reactions) -> pass B (force-controlled
release [+ compression]).  Calibration mode iterates the per-suture closure so that the anchor-pair gaps
after release (t=2) match a morph level's per-pair gaps.

  python3 revision_tighten_driver.py calibrate SCEN LEVEL      (LEVEL in tighten/targets.json, e.g. app12)
  python3 revision_tighten_driver.py full SCEN TAG c1 c2 c3     (full model to 2000 N with given closures)
"""
import sys, os, json, subprocess, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, "tighten")
FEBIO = "/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4"
PAIRS = json.load(open(os.path.join(HERE, "suture_config.json")))["suture3"]
ENV = dict(os.environ, OMP_NUM_THREADS="2")

def febio(stem):
    for ext in (".log", ".xplt"):
        if os.path.exists(stem + ext): os.remove(stem + ext)
    subprocess.run([FEBIO, "-i", stem + ".feb", "-silent"], cwd=T, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV)
    txt = open(stem + ".log", errors="ignore").read()
    return "N O R M A L" in txt[-5000:]

def gaps(stem, t):
    best = {}
    for blk in re.split(r'Data Record #\d+', open(stem + ".log", encoding="ISO-8859-1", errors="ignore").read()):
        mt = re.search(r'Time\s*=\s*([\d.eE+-]+)', blk); md = re.search(r'Data\s*=\s*([a-zA-Z;]+)', blk)
        if mt and md and md.group(1) == "x;y;z" and abs(float(mt.group(1)) - t) < 1e-6:
            for m in re.finditer(r'^\s*(\d+)\s+([-\d.eE+ ]+)$', blk, re.M):
                best[int(m.group(1))] = np.array([float(v) for v in m.group(2).split()])
    return [float(np.linalg.norm(best[a] - best[b])) for a, b in PAIRS]

def run(scen, tag, cl, full):
    gen = os.path.join(HERE, "revision_gen_tighten.py"); cs = [f"{c:.4f}" for c in cl]
    a = os.path.join(T, f"{tag}_A"); b = os.path.join(T, f"{tag}")
    subprocess.run(["python3", gen, scen, a + ".feb", *cs, "--passA"], cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    if not febio(a): return None, "passA failed"
    subprocess.run(["python3", os.path.join(HERE, "revision_tighten_react.py"), a + ".log", a + "_react.json"], cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    args = ["python3", gen, scen, b + ".feb", *cs, "--react", a + "_react.json"] + ([] if full else ["--short"])
    subprocess.run(args, cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    ok = febio(b)
    if not ok and not full: return None, "passB failed"
    return gaps(b, 2.0), ("ok" if ok else "passB incomplete")

if __name__ == "__main__":
    mode, scen = sys.argv[1], sys.argv[2]
    if mode == "calibrate":
        level = sys.argv[3]; tgt = np.array(json.load(open(os.path.join(T, "targets.json")))[level])
        g_open = np.array([2.538, 2.562, 3.267])
        cl = np.clip(g_open - tgt, 0, None)            # first guess: closure = gap reduction
        if len(sys.argv) > 6: cl = np.array([float(x) for x in sys.argv[4:7]])   # optional starting closures
        hist = []
        for it in range(6):
            g, st = run(scen, f"cal_{level}_{scen}_it{it}", cl, False)
            print(f"it{it} closures {np.round(cl,3)} -> gaps(t=2) {None if g is None else np.round(g,3)} [{st}]", flush=True)
            if g is None: break
            g = np.array(g); hist.append((cl.copy(), g))
            if np.max(np.abs(g - tgt)) < 0.03: break
            cl = np.clip(cl + (g - tgt) * 1.0, 0, None)       # spring-back correction
        json.dump([(c.tolist(), g.tolist()) for c, g in hist], open(os.path.join(T, f"cal_{level}_{scen}.json"), "w"))
        if hist: print("FINAL", level, scen, "closures", np.round(hist[-1][0], 4).tolist(), "gaps", np.round(hist[-1][1], 3).tolist(), "target", tgt.tolist())
    elif mode == "full":
        tag = sys.argv[3]; cl = [float(x) for x in sys.argv[4:7]]
        g, st = run(scen, tag, cl, True); print(tag, scen, "gaps(t=2)", g, st)
