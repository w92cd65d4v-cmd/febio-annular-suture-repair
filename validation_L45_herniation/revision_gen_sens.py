#!/usr/bin/env python3
"""Revision R1 sensitivity models (coarse, original mesh) -> sens/

A. Preload sensitivity (R1-5, R2-2): expansive NP preload 0.1 and 0.5 MPa (baseline 1.0 MPa) for
   intact, defect, suture1, suture3 (open gap) and app12; pure and flex.  Intact also at 1.0 MPa
   (re-run only to add disc-deformation logging).
B. Independent suture placement (R1-3, R2-6): one 50 N/mm suture at EACH of the 20 candidate anchor
   pairs in turn, plus geometry-only outer-surface layouts with 1, 2, 3 sutures (chosen from position
   alone, without the opening data).

Every model logs: region nodes x;y;z (as *_map), rim/region element stresses + NP core elements (IDP),
and the L4 rigid body position (axial FSU displacement = disc deformation). Plot output = displacement only.
"""
import re, os, json
from make_suture import build

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "sens"); os.makedirs(OUT, exist_ok=True)
S = json.load(open(os.path.join(HERE, "metric_sets.json")))
PAIRS = json.load(open(os.path.join(HERE, "suture_pairs.json")))["pairs"]
BASE = {"intact": "comp_intact_%s_map.feb", "defect": "comp_defect_%s_map.feb", "suture1": "comp_suture1_%s_map.feb",
        "suture3": "comp_suture3_%s_map.feb", "app12": "comp_suture3_app12_%s_map.feb"}
# geometry-only outer-surface layouts (radial >= +0.8 mm, i.e. outer half of the slot; spread along the slot height)
OUTER = {"outer1": [11], "outer2": [6, 16], "outer3": [0, 11, 16]}

def instrument(txt, preload=None):
    m = re.search(r'<element_data data="sx;sy;sz;sxy;syz;sxz">([^<]+)</element_data>', txt)
    ids = [int(x) for x in m.group(1).split(",")]
    ids += [e for e in S["core"] + S["rim"] if e not in ids]
    txt = txt[:m.start(1)] + ",".join(map(str, ids)) + txt[m.end(1):]
    txt, n = re.subn(r'(<rigid_body_data data="Mx;My;Mz">11</rigid_body_data>)',
                     r'\1\n\t\t\t<rigid_body_data data="x;y;z">11</rigid_body_data>', txt, count=1)
    assert n == 1
    txt, n = re.subn(r'<plotfile type="febio">.*?</plotfile>',
                     '<plotfile type="febio">\n\t\t\t<var type="displacement"/>\n\t\t</plotfile>', txt, count=1, flags=re.S)
    assert n == 1
    if preload is not None:
        txt, n1 = re.subn(r'<pressure lc="1">-1</pressure>', f'<pressure lc="1">-{preload:g}</pressure>', txt)
        txt, n2 = re.subn(r'<pressure lc="1">1</pressure>', f'<pressure lc="1">{preload:g}</pressure>', txt)
        assert n1 == 1 and n2 == 1
    return txt

def write(name, txt):
    open(os.path.join(OUT, name + ".feb"), "w", encoding="ISO-8859-1").write(txt)

if __name__ == "__main__":
    names = []
    for scen in ("pure", "flex"):
        for st, f in BASE.items():
            base = open(os.path.join(HERE, f % scen), encoding="ISO-8859-1").read()
            for p in (0.1, 0.5):
                n = f"pl{int(p*100):03d}_{st}_{scen}"; write(n, instrument(base, p)); names.append(n)
            if st == "intact":
                n = f"pl100_{st}_{scen}"; write(n, instrument(base, 1.0)); names.append(n)
        defect = open(os.path.join(HERE, BASE["defect"] % scen), encoding="ISO-8859-1").read()
        for i, pr in enumerate(PAIRS):
            n = f"single_p{i:02d}_{scen}"; write(n, instrument(build(defect, [pr], 50.0))); names.append(n)
        for tag, idx in OUTER.items():
            n = f"{tag}_{scen}"; write(n, instrument(build(defect, [PAIRS[i] for i in idx], 50.0))); names.append(n)
    json.dump({"outer_layouts": {k: [PAIRS[i] for i in v] for k, v in OUTER.items()}, "models": names},
              open(os.path.join(OUT, "sens_manifest.json"), "w"), indent=1)
    print(len(names), "models written to", OUT)
