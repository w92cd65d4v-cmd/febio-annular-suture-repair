#!/usr/bin/env python3
"""RECONSTRUCTED suture-repair generator (issue #2: original make_suture.py was lost).

Reverse-engineered from the existing, validated comp_suture{1,3}_{scen}.feb files by
diffing them against the unrepaired comp_defect_{scen}.feb base. Verified to reproduce
the suture material/discrete-set/load-curve blocks of those files (see verify_repro.py).

Transformation applied to comp_defect_{scen}.feb:
  1. <DiscreteSet name="SUTURE"> with one <delem>n1,n2</delem> per suture pair
     (pairs from suture_config.json: 'suture1' = 1 pair, 'suture3' = 3 pairs);
  2. <discrete_material id="8" name="SUTURE" type="nonlinear spring"> driven by
     <force lc="10">1</force>  (tension-only spring; force from load curve 10);
  3. <discrete dmat="8" discrete_set="SUTURE"/> binding material to the set;
  4. <loadcurve id="10" type="linear"> with points (-100,0)(0,0)(100, E*100) giving a
     tension-only stiffness of E N/mm (compression carries zero force). Base E = 50 N/mm.

The suture pairs are anchored only on slot-edge node pairs that OPEN under load
(determined upstream and stored in suture_config.json / suture_pairs.json).

Usage:  python3 make_suture.py [output_dir]      (default output_dir = '.')
"""
import re, os, sys, json

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
CFG = json.load(open("suture_config.json"))     # {"suture1":[[a,b]], "suture3":[[..],[..],[..]]}
SCENS = ["pure", "flex"]
E_BASE = 50.0                                    # N/mm (default suture stiffness)

def suture_loadcurve(E):
    f = E * 100.0                                # point (100, E*100) -> slope = E N/mm
    return ('\t\t\t<loadcurve id="10" type="linear">\n'
            '\t\t\t<point>-100,0</point>\n'
            '\t\t\t<point>0,0</point>\n'
            f'\t\t\t<point>100,{f:g}</point>\n'
            '\t\t</loadcurve>\n')

def discrete_set(pairs):
    rows = ''.join(f'\t\t\t<delem>{a},{b}</delem>\n' for a, b in pairs)
    return f'\t\t\t<DiscreteSet name="SUTURE">\n{rows}\t\t</DiscreteSet>\n'

SUT_MAT = ('\t\t<discrete_material id="8" name="SUTURE" type="nonlinear spring">\n'
           '\t\t\t<force lc="10">1</force>\n'
           '\t\t</discrete_material>\n')
SUT_BIND = '\t\t\t<discrete dmat="8" discrete_set="SUTURE"/>\n'

def build(defect_txt, pairs, E):
    txt = defect_txt
    # 1. DiscreteSet after the LF_L4S set (last ligament set, before discrete_materials)
    txt, n = re.subn(r'(<DiscreteSet name="LF_L4S">.*?</DiscreteSet>\n)',
                     lambda m: m.group(1) + discrete_set(pairs), txt, count=1, flags=re.S)
    assert n == 1, "LF_L4S DiscreteSet anchor not found"
    # 2. suture material before the first <discrete dmat="1"
    txt, n = re.subn(r'(\n)(\t*<discrete dmat="1")', lambda m: '\n' + SUT_MAT + m.group(2), txt, count=1)
    assert n == 1, "discrete-material anchor not found"
    # 3. bind suture after the dmat=7 (LF) entry
    txt, n = re.subn(r'(<discrete dmat="7" discrete_set="LF_L4S"/>\n)',
                     lambda m: m.group(1) + SUT_BIND, txt, count=1)
    assert n == 1, "dmat=7 anchor not found"
    # 4. suture load curve before </LoadData>
    txt, n = re.subn(r'(</LoadData>)', lambda m: suture_loadcurve(E) + '\t' + m.group(1), txt, count=1)
    assert n == 1, "</LoadData> anchor not found"
    return txt

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for scen in SCENS:
        base = open(f"comp_defect_{scen}.feb", encoding="ISO-8859-1").read()
        for cfg, pairs in CFG.items():            # suture1, suture3
            out = build(base, pairs, E_BASE)
            path = os.path.join(OUT, f"comp_{cfg}_{scen}.feb")
            open(path, "w", encoding="ISO-8859-1").write(out)
            print(f"  wrote {path}  ({len(pairs)} suture(s), E={E_BASE:g} N/mm)")
    print("done")
