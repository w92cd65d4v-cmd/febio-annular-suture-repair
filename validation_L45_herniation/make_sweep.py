#!/usr/bin/env python3
"""RECONSTRUCTED suture-stiffness sweep generator (issue #2: original make_sweep.py lost).

Reverse-engineered from the existing comp_s{1,3}_E{10,200}_{scen}.feb files: each is the
corresponding base suture model (comp_suture{1,3}_{scen}.feb, E=50 N/mm) with ONLY the
suture load curve (id=10) endpoint rescaled to the swept stiffness. Verified to reproduce
those files (see verify_repro.py).

Sweep: E in {10, 200} N/mm (the 50 N/mm case IS the base comp_suture* model).
Mapping cfg->tag: suture1->s1, suture3->s3.

Usage:  python3 make_sweep.py [output_dir]      (default '.')
"""
import re, os, sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
SCENS = ["pure", "flex"]
E_SWEEP = [10, 200]                 # 50 N/mm is the base comp_suture* model
TAG = {"suture1": "s1", "suture3": "s3"}

def reslope(txt, E):
    """Rewrite the suture load curve (id=10) tension endpoint to (100, E*100)."""
    f = E * 100.0
    new, n = re.subn(r'(<loadcurve id="10" type="linear">\s*<point>-100,0</point>\s*'
                     r'<point>0,0</point>\s*<point>100,)[0-9.eE+]+(</point>)',
                     lambda m: m.group(1) + f'{f:g}' + m.group(2), txt, count=1, flags=re.S)
    assert n == 1, "suture loadcurve id=10 endpoint not found"
    return new

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for scen in SCENS:
        for cfg, tag in TAG.items():
            base = open(f"comp_{cfg}_{scen}.feb", encoding="ISO-8859-1").read()
            for E in E_SWEEP:
                out = reslope(base, E)
                path = os.path.join(OUT, f"comp_{tag}_E{E}_{scen}.feb")
                open(path, "w", encoding="ISO-8859-1").write(out)
                print(f"  wrote {path}  (E={E} N/mm)")
    print("done")
