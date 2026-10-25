#!/usr/bin/env python3
"""Build axial-compression herniation simulations (intact vs defect) for the L4-L5 model.

Replaces the pure-moment loading with an axial compressive force ramped 0 -> 2000 N on
the superior rigid endplate (inferior endplate fixed), in two scenarios:
  pure  : pure axial compression (endplate kept level: x,y,Rx,Ry,Rz fixed; z free, loaded)
  flex  : combined flexion (7.5 N.m) + axial compression (x,y,Ry,Rz fixed; z & Rx free, loaded)
Load curve 2 is made linear so applied force = 2000*(t-1)/9 for t in [1,10] (exact).
Instrumentation added to the logfile (explicit id lists):
  - nucleus 'mouth' node positions x;y;z  -> posterolateral protrusion through the defect
  - annulus rim element stress             -> stress concentration at the defect edge
  - nucleus core element stress            -> intradiscal pressure = -(sx+sy+sz)/3
Outputs 4 input files: comp_{intact,defect}_{pure,flex}.feb
"""
import re, json

S = json.load(open("../validation_L45_defect/metric_sets.json"))
mouth = ",".join(map(str, S["mouth"]))
rim   = ",".join(map(str, S["rim"]))
core  = ",".join(map(str, S["core"]))

BC = {
 "pure": ('\t\t<rigid_body mat="11">\n'
          '\t\t\t<fixed bc="x"/>\n\t\t\t<fixed bc="y"/>\n'
          '\t\t\t<fixed bc="Rx"/>\n\t\t\t<fixed bc="Ry"/>\n\t\t\t<fixed bc="Rz"/>\n'
          '\t\t\t<force bc="z" lc="2">-2000</force>\n'
          '\t\t</rigid_body>'),
 "flex": ('\t\t<rigid_body mat="11">\n'
          '\t\t\t<fixed bc="x"/>\n\t\t\t<fixed bc="y"/>\n'
          '\t\t\t<fixed bc="Ry"/>\n\t\t\t<fixed bc="Rz"/>\n'
          '\t\t\t<force bc="z" lc="2">-2000</force>\n'
          '\t\t\t<force bc="Rx" lc="2">7500</force>\n'
          '\t\t</rigid_body>'),
}
LOG_ADD = (f'\t\t\t<node_data data="x;y;z">{mouth}</node_data>\n'
           f'\t\t\t<element_data data="sx;sy;sz;sxy;syz;sxz">{rim}</element_data>\n'
           f'\t\t\t<element_data data="sx;sy;sz">{core}</element_data>\n')

BASES = {"intact": "../validation_L45/L45_flexion.feb",
         "defect": "../validation_L45_defect/L45d_flexion.feb"}
PH = "@@RIGID11@@"

for model, basef in BASES.items():
    txt = open(basef, encoding="ISO-8859-1").read()
    # 1. swap the superior-endplate rigid-body block for a placeholder
    txt, n = re.subn(r'\t*<rigid_body mat="11">.*?</rigid_body>', PH, txt, count=1, flags=re.S)
    assert n == 1, (model, n)
    # 2. linear load curve 2 (exact force scheduling)
    txt = txt.replace('<loadcurve id="2" type="smooth">', '<loadcurve id="2" type="linear">')
    # 3. add instrumentation before </logfile>
    txt = txt.replace('\t\t</logfile>', LOG_ADD + '\t\t</logfile>')
    assert LOG_ADD in txt
    for scen, block in BC.items():
        out = txt.replace(PH, block)
        assert PH not in out
        fn = f"comp_{model}_{scen}.feb"
        open(fn, "w", encoding="ISO-8859-1").write(out)
        print(f"  wrote {fn}")
print("done")
