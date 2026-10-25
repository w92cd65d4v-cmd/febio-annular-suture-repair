#!/usr/bin/env python3
"""Generate the SENSITIVITY-STUDY models (issue #9), reusing the validated pipeline.

(A) Defect-size sensitivity (vary annulotomy height about the 18-element baseline):
      defectS  = the 10 baseline elements nearest the defect centre in z (shorter slot)
      defectL  = baseline 18 + face-adjacent posterior elements (taller slot)
    Built by deleting the chosen element set from the instrumented comp_intact_{scen}.feb
    (validated: deleting the baseline 18 reproduces comp_defect). The logfile node/element
    ID lists are pruned of any deleted entities so the run stays valid; the protrusion/rim/
    core metrics are then measured over the surviving instrumented nodes/elements.

(B) Suture number/placement sensitivity (defect unchanged -> existing metric sets stay valid):
      suture2  = 2 sutures (the two strongest-opening pairs)
      suture3b = 3 sutures at an ALTERNATIVE placement (different opening pairs)
    Built by inserting the tension-only suture spring (base 50 N/mm) into comp_defect_{scen}.feb,
    exactly as make_suture.py.

Writes comp_{defectS,defectL,suture2,suture3b}_{pure,flex}.feb. Does NOT run FEBio
(see run_sensitivity.sh). Does NOT touch any existing file.
"""
import re, json, numpy as np

SCENS = ["pure", "flex"]
BASE_IDS = json.load(open("../validation_L45/defect_spec.json"))["ids"]
PAIRS = json.load(open("suture_pairs.json"))["pairs"]           # candidate opening node-pairs
CFG3 = json.load(open("suture_config.json"))["suture3"]

# ---------- geometry from the intact mesh (for defect-size sets) ----------
itxt = open("comp_intact_pure.feb", encoding="ISO-8859-1").read()
nsec = itxt[itxt.index("<Nodes"):itxt.index("</Nodes>")]
NODE = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
        for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', nsec)}
POST = {}   # posterior-lamella elem id -> set(nodes)
for m in re.finditer(r'<Elements type="hex8" mat="\d+" name="(L45_POST_LAM\d)">(.*?)</Elements>', itxt, re.S):
    for em in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', m.group(2)):
        POST[int(em.group(1))] = set(int(v) for v in em.group(2).split(','))
cen = np.mean([np.mean([NODE[n] for n in POST[e]], axis=0) for e in BASE_IDS], axis=0)

def zc(eid):  # element centroid z
    return np.mean([NODE[n][2] for n in POST[eid]], axis=0)

# small: 10 baseline elements closest to centre height (shorter slot)
small = set(sorted(BASE_IDS, key=lambda e: abs(zc(e) - cen[2]))[:10])
# large: baseline + face-adjacent posterior elements that EXTEND THE HEIGHT only
# (centroid above the top / below the bottom of the baseline z-range), so the slot
# grows taller while preserving the lateral (circumferential) instrumented rim.
base_nodes = set().union(*(POST[e] for e in BASE_IDS))
zlo = min(zc(e) for e in BASE_IDS); zhi = max(zc(e) for e in BASE_IDS)
large = set(BASE_IDS)
for e, ns in POST.items():
    if e not in large and len(ns & base_nodes) >= 4 and not (zlo - 1e-6 <= zc(e) <= zhi + 1e-6):
        large.add(e)
print(f"defect-size sets: small={len(small)} elems, baseline=18, large={len(large)} elems (height-extended)")

# ---------- cleaning (delete elem set, orphans, pressure faces, prune loggers) ----------
def clean(txt, del_ids):
    del_ids = set(del_ids)
    delnodes = set()
    for eid in del_ids:
        m = re.search(r'<elem id="%d">\s*([^<]+)</elem>' % eid, txt)
        delnodes.update(int(v) for v in m.group(1).split(','))
    for eid in del_ids:
        txt, n = re.subn(r'[ \t]*<elem id="%d">[^<]*</elem>\n' % eid, '', txt); assert n == 1, eid
    used = set(int(v) for blk in re.finditer(r'<Elements [^>]*>(.*?)</Elements>', txt, re.S)
               for e in re.finditer(r'<elem id="\d+">\s*([^<]+)</elem>', blk.group(1))
               for v in e.group(1).split(','))
    orphans = delnodes - used
    for nid in orphans:
        txt = re.sub(r'[ \t]*<node id="%d">[^<]*</node>\n' % nid, '', txt)
    # pressure faces with all 4 nodes deleted
    def strip(m):
        body = re.sub(r'[ \t]*<quad4 id="\d+">\s*([^<]+)</quad4>\n',
                      lambda fm: '' if all(int(v) in delnodes for v in fm.group(1).split(',')) else fm.group(0),
                      m.group(2))
        return m.group(1) + body + m.group(3)
    txt = re.sub(r'(<Surface name="PressureLoad0[12]">)(.*?)(</Surface>)', strip, txt, flags=re.S)
    # prune logfile ID lists of deleted nodes/elements
    surv_n = set(int(m.group(1)) for m in re.finditer(r'<node id="(\d+)">', txt))
    surv_e = set(int(m.group(1)) for m in re.finditer(r'<elem id="(\d+)">', txt))
    def prune(tag, keep):
        def rep(m):
            ids = [int(x) for x in m.group(2).split(',') if int(x) in keep]
            return m.group(1) + ','.join(map(str, ids)) + m.group(3)
        return rep
    txt = re.sub(r'(<node_data[^>]*>)([0-9,]+)(</node_data>)', prune('n', surv_n), txt)
    txt = re.sub(r'(<element_data[^>]*>)([0-9,]+)(</element_data>)', prune('e', surv_e), txt)
    return txt

# ---------- suture insertion (mirrors make_suture.py) ----------
def add_suture(defect_txt, pairs, E=50.0):
    txt = defect_txt
    ds = '\t\t\t<DiscreteSet name="SUTURE">\n' + ''.join(f'\t\t\t<delem>{a},{b}</delem>\n' for a, b in pairs) + '\t\t</DiscreteSet>\n'
    txt, n = re.subn(r'(<DiscreteSet name="LF_L4S">.*?</DiscreteSet>\n)', lambda m: m.group(1) + ds, txt, 1, flags=re.S); assert n == 1
    mat = '\t\t<discrete_material id="8" name="SUTURE" type="nonlinear spring">\n\t\t\t<force lc="10">1</force>\n\t\t</discrete_material>\n'
    txt, n = re.subn(r'(\n)(\t*<discrete dmat="1")', lambda m: '\n' + mat + m.group(2), txt, 1); assert n == 1
    txt, n = re.subn(r'(<discrete dmat="7" discrete_set="LF_L4S"/>\n)', lambda m: m.group(1) + '\t\t\t<discrete dmat="8" discrete_set="SUTURE"/>\n', txt, 1); assert n == 1
    lc = f'\t\t\t<loadcurve id="10" type="linear">\n\t\t\t<point>-100,0</point>\n\t\t\t<point>0,0</point>\n\t\t\t<point>100,{E*100:g}</point>\n\t\t</loadcurve>\n'
    txt, n = re.subn(r'(</LoadData>)', lambda m: lc + '\t' + m.group(1), txt, 1); assert n == 1
    return txt

# alternative 3-suture placement: 3 well-separated opening pairs distinct from CFG3
cfg3_set = {tuple(p) for p in CFG3}
alt = [p for p in PAIRS if tuple(p) not in cfg3_set]
suture3b = [alt[0], alt[len(alt)//2], alt[-1]]
suture2 = CFG3[:2]
print(f"suture2 pairs = {suture2}")
print(f"suture3b (alt placement) pairs = {suture3b}")

if __name__ == "__main__":
    # self-check: deleting the baseline 18 from comp_intact must drop exactly 18 elems
    chk = clean(itxt, BASE_IDS)
    assert len(re.findall(r'<elem id=', itxt)) - len(re.findall(r'<elem id=', chk)) == 18
    print("self-check OK: clean(baseline 18) removes exactly 18 elements")
    for scen in SCENS:
        intact = open(f"comp_intact_{scen}.feb", encoding="ISO-8859-1").read()
        defect = open(f"comp_defect_{scen}.feb", encoding="ISO-8859-1").read()
        for tag, ids in (("defectS", small), ("defectL", large)):
            open(f"comp_{tag}_{scen}.feb", "w", encoding="ISO-8859-1").write(clean(intact, ids))
            print(f"  wrote comp_{tag}_{scen}.feb  ({len(ids)} elems deleted)")
        for tag, pairs in (("suture2", suture2), ("suture3b", suture3b)):
            open(f"comp_{tag}_{scen}.feb", "w", encoding="ISO-8859-1").write(add_suture(defect, pairs))
            print(f"  wrote comp_{tag}_{scen}.feb  ({len(pairs)} sutures)")
    print("done")
