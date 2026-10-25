#!/usr/bin/env python3
"""Reviewer R2-5: force-based apposition ("tighten - knot - release" + wall contact) instead of the
stress-free morph.

Start: open-gap three-suture model (comp_suture3_<scen>_map.feb; unmorphed defect, three load-targeted
anchor pairs).  Three analysis steps (time is cumulative):
  Step 1  t 0-1   preload ramps (as before) while each suture's tissue BITE (wall nodes within 1.0 mm of the
                  anchor, one patch on each slot wall) is pulled towards the opposite bite by a prescribed,
                  ramped closure c_k (half on each side, along the anchor-pair direction).
  Step 2  t 1-2   prescribed closure released ("knot"): every bite-node pair of suture k is joined by a
                  tension-only spring whose rest length is its tightened length (shifted curve), total
                  stiffness 50 N/mm per suture.  Tissue springs back against the knotted sutures ->
                  residual stress; no axial load.
  Step 3  t 2-11  axial compression 0 -> 2000 N (flex model: as before).  Reference state = t = 2.
The two opposing slot walls (18 + 18 annulus faces) form a frictionless sliding-elastic contact pair.
Logging as revision_gen_sens.instrument.

Usage: python3 revision_gen_tighten.py SCEN OUT.feb c1 [c2 ...] [--short] [--passA] [--react R.json] [--pairs i,j,k]
(--short: steps 1-2 only; --pairs: anchor-pair indices into suture_pairs.json, default = the three load-targeted pairs)
"""
import sys, re, json, os, collections, itertools, numpy as np
from revision_mesh_refine import parse
from revision_gen_sens import instrument
HERE = os.path.dirname(os.path.abspath(__file__))
K_SUTURE, BITE_R, KNOT_EPS = 50.0, 1.0, 0.1
FACES = [(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(0,3,2,1),(4,5,6,7)]

def slot_walls():
    Ni, Bi = parse(open(os.path.join(HERE, "comp_intact_pure_map.feb"), encoding="ISO-8859-1").read())
    Ei = {e: c for b in Bi for e, c in b["els"]}
    D = set(json.load(open(os.path.join(HERE, "region_sets.json")))["defect"])
    tg = np.r_[json.load(open(os.path.join(HERE, "suture_pairs.json")))["tang"], 0]
    fm = collections.defaultdict(list)
    for e, c in Ei.items():
        for f in FACES:
            q = tuple(c[i] for i in f); fm[frozenset(q)].append((e, q))
    A, B = [], []
    for lst in fm.values():
        if len(lst) == 2 and ((lst[0][0] in D) != (lst[1][0] in D)):
            e, q = [x for x in lst if x[0] not in D][0]
            X = np.array([Ni[n] for n in q]); nrm = np.cross(X[2]-X[0], X[3]-X[1]); nrm /= np.linalg.norm(nrm)
            if nrm @ tg > 0.6: A.append(q)
            elif nrm @ tg < -0.6: B.append(q)
    return A, B

def bites(N, pairs, A, B):
    An = {n for q in A for n in q}; Bn = {n for q in B for n in q}
    anc = [(a, b) if a in An else (b, a) for a, b in pairs]
    out = []
    for k, (a, b) in enumerate(anc):
        near = lambda n, side: min(range(len(anc)), key=lambda j: np.linalg.norm(N[n] - N[anc[j][side]]))
        ba = sorted(n for n in An if np.linalg.norm(N[n]-N[a]) <= BITE_R and near(n, 0) == k)
        bb = sorted(n for n in Bn if np.linalg.norm(N[n]-N[b]) <= BITE_R and near(n, 1) == k)
        u = (N[b] - N[a]) / np.linalg.norm(N[b] - N[a])
        out.append(dict(a=a, b=b, A=ba, B=bb, u=u))
    return out

def make(scen, closures, short=False, passA=False, reactions=None, pair_idx=None):
    txt = open(os.path.join(HERE, f"comp_suture3_{scen}_map.feb"), encoding="ISO-8859-1").read()
    N, _ = parse(txt)
    if pair_idx is None:   # default: the three load-targeted anchor pairs
        pairs = json.load(open(os.path.join(HERE, "suture_config.json")))["suture3"]
    else:                  # any candidate anchor pairs (indices into suture_pairs.json), e.g. outer-surface layouts
        allp = json.load(open(os.path.join(HERE, "suture_pairs.json")))["pairs"]; pairs = [allp[i] for i in pair_idx]
    assert len(pairs) == len(closures), (pairs, closures)
    A, B = slot_walls(); BT = bites(N, pairs, A, B)
    # --- remove the original single-spring sutures
    txt, n1 = re.subn(r'\t*<DiscreteSet name="SUTURE">.*?</DiscreteSet>\n', '', txt, flags=re.S)
    txt, n2 = re.subn(r'\t*<discrete_material id="8" name="SUTURE".*?</discrete_material>\n', '', txt, flags=re.S)
    txt, n3 = re.subn(r'\t*<discrete dmat="8" discrete_set="SUTURE"/>\n', '', txt)
    txt, n4 = re.subn(r'\t*<loadcurve id="10" type="linear">.*?</loadcurve>\n', '', txt, flags=re.S)
    assert (n1, n2, n3, n4) == (1, 1, 1, 1)
    # --- bite springs (fan: every A-bite node to every B-bite node), shifted to the tightened length
    sets, mats, binds, lcs, nsets = [], [], [], [], []
    sid, info = 0, []
    for k, (bt, c) in enumerate(zip(BT, closures)):
        springs = list(itertools.product(bt["A"], bt["B"])); kk = K_SUTURE / len(springs)
        for na, nb in springs:
            sid += 1
            L0 = np.linalg.norm(N[nb] - N[na])
            Lt = np.linalg.norm((N[nb] - c/2*bt["u"]) - (N[na] + c/2*bt["u"]))
            dl = max(L0 - Lt, 0.0) + KNOT_EPS   # knot pretension: rest length KNOT_EPS shorter than the tightened length
            name, mid, lid = f"SUT{sid:02d}", 7 + sid, 9 + sid
            sets.append(f'\t\t<DiscreteSet name="{name}">\n\t\t\t<delem>{na},{nb}</delem>\n\t\t</DiscreteSet>\n')
            mats.append(f'\t\t<discrete_material id="{mid}" name="{name}" type="nonlinear spring">\n\t\t\t<force lc="{lid}">1</force>\n\t\t</discrete_material>\n')
            binds.append(f'\t\t\t<discrete dmat="{mid}" discrete_set="{name}"/>\n')
            lcs.append(f'\t\t<loadcurve id="{lid}" type="linear">\n\t\t\t<point>-100,0</point>\n\t\t\t<point>{-dl:.6g},0</point>\n'
                       f'\t\t\t<point>0,{kk*dl:.6g}</point>\n\t\t\t<point>100,{kk*(100+dl):.6g}</point>\n\t\t</loadcurve>\n')
        nsets.append(f'\t\t<NodeSet name="BITE_A{k+1}">\n' + "".join(f'\t\t\t<node id="{n}"/>\n' for n in bt["A"]) + '\t\t</NodeSet>\n')
        nsets.append(f'\t\t<NodeSet name="BITE_B{k+1}">\n' + "".join(f'\t\t\t<node id="{n}"/>\n' for n in bt["B"]) + '\t\t</NodeSet>\n')
        info.append(dict(anchor=(bt["a"], bt["b"]), biteA=bt["A"], biteB=bt["B"], closure=c, n_springs=len(springs)))
    surf = lambda name, F: f'\t\t<Surface name="{name}">\n' + "".join(
        f'\t\t\t<quad4 id="{i+1}">{",".join(map(str, q))}</quad4>\n' for i, q in enumerate(F)) + '\t\t</Surface>\n'
    pair = '\t\t<SurfacePair name="SLOT">\n\t\t\t<master surface="SLOT_A"/>\n\t\t\t<slave surface="SLOT_B"/>\n\t\t</SurfacePair>\n'
    txt, n = re.subn(r'(\t*<DiscreteSet name="ITL_L4S">)', lambda m: "".join(nsets) + surf("SLOT_A", A) + surf("SLOT_B", B) + pair + m.group(1), txt, count=1); assert n == 1
    txt, n = re.subn(r'(</Geometry>)', lambda m: "".join(sets) + "\t" + m.group(1), txt, count=1); assert n == 1
    txt, n = re.subn(r'(\n)(\t*<discrete dmat="1")', lambda m: "\n" + "".join(mats) + m.group(2), txt, count=1); assert n == 1
    txt, n = re.subn(r'(<discrete dmat="7" discrete_set="LF_L4S"/>\n)', lambda m: m.group(1) + "".join(binds), txt, count=1); assert n == 1
    con = ('\t\t<contact type="sliding-elastic" name="slotwall" surface_pair="SLOT">\n'
           '\t\t\t<laugon>0</laugon>\n\t\t\t<penalty>1</penalty>\n\t\t\t<auto_penalty>1</auto_penalty>\n'
           '\t\t\t<two_pass>1</two_pass>\n\t\t\t<symmetric_stiffness>0</symmetric_stiffness>\n'
           '\t\t\t<search_radius>3</search_radius>\n\t\t\t<fric_coeff>0</fric_coeff>\n\t\t</contact>\n')
    txt, n = re.subn(r'(\t*</Contact>)', lambda m: con + m.group(1), txt, count=1); assert n == 1
    # --- load curves: compression starts at t=2; prescribed closure ramp (lc 11)
    txt, n = re.subn(r'(<loadcurve id="2" type="linear">\s*<point>0,0</point>\s*)<point>1,0</point>(\s*)<point>10(\.0)?,1</point>',
                     r'\1<point>2,0</point>\2<point>11,1</point>', txt, count=1); assert n == 1
    LC_RAMP = 10 + sid
    txt, n = re.subn(r'(</LoadData>)', lambda m: "".join(lcs) + f'\t\t<loadcurve id="{LC_RAMP}" type="linear" extend="constant">\n\t\t\t<point>0,0</point>\n\t\t\t<point>1,1</point>\n\t\t</loadcurve>\n\t' + m.group(1), txt, count=1); assert n == 1
    # --- steps
    m = re.search(r'\t*<Control>.*?</Control>\n', txt, re.S); ctrl = m.group(0); txt = txt[:m.start()] + txt[m.end():]
    ctl = lambda ns: re.sub(r'<time_steps>\d+</time_steps>', f'<time_steps>{ns}</time_steps>', ctrl)
    pres = ""
    for k, (bt, c) in enumerate(zip(BT, closures)):
        for side, sgn in (("A", +1), ("B", -1)):
            for j, ax in enumerate("xyz"):
                pres += f'\t\t\t<prescribe bc="{ax}" node_set="BITE_{side}{k+1}"><scale lc="{LC_RAMP}">{sgn*c/2*bt["u"][j]:.6g}</scale><relative>0</relative></prescribe>\n'
    steps = f'\t<Step name="tighten">\n{ctl(10)}\t\t<Boundary>\n{pres}\t\t</Boundary>\n\t</Step>\n'
    if passA:
        bnodes = sorted({n for bt in BT for n in bt["A"] + bt["B"]})
        txt, n = re.subn(r'(\t*</logfile>)', lambda m: f'\t\t\t<node_data data="Rx;Ry;Rz" name="react">{",".join(map(str, bnodes))}</node_data>\n' + m.group(1), txt, count=1); assert n == 1
        txt, n = re.subn(r'(\t*<Output>)', lambda m: steps + m.group(1), txt, count=1); assert n == 1
        return instrument(txt), info
    # force-controlled release: the bite reaction forces at t=1 are applied as nodal loads ramped 1 -> 0 over t=1-2
    LC_REL = LC_RAMP + 1
    txt, n = re.subn(r'(</LoadData>)', lambda m: f'\t\t<loadcurve id="{LC_REL}" type="linear" extend="constant">\n\t\t\t<point>1,1</point>\n\t\t\t<point>2,0</point>\n\t\t</loadcurve>\n\t' + m.group(1), txt, count=1); assert n == 1
    nl, rsets = "", ""
    for nd, r in sorted(reactions.items()):
        rsets += f'\t\t<NodeSet name="RN{nd}">\n\t\t\t<node id="{nd}"/>\n\t\t</NodeSet>\n'
        for j, ax in enumerate("xyz"):
            nl += f'\t\t\t<nodal_load bc="{ax}" node_set="RN{nd}"><scale lc="{LC_REL}">{r[j]:.8g}</scale></nodal_load>\n'
    txt, n = re.subn(r'(\t*<Surface name="SLOT_A">)', lambda m: rsets + m.group(1), txt, count=1); assert n == 1
    ctl_rel = re.sub(r'<step_size>[\d.]+</step_size>', '<step_size>0.025</step_size>', ctl(40))
    ctl_rel = re.sub(r'<dtmax>[\d.]+</dtmax>', '<dtmax>0.025</dtmax>', ctl_rel)
    steps += f'\t<Step name="release">\n{ctl_rel}\t\t<Loads>\n{nl}\t\t</Loads>\n\t</Step>\n'
    if not short:
        steps += f'\t<Step name="compress">\n{ctl(90)}\t</Step>\n'
    txt, n = re.subn(r'(\t*<Output>)', lambda m: steps + m.group(1), txt, count=1); assert n == 1
    return instrument(txt), info

if __name__ == "__main__":
    scen, out = sys.argv[1], sys.argv[2]
    cl = []
    for x in sys.argv[3:]:
        if x.startswith("--"): break
        cl.append(float(x))
    short = "--short" in sys.argv; passA = "--passA" in sys.argv
    pair_idx = [int(i) for i in sys.argv[sys.argv.index("--pairs") + 1].split(",")] if "--pairs" in sys.argv else None
    reactions = None
    if "--react" in sys.argv:
        reactions = {int(k): v for k, v in json.load(open(sys.argv[sys.argv.index("--react") + 1])).items()}
    txt, info = make(scen, cl, short, passA, reactions, pair_idx)
    open(out, "w", encoding="ISO-8859-1").write(txt)
    json.dump(info, open(out.replace(".feb", "_tighten.json"), "w"), default=lambda o: list(o) if isinstance(o, tuple) else str(o))
    print(out, [(i["anchor"], len(i["biteA"]), len(i["biteB"]), i["closure"]) for i in info])
