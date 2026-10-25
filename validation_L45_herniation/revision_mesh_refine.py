#!/usr/bin/env python3
"""Reviewer R1-10 / R2-1 / R2-7: local mesh-sensitivity models.

Pure h-refinement of a fixed zone around the posterolateral defect: every hex in the zone is split
n x n x n by its own trilinear map, so the geometry (including any apposition morph) is unchanged.
Nodes left hanging on the fine/coarse transition are tied to their coarse edge/face by augmented-
Lagrangian linear constraints whose coefficients are the exact trilinear weights (watertight mesh).
Children inherit the parent's material block and node orientation, so the local mat_axis (fiber
frame) is reproduced. Pressure-preload faces in the zone are subdivided; springs/sutures, ligaments,
rigid bodies and loading are untouched.

Zone (fixed across all states, defined once on the intact model -> REFINE_ZONE.json):
  deformable elements with centroid <= 8 mm from the defect centre, plus the one-element ring around
  the 57 rim elements and the NP-defect interface (so no rim/interface node lies on the transition).

Usage: python3 revision_mesh_refine.py BASE.feb OUT.feb N
"""
import re, sys, json, os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
XI = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
RIGID = {"L4UVEND", "L5LVEND"}
PENALTY, TOL, MAXAUG = 10.0, 1e-2, 50

def parse(txt):
    ns = txt[txt.index("<Nodes"):txt.index("</Nodes>")]
    N = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(",")])
         for m in re.finditer(r'<node id="(\d+)">\s*([^<]+)</node>', ns)}
    blocks = []
    for b in re.finditer(r'<Elements type="hex8" mat="(\d+)" name="([^"]+)">(.*?)</Elements>', txt, re.S):
        els = [(int(e.group(1)), [int(v) for v in e.group(2).split(",")])
               for e in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', b.group(3))]
        blocks.append(dict(mat=b.group(1), name=b.group(2), els=els, span=b.span()))
    return N, blocks

def make_zone():
    txt = open(os.path.join(HERE, "comp_intact_pure_map.feb"), encoding="ISO-8859-1").read()
    N, blocks = parse(txt)
    S = json.load(open(os.path.join(HERE, "metric_sets.json")))
    C = np.array(json.load(open(os.path.join(HERE, "region_sets.json")))["C"])
    E = {e: (c, b["name"]) for b in blocks for e, c in b["els"]}
    zone = {e for e, (c, nm) in E.items() if nm not in RIGID and np.linalg.norm(np.mean([N[n] for n in c], 0) - C) <= 8.0}
    key = set(S["rim"]) | {e for e, (c, nm) in E.items() if set(c) & set(S["mouth"])}
    kn = {n for e in key for n in E[e][0]}
    zone |= {e for e, (c, nm) in E.items() if nm not in RIGID and set(c) & kn}
    json.dump(sorted(zone), open(os.path.join(HERE, "REFINE_ZONE.json"), "w"))
    return zone

def refine(base, out, n):
    zf = os.path.join(HERE, "REFINE_ZONE.json")
    ZONE = set(json.load(open(zf))) if os.path.exists(zf) else make_zone()
    txt = open(base, encoding="ISO-8859-1").read()
    N, blocks = parse(txt)
    E = {e: c for b in blocks for e, c in b["els"]}
    NAME = {e: b["name"] for b in blocks for e, _ in b["els"]}
    ref = {e for e in E if e in ZONE}
    unref = [e for e in E if e not in ref]
    # node -> unrefined elements (for hanging detection)
    n2u = {}
    for e in unref:
        for nd in E[e]:
            n2u.setdefault(nd, set()).add(e)
    nid = max(N) + 1
    key2id, newpos, support = {}, {}, {}
    def node_for(W):                      # W: {corner node: integer weight}, sum = n^3
        nonlocal nid
        k = tuple(sorted((a, w) for a, w in W.items() if w))
        if len(k) == 1:
            return k[0][0]
        if k not in key2id:
            key2id[k] = nid
            newpos[nid] = sum(w * N[a] for a, w in k) / n**3
            support[nid] = k
            nid += 1
        return key2id[k]
    children, next_e = {}, max(E) + 1
    for e in sorted(ref):
        c = E[e]
        grid = {}
        for i, j, k in itertools.product(range(n + 1), repeat=3):
            W = {}
            for a, (xa, ya, za) in enumerate(XI):
                W[c[a]] = W.get(c[a], 0) + (i if xa > 0 else n - i) * (j if ya > 0 else n - j) * (k if za > 0 else n - k)
            grid[i, j, k] = node_for(W)
        ch = []
        for i, j, k in itertools.product(range(n), repeat=3):
            conn = [grid[i + (xa > 0), j + (ya > 0), k + (za > 0)] for xa, ya, za in XI]
            ch.append((next_e, conn)); next_e += 1
        children[e] = ch
    # hanging nodes: support edge/face shared with an unrefined element
    hang = {}
    for h, k in support.items():
        nodes = [a for a, _ in k]
        if len(nodes) == 8:
            continue
        cand = set.intersection(*[n2u.get(a, set()) for a in nodes])
        if cand:
            if any(NAME[x] in RIGID for x in cand):
                raise RuntimeError(f"hanging node {h} on a rigid element")
            tot = sum(w for _, w in k)
            hang[h] = [(a, w / tot) for a, w in k]
    n_written = len(newpos)
    # ---- rebuild text ----
    out_txt = txt
    # nodes
    nblk = "".join(f'\t\t\t<node id="{i}">{p[0]:.9g},{p[1]:.9g},{p[2]:.9g}</node>\n' for i, p in sorted(newpos.items()))
    ie = out_txt.index("</Nodes>")
    out_txt = out_txt[:ie] + nblk + "\t\t" + out_txt[ie:]
    # elements (rebuild each affected block)
    def rebuild(m):
        mat, name, body = m.group(1), m.group(2), m.group(3)
        els = [(int(x.group(1)), [int(v) for v in x.group(2).split(",")])
               for x in re.finditer(r'<elem id="(\d+)">\s*([^<]+)</elem>', body)]
        if not any(e in ref for e, _ in els):
            return m.group(0)
        lines = []
        for e, c in els:
            for ce, cc in (children[e] if e in ref else [(e, c)]):
                lines.append(f'\t\t\t<elem id="{ce}">{",".join(map(str, cc))}</elem>\n')
        return f'<Elements type="hex8" mat="{mat}" name="{name}">\n' + "".join(lines) + "\t\t</Elements>"
    out_txt = re.sub(r'<Elements type="hex8" mat="(\d+)" name="([^"]+)">(.*?)</Elements>', rebuild, out_txt, flags=re.S)
    # surfaces
    refset_nodes = {e: set(E[e]) for e in ref}
    node2ref = {}
    for e in ref:
        for a in E[e]:
            node2ref.setdefault(a, set()).add(e)
    added_nodes = {}
    def sub_surface(m):
        name, body = m.group(1), m.group(2)
        faces = [[int(v) for v in f.split(",")] for f in re.findall(r'<quad4 id="\d+">\s*([^<]+)</quad4>', body)]
        outf, touched = [], 0
        for f in faces:
            owners = set.intersection(*[node2ref.get(a, set()) for a in f])
            if not owners:
                outf.append(f); continue
            if "slave" in name or "master" in name:
                raise RuntimeError(f"contact surface {name} touches refinement zone")
            touched += 1
            a, b, c, d = f
            g = {}
            for s, t in itertools.product(range(n + 1), repeat=2):
                W = {a: (n - s) * (n - t) * n, b: s * (n - t) * n, c: s * t * n, d: (n - s) * t * n}
                g[s, t] = node_for(W)
                added_nodes.setdefault(name, set()).add(g[s, t])
            for s, t in itertools.product(range(n), repeat=2):
                outf.append([g[s, t], g[s + 1, t], g[s + 1, t + 1], g[s, t + 1]])
        body2 = "".join(f'\t\t\t<quad4 id="{i+1}">{",".join(map(str, f))}</quad4>\n' for i, f in enumerate(outf))
        return f'<Surface name="{name}">\n{body2}\t\t</Surface>'
    out_txt = re.sub(r'<Surface name="([^"]+)">(.*?)</Surface>', sub_surface, out_txt, flags=re.S)
    if len(newpos) != n_written:
        raise RuntimeError("surface subdivision created nodes not present in the refined elements")
    # node sets named like surfaces: add new nodes
    for name, nodes in added_nodes.items():
        m = re.search(rf'(<NodeSet name="{name}">)(.*?)(</NodeSet>)', out_txt, re.S)
        if m:
            have = {int(x) for x in re.findall(r'<node id="(\d+)"/>', m.group(2))}
            extra = "".join(f'\t\t\t<node id="{x}"/>\n' for x in sorted(nodes - have))
            out_txt = out_txt[:m.end(2)] + extra + "\t\t" + out_txt[m.end(2):]
    # linear constraints
    lc = [f'\t\t<constraint type="linear constraint">\n\t\t\t<tol>{TOL}</tol>\n\t\t\t<penalty>{PENALTY}</penalty>\n\t\t\t<maxaug>{MAXAUG}</maxaug>\n']
    for h, par in sorted(hang.items()):
        for dof in "xyz":
            lc.append(f'\t\t\t<linear_constraint><node id="{h}" bc="{dof}">1</node>' +
                      "".join(f'<node id="{a}" bc="{dof}">{-w:.12g}</node>' for a, w in par) + "</linear_constraint>\n")
    lc.append("\t\t</constraint>\n")
    out_txt = out_txt.replace("<Constraints>\n\t</Constraints>", "<Constraints>\n" + "".join(lc) + "\t</Constraints>", 1)
    assert "linear constraint" in out_txt
    # logging: element stresses of children of logged parents
    m = re.search(r'<element_data data="sx;sy;sz;sxy;syz;sxz">([^<]+)</element_data>', out_txt)
    logged = [int(x) for x in m.group(1).split(",")]
    S = json.load(open(os.path.join(HERE, "metric_sets.json")))
    logged = sorted(set(logged) | set(S["core"]) | set(S["rim"]))
    newlog = [ce for e in logged if e in E for ce, _ in (children[e] if e in ref else [(e, None)])]
    out_txt = out_txt[:m.start(1)] + ",".join(map(str, newlog)) + out_txt[m.end(1):]
    # logging: nodes (original list + every new node on the NP-defect interface patch)
    mouth = set(S["mouth"])
    iface_new = [h for h, k in support.items() if {a for a, _ in k} <= mouth]
    m = re.search(r'<node_data data="x;y;z">([^<]+)</node_data>', out_txt)
    zone_nodes = {a for ch in children.values() for _, cc in ch for a in cc}
    base_list = [int(x) for x in m.group(1).split(",")]
    lst = base_list + sorted((set(iface_new) | zone_nodes) - set(base_list))   # all refined-zone nodes (diagnostics)
    out_txt = out_txt[:m.start(1)] + ",".join(map(str, lst)) + out_txt[m.end(1):]
    # solver control: allow finer automatic step cuts in the refined models (static, path-independent)
    out_txt = out_txt.replace("<dtmin>0.01</dtmin>", "<dtmin>0.001</dtmin>", 1).replace("<max_retries>10</max_retries>", "<max_retries>30</max_retries>", 1)
    out_txt = out_txt.replace("<max_ups>20</max_ups>", "<max_ups>0</max_ups>", 1)   # full Newton: BFGS updates diverge with the constraint penalty terms
    assert "<max_ups>0</max_ups>" in out_txt
    assert "<dtmin>0.001</dtmin>" in out_txt and "<max_retries>30</max_retries>" in out_txt
    open(out, "w", encoding="ISO-8859-1").write(out_txt)
    meta = dict(base=os.path.basename(base), n=n, refined_parents=len(ref), n_children=len(ref) * n**3,
                new_nodes=len(newpos), hanging_nodes=len(hang), constraints=3 * len(hang),
                elements_total=len(E) - len(ref) + len(ref) * n**3, nodes_total=len(N) + len(newpos),
                children={str(e): [ce for ce, _ in children[e]] for e in ref},
                iface_new_nodes=sorted(iface_new),
                surfaces_subdivided={k: len(v) for k, v in added_nodes.items()})
    json.dump(meta, open(out.replace(".feb", "_refine.json"), "w"))
    print({k: v for k, v in meta.items() if k not in ("children", "iface_new_nodes")})

if __name__ == "__main__":
    refine(sys.argv[1], sys.argv[2], int(sys.argv[3]))
