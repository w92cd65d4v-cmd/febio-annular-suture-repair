#!/usr/bin/env python3
"""Create the posterolateral-left full-thickness box annulotomy in the L4-L5 model.

Deletes 18 posterior-annulus hex elements forming a transannular slot
(~4.0 mm high x ~3.2 mm wide x full radial thickness, all 5 lamellae), at the
paracentral posterolateral-left corner (x~-5.2, posterior, mid-height z~15.1).
Nucleus left intact. Then cleans up: orphan nodes, the affected pressure-load
faces, and any node-set references. Applied to each of the 6 converted (FEBio
4.12-ready, fiber-faithful) L4-L5 ROM load cases.
"""
import re, os, json

SRC_DIR = "../validation_L45"
OUT_DIR = "."
DEFECT_IDS = [32480, 32488, 32520, 32528, 32536, 32544, 32568, 32576, 32584,
              32592, 32616, 32624, 32632, 32640, 32664, 32672, 32680, 32688]
LOADS = ["flexion", "extension", "latbend_pos", "latbend_neg", "axrot_pos", "axrot_neg"]

def apply_defect(txt):
    dset = set(DEFECT_IDS)
    # 1. nodes belonging to the deleted elements
    delnodes = set()
    for eid in DEFECT_IDS:
        m = re.search(r'<elem id="%d">\s*([^<]+)</elem>' % eid, txt)
        delnodes.update(int(v) for v in m.group(1).split(','))
    # 2. delete the 18 element entries
    for eid in DEFECT_IDS:
        txt, n = re.subn(r'[ \t]*<elem id="%d">[^<]*</elem>\n' % eid, '', txt)
        assert n == 1, (eid, n)
    # 3. orphan nodes = del-element nodes no longer used by any remaining element
    used = set()
    for blk in re.finditer(r'<Elements [^>]*>(.*?)</Elements>', txt, re.S):
        for e in re.finditer(r'<elem id="\d+">\s*([^<]+)</elem>', blk.group(1)):
            used.update(int(v) for v in e.group(1).split(','))
    orphans = delnodes - used
    # 4. remove orphan <node> defs and node-set refs
    for nid in orphans:
        txt = re.sub(r'[ \t]*<node id="%d">[^<]*</node>\n' % nid, '', txt)
        txt = re.sub(r'[ \t]*<node id="%d"/>\n' % nid, '', txt)
    # 5. remove pressure-load faces lying on deleted elements (all 4 nodes deleted)
    removed_faces = 0
    def strip_faces(m):
        nonlocal removed_faces
        body = m.group(2)
        def keep(fm):
            nonlocal removed_faces
            ns = [int(v) for v in fm.group(1).split(',')]
            if all(x in delnodes for x in ns):
                removed_faces += 1
                return ''
            return fm.group(0)
        body = re.sub(r'[ \t]*<quad4 id="\d+">\s*([^<]+)</quad4>\n', keep, body)
        return m.group(1) + body + m.group(3)
    txt = re.sub(r'(<Surface name="PressureLoad0[12]">)(.*?)(</Surface>)', strip_faces, txt, flags=re.S)
    return txt, len(orphans), removed_faces

os.makedirs(OUT_DIR, exist_ok=True)
for ld in LOADS:
    src = os.path.join(SRC_DIR, f"L45_{ld}.feb")
    txt = open(src, encoding="ISO-8859-1").read()
    txt, norph, nfaces = apply_defect(txt)
    # sanity: all 18 gone, no dangling refs to deleted elements
    for eid in DEFECT_IDS:
        assert f'<elem id="{eid}">' not in txt
    out = os.path.join(OUT_DIR, f"L45d_{ld}.feb")
    open(out, "w", encoding="ISO-8859-1").write(txt)
    print(f"  L45d_{ld}: removed 18 elems, {norph} orphan nodes, {nfaces} pressure faces -> {out}")
print("done")
