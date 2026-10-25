#!/usr/bin/env python3
"""Convert the published FEBio 2.5 L1-L2 lumbar model to FEBio 4.12-compatible
input files for a 3-plane (+/-) range-of-motion validation, preserving the
original annulus fibrosus fiber field exactly.

Fixes applied (2.5 -> 4.12):
  1. fiber-exp-pow orientation: bare <theta>/<phi>  ->  <fiber type="angles">..</fiber>
  2. local material axis hoisted from each fiber solid up to the solid-mixture
     material level as <mat_axis type="local">; per-fiber mat_axis/fiber-local removed
  3. contact type "sliding2"  ->  "sliding-facet-on-facet"  (direct F2F successor)
  4. logfile rigid_body_data target material 1 (non-rigid) -> 11 (loaded superior endplate)
Each variant differs ONLY in the moment applied to the superior rigid body (mat 11).
"""
import re, os, sys

# usage: convert.py [SOURCE_FEB] [PREFIX]   (defaults: 49F L1-L2)
# Works for any segment in this dataset: rigid bodies (mat 11 superior / 14 inferior),
# the Rx moment line, the logfile bug, and the annulus fiber spec are identical across levels.
REPO = "../FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine"
SRC = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/49F_L12_hybrid.feb"
PREFIX = sys.argv[2] if len(sys.argv) > 2 else "L12"
OUT = "."

txt = open(SRC, encoding="ISO-8859-1").read()
a, b = txt.index("<Material>"), txt.index("</Material>")
head, mat_sec, tail = txt[:a], txt[a:b], txt[b:]

def fix_material(m):
    mid, name, typ, body = m.groups()
    if typ != "solid mixture" or "fiber-exp-pow" not in body:
        return m.group(0)
    triples = set(x.strip() for x in
                  re.findall(r'<mat_axis type="local">([^<]+)</mat_axis>', body)
                + re.findall(r'<fiber type="local">([^<]+)</fiber>', body))
    assert len(triples) == 1, (mid, triples)
    triple = triples.pop()
    body = re.sub(r'[ \t]*<mat_axis type="local">[^<]+</mat_axis>\n?', '', body)
    body = re.sub(r'[ \t]*<fiber type="local">[^<]+</fiber>\n?', '', body)
    body = re.sub(r'(?P<ind>[ \t]*)<theta>([^<]*)</theta>\s*<phi>([^<]*)</phi>',
                  lambda x: (f'{x.group("ind")}<fiber type="angles">'
                             f'<theta>{x.group(2)}</theta><phi>{x.group(3)}</phi></fiber>'),
                  body)
    return (f'<material id="{mid}" name="{name}" type="{typ}">\n'
            f'\t\t\t<mat_axis type="local">{triple}</mat_axis>{body}</material>')

mat_sec = re.sub(r'<material id="(\d+)" name="([^"]+)" type="([^"]+)">(.*?)</material>',
                 fix_material, mat_sec, flags=re.S)
base = head + mat_sec + tail

# global (non-material) fixes
base = base.replace('type="sliding2"', 'type="sliding-facet-on-facet"')
# symmetric_stiffness is not a valid sliding-facet-on-facet tag in v4 (solver-level now);
# the F2F non-symmetric form converges better per the manual. Strip the 2 active contact tags.
assert base.count('<symmetric_stiffness>1</symmetric_stiffness>') == 2
base = re.sub(r'[ \t]*<symmetric_stiffness>1</symmetric_stiffness>\n?', '', base)
base = base.replace('<rigid_body_data data="qx;qy;qz;qw">1</rigid_body_data>',
                    '<rigid_body_data data="qx;qy;qz;qw">11</rigid_body_data>')
base = base.replace('<rigid_body_data data="Mx;My;Mz">1</rigid_body_data>',
                    '<rigid_body_data data="Mx;My;Mz">11</rigid_body_data>')

# integrity checks
assert base.count('<fiber type="angles">') == 40
assert base.count('<mat_axis type="local">') == 20
assert '<fiber type="local">' not in base
assert 'sliding2' not in base
assert base.count('type="sliding-facet-on-facet"') == 2

FORCE_OLD = '<force bc="Rx" lc="2">7500</force>'
assert base.count(FORCE_OLD) == 1
VARIANTS = {
    f"{PREFIX}_flexion":     ("Rx",  7500),
    f"{PREFIX}_extension":   ("Rx", -7500),
    f"{PREFIX}_latbend_pos": ("Ry",  7500),
    f"{PREFIX}_latbend_neg": ("Ry", -7500),
    f"{PREFIX}_axrot_pos":   ("Rz",  7500),
    f"{PREFIX}_axrot_neg":   ("Rz", -7500),
}
for nm, (bc, val) in VARIANTS.items():
    new = base.replace(FORCE_OLD, f'<force bc="{bc}" lc="2">{val}</force>')
    open(os.path.join(OUT, nm + ".feb"), "w", encoding="ISO-8859-1").write(new)
    print(f"  {nm:18s} bc={bc} M={val:+5d} N.mm")
print("wrote", len(VARIANTS), "FEBio 4.12 variants")
