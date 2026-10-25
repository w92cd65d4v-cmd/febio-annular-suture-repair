#!/usr/bin/env python3
"""Create logging-only intact reruns that include the planned annulotomy elements.

The earlier map files logged the local region but omitted the 18 elements later
deleted to create the annulotomy. These files preserve the mechanical model and
loading exactly, add those IDs to the stress logger, and suppress the plot file.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
region = json.loads((HERE / "region_sets.json").read_text())
stress_ids = sorted(set(region["elems"]) | set(region["defect"]))

for scenario in ("pure", "flex"):
    source = HERE / f"comp_intact_{scenario}_map.feb"
    target = HERE / f"comp_intact_{scenario}_fullstress.feb"
    text = source.read_text(encoding="ISO-8859-1")
    text, removed = re.subn(r"\s*<plotfile\b.*?</plotfile>", "", text, count=1, flags=re.S)
    assert removed == 1
    logger = (
        "\n\t\t<logfile>\n"
        "\t\t\t<element_data data=\"sx;sy;sz;sxy;syz;sxz\">"
        + ",".join(map(str, stress_ids))
        + "</element_data>\n\t\t</logfile>"
    )
    text, replaced = re.subn(r"\s*<logfile>.*?</logfile>", logger, text, count=1, flags=re.S)
    assert replaced == 1
    target.write_text(text, encoding="ISO-8859-1")
    print(f"wrote {target.name}: {len(stress_ids)} stress elements")
