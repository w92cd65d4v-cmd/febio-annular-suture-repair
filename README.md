# febio-annular-suture-repair

Scripts and result tables for:

> Li G-S, Chen G-C, Tu W-T, Zhang Z-P, Song E, Zhou J-J. *Biomechanical determinants of annular suture repair in an idealized post-discectomy annular defect: a finite element study.* Frontiers in Bioengineering and Biotechnology (under review).

The study modifies a public FEBio model of the human L4–L5 functional spinal unit. It creates a posterolateral full-thickness annular defect, adds suture repairs, and compares contained nucleus pulposus (NP) protrusion and annular-rim stress under axial compression with and without flexion.

**To avoid copyright concerns, the original source model is not distributed in this repository.** Users may obtain it from the external download link given in section 1. This repository contains only our own material:

- the Python scripts that modify the source model and analyze the results;
- small data files that list node and element ID numbers (no model geometry);
- the numerical result tables behind the figures and tables of the paper.

No finite element input files (`.feb`), solver logs (`.log`), or result files (`.xplt`) are included, because they contain the geometry of the source model.

## 1. Source model

Finley SM, Brodke DS, Spina NT, DeDen CA, Ellis BJ (2018). FEBio finite element models of the human lumbar spine. *Computer Methods in Biomechanics and Biomedical Engineering* 21, 444–452. doi:10.1080/10255842.2018.1478967

- Download: FEBio model repository, <https://repo.febio.org/modelRepo/api/v1.05/files/0/57> (file `FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine.zip`). The same models are described on the authors' SimTK project page, <https://simtk.org/projects/febiolumbar>. Please follow the terms of use and the citation request of the source.
- FEBio: <https://febio.org>

We used the L4–L5 segment file `49F_L45_hybrid.feb` (FEBio specification 2.5) of the 49-year-old female specimen. The zip file contains the model files without an enclosing folder. Extract it into the folder

```
FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine/
```

so that the scripts find `FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine/49F_L45_hybrid.feb`.

The file we used has this SHA-256 checksum:

```
04d44d801559cecfb6f52e14277903f1c035047c6ea1e3e2c4c011da379e1f49
```

Check it with `shasum -a 256 "FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine/49F_L45_hybrid.feb"`. A different checksum means a different version of the source file, and the element and node IDs used by our scripts may not match.

## 2. Requirements

- Python 3.9 or later with `numpy` and `matplotlib` (see `requirements.txt`).
- FEBio 4.12 command-line solver (`febio4`, part of FEBio Studio). The scripts that run FEBio themselves (`revision_tighten_driver.py`, `revision_tighten_outer.py`) use the macOS path `/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4`. Change the `FEBIO` variable at the top of these scripts for another system.
- All other models are solved with `febio4 -i <model>.feb`.

These are the research scripts as they were used for the paper. They are not a packaged software tool. Run each script from its own folder, because the scripts use relative paths.

## 3. Folder layout

```
validation_L45/              conversion to FEBio 4.12; range-of-motion benchmark of the intact model
validation_L45_defect/       creation of the 18-element posterolateral annular defect
validation_L45_herniation/   compression models, suture repairs, apposition, sensitivity analyses, mesh refinement,
                             analysis and figure scripts
results/                     numerical result tables used in the paper
```

## 4. Building the models (in order)

All steps start from the source file and write new `.feb` files next to the scripts. Element and node IDs refer to the source mesh.

| Step | Folder and command | Output models | Purpose |
|---|---|---|---|
| 1 | `validation_L45`: `python3 convert.py "../FBSRepo_User Projects_FEBio FEMs Human Lumbar Spine/49F_L45_hybrid.feb" L45` | `L45_{flexion,extension,latbend_pos,latbend_neg,axrot_pos,axrot_neg}.feb` | Conversion from FEBio 2.5 to 4.12; intact model under 7.5 N m pure moments (benchmark). `postprocess.py` extracts the range of motion. |
| 2 | `validation_L45_defect`: `python3 make_defect.py` | `L45d_*.feb` | Removes the 18 posterior annular elements (`defect_spec.json`) and cleans up the affected nodes and pressure faces. |
| 3 | `validation_L45_herniation`: `python3 make_compression.py` | `comp_{intact,defect}_{pure,flex}.feb` | Axial compression 0–2,000 N, with or without a 7.5 N m flexion moment, after the 1-MPa expansive nucleus-boundary pressure preload. |
| 4 | `python3 make_suture.py` | `comp_suture{1,3}_{pure,flex}.feb` | Open-gap repairs: tension-only springs (50 N/mm) at the load-targeted anchor pairs (`suture_config.json`). |
| 5 | `python3 make_sensitivity.py`; `python3 make_sweep.py` | `comp_defect{S,L}_*`, `comp_suture2_*`, `comp_suture3b_*`; `comp_s{1,3}_E{10,200}_*` | Defect size, suture number, off-target placement, and axial spring stiffness (10 and 200 N/mm). |
| 6 | `python3 make_apposed.py`; `python3 make_apposed3_levels.py` | `comp_suture3_apposed_*`; `comp_suture3_app{17,15,12,10b,08,06}_*` | Morphed partial edge apposition: the rest geometry is morphed to a target mean residual anchor-pair gap and reset as stress-free. |
| 7 | `python3 revision_gen_sens.py` | `sens/pl010_*`, `sens/pl050_*`, `sens/pl100_intact_*`, `sens/single_p00`–`p19_*`, `sens/outer{1,2,3}_*` | Preload sensitivity (0.1 and 0.5 MPa), independent single-suture placement at each of the 20 candidate anchor pairs, and outer-half layouts. Uses the logging variants of step 9 as bases. |
| 8 | `python3 revision_tighten_driver.py calibrate <pure\|flex> <app17\|app12\|app08>`, then `python3 revision_tighten_driver.py full <scen> <tag> c1 c2 c3`; `python3 revision_tighten_outer.py calibrate <outer1\|outer3>`, then `full` | `tighten/*`, `tighten_outer/*` | Force-based apposition: sutures tightened, knotted, and released with frictionless wall contact (`revision_gen_tighten.py`). The calibrated closures are stored in `tighten/*.json` and `tighten_outer/*.json`. These scripts run FEBio themselves. |
| 9 | see section 5 | `comp_*_{pure,flex}_map.feb`; `comp_intact_{pure,flex}_fullstress.feb` (`make_intact_fullstress_models.py`) | Logging variants: identical models that record the region around the defect. |
| 10 | `python3 revision_mesh_refine.py <BASE>_map.feb <OUT>.feb <2\|3>` | `*_h2.feb`, `*_h3.feb` | Local mesh refinement (h/2, h/3) of the zone in `REFINE_ZONE.json`, with hanging nodes tied by linear constraints. Bases: intact, defect, open-gap three-suture, and app12 logging variants. |

## 5. Steps without a dedicated script

Three small edits were made directly in the model files. They are described here so that they can be repeated.

- **Logging variants (`*_map.feb`).** Same model as the base file. In the `<logfile>` section, the three original entries are replaced by two: `<node_data data="x;y;z">` with the node IDs in `region_sets.json` (`"nodes"`, 591 nodes), and `<element_data data="sx;sy;sz;sxy;syz;sxz">` with the element IDs in `region_sets.json` (`"elems"`, 374 elements). This was done for the intact, defect, one-suture, three-suture, and app12 models under both loading conditions.
- **0.1-MPa preload benchmark.** Copies of the step-1 files in which the two nucleus-boundary pressure magnitudes are changed from `-1` and `1` to `-0.1` and `0.1` (`<pressure lc="1">`).
- **Three-suture pretension (`comp_s3pre_*`).** Copy of `comp_suture3_*` in which the suture load curve (`<loadcurve id="10">`) is shifted to give a 0.2-mm shortened rest length: points (-100, 0), (-0.2, 0), (0, 10), (100, 5010).

## 6. Analysis and result tables

The analysis scripts read the FEBio logs. The protrusion metric is the mean outward displacement of the eight non-anchor nodes of the NP–defect interface (`metric_sets.json`, `"mouth"`), referenced to the start of compression. Restoration is 100 × (defect − repair)/(defect − intact).

| Table in `results/` | Produced by | Content |
|---|---|---|
| `RESULTS_MASTER_TABLE.csv` | compiled from the FEBio logs; apposition-series values identical to `apposed3_doseresponse.json` from `analyze_levels3.py`; read by `make_results_figs.py` and `make_results_figs_R1.py` | Primary outcomes of all primary states and the apposition series (Table 2 of the paper). |
| `RESULTS_SECONDARY_TABLE.csv` | compiled from the FEBio logs; read by `make_results_figs.py` | Secondary outcomes. |
| `REVISION_SENS_RESULTS.csv` | `revision_sens_analyze.py` | Preload sensitivity, single-suture placement, outer-half layouts. |
| `REVISION_PLACEMENT.csv` | `revision_placement_analyze.py` | Independent suture placement summary. |
| `REVISION_METRIC_ROBUSTNESS.csv`, `REVISION_NODE_LAYOUT.csv` | `revision_metric_robustness.py` | Alternative protrusion definitions, leave-one-node-out, rim-stress measures; interface node layout. |
| `REVISION_MORPH_DIAGNOSTICS*.csv` | `revision_morph_diagnostics.py` | Effect of the apposition morph on the mesh and fiber architecture. |
| `REVISION_MESHSENS_PARTIAL_t4.2.csv`, `_t4.4.csv` | `revision_meshsens_partial.py 4.2` / `4.4` | Mesh sensitivity up to 711 N and 756 N. |
| `METHODS_PARAMETERS.csv` | compiled by hand from the model files | Model parameters and settings. |

Metrics of the force-based apposition runs are computed with `revision_tighten_metrics.py`. Figures: `make_fig1_R2.py` (Figure 1; uses `make_model_construction_figs.py`), `make_results_figs_R1.py` (Figures 2, 3, 5, 6), `make_stress_maps.py` (Figure 4).

## 7. Notes

- All values are numerical outputs of one deterministic model. Absolute protrusion, restoration percentages, and stresses at 2,000 N were not mesh-converged; local mesh refinement could be performed only up to about 750 N (see the paper).
- `make_suture.py` and `make_sweep.py` were reconstructed during the project from the model files they had produced; their headers say so.
- `make_manuscript_figs.py`, `make_results_figs.py`, and `make_whole_fsu_stress_view.py` first look for the apposition-series logs in a network folder used during the project and otherwise read them from `nas_copy/`.

## 8. Citation and license

If you use this material, please cite the paper above, the source model (Finley et al., 2018), and FEBio:

Maas SA, Ellis BJ, Ateshian GA, Weiss JA (2012). FEBio: finite elements for biomechanics. *Journal of Biomechanical Engineering* 134, 011005. doi:10.1115/1.4005694

- Code: MIT License (`LICENSE`).
- Result tables in `results/` and the data files: Creative Commons Attribution 4.0 (`LICENSE-DATA.md`).
- These licenses cover only the material in this repository. The source model remains under the terms set by its authors.
