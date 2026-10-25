#!/usr/bin/env python3
"""Revised main-text Figures 2, 3, 5 and 6 for the round-1 revision (Reviewer 2 figure comments).
Fig 2 = R1 (larger type); Fig 3 = dose-response emphasised, mesh-valid vs distorted separated, force-based apposition added;
Fig 5 = simplified single panel; Fig 6 = summary of containment recovery versus rim stress (replaces the rim-stress bars).
The submitted figures are left untouched; outputs are results_figures/fig{2,3,5,6}_rev1.{png,pdf,tiff}."""
import os, csv, json, re, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from make_results_figs import series, logp, GAP, SJ           # original loaders (no plotting at import)
from revision_sens_analyze import metrics, parse as sparse

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "results_figures")
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GREY, INK, MUTED = "#8a8a85", "#1f1f1e", "#6b6b66"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12.5,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10.5, "axes.linewidth": 0.8,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e2",
                     "grid.linewidth": 0.6, "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
def save(fig, name):
    for ext in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{name}.tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig); print("wrote", name)
MASTER = {(r["loading"], r["state"]): r for r in csv.DictReader(open(os.path.join(HERE, "..", "RESULTS_MASTER_TABLE.csv")))}
SC = {"pure": "Pure axial compression", "flex": "Flexion + axial compression"}
TI = os.path.join(HERE, "tighten")
PAIRS = json.load(open(os.path.join(HERE, "suture_config.json")))["suture3"]
def tight(level, scen):
    stem = os.path.join(TI, f"tight_{level}_{scen}")
    m = metrics(stem + ".log", t0=2.0)
    X = dict(sparse(stem + ".log")); x2 = X[min(X, key=lambda t: abs(t - 2.0))]["x"]
    m["gap"] = float(np.mean([np.linalg.norm(x2[a] - x2[b]) for a, b in PAIRS])); return m
def rest(scen, auc):
    d, i = float(MASTER[(scen, "defect")]["positive_AUC_mmN"]), float(MASTER[(scen, "intact")]["positive_AUC_mmN"])
    return 100 * (d - auc) / (d - i)

# ---------------- Figure 2 ----------------
def fig2():
    states = [("intact", "Intact", GREY, ":"), ("defect", "Unrepaired defect", INK, "-"), ("suture1", "Open gap, 1 load-targeted suture", C4, "-"),
              ("bridge", "Open gap, 3 load-targeted sutures", C2, "--"), ("app12", "3 sutures, partial apposition (app12)", C1, "-")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), sharey=True)
    for ax, scen in zip(axes, ("pure", "flex")):
        for st, lab, col, ls in states:
            r = series(logp(st, scen)); F = [x[0] for x in r]; P = [x[1] for x in r]
            ax.plot(F, P, ls, color=col, lw=2.2, label=lab, zorder=3)
            k = int(np.argmax(P)); ax.plot(F[k], P[k], "o", color=col, ms=7, mec="white", mew=1, zorder=4)
            ax.plot(F[-1], P[-1], "s", color=col, ms=6, mec="white", mew=1, zorder=4)
        ax.axhline(0, color=GREY, lw=0.8); ax.set_title(SC[scen]); ax.set_xlabel("Axial compressive force (N)")
    axes[0].set_ylabel("Contained protrusion (mm)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    axes[1].text(0.98, 0.03, "● path maximum   ■ 2,000-N endpoint", transform=axes[1].transAxes, ha="right", fontsize=10, color=MUTED)
    fig.tight_layout(); save(fig, "fig2_rev1")

# ---------------- Figure 3 ----------------
def fig3():
    valid = ["app17", "app15", "app12", "app10b", "app08"]; ub = ["app06", "near"]
    MT = {}
    for scen in ("pure", "flex"):
        for st in valid + ub:
            r = series(logp(st, scen)); P = np.array([x[1] for x in r]); F = np.array([x[0] for x in r])
            MT[(scen, st)] = (P.max(), np.trapz(np.maximum(P, 0), F))
    TT = {(lv, scen): tight(lv, scen) for lv in ("app17", "app12") for scen in ("pure", "flex")}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.9), gridspec_kw=dict(width_ratios=[1, 1, 0.7]))
    sty = {"pure": dict(ls="-", marker="o", color=C1), "flex": dict(ls="--", marker="s", color=C3)}
    for k, (ax, ylab) in enumerate(zip(axes[:2], ("Peak contained protrusion (mm)", "Positive-only AUC (mm·N)"))):
        ax.axvspan(0.5, 0.70, color="#f1e3dc", lw=0, zorder=0)
        for scen in ("pure", "flex"):
            s = sty[scen]
            ax.plot([GAP[x] for x in valid], [MT[(scen, x)][k] for x in valid], ls=s["ls"], marker=s["marker"], color=s["color"], ms=7, lw=2,
                    label=f"Morphed apposition — {SC[scen].lower()}")
            ax.plot([GAP[x] for x in ub], [MT[(scen, x)][k] for x in ub], s["marker"], ms=7, mfc="white", mec=s["color"], mew=1.5, ls="none")
            ref = float(MASTER[(scen, "intact")]["peak_prot_mm" if k == 0 else "positive_AUC_mmN"])
            ax.axhline(ref, color=s["color"], lw=1, ls=":", alpha=0.8)
            tx = [TT[(lv, scen)]["gap"] for lv in ("app17", "app12")]; ty = [TT[(lv, scen)]["peak" if k == 0 else "posAUC"] for lv in ("app17", "app12")]
            ax.plot(tx, ty, "D", ms=8, color=C2 if scen == "pure" else "white", mec=C2, mew=1.6, ls="none",
                    label=f"Tightened + contact — {SC[scen].lower()}")
        ax.set_xlim(1.8, 0.5); ax.set_xlabel("Mean residual anchor-pair gap (mm)"); ax.set_ylabel(ylab)
        ax.text(0.6, ax.get_ylim()[1], "distorted\n(upper bound)", ha="center", va="top", fontsize=9.5, color="#a5461c")
    axes[0].text(1.78, float(MASTER[("pure", "intact")]["peak_prot_mm"]), "intact (dotted)", fontsize=9.5, color=MUTED, va="bottom")
    ax = axes[2]; ax.axvspan(0.5, 0.70, color="#f1e3dc", lw=0, zorder=0)
    ax.plot([GAP[x] for x in valid], [SJ[x] for x in valid], "o", color=INK, ms=7)
    ax.plot([GAP[x] for x in ub], [SJ[x] for x in ub], "o", ms=7, mfc="white", mec=INK, mew=1.5)
    ax.axhline(0.2, color="#a5461c", ls="--", lw=1.2); ax.text(1.78, 0.205, "0.2 screening threshold", fontsize=9.5, color="#a5461c", va="bottom")
    ax.set_xlim(1.8, 0.5); ax.set_xlabel("Mean residual anchor-pair gap (mm)"); ax.set_ylabel("Minimum scaled Jacobian")
    for a, lab in zip(axes, "ABC"): a.text(-0.16, 1.03, lab, transform=a.transAxes, fontsize=15, fontweight="bold")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.42, -0.13))
    fig.tight_layout(); save(fig, "fig3_rev1")
    return TT

# ---------------- Figure 5 ----------------
def fig5():
    rows = list(csv.DictReader(open(os.path.join(HERE, "..", "RESULTS_SECONDARY_TABLE.csv"))))
    def get(scen, var): return float(next(r for r in rows if r["loading"] == scen and r["variant"] == var)["positive_AUC_mmN"])
    groups = [("Suture number\n(50 N/mm)", [("1", "1 suture (E50)"), ("2", "2 sutures (E50)"), ("3", "3 sutures (E50)=bridge")]),
              ("Axial spring stiffness,\n1 suture (N/mm)", [("10", "1 suture E10"), ("50", "1 suture (E50)"), ("200", "1 suture E200")]),
              ("Axial spring stiffness,\n3 sutures (N/mm)", [("10", "3 sutures E10"), ("50", "3 sutures (E50)=bridge"), ("200", "3 sutures E200")])]
    fig, ax = plt.subplots(figsize=(10, 4.4)); x0 = 0; ticks, tl = [], []
    for gi, (gname, items) in enumerate(groups):
        for j, (lab, var) in enumerate(items):
            x = x0 + j
            for scen, mk, col, dx in (("pure", "o", C1, -0.12), ("flex", "s", C3, 0.12)):
                ax.plot(x + dx, rest(scen, get(scen, var)), mk, color=col, ms=9, mec="white", mew=1, label=SC[scen] if (gi, j) == (0, 0) else None)
            ticks.append(x); tl.append(lab)
        ax.text(x0 + 1, -24, gname, ha="center", va="top", fontsize=11, color=INK)
        x0 += 4
    ax.axhline(100, color=GREY, ls=":", lw=1); ax.text(10.6, 101, "intact", color=MUTED, fontsize=10, va="bottom", ha="right")
    ax.axhline(0, color=INK, lw=0.8); ax.text(10.6, 1, "unrepaired defect", color=MUTED, fontsize=10, va="bottom", ha="right")
    ax.set_xticks(ticks); ax.set_xticklabels(tl); ax.set_ylim(-5, 110); ax.set_xlim(-0.6, 10.8)
    ax.set_ylabel("Positive-only AUC restoration (%)"); ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62))
    fig.subplots_adjust(bottom=0.25); save(fig, "fig5_rev1")

# ---------------- Figure 6 ----------------
def fig6(TT):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), sharey=True)
    valid = ["app17", "app15", "app12", "app10b", "app08"]
    for ax, scen in zip(axes, ("pure", "flex")):
        M = lambda st, k: float(MASTER[(scen, st)][k])
        ri = M("intact", "rim_VM_endpoint_MPa")
        ax.axhline(ri, color=GREY, ls=":", lw=1.2); ax.axvline(100, color=GREY, ls=":", lw=1.2)
        ax.plot(100, ri, "*", ms=15, color=GREY, mec="white", zorder=5, label="Intact")
        ax.plot(0, M("defect", "rim_VM_endpoint_MPa"), "X", ms=11, color=INK, mec="white", zorder=5, label="Unrepaired defect")
        for st, mk, lab in (("suture1", "o", "Open gap, 1 suture"), ("bridge", "s", "Open gap, 3 sutures")):
            ax.plot(M(st, "posAUC_restoration_pct"), M(st, "rim_VM_endpoint_MPa"), mk, ms=9, color=C4, mec="white", zorder=5, label=lab)
        ax.plot([M(s, "posAUC_restoration_pct") for s in valid], [M(s, "rim_VM_endpoint_MPa") for s in valid], "-o", color=C1, ms=7, lw=1.6,
                mec="white", zorder=4, label="Morphed apposition (app17→app08)")
        ax.plot([M(s, "posAUC_restoration_pct") for s in ("app06", "near")], [M(s, "rim_VM_endpoint_MPa") for s in ("app06", "near")], "o",
                ms=7, mfc="white", mec=C1, mew=1.5, ls="none", label="Morphed, distorted upper bound")
        tx = [rest(scen, TT[(lv, scen)]["posAUC"]) for lv in ("app17", "app12")]; ty = [TT[(lv, scen)]["rim_max"] for lv in ("app17", "app12")]
        ax.plot(tx, ty, "--D", color=C2, ms=9, mec="white", lw=1.4, zorder=5, label="Tightened + contact (app17, app12 gaps)")
        for lv, x, y in zip(("app17", "app12"), tx, ty): ax.annotate(lv, (x, y), xytext=(-4, 14) if lv == "app17" else (8, 0), textcoords="offset points", fontsize=10, color="#a5461c", va="center", ha="right" if lv == "app17" else "left")
        ax.text(2, ri, "intact rim stress", fontsize=10, color=MUTED, va="bottom")
        ax.set_title(SC[scen]); ax.set_xlabel("Containment recovery: positive-only AUC restoration (%)"); ax.set_xlim(-10, 150)
    axes[0].set_ylabel("Peak annular-rim von Mises stress (MPa)"); axes[0].set_ylim(0, 12.5)
    h, l = axes[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.2))
    fig.tight_layout(); save(fig, "fig6_rev1")

if __name__ == "__main__":
    fig2(); TT = fig3(); fig5(); fig6(TT)
