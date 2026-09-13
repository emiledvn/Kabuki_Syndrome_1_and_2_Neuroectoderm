#!/usr/bin/env python3
"""Combined figure for the 36 TFs significant in BOTH contrasts (A07b):
  A = per-TF tornado of those 36 TFs (grouped by family), taller boxes
  B = footprint vs accessibility-matched background (the binding map), as in
      tornado_perTF_and_background_combined.
Shared 6 tracks, +/-3kb, one shared per-column colour scale."""
import gzip, json, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SCR = os.path.join(REPO_ROOT, "results", "atac", "tornado_intermediate")
OUT = os.path.join(REPO_ROOT, "results", "figures", "tf_footprinting_network",
                    "tornado_bothSigTF_and_background_combined.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

TRACK_NAMES = ["KMT2D", "KDM6A", "FLAG", "H3K4me1", "H3K27ac", "H3K27me3"]
TRACK_COLORS = {"KMT2D": "#7570B3", "KDM6A": "#1B9E77", "FLAG": "#4D4D4D",
                "H3K4me1": "#1baf7a", "H3K27ac": "#eb6834", "H3K27me3": "#4a3aa7"}


def load(path):
    with gzip.open(path, "rt") as fh:
        hdr = json.loads(fh.readline().lstrip("@"))
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    return hdr, np.array([[float(x) if x not in ("nan", "") else 0.0 for x in r[6:]] for r in rows])


hA, dA = load(f"{SCR}/matrix_bothsig36.gz")
hB, dB = load(f"{SCR}/matrix_edu04_6col.gz")
sbA, gbA = hA["sample_boundaries"], hA["group_boundaries"]
sbB, gbB = hB["sample_boundaries"], hB["group_boundaries"]

# Per-track colour-scale ceiling: 99th percentile of the pooled signal across
# both panels' matrices for that track. Originally read from a static
# shared_vmax.json produced in an untraced scratch session; computed here
# instead so the figure has no external artefact dependency. FLAG shares
# KMT2D's ceiling by design (same MLL4 antibody y-scale, for visual
# comparability against the KMT2D-specific signal).
def track_vmax(ci, pct=99):
    pooled = np.concatenate([dA[:, sbA[ci]:sbA[ci + 1]].ravel(), dB[:, sbB[ci]:sbB[ci + 1]].ravel()])
    return float(np.percentile(pooled, pct))

vmax_by_track = {name: track_vmax(ci) for ci, name in enumerate(TRACK_NAMES)}
vmax_by_track["FLAG"] = vmax_by_track["KMT2D"]

COLS = [(name, TRACK_COLORS[name], vmax_by_track[name]) for name in TRACK_NAMES]
NCOL = len(COLS)
CMAPS = [LinearSegmentedColormap.from_list("c", ["white", c]) for _, c, _ in COLS]
up, down, bs = hA["upstream"][0], hA["downstream"][0], hA["bin size"][0]
nb = sbA[1] - sbA[0]
xt = [0, nb // 2, nb - 1]
xtl = [f"-{up/1000:g}", "0", f"+{down/1000:g}"]

# TF order + families
fam_rows, syms, fams = [], [], []
with open(f"{SCR}/bothsig36_order.tsv") as fh:
    for fam, sym, motif, n in csv.reader(fh, delimiter="\t"):
        fams.append(fam); syms.append(sym)
n_gA = len(syms)
assert n_gA == len(gbA) - 1
# family blocks (contiguous)
fam_blocks = []
i = 0
while i < n_gA:
    j = i
    while j < n_gA and fams[j] == fams[i]:
        j += 1
    fam_blocks.append((fams[i], i, j))
    i = j

# sort rows within each TF by KMT2D+KDM6A mean
for g in range(n_gA):
    g0, g1 = gbA[g], gbA[g + 1]
    k = dA[g0:g1, sbA[0]:sbA[1]].mean(1) + dA[g0:g1, sbA[1]:sbA[2]].mean(1)
    dA[g0:g1] = dA[g0:g1][np.argsort(-k)]
totalA = gbA[-1]

fig = plt.figure(figsize=(18.5, 8.6))

# ===================== PANEL A : 36 both-sig TFs =====================
FAM_SHORT = {"AP-1 (Fos/Jun)": "AP-1", "bZIP (BACH/Maf)": "bZIP", "FOX (Forkhead)": "FOX",
             "Nuclear receptor": "NR", "RFX": "RFX", "SOX (HMG-box)": "SOX",
             "Zinc finger": "ZnF", "Homeodomain": "HD"}
A_TOP, A_BOT = 0.925, 0.130
gsA = fig.add_gridspec(1, NCOL + 1, left=0.575, right=0.985, top=A_TOP, bottom=A_BOT,
                       width_ratios=[0.60] + [1] * NCOL, wspace=0.34)
axm = fig.add_subplot(gsA[0, 0])
axm.set_xlim(0, 1); axm.set_ylim(totalA, 0); axm.set_xticks([]); axm.set_yticks([])
for sp in axm.spines.values():
    sp.set_visible(False)
axm.axvline(1.0, color="#333333", lw=1.0, clip_on=False)
for g in range(n_gA):
    mid = (gbA[g] + gbA[g + 1]) / 2
    axm.plot([0.90, 1.0], [mid, mid], color="#333333", lw=1.0, clip_on=False)
    axm.text(0.85, mid, syms[g], fontsize=8, ha="right", va="center")

# family brackets + short codes, in the figure margin left of panel A
from matplotlib.lines import Line2D
def yfig(row):
    return A_TOP - (row / totalA) * (A_TOP - A_BOT)
for fam, i0, i1 in fam_blocks:
    yt, yb = yfig(gbA[i0]), yfig(gbA[i1])
    p = (yt - yb) * 0.06
    fig.add_artist(Line2D([0.560, 0.560], [yb + p, yt - p], transform=fig.transFigure,
                          color="#666666", lw=1.7, solid_capstyle="round"))
    fig.text(0.546, (yt + yb) / 2, FAM_SHORT[fam], rotation=90, fontsize=9,
             fontweight="bold", ha="center", va="center", color="#333333")

A_axes = []
for ci, (label, color, vmax) in enumerate(COLS):
    ax = fig.add_subplot(gsA[0, ci + 1], sharey=axm)
    ax.imshow(dA[:, sbA[ci]:sbA[ci + 1]], aspect="auto", cmap=CMAPS[ci],
              vmin=0, vmax=vmax, interpolation="nearest")
    ax.set_title(label, fontsize=12, fontweight="bold")
    ax.set_xticks(xt); ax.set_xticklabels(xtl, fontsize=9)
    ax.set_yticks([])
    for g in range(1, n_gA):
        ax.axhline(gbA[g], color="#e6e6e6", lw=0.4)
    for fam, i0, i1 in fam_blocks[1:]:
        ax.axhline(gbA[i0], color="#8a8a8a", lw=0.9)
    for sp in ax.spines.values():
        sp.set_linewidth(0.7); sp.set_color("#9a9a9a")
    A_axes.append(ax)

# ===================== PANEL B : footprint vs background (binding map) =====================
GROUP_LABELS = ["TF-bound\n(high-conf footprint)",
                "non-footprinted consensus\nATAC peak (Background)"]
grpsB = [0, 1]
heightsB = [gbB[g + 1] - gbB[g] for g in grpsB]


def sort_within_B(g0, g1):
    k = dB[g0:g1, sbB[0]:sbB[1]].mean(1) + dB[g0:g1, sbB[1]:sbB[2]].mean(1)
    return g0 + np.argsort(-k)


ordersB = [sort_within_B(gbB[g], gbB[g + 1]) for g in grpsB]
grp_c = ["#111111", "#c0392b"]
# KMT2D profile ceiling -- FLAG profile is drawn on this same (MLL4) y-scale
kmt2d_prof_top = max(dB[gbB[g]:gbB[g + 1], sbB[0]:sbB[1]].mean(0).max() for g in grpsB) * 1.08
gsB = fig.add_gridspec(2, NCOL, left=0.120, right=0.470, top=0.925, bottom=0.130,
                       height_ratios=[1.0, 6.2], wspace=0.46, hspace=0.07)
B_axes = []
for ci, (label, color, vmax) in enumerate(COLS):
    axp = fig.add_subplot(gsB[0, ci])
    for gi, g in enumerate(grpsB):
        m = dB[gbB[g]:gbB[g + 1], sbB[ci]:sbB[ci + 1]].mean(0)
        axp.plot(m, lw=1.7, color=grp_c[gi], ls="-" if gi == 0 else (0, (4, 2)))
    axp.set_title(label, fontsize=12.5, fontweight="bold")
    axp.set_xticks(xt); axp.set_xticklabels([])
    axp.tick_params(labelsize=9, pad=1.5)
    axp.yaxis.set_major_locator(plt.MaxNLocator(2))
    axp.margins(x=0)
    axp.set_ylim(0, kmt2d_prof_top if label == "FLAG" else None)
    if ci == 0:
        axp.set_ylabel("mean signal", fontsize=10.5)
    axh = fig.add_subplot(gsB[1, ci])
    stacked = np.vstack([dB[gbB[grpsB[gi]]:gbB[grpsB[gi] + 1], sbB[ci]:sbB[ci + 1]][ordersB[gi] - gbB[grpsB[gi]]]
                         for gi in range(2)])
    axh.imshow(stacked, aspect="auto", cmap=CMAPS[ci], vmin=0, vmax=vmax, interpolation="nearest")
    axh.set_xticks(xt); axh.set_xticklabels(xtl, fontsize=9)
    axh.set_yticks([])
    axh.axhline(heightsB[0], color="k", lw=0.9)
    if ci == 0:
        yc, y0 = [], 0
        for h in heightsB:
            yc.append(y0 + h / 2); y0 += h
        axh.set_yticks(yc); axh.set_yticklabels(GROUP_LABELS, fontsize=9.5)
    for sp in axh.spines.values():
        sp.set_linewidth(0.7)
    B_axes.append(axh)

fig.legend([plt.Line2D([], [], color="#111111", ls="-", lw=1.8),
            plt.Line2D([], [], color="#c0392b", ls=(0, (4, 2)), lw=1.8)],
           ["footprint", "background"],
           loc="upper left", bbox_to_anchor=(0.120, 0.978), ncol=2, frameon=False,
           fontsize=10, columnspacing=3, handlelength=1.9)

fig.canvas.draw()
for axes in (A_axes, B_axes):
    for ci, ax in enumerate(axes):
        pos = ax.get_position()
        cax = fig.add_axes([pos.x0, 0.082, pos.width, 0.012])
        cb = fig.colorbar(ax.images[0], cax=cax, orientation="horizontal")
        cb.set_ticks([0, COLS[ci][2]]); cb.ax.set_xticklabels(["0", f"{COLS[ci][2]:.1f}"], fontsize=8)
        cb.ax.tick_params(length=2, pad=1.5); cb.outline.set_linewidth(0.4)

fig.text(0.012, 0.990, "A", fontsize=24, fontweight="bold", va="top")
fig.text(0.520, 0.990, "B", fontsize=24, fontweight="bold", va="top")
fig.text((0.575 + 0.985) / 2, 0.040, "distance from motif centre (kb)", ha="center", fontsize=11)
fig.text((0.120 + 0.470) / 2, 0.040, "distance from motif centre (kb)", ha="center", fontsize=11)

fig.savefig(OUT, dpi=300)
fig.savefig(OUT.replace(".png", ".pdf"))
print("wrote", OUT)
