# Tornado plot chain (Figure4: TF-binding + histone-mark tornado)

Reproduces `tornado_bothSigTF_and_background_combined.pdf`/`.png`: Panel A
(36 TFs significant in both KDM6A_ko_vs_WT and KMT2D_Het_vs_WT, tornado by
TF family) + Panel B (TOBIAS-bound footprints vs. non-footprinted background,
binding map), both over 6 tracks (KMT2D/MLL4 ChIP, KDM6A/UTX ChIP,
FLAG-control ChIP, H3K4me1, H3K27ac, H3K27me3).

## Run order

1. `00_fetch_encode_tracks.sh` -- fully automated. Downloads 3 public ENCODE
   H1-hESC bigWigs into `data/external/encode_h1hesc/`.
2. `01_fetch_akiyama_chip.md` -- **not automated, manual/external step.**
   Documents the GEO/SRA accessions for the MLL4/UTX/FLAG-control ChIP-seq
   (Akiyama et al., GSE301295) and what to align and where to place the
   resulting 3 bigWigs (`data/external/akiyama_chip/`).
3. `02_define_bound_highconf_canon.sh` -- canonical-chromosome filter over
   the merged high-confidence WT-bound TOBIAS footprint set. **Has its own
   documented upstream gap** -- see the script's header: the merged input
   file itself has no generating script yet.
4. `03_define_background_regions.R` -- defines the Panel B background region
   set (consensus ATAC peaks minus footprinted regions, downsampled,
   seed=42). **Not a byte-identical reproduction of the original figure** --
   the original background set's exact derivation is lost; this is a fresh,
   documented, scientifically equivalent replacement. Say so wherever this
   panel is described.
5. `bothsig36_beds.sh` -- per-TF top-200 high-confidence bound sites for the
   36 both-significant TFs (Panel A).
6. `bothsig36_matrix.sh` + `edu04_matrix2.sh` -- deepTools `computeMatrix`
   for Panel A and Panel B respectively.
7. `bothsig36_combined_figure.py` -- renders the combined two-panel figure.
   Color-scale ceilings are now computed on the fly from the matrices
   (99th percentile per track) rather than read from a separate
   `shared_vmax.json`, which had the same untraced-scratch-session problem
   as the files above and has been removed as a dependency entirely.

Steps 5-7 write into `results/atac/tornado_intermediate/` (gitignored,
regenerable) and the final figure into
`results/figures/tf_footprinting_network/`.

## What's fully automated vs. not

Fully automated and deterministic: steps 00, 02 (given its input), 03
(given its inputs, and given the seed), 5, 6, 7.

Requires external, one-time manual work: step 01 (ChIP-seq realignment from
raw reads -- a full companion pipeline, cited not vendored) and the
`all_WT_bound_highconf.bed` merge referenced in step 02's header (no script
exists yet; needs writing, or the file supplied).
