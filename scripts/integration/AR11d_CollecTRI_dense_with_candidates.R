#!/usr/bin/env Rscript
# AR11d_CollecTRI_dense_with_candidates.R -- denser per-contrast network
# figures: the orthogonally supported ("agreeing only", AR11b) curated
# edges, PLUS a capped set of high-confidence AR06 edges whose TF is
# well-curated in CollecTRI generally but this specific (TF, gene) pair is
# not (AR11c) -- rendered as a visually distinct "candidate" layer (gray,
# dashed, lower alpha), never colored/treated as if independently
# confirmed.
#
# Rule (SAME for both contrasts, for methods consistency): fill each
# already-present TF's total shown edges (curated + candidate) up to
# MAX_TOTAL_PER_TF, prioritizing curated edges first, then topping up with
# the highest-effect-size candidates. This naturally gives sparse hubs a
# meaningful boost while adding nothing to already-dense hubs (e.g.
# KMT2D_Het_vs_WT's SP1/JUN/FOS, all already over the cap) -- avoiding the
# overload a flat candidates-per-TF cap would cause on the denser network.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR11d_CollecTRI_dense_with_candidates.R
# Requires: AR11b_CollecTRI_curated_network.R, AR11c_CollecTRI_candidate_uncurated_links.R
# Self-checkpointing: skips entirely if the last contrast's output PDF exists.

suppressPackageStartupMessages({ library(igraph); library(dplyr); library(readr); library(yaml) })
source("scripts/AR11_network_plot_helpers.R")

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR11d] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)
CONTRASTS <- sapply(cfg$contrasts, function(ct) ct$name)
MAX_TOTAL_PER_TF <- 8
# Per-TF cap on GREAT peak-gene link confidence, applied AFTER
# MAX_TOTAL_PER_TF (2026-08-19, per user request). GREAT itself has no
# numeric "confidence" for a peak-gene assignment -- what it DOES give is
# distanceToTSS, already computed per (TF,gene) edge in
# AR06_direct_site_edges_{contrast}.csv (annotation_class buckets it into
# Promoter <=2kb / Proximal <=20kb / Distal <=100kb / Intergenic). An
# earlier attempt at a hard Promoter/Proximal-only FILTER (tmp scratch, not
# committed) worked well for KMT2D (217->71 edges) but nearly destroyed
# KDM6A (121->15 edges, 1 TF survived a >=3-target follow-up filter) --
# most of this pipeline's edges are distal-supported, so an absolute
# distance threshold isn't contrast-robust. A RANK + CAP is: within each
# TF's (already MAX_TOTAL_PER_TF-capped) edge set, sort by ascending
# |distanceToTSS| (most TSS-proximal = most confident first) and keep only
# the top PROXIMITY_CAP_PER_TF -- this can never zero out a sparse
# contrast (it keeps a TF's N most-confident links whatever N looks like
# for that TF), so the same constant is meant to be reused across
# contrasts without needing per-contrast tuning. Inf = off (no additional
# trimming beyond MAX_TOTAL_PER_TF) -- the default here, so this change
# alone does not alter any existing figure's edge set; lower it to actually
# apply the proximity cap.
PROXIMITY_CAP_PER_TF <- Inf
# Symmetric cap on the GENE side (2026-08-19, per user request): when a
# gene already has more than GENE_CAP_PER_TARGET incoming TF links, keep
# only its GENE_CAP_PER_TARGET most confident (closest-to-TSS) ones and
# drop the rest -- same distanceToTSS confidence metric as
# PROXIMITY_CAP_PER_TF above, just ranked within each TARGET gene instead
# of within each TF. Unlike a gene in-degree FILTER (drop genes with <=N
# links entirely -- tried as a scratch experiment; at >4 links it stranded
# 19/29 KDM6A TFs and 23/41 KMT2D TFs whose only targets didn't clear the
# bar), this only trims the excess edges on already-oversubscribed genes; a
# gene with few links, or a TF whose targets are all low-degree, is
# untouched -- can't zero out a TF the way the in-degree filter did.
# Inf = off (default, no behaviour change).
# Tried cap=4, then cap=3 (2026-08-19) -- reverted back to off (Inf) at
# user's request: it's a GREAT peak->gene ASSIGNMENT-confidence ranking
# (is this differentially-bound site correctly attributed to this gene, vs.
# a nearer neighbour), NOT a regulatory-strength ranking -- distal/
# enhancer-mediated regulation is common and often functionally MORE
# important than promoter-proximal binding in development, so this cap
# would have systematically favoured generic promoter-proximal regulators
# (e.g. Sp/KLF GC-box binders) over real distal drivers. Machinery (the
# gene_proximity_rank column, always exported) kept available/adaptable in
# case a properly-caveated use for it comes up later; just switched off.
GENE_CAP_PER_TARGET <- Inf

PLOTS_DIR  <- "results/integration"
TABLES_DIR <- "results/tables"
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
# Bump VERSION on every aesthetic tweak from here on (2026-08-18, per user
# request): each iteration writes a NEW file instead of overwriting the
# previous one, so earlier versions stay around for comparison.
VERSION <- "v13"
DONE_MARKER <- file.path(PLOTS_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_%s.pdf", tail(CONTRASTS, 1), VERSION))

if (file.exists(DONE_MARKER)) {
  cat("[AR11d] Already complete (", DONE_MARKER, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

curated_path <- file.path(TABLES_DIR, "AR11_CollecTRI_orthogonally_supported_edges.csv")
candidates_path <- file.path(TABLES_DIR, "AR11_candidate_uncurated_links.csv")
if (!file.exists(curated_path)) stop("[AR11d] ERROR: ", curated_path, " not found -- run AR11b_CollecTRI_curated_network.R first.")
if (!file.exists(candidates_path)) stop("[AR11d] ERROR: ", candidates_path, " not found -- run AR11c_CollecTRI_candidate_uncurated_links.R first.")

UP_COL <- "#D55E00"; DOWN_COL <- "#0072B2"; TF_NEUTRAL_COL <- "gray90"
ACT_COL <- "#009E73"; REP_COL <- "#CC79A7"; CAND_COL <- "gray55"
LFC_PAL <- colorRampPalette(c(DOWN_COL, "white", UP_COL))(201)
# Fixed colour-scale cap, not the per-figure max -- see AR11b for rationale.
LFC_CAP <- 4
lfc_to_col <- function(lfc, lfc_range) {
  idx <- round((lfc / lfc_range + 1) / 2 * 200) + 1
  LFC_PAL[pmin(pmax(idx, 1), 201)]
}

for (CT in CONTRASTS) {
  curated <- read_csv(curated_path, show_col_types = FALSE) %>%
    filter(contrast == CT, agrees_with_diffbind_call) %>%
    mutate(edge_class = "curated")

  existing_tfs <- unique(curated$from)
  n_curated_per_tf <- curated %>% count(from, name = "n_curated")
  candidates <- read_csv(candidates_path, show_col_types = FALSE) %>%
    filter(contrast == CT, from %in% existing_tfs) %>%
    left_join(n_curated_per_tf, by = "from") %>%
    mutate(slots_left = pmax(MAX_TOTAL_PER_TF - n_curated, 0)) %>%
    filter(slots_left > 0) %>%
    group_by(from) %>% arrange(desc(abs(site_log2fc_corrected)), .by_group = TRUE) %>%
    mutate(rk = row_number()) %>% filter(rk <= slots_left) %>% ungroup() %>%
    mutate(curated_direction = regulatory_mode, edge_class = "candidate", agrees_with_diffbind_call = NA,
           curation_effort = NA_real_, curated_weight = NA_real_, diffbind_dir = NA_character_, deg_dir = NA_character_)

  common_cols <- intersect(colnames(curated), colnames(candidates))
  edges_ct <- bind_rows(curated %>% select(all_of(common_cols)), candidates %>% select(all_of(common_cols)))
  cat(sprintf("[AR11d] %s: %d curated + %d candidate (fill to %d/TF total, %d TFs topped up) = %d total edges\n",
              CT, nrow(curated), nrow(candidates), MAX_TOTAL_PER_TF, n_distinct(candidates$from), nrow(edges_ct)))

  # GREAT peak-gene link confidence (distanceToTSS) -- join in regardless of
  # whether PROXIMITY_CAP_PER_TF is active, so it's always available as
  # exported metadata; annotation_class/distanceToTSS covers every edge here
  # (verified: 100% join match for both contrasts, since curated/candidate
  # both trace back to the same AR06 direct-site High/Medium-tier table).
  site_info_path <- file.path(TABLES_DIR, sprintf("AR06_direct_site_edges_%s.csv", CT))
  site_info <- read_csv(site_info_path, show_col_types = FALSE) %>%
    dplyr::select(TF, SYMBOL, annotation_class, distanceToTSS)
  edges_ct <- edges_ct %>%
    left_join(site_info, by = c("from" = "TF", "to" = "SYMBOL")) %>%
    group_by(from) %>% arrange(abs(distanceToTSS), .by_group = TRUE) %>%
    mutate(tf_proximity_rank = row_number()) %>% ungroup()
  if (is.finite(PROXIMITY_CAP_PER_TF)) {
    n_before <- nrow(edges_ct)
    edges_ct <- edges_ct %>% filter(tf_proximity_rank <= PROXIMITY_CAP_PER_TF)
    cat(sprintf("[AR11d] %s: proximity cap (top %d most TSS-proximal per TF) -- %d -> %d edges\n",
                CT, PROXIMITY_CAP_PER_TF, n_before, nrow(edges_ct)))
  }

  # Symmetric cap on the GENE side (see GENE_CAP_PER_TARGET definition near
  # the top of the script for the full rationale): when a gene already has
  # more than GENE_CAP_PER_TARGET incoming TF links, keep only its most
  # confident (closest-to-TSS) ones.
  edges_ct <- edges_ct %>%
    group_by(to) %>% arrange(abs(distanceToTSS), .by_group = TRUE) %>%
    mutate(gene_proximity_rank = row_number()) %>% ungroup()
  if (is.finite(GENE_CAP_PER_TARGET)) {
    n_before <- nrow(edges_ct)
    edges_ct <- edges_ct %>% filter(gene_proximity_rank <= GENE_CAP_PER_TARGET)
    cat(sprintf("[AR11d] %s: gene-side cap (top %d most TSS-proximal TFs per gene) -- %d -> %d edges\n",
                CT, GENE_CAP_PER_TARGET, n_before, nrow(edges_ct)))
  }

  all_nodes <- union(edges_ct$from, edges_ct$to)
  tf_panel_ct <- unique(edges_ct$from)
  tf_diffbind <- edges_ct %>% group_by(from) %>%
    summarise(tf_binding_change = mean(tf_binding_change), .groups = "drop") %>% rename(name = from)
  gene_lfc <- edges_ct %>% distinct(to, gene_log2FoldChange) %>% rename(name = to)
  vertices_df <- tibble(name = all_nodes) %>% left_join(tf_diffbind, by = "name") %>% left_join(gene_lfc, by = "name")

  g <- graph_from_data_frame(edges_ct, directed = TRUE, vertices = vertices_df)
  V(g)$is_TF <- V(g)$name %in% tf_panel_ct
  lfc_range <- LFC_CAP

  # Machine-readable export of exactly what's drawn (post per-TF capping),
  # not the full upstream tables -- one edges + one nodes TSV per contrast.
  write_tsv(edges_ct, file.path(TABLES_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_edges.tsv", CT)))
  nodes_export <- tibble(name = V(g)$name, is_TF = V(g)$is_TF,
                          diffbind_direction = ifelse(V(g)$is_TF, ifelse(V(g)$tf_binding_change >= 0, "UP", "DOWN"), NA_character_),
                          tf_binding_change = V(g)$tf_binding_change, gene_log2FoldChange = V(g)$gene_log2FoldChange)
  write_tsv(nodes_export, file.path(TABLES_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_nodes.tsv", CT)))
  lbl_cex <- 0.8
  vcol <- ifelse(V(g)$is_TF, TF_NEUTRAL_COL, lfc_to_col(V(g)$gene_log2FoldChange, lfc_range))
  vshape <- "circle"
  # TF nodes are ONE uniform size (not degree-scaled) -- sized to comfortably
  # fit this graph's longest TF name at this cex (bold, *1.08 like
  # tf_label_vsize used to apply), so every TF label sits INSIDE its own
  # circle instead of being repelled outside it. Gene/DEG nodes stay
  # uniform too, sized as a fixed fraction of the TF circle (gene_ratio =
  # 0.47, matching the geneRadius:tfRadius = 8:17 ratio the user landed on
  # via the interactive tuner at tmp/AR11d_param_sweep/regulon_tuner.html)
  # -- labels placed OUTSIDE via repel_labels() since a dot this small can't
  # contain its own gene name.
  tf_names <- V(g)$name[V(g)$is_TF]
  tf_vsize_uniform <- max(gene_label_vsize(tf_names, lbl_cex * 1.08))
  gene_vsize <- tf_vsize_uniform * 0.32
  vsize <- ifelse(V(g)$is_TF, tf_vsize_uniform, gene_vsize)
  vframe <- ifelse(V(g)$is_TF, "black", NA)
  vframe_width <- ifelse(V(g)$is_TF, 2.4, 1)

  is_candidate <- E(g)$edge_class == "candidate"
  # Candidate edges get curated_direction from regulatory_mode (site-level
  # Activating/Repressive, see the mutate() building `candidates` above) --
  # colour them by that direction too, same as curated edges, instead of a
  # flat CAND_COL that threw the direction information away. Dashing (elty)
  # remains the only "not yet literature-confirmed" cue.
  ecol_base <- ifelse(E(g)$curated_direction == "Activating", ACT_COL, REP_COL)
  ealpha <- ifelse(is_candidate, 0.6, 0.9)
  ewidth <- ifelse(is_candidate, 1.6, 1.3)
  ecol <- mapply(function(c, a) adjustcolor(c, alpha.f = a), ecol_base, ealpha)
  elty <- ifelse(is_candidate, 2, 1)

  # Plain layout_with_fr(), no community-weighting: that experiment (also
  # 2026-08-18, briefly) pulled communities apart into visually distinct
  # blobs but broke the hub-centred organisation this network reads best
  # with -- AP-1 (JUN/JUND/FOS/FOSL1/JUNB) really is the shared hub of most
  # of this graph, and plain FR already places it centrally with everything
  # else arranged by how strongly it's pulled toward that hub vs. others.
  # declutter_extra=1.6 (not the earlier 1.8 or 3): once TF nodes became
  # uniformly (and larger) sized, the SAME declutter_extra as before gave
  # them much more collision clearance than intended and blew the whole
  # layout apart -- TFs got pushed to the rim with all gene nodes collapsing
  # into the middle. Re-validated with a small-multiple param sweep
  # (tmp/AR11d_param_sweep/sheet5_uniform_tf.png) against the NEW uniform
  # TF/gene sizing above: 1.6 was the largest value that still kept the
  # graph filled and organized around the AP-1 hub; 2.2 already showed the
  # same rim/collapse failure re-emerging. The layout is then NOT
  # force-recompressed back into a fixed +-1 box afterward (that
  # recompression is exactly what caused this file's original "collapsed"
  # complaint): `lim` grows only as much as this extra padding actually
  # needs, while the PDF canvas itself stays fixed at 16x14in.
  set.seed(1)
  lay <- layout_with_fr(g, niter = 15000)
  lay <- norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  lay <- declutter_nodes(lay, declutter_radius_pad(vsize, extra = 2.1))
  rng <- max(1, max(abs(lay)))
  lim <- c(-(rng + 0.1), rng + 0.1)
  # Extra empty margin below the lowest possible node position, reserved
  # only for the bottom legends -- "bottomleft"/"bottomright" anchor at
  # ylim's actual lower edge, and since `lay` can place any node anywhere
  # up to that same edge, a node occasionally landed exactly where the
  # legend text was (e.g. TCF3 under the gradient-legend title). Nodes
  # never extend past `lim` regardless of what ylim we draw with, so
  # widening ylim's bottom alone creates real clear space for the legend
  # without moving anything else.
  ylim_plot <- c(lim[1] - 0.4, lim[2])
  is_tf_v <- V(g)$is_TF

  fpath <- file.path(PLOTS_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_%s.pdf", CT, VERSION))
  pdf(fpath, width = 16, height = 14)
  plot(g, layout = lay, rescale = FALSE, xlim = lim, ylim = ylim_plot,
       vertex.color = vcol, vertex.frame.color = vframe, vertex.frame.width = vframe_width, vertex.shape = vshape,
       vertex.size = vsize, vertex.label = NA,
       edge.color = ecol, edge.width = ewidth, edge.lty = elty, edge.arrow.mode = 0, edge.curved = 0, main = "")
  # Fixed arrow/T-bar size for ALL edges, matched to the interactive tuner's
  # markerSize choice -- draw_edge_direction_markers() still uses the REAL
  # `vsize` for POSITIONING each marker at its target's actual border
  # (critical: passing a fake uniform vsize for that too, as an earlier
  # version did, put arrows on big TF nodes' borders in the wrong place,
  # landing inside the circle), but sizes the glyph itself from the fixed
  # `marker_vsize = gene_vsize` so an edge into a big TF and one into a
  # small gene dot both get the same-size arrowhead/T-bar.
  draw_edge_direction_markers(g, lay, vsize,
                               is_activating = E(g)$curated_direction == "Activating",
                               col = ecol, lwd = ewidth, arrow_len_f = 0.9, arrow_width_f = 0.25,
                               marker_vsize = gene_vsize)
  text(lay[is_tf_v, 1], lay[is_tf_v, 2], V(g)$name[is_tf_v], cex = lbl_cex, col = "black", family = "sans", font = 2)
  # Gene labels (v5): offset just ABOVE the dot (dot radius + small gap),
  # like the interactive tuner artifact, not centred through the middle of
  # it (v4's attempt at "on top of" the dot, which read as text crossing
  # through the dot rather than sitting above it) and not repelled away
  # with a leader line either.
  gene_r <- gene_vsize / 200
  text(lay[!is_tf_v, 1], lay[!is_tf_v, 2] + gene_r + 0.012, V(g)$name[!is_tf_v], cex = lbl_cex, col = "black", family = "sans")

  title(main = sprintf("Orthogonally supported edges + candidate uncurated links -- %s", CT), cex.main = 1.8, line = 1.8)
  mtext(sprintf("Solid = curated (CollecTRI). Dashed = high-confidence AR06 edge, TF is CollecTRI-covered but this pair is not curated (candidate, not confirmed) -- colour still shows the inferred direction. Each TF filled to max %d total edges shown.", MAX_TOTAL_PER_TF),
        side = 3, line = 0.6, cex = 0.85, col = "gray30")

  # Position the gradient legend to the RIGHT of the Nodes legend's ACTUAL
  # measured rectangle (legend() returns it invisibly), at the same
  # vertical band -- instead of a hardcoded coordinate tuned for the old,
  # shorter legend text (this file's Nodes legend text is longer/bigger
  # now, so a fixed position collided with it), and instead of stacking it
  # BELOW, which pushed it far enough down to fall outside the page's
  # bottom margin and get silently clipped.
  nodes_leg <- legend("bottomleft", bty = "n", cex = 1.0, title = "Nodes", title.adj = 0,
         legend = c("TF (gray, outlined, size = out-degree)", "Gene/DEG (colour = own log2FC, uniform size)"),
         pt.bg = c(TF_NEUTRAL_COL, "white"), pch = 21, pt.lwd = c(2, 1), col = "black", pt.cex = c(2, 1.2))
  gx0 <- nodes_leg$rect$left + nodes_leg$rect$w + 0.05
  gx1 <- gx0 + 0.5
  gy1 <- nodes_leg$rect$top - 0.02
  gy0 <- gy1 - 0.08
  draw_gradient_legend(gx0, gx1, gy0, gy1, sprintf("<= -%.0f", lfc_range), sprintf(">= +%.0f", lfc_range),
                        "Gene nodes: colour = own RNA log2FC (capped)", pal = LFC_PAL, cex_title = 1.0, cex_lab = 0.9)
  linestyle_leg <- legend("bottomright", bty = "n", cex = 1.0, title = "Edges: line style", title.adj = 0,
         legend = c("Curated (CollecTRI, solid)", "Candidate, not in CollecTRI (dashed, not confirmed)"),
         col = "gray20", lty = c(1, 2), lwd = c(1.3, 1.6))
  # Explicit direction legend (v4, per user request) drawn ABOVE the line-
  # style legend, using the exact same triangle/T-bar drawing code as the
  # markers on the graph (not an approximated pch symbol), so what's shown
  # here matches what's actually drawn -- positioned from the line-style
  # legend's own measured rect, same pattern as the Nodes/gradient stack.
  dir_leg <- legend(x = linestyle_leg$rect$left, y = linestyle_leg$rect$top + linestyle_leg$rect$h + 0.12,
         bty = "n", cex = 1.0, title = "Edges: direction (fixed marker size)", title.adj = 0,
         legend = c("Activating", "Repressive"), lty = 0, pch = NA)
  sym_x <- dir_leg$rect$left + 0.06
  r_fixed <- gene_vsize / 200
  arrow_len <- r_fixed * 1.3; arrow_width <- r_fixed * 0.7; bar_half <- r_fixed * 0.65
  y1 <- dir_leg$text$y[1]; y2 <- dir_leg$text$y[2]
  segments(sym_x - 0.065, y1, sym_x, y1, col = ACT_COL, lwd = 1.8)
  polygon(c(sym_x + arrow_len, sym_x, sym_x), c(y1, y1 + arrow_width, y1 - arrow_width), col = ACT_COL, border = NA)
  segments(sym_x - 0.065, y2, sym_x, y2, col = REP_COL, lwd = 1.8)
  segments(sym_x, y2 - bar_half, sym_x, y2 + bar_half, col = REP_COL, lwd = 1.8 * 1.5)
  mtext("Not validated TF-target relationships in neuroectoderm specifically. CollecTRI support reflects a documented relationship in human cells generally, not tissue-specific validation. Candidate edges are hypothesis-generating only.",
        side = 1, line = 3.5, cex = 0.65, col = "gray40", adj = 0)
  dev.off()
  cat(sprintf("[DONE] %s -- %d nodes, %d edges (%d curated, %d candidate)\n", fpath, vcount(g), ecount(g), nrow(curated), nrow(candidates)))
}

cat("\n[DONE] AR11d_CollecTRI_dense_with_candidates complete\n")
