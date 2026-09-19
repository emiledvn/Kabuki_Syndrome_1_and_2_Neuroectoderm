#!/usr/bin/env Rscript
# AR11b_CollecTRI_curated_network.R -- intersects AR06's chromatin-
# accessibility-inferred regulatory links (direct-site network) with
# CollecTRI (literature-curated regulatory interactions, AR11) to produce
# the "orthogonally supported edges" network: only TF-gene pairs with
# independent literature evidence are kept, sidestepping the ambiguity
# raised by cases like SOX10->WNT3A or RARB->CYP26A1 (genome-wide DiffBind
# trend vs. locus-specific footprint-depth trend disagreeing in sign).
# CollecTRI covers 124/148 (83.8%) of the full TOBIAS TF panel (union of
# both contrasts) -- above a 60% coverage threshold, so no DoRothEA fallback
# is used (confirmed empirically to add zero additional TFs anyway, since
# DoRothEA is itself one of CollecTRI's input sources).
#
# TERMINOLOGY (enforced in every label/legend/comment below):
#   TOBIAS footprint score      -> "inferred occupancy" / "footprint depth",
#                                   never "binding"
#   GREAT peak->gene assignment -> flagged explicitly as proximity-based
#   CollecTRI edges              -> "literature-curated regulatory interactions"
#   AR06 full network            -> "chromatin-accessibility-inferred regulatory links"
#   this curated-intersection    -> "orthogonally supported edges"
# Disclaimers (stated once, in the figure caption, not repeated per-element):
#   - not validated TF-target relationships in neuroectoderm specifically
#   - CollecTRI support = documented relationship in human cells generally,
#     not a tissue-specific validation
#   - GREAT proximity-based assignment may attribute a chromatin-accessibility
#     change to the wrong nearby gene
#   - hypothesis-generating, not mechanistic proof of regulation
#
# NOTE on intersection semantics (added when restoring this script into the
# publication repo, 2026-09-19): "orthogonally supported" below means the
# (TF, gene) PAIR exists in CollecTRI at all -- the inner_join has no sign
# check of its own. `agrees_with_diffbind_call` is a SEPARATE, additional
# flag comparing CollecTRI's curated direction against this TF's genome-wide
# DiffBind/BINDetect binding-change direction paired with the target gene's
# own DEG direction (not the site-level regulatory_mode AR06 used to build
# the edge). AR11d_CollecTRI_dense_with_candidates.R (and therefore the
# force-simulation figure downstream of it) draws an edge solid ("curated")
# only when BOTH hold: pair-presence in CollecTRI, AND
# agrees_with_diffbind_call == TRUE. Do not describe this as a single
# "consistent sign" check in the methods text -- it is two distinct checks.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/integration/AR11b_CollecTRI_curated_network.R
# Requires: AR06b_Direct_Site_Edges.R (AR06_graph_direct_site_<contrast>.rds), AR11_CollecTRI_fetch_regulons.R
# Self-checkpointing: skips entirely if results/tables/AR11_CollecTRI_orthogonally_supported_edges.csv exists.

suppressPackageStartupMessages({ library(igraph); library(dplyr); library(readr); library(yaml) })
source("scripts/integration/AR11_network_plot_helpers.R")

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR11b] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)
CONTRASTS <- sapply(cfg$contrasts, function(ct) ct$name)

RDS_DIR    <- "results/RDS"
PLOTS_DIR  <- "results/integration"
TABLES_DIR <- "results/tables"
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
CSV_MARKER <- file.path(TABLES_DIR, "AR11_CollecTRI_orthogonally_supported_edges.csv")
# Aesthetic pass (2026-08-03): checkpoint on the restyled figure, not the
# (unchanged, deterministic) edges CSV -- otherwise this script would skip
# before ever reaching the new plotting code, since CSV_MARKER already
# exists from the original run. Old "_v1" PDFs are left untouched; this run
# only adds new "_v2"-suffixed figures.
DONE_MARKER <- file.path(PLOTS_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_v3.pdf", CONTRASTS[1]))

if (file.exists(DONE_MARKER)) {
  cat("[AR11b] Already complete (", DONE_MARKER, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

collectri_path <- file.path(RDS_DIR, "AR11_collectri_signed.rds")
if (!file.exists(collectri_path)) stop("[AR11b] ERROR: ", collectri_path, " not found -- run AR11_CollecTRI_fetch_regulons.R first.")
graph_paths_exist <- all(sapply(CONTRASTS, function(ct) file.exists(file.path(RDS_DIR, sprintf("AR06_graph_direct_site_%s.rds", ct)))))
if (!graph_paths_exist) stop("[AR11b] ERROR: AR06_graph_direct_site_*.rds missing -- run AR06b_Direct_Site_Edges.R first.")

UP_COL <- "#D55E00"; DOWN_COL <- "#0072B2"; TF_NEUTRAL_COL <- "gray90"
ACT_COL <- "#009E73"; REP_COL <- "#CC79A7"
LFC_PAL <- colorRampPalette(c(DOWN_COL, "white", UP_COL))(201)
# Fixed colour-scale cap, not the per-figure max: >95% of DEG nodes across
# this project's networks have |log2FC| < 4, but a handful of outliers
# (e.g. some genes reach |log2FC| > 12) would otherwise stretch the whole
# gradient and crush the typical range into a narrow near-white band.
# Outliers beyond the cap saturate to the endpoint colour instead.
LFC_CAP <- 4
lfc_to_col <- function(lfc, lfc_range) {
  idx <- round((lfc / lfc_range + 1) / 2 * 200) + 1
  LFC_PAL[pmin(pmax(idx, 1), 201)]
}

## ---- 1. Intersect each contrast's chromatin-accessibility-inferred
## regulatory links with the literature-curated regulatory interactions set.
collectri <- readRDS(collectri_path) %>%
  rename(curated_weight = weight)  # AR06 edges already have their own `weight` column

curated_edges <- list()
for (ct in CONTRASTS) {
  g <- readRDS(file.path(RDS_DIR, sprintf("AR06_graph_direct_site_%s.rds", ct)))
  el <- igraph::as_data_frame(g, what = "edges")
  joined <- el %>%
    inner_join(collectri, by = c("from" = "source", "to" = "target")) %>%
    mutate(curated_direction = ifelse(curated_weight > 0, "Activating", "Repressive"),
           # DiffBind call: the TF's own genome-wide BINDetect differential
           # inferred-occupancy call (tf_binding_change -- constant per TF,
           # NOT the site-specific value) vs. the target gene's own DEG
           # direction. Same-direction movement agrees with a curated
           # activator; opposite-direction movement agrees with a curated
           # repressor (de-repression logic). Using the TF-level DiffBind
           # call (not a locus-specific one) keeps this consistent with how
           # the TF node itself is styled below -- one statistic, not two
           # conflated ones (see the SOX10->WNT3A discussion).
           diffbind_dir = ifelse(tf_binding_change > 0, "UP", "DOWN"),
           deg_dir = ifelse(gene_log2FoldChange > 0, "UP", "DOWN"),
           agrees_with_diffbind_call = (curated_direction == "Activating" & diffbind_dir == deg_dir) |
                                        (curated_direction == "Repressive" & diffbind_dir != deg_dir))
  cat(sprintf("[AR11b] %s: %d chromatin-accessibility-inferred edges -> %d orthogonally supported (in CollecTRI) -- %.1f%% agree with DiffBind call\n",
              ct, nrow(el), nrow(joined), 100 * mean(joined$agrees_with_diffbind_call)))
  curated_edges[[ct]] <- joined
}
curated_all <- bind_rows(curated_edges)
write_csv(curated_all, CSV_MARKER)

## ---- 2. Cross-contrast annotation only (no merging): for each contrast's
## edges, flag whether the same (TF, gene, curated_direction) triple is also
## orthogonally supported in the *other* contrast -- shown as a bold border
## within each contrast's own, separate figure, not combined into one graph.
other_edge_key <- function(ct) {
  curated_all %>% filter(contrast != ct) %>%
    transmute(key = paste(from, to, curated_direction)) %>% pull(key) %>% unique()
}
for (ct in CONTRASTS) {
  curated_edges[[ct]] <- curated_edges[[ct]] %>%
    mutate(also_in_other_contrast = paste(from, to, curated_direction) %in% other_edge_key(ct))
}

## ---- 4. Build and render one curated-intersection network per contrast.
## TF nodes: all circles now (see AR11_network_plot_helpers.R header for why
## the old shape=DiffBind-direction encoding was dropped) -- a thicker
## outline is the only TF-vs-gene visual distinction. Gene (DEG) nodes:
## filled on a continuous log2FC colour scale, sized to fit their own label,
## no outline. Edge colour = curated direction; edge endpoint marker (arrow/
## T-bar) = Activating/Repressive, redundant with colour on purpose.
render_curated_network <- function(ct, agreeing_only = FALSE) {
  edges_ct <- curated_edges[[ct]]
  if (agreeing_only) edges_ct <- edges_ct %>% filter(agrees_with_diffbind_call)
  all_nodes <- union(edges_ct$from, edges_ct$to)
  tf_panel_ct <- unique(edges_ct$from)

  # A TF gene symbol can map to >1 motif/motif-cluster (e.g. RXRB matches
  # both C_NR2C2 and C_RXRG JASPAR motifs) with distinct per-motif DiffBind
  # scores -- average them to one representative value per TF symbol.
  tf_diffbind <- edges_ct %>% group_by(from) %>%
    summarise(tf_binding_change = mean(tf_binding_change), .groups = "drop") %>% rename(name = from)
  gene_lfc <- edges_ct %>% distinct(to, gene_log2FoldChange) %>% rename(name = to)
  vertices_df <- tibble(name = all_nodes) %>%
    left_join(tf_diffbind, by = "name") %>%
    left_join(gene_lfc, by = "name")

  g_curated <- graph_from_data_frame(edges_ct, directed = TRUE, vertices = vertices_df)
  V(g_curated)$is_TF <- V(g_curated)$name %in% tf_panel_ct

  lfc_range <- LFC_CAP
  lbl_cex <- 1.05
  vcol <- ifelse(V(g_curated)$is_TF, TF_NEUTRAL_COL, lfc_to_col(V(g_curated)$gene_log2FoldChange, lfc_range))
  vshape <- "circle"
  tf_degree_vsize <- pmin(5 + degree(g_curated, mode = "out"), 15)
  vsize <- ifelse(V(g_curated)$is_TF,
                   tf_label_vsize(V(g_curated)$name, lbl_cex, tf_degree_vsize),
                   gene_label_vsize(V(g_curated)$name, lbl_cex))
  vframe <- ifelse(V(g_curated)$is_TF, "black", NA)  # no outline on gene dots
  vframe_width <- ifelse(V(g_curated)$is_TF, 3.6, 1)

  ecol_base <- ifelse(E(g_curated)$curated_direction == "Activating", ACT_COL, REP_COL)
  ealpha <- ifelse(E(g_curated)$also_in_other_contrast, 0.95, 0.7)
  is_disagree <- if (agreeing_only) rep(FALSE, ecount(g_curated)) else !E(g_curated)$agrees_with_diffbind_call
  # Base width up across the board (fix #1); dashed/disagreeing edges get a
  # much bigger extra bump on top (fix #2) -- R's dashed lty pattern has a
  # fixed dash/gap length regardless of lwd, so a thin dashed line reads as
  # much fainter than a solid line of the same width; matching that
  # perceived weight needs a substantial width bump, not a token one.
  ewidth <- ifelse(E(g_curated)$also_in_other_contrast, 3.6, 2.2) + ifelse(is_disagree, 1.6, 0)
  ecol <- mapply(function(c, a) adjustcolor(c, alpha.f = a), ecol_base, ealpha)
  elty <- if (agreeing_only) rep(1, ecount(g_curated)) else ifelse(E(g_curated)$agrees_with_diffbind_call, 1, 2)

  set.seed(1)
  lay <- layout_with_fr(g_curated, niter = 5000)
  lay <- norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  lay <- declutter_nodes(lay, declutter_radius_pad(vsize))
  lay <- norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  is_tf_v <- V(g_curated)$is_TF
  tf_rep_pos <- repel_labels(lay[is_tf_v, 1], lay[is_tf_v, 2], V(g_curated)$name[is_tf_v], lbl_cex)

  suffix <- if (agreeing_only) "_agreeingOnly" else ""
  fpath <- file.path(PLOTS_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s%s_v3.pdf", ct, suffix))
  pdf(fpath, width = 16, height = 14)
  plot(g_curated, layout = lay, rescale = FALSE, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
       vertex.color = vcol, vertex.frame.color = vframe, vertex.frame.width = vframe_width, vertex.shape = vshape,
       vertex.size = vsize, vertex.label = NA,
       edge.color = ecol, edge.width = ewidth, edge.lty = elty, edge.arrow.mode = 0, edge.curved = 0.05,
       main = "")
  draw_edge_direction_markers(g_curated, lay, vsize,
                               is_activating = E(g_curated)$curated_direction == "Activating",
                               col = ecol, lwd = ewidth)
  moved <- tf_rep_pos$moved
  if (any(moved)) segments(lay[is_tf_v, 1][moved], lay[is_tf_v, 2][moved], tf_rep_pos$x[moved], tf_rep_pos$y[moved],
                            col = adjustcolor("gray40", 0.4), lwd = 0.4)
  text(tf_rep_pos$x, tf_rep_pos$y, V(g_curated)$name[is_tf_v], cex = lbl_cex, col = "black", family = "sans", font = 2)
  text(lay[!is_tf_v, 1], lay[!is_tf_v, 2], V(g_curated)$name[!is_tf_v], cex = lbl_cex, col = "black", family = "sans")

  title(main = sprintf("Orthogonally supported edges -- %s%s", ct, if (agreeing_only) " (agreeing only)" else ""),
        cex.main = 1.4, line = 1.8)
  mtext(if (agreeing_only)
          "chromatin-accessibility-inferred regulatory links intersected with literature-curated regulatory interactions (CollecTRI), restricted to edges where the DiffBind call agrees with CollecTRI's curated direction"
        else
          "chromatin-accessibility-inferred regulatory links intersected with literature-curated regulatory interactions (CollecTRI)",
        side = 3, line = 0.6, cex = 0.9, col = "gray30")

  legend("bottomleft", bty = "n", cex = 1.05, title = "Nodes", title.adj = 0,
         legend = c("TF (thick outline)", "Gene/DEG (colour = own log2FC)"),
         pt.bg = c(TF_NEUTRAL_COL, "white"), pch = 21, pt.lwd = c(3, 1), col = "black", pt.cex = 1.8)
  draw_gradient_legend(-1.05, -0.48, -1.32, -1.15, sprintf("<= -%.0f", lfc_range), sprintf(">= +%.0f", lfc_range),
                        "Gene nodes: colour = own RNA log2FC (no outline, capped)", pal = LFC_PAL,
                        cex_title = 1.0, cex_lab = 0.9)

  edge_legend <- c("Activating (arrowhead)", "Repressive (T-bar)", "Also supported in the other contrast (bold)")
  edge_col_l <- c(ACT_COL, REP_COL, "gray20"); edge_lty_l <- c(1, 1, 1); edge_lwd_l <- c(2.8, 2.8, 3.6)
  if (!agreeing_only) {
    edge_legend <- c(edge_legend, "Solid = agrees with DiffBind call", "Dashed = disagrees (colour unchanged)")
    edge_col_l <- c(edge_col_l, "gray20", "gray20"); edge_lty_l <- c(edge_lty_l, 1, 2); edge_lwd_l <- c(edge_lwd_l, 1.2, 2.8)
  }
  legend("bottomright", bty = "n", cex = 1.05, title = "Edges: colour = CollecTRI-curated direction", title.adj = 0,
         legend = edge_legend, col = edge_col_l, lty = edge_lty_l, lwd = edge_lwd_l)

  mtext("Not validated TF-target relationships in neuroectoderm specifically. CollecTRI support reflects a documented relationship in human cells generally, not tissue-specific validation.\nGREAT peak-to-gene assignment is proximity-based and may attribute a chromatin-accessibility change to the wrong nearby gene. Hypothesis-generating, not mechanistic proof of regulation.",
        side = 1, line = 4, cex = 0.68, col = "gray40", adj = 0)
  dev.off()
  cat(sprintf("[DONE] %s -- %d nodes, %d edges (%d also supported in the other contrast)\n",
              fpath, vcount(g_curated), ecount(g_curated), sum(E(g_curated)$also_in_other_contrast)))
}

for (ct in CONTRASTS) {
  render_curated_network(ct, agreeing_only = FALSE)
  render_curated_network(ct, agreeing_only = TRUE)
}

cat("\n[DONE] AR11b_CollecTRI_curated_network complete\n")
