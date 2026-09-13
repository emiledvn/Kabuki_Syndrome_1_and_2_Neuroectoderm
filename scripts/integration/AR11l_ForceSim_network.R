#!/usr/bin/env Rscript
# AR11l_ForceSim_network.R -- denseWithCandidates network (AR11d), laid out
# with a port of the interactive tuner's ACTUAL force simulation
# (tmp/AR11d_param_sweep/regulon_tuner.html tick(), lines 289-348:
# repulsion + springs + collision + centering, all applied together every
# tick) instead of layout_with_fr()+declutter_nodes() (a fundamentally
# different algorithm: one-shot static layout, then a SEPARATE collision
# pass after the fact). This is the root-cause fix for "it's very different
# from the previous interactive one" -- not more parameter tuning of the FR
# layout, an actually different layout algorithm.
#
# Parameter values (charge/linkDist/linkStr/collide/center/tfR/geneR) are
# the ones the user found via the interactive tool for KMT2D; reused as-is
# for KDM6A too since these are generic "what looks organic" physics
# constants, not contrast-specific.
#
# Colour scheme deliberately does NOT match the artifact 100%: the artifact
# is dark-themed (bg #10141a) with TF fill/stroke calibrated for a dark
# canvas (#d7dce2 fill / #f4f6f8 stroke -- a near-invisible outline against
# white). Every other AR11 figure is a white-background PDF for
# print/journal-figure conventions, so this keeps that: white background,
# gene log2FC gradient reuses the pipeline's existing LFC_PAL (not the
# artifact's cream-tinted version) for consistency with every other AR11
# figure. Only the LAYOUT algorithm and edge/marker geometry are ported
# faithfully; TF border colour is darkened for visibility on white.
#
# Env: ks_1_2_r. Run from repo root.
# Requires: AR11d_CollecTRI_dense_with_candidates.R (for its edges/nodes TSV exports)
# Bump VERSION on every change -- never overwrites a previous version's PDF.

suppressPackageStartupMessages({ library(igraph); library(dplyr); library(readr); library(yaml) })

PLOTS_DIR  <- "results/integration"
TABLES_DIR <- "results/tables"
VERSION <- "v14"
# Seed picked by the user from a 6-way visual sweep (see
# tmp/AR11d_param_sweep/sheet6_forceSim_seed_sweep_KMT2D.png) -- the
# simulation is deterministic given the same random initial placement, but
# a nonlinear multi-body system like this converges to a genuinely
# different final arrangement per seed, not just a rotation/reflection of
# the same one.
SEED <- 4
# GENE_CAP_PER_TARGET tried (v3: cap=4, v4: cap=3, filenames tagged
# "TSSproxCap") then reverted back to off at user's request (see AR11d for
# why) -- back to the plain filename, no cap tag, matching v1/v2.

CONFIG <- "config/pipeline_config.yaml"
cfg <- yaml::read_yaml(CONFIG)
CONTRASTS <- sapply(cfg$contrasts, function(ct) ct$name)

## ---- physics params (from the interactive tuner, KMT2D-derived) ----
CHARGE <- 330; LINK_DIST <- 82; LINK_STR <- 34; COLLIDE <- 270; CENTER <- 40
TF_R <- 17; GENE_R <- 8; MARKER <- 4
WIDTH <- 1200; HEIGHT <- 900
NTICKS <- 900

UP_COL <- "#D55E00"; DOWN_COL <- "#0072B2"; TF_FILL <- "#d7dce2"; TF_BORDER <- "#8a93a3"
ACT_COL <- "#009E73"; REP_COL <- "#CC79A7"
LFC_PAL <- colorRampPalette(c(DOWN_COL, "white", UP_COL))(201)
LFC_CAP <- 4
lfc_to_col <- function(lfc, r) LFC_PAL[pmin(pmax(round((lfc / r + 1) / 2 * 200) + 1, 1), 201)]

# Jaccard-scaled TF-TF repulsion (user idea, 2026-08-19): TFs with more
# similar target-gene sets (same Jaccard metric as AR11h/AR11i) repel each
# other less, so family-like clusters emerge organically from the actual
# data instead of requiring an explicit collapse step. Only applied to
# TF-TF pairs -- gene-gene and TF-gene repulsion untouched (scale=1).
# REPULSION_JACCARD_STRENGTH=0.95 with MIN_SCALE=0.08 (pushed stronger than
# the first trial, which used 0.85/0.15, per user request): a jaccard=1 TF
# pair's repulsion drops to 5% of baseline, a jaccard=0 pair is unaffected.
# Floor is NOT zero -- some residual repulsion is kept even for
# near-identical TFs so the collision pass (below) still has a pairwise
# force to resolve against instead of nodes drifting into a true fixed point.
# Switched back off (2026-08-19): user compared all saved versions and
# preferred v6 (plain uniform repulsion, no Jaccard scaling, no TF
# outline) over v8-v10's Jaccard-scaled-repulsion look -- kept the
# machinery here, adaptable, in case it's wanted again later.
USE_JACCARD_REPULSION <- FALSE
REPULSION_JACCARD_STRENGTH <- 0.95
MIN_SCALE <- 0.08

compute_tf_repulsion_scale <- function(g, is_tf_v) {
  nn <- vcount(g)
  scale_mat <- matrix(1, nn, nn)
  if (!USE_JACCARD_REPULSION) return(scale_mat)
  tf_idx <- which(is_tf_v)
  targets_by_tf <- lapply(tf_idx, function(i) igraph::neighbors(g, i, mode = "out"))
  names(targets_by_tf) <- tf_idx
  for (a in seq_along(tf_idx)) for (b in seq_along(tf_idx)) {
    if (a >= b) next
    i <- tf_idx[a]; j <- tf_idx[b]
    ti <- targets_by_tf[[as.character(i)]]; tj <- targets_by_tf[[as.character(j)]]
    u <- length(union(ti, tj))
    jac <- if (u == 0) 0 else length(intersect(ti, tj)) / u
    s <- max(1 - REPULSION_JACCARD_STRENGTH * jac, MIN_SCALE)
    scale_mat[i, j] <- s; scale_mat[j, i] <- s
  }
  scale_mat
}

# White halo behind label text (2026-08-19, per user request): edges are
# two different colours (green/pink) crossing at all angles, so no single
# edge colour choice keeps black text legible everywhere it might cross
# behind a label -- a white halo (draw the string offset in several
# directions in white first, then black on top) guarantees contrast
# against ANY background behind it, the standard fix for label-over-line
# legibility (no native text-outline in base R graphics).
halo_text <- function(x, y, labels, cex, font = 1, col = "black", halo_col = "white", halo_r = 0.0016 * WIDTH) {
  for (ang in seq(0, 2 * pi, length.out = 8)) {
    text(x + cos(ang) * halo_r, y + sin(ang) * halo_r, labels, cex = cex, font = font, col = halo_col)
  }
  text(x, y, labels, cex = cex, font = font, col = col)
}

# Gene-label decluttering (2026-08-19, per user request: bigger gene text
# must not start overlapping in the dense JUN/FOS core). Same overlap-
# repulsion idea as AR11_network_plot_helpers.R's repel_labels(), but that
# helper's width/height formulas assume the -1..1 igraph coordinate range
# used elsewhere in AR11 -- this script's custom force-sim runs in its own
# ~1200x900 pixel-unit space (see WIDTH/HEIGHT), so a version calibrated to
# that scale instead of reusing the helper as-is.
declutter_labels_px <- function(x, y, labels, cex, iter = 400) {
  n <- length(x)
  if (n <= 1) return(list(x = x, y = y))
  w <- nchar(labels) * cex * 12 + 6
  h <- cex * 22
  px <- x; py <- y
  for (it in seq_len(iter)) {
    moved_any <- FALSE
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i >= j) next
      dx <- px[i] - px[j]; dy <- py[i] - py[j]
      ox <- (w[i] + w[j]) / 2 - abs(dx); oy <- (h[i] + h[j]) / 2 - abs(dy)
      if (is.na(ox) || is.na(oy)) next
      if (ox > 0 && oy > 0) {
        moved_any <- TRUE
        if (ox < oy) { s <- ox * 0.5 + 0.3; sg <- if (dx >= 0) 1 else -1; px[i] <- px[i] + sg * s; px[j] <- px[j] - sg * s }
        else         { s <- oy * 0.5 + 0.3; sg <- if (dy >= 0) 1 else -1; py[i] <- py[i] + sg * s; py[j] <- py[j] - sg * s }
      }
    }
    if (!moved_any) break
  }
  list(x = px, y = py)
}

run_force_sim <- function(g, is_tf_v) {
  nn <- vcount(g)
  r <- ifelse(is_tf_v, TF_R, GENE_R)
  repulsion_scale <- compute_tf_repulsion_scale(g, is_tf_v)
  set.seed(SEED)
  x <- runif(nn, WIDTH * 0.3, WIDTH * 0.7)
  y <- runif(nn, HEIGHT * 0.3, HEIGHT * 0.7)
  vx <- rep(0, nn); vy <- rep(0, nn)
  el <- igraph::ends(g, igraph::E(g), names = FALSE)
  ei <- el[, 1]; ej <- el[, 2]
  lstr <- LINK_STR / 1000; cstr <- CENTER / 8000; pad <- COLLIDE / 100

  for (tick in seq_len(NTICKS)) {
    dx <- outer(x, x, "-"); dy <- outer(y, y, "-")
    d2 <- dx^2 + dy^2; diag(d2) <- Inf; d2 <- pmax(d2, 0.01); d <- sqrt(d2)
    force <- (CHARGE * 4 * repulsion_scale) / d2
    fx <- (dx / d) * force; fy <- (dy / d) * force
    vx <- vx + rowSums(fx); vy <- vy + rowSums(fy)

    dxl <- x[ej] - x[ei]; dyl <- y[ej] - y[ei]
    dl <- sqrt(dxl^2 + dyl^2); dl[dl < 0.01] <- 0.01
    diff <- (dl - LINK_DIST) * lstr
    fxl <- (dxl / dl) * diff; fyl <- (dyl / dl) * diff
    for (k in seq_along(ei)) {
      vx[ei[k]] <- vx[ei[k]] + fxl[k]; vy[ei[k]] <- vy[ei[k]] + fyl[k]
      vx[ej[k]] <- vx[ej[k]] - fxl[k]; vy[ej[k]] <- vy[ej[k]] - fyl[k]
    }

    vx <- vx + (WIDTH / 2 - x) * cstr
    vy <- vy + (HEIGHT / 2 - y) * cstr

    Rsum <- outer(r, r, "+"); minD <- Rsum * (1 + pad * 0.5) + 2
    dxc <- outer(x, x, "-"); dyc <- outer(y, y, "-")
    dc <- sqrt(dxc^2 + dyc^2); diag(dc) <- Inf
    overlap <- pmax((minD - dc) / 2, 0)
    ux <- dxc / dc; uy <- dyc / dc
    ux[!is.finite(ux)] <- 0; uy[!is.finite(uy)] <- 0
    x <- x + rowSums(ux * overlap); y <- y + rowSums(uy * overlap)

    vx <- vx * 0.82; vy <- vy * 0.82
    x <- x + vx * 0.9; y <- y + vy * 0.9
    x <- pmax(r + 4, pmin(WIDTH - r - 4, x)); y <- pmax(r + 4, pmin(HEIGHT - r - 4, y))
  }
  list(x = x, y = y, r = r)
}

for (CT in CONTRASTS) {
  edges_path <- file.path(TABLES_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_edges.tsv", CT))
  nodes_path <- file.path(TABLES_DIR, sprintf("AR11_CollecTRI_orthogonally_supported_network_%s_denseWithCandidates_nodes.tsv", CT))
  if (!file.exists(edges_path) || !file.exists(nodes_path)) stop("[AR11l] ERROR: run AR11d_CollecTRI_dense_with_candidates.R first.")
  e <- read_tsv(edges_path, show_col_types = FALSE)
  n <- read_tsv(nodes_path, show_col_types = FALSE)

  all_nodes <- union(e$from, e$to)
  vdf <- n %>% filter(name %in% all_nodes)
  g <- graph_from_data_frame(e, directed = TRUE, vertices = vdf)
  V(g)$is_TF <- V(g)$name %in% unique(e$from)
  is_tf_v <- V(g)$is_TF

  sim <- run_force_sim(g, is_tf_v)
  x <- sim$x; y <- sim$y; r <- sim$r
  el <- igraph::ends(g, igraph::E(g), names = FALSE)
  ei <- el[, 1]; ej <- el[, 2]

  fpath <- file.path(PLOTS_DIR, sprintf("AR11_CollecTRI_network_%s_denseWithCandidates_forceSim_%s.pdf", CT, VERSION))
  pdf(fpath, width = 16, height = 13)
  par(mar = c(2, 1, 4, 1), bg = "white")
  plot(NA, xlim = c(0, WIDTH), ylim = c(-HEIGHT * 0.22, HEIGHT), xlab = "", ylab = "", axes = FALSE, asp = 1)

  is_candidate <- E(g)$edge_class == "candidate"
  ecol_base <- ifelse(E(g)$curated_direction == "Activating", ACT_COL, REP_COL)
  # Fainter than before (was 0.55/0.85) so edges crossing near a label don't
  # compete with the text for attention.
  ecol <- mapply(function(c, a) adjustcolor(c, alpha.f = a), ecol_base, ifelse(is_candidate, 0.38, 0.6))
  for (k in seq_along(ei)) {
    a <- ei[k]; b <- ej[k]
    ddx <- x[b] - x[a]; ddy <- y[b] - y[a]
    dd <- sqrt(ddx^2 + ddy^2); if (dd < 0.01) dd <- 0.01
    ux0 <- ddx / dd; uy0 <- ddy / dd
    sx <- x[a] + ux0 * r[a]; sy <- y[a] + uy0 * r[a]
    ex <- x[b] - ux0 * (r[b] + MARKER * 0.9); ey <- y[b] - uy0 * (r[b] + MARKER * 0.9)
    segments(sx, sy, ex, ey, col = ecol[k], lwd = if (is_candidate[k]) 1.4 else 1.6, lty = if (is_candidate[k]) 2 else 1)
    m <- MARKER
    if (E(g)$curated_direction[k] == "Activating") {
      bx <- ex - ux0 * m; by <- ey - uy0 * m
      px <- -uy0; py <- ux0
      polygon(c(ex, bx + px * m * 0.55, bx - px * m * 0.55), c(ey, by + py * m * 0.55, by - py * m * 0.55), col = ecol[k], border = NA)
    } else {
      px <- -uy0; py <- ux0
      segments(ex + px * m * 0.5, ey + py * m * 0.5, ex - px * m * 0.5, ey - py * m * 0.5, col = ecol[k], lwd = 2)
    }
  }

  gene_cols <- lfc_to_col(V(g)$gene_log2FoldChange, LFC_CAP)
  for (i in seq_len(vcount(g))) {
    if (is_tf_v[i]) symbols(x[i], y[i], circles = r[i], inches = FALSE, add = TRUE, bg = TF_FILL, fg = TF_BORDER, lwd = 1.4)
    else symbols(x[i], y[i], circles = r[i], inches = FALSE, add = TRUE, bg = gene_cols[i], fg = NA)
  }
  tf_idx <- which(is_tf_v); gene_idx <- which(!is_tf_v)
  GENE_LABEL_CEX <- 0.82
  gene_label_x0 <- x[gene_idx]; gene_label_y0 <- y[gene_idx] + r[gene_idx] + 5
  gene_labels <- V(g)$name[gene_idx]
  gl <- declutter_labels_px(gene_label_x0, gene_label_y0, gene_labels, GENE_LABEL_CEX)
  text(x[tf_idx], y[tf_idx], V(g)$name[tf_idx], cex = 0.78, font = 2, col = "black")
  halo_text(gl$x, gl$y, gene_labels, cex = GENE_LABEL_CEX, col = "black")

  title(main = sprintf("Orthogonally supported edges + candidate uncurated links -- %s (force-simulation layout)", CT), cex.main = 1.6, line = 1.8)
  mtext("Layout: live repulsion+spring+collision+centering simulation (ported from the interactive tuner), not a static FR layout. Solid = curated (CollecTRI), dashed = candidate (not confirmed).",
        side = 3, line = 0.5, cex = 0.78, col = "gray30")

  ## legends, in this script's own pixel coordinate system
  legend(0, -HEIGHT * 0.02, bty = "n", cex = 1.2, title = "Nodes", title.adj = 0,
         legend = c("TF", "Gene/DEG (colour = own log2FC)"),
         pt.bg = c(TF_FILL, "white"), pch = 21, pt.lwd = c(1.4, 1), col = c(TF_BORDER, "black"), pt.cex = c(2.2, 1.4))
  legend(WIDTH * 0.28, -HEIGHT * 0.02, bty = "n", cex = 1.2, title = "Edges", title.adj = 0,
         legend = c("Activating", "Repressive", "Curated (solid)", "Candidate (dashed)"),
         col = c(ACT_COL, REP_COL, "gray20", "gray20"), lty = c(1, 1, 1, 2), lwd = 1.6)
  gx0 <- WIDTH * 0.58; gx1 <- gx0 + WIDTH * 0.22; gy0 <- -HEIGHT * 0.10; gy1 <- gy0 + HEIGHT * 0.035
  xs <- seq(gx0, gx1, length.out = length(LFC_PAL) + 1)
  rect(xs[-length(xs)], gy0, xs[-1], gy1, col = LFC_PAL, border = NA)
  rect(gx0, gy0, gx1, gy1, border = "black", lwd = 1)
  text(gx0, gy1 + HEIGHT * 0.02, "Gene nodes: colour = own RNA log2FC (capped)", cex = 1.05, adj = c(0, 0), font = 2)
  text(gx0, gy0 - HEIGHT * 0.015, sprintf("<= -%.0f", LFC_CAP), cex = 0.95, adj = c(0, 1))
  text(gx1, gy0 - HEIGHT * 0.015, sprintf(">= +%.0f", LFC_CAP), cex = 0.95, adj = c(1, 1))
  mtext("Not validated TF-target relationships in neuroectoderm specifically. CollecTRI support reflects a documented relationship in human cells generally, not tissue-specific validation. Candidate edges are hypothesis-generating only.",
        side = 1, line = 0.5, cex = 0.6, col = "gray40", adj = 0)
  dev.off()
  cat(sprintf("[DONE] %s -- %d nodes, %d edges\n", fpath, vcount(g), ecount(g)))
}

cat("\n[DONE] AR11l_ForceSim_network complete\n")
