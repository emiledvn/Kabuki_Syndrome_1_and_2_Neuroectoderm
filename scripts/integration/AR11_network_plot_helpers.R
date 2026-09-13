# AR11_network_plot_helpers.R -- shared plotting helpers for the AR11
# CollecTRI network figures (AR11b/d/e/f). Sourced (source(), not run
# standalone) from each of those scripts, which must already have LFC_PAL,
# ACT_COL, REP_COL etc. defined before calling draw_gradient_legend()/
# draw_edge_direction_markers() below (passed explicitly as arguments, not
# looked up from the caller's environment, so there's no ordering footgun).
#
# Aesthetic pass (2026-08-03), replacing the per-script duplicated versions
# of these functions with one shared copy:
#   - TF node SHAPE no longer encodes the genome-wide DiffBind direction
#     (circle=UP/square=DOWN). Dropped because AR06_sigTF_site_direction_
#     tally_*.csv showed many TFs' genome-wide call is a near-coinflip
#     across their own bound sites (minority_frac up to ~0.50) -- the shape
#     was implying a confidence the aggregate number doesn't always have.
#     TF identity is now shown by a thicker circle outline only (all nodes
#     are circles); genome-wide DiffBind direction remains available in the
#     underlying tables/other panels, just not as a node shape here.
#   - Gene (DEG) nodes are now sized to fit their own label text, so the
#     label sits centered inside the node instead of being repelled outside
#     it -- declutter_nodes() (unchanged) then pushes same-sized-as-ever TF
#     circles and these now-bigger gene circles apart so they don't overlap.
#     TF labels still use repel_labels() (external placement + leader line),
#     since TF node size is degree-scaled, not label-fit.
#   - Edges get an explicit arrow (Activating) vs T-bar (Repressive) marker
#     at the target end, on top of the existing colour coding, since colour
#     alone doesn't survive grayscale printing or colourblind viewing.
#     igraph has no built-in flat/T-bar arrow end, so these are drawn
#     manually in plot-coordinate space after the base igraph plot() call
#     (same approach this codebase already uses for repel_labels/
#     declutter_nodes -- manual geometry over an unsupported native option).

repel_labels <- function(x, y, labels, cex, iter = 600) {
  n <- length(x)
  if (n <= 1) return(list(x = x, y = y, moved = rep(FALSE, n)))
  # Iteration cap raised 2026-08-18 (was 60/150 for n>80/n>40): that cap was
  # tuned for performance on graphs up to ~150 nodes where it left visible
  # unresolved label overlap at larger cex. O(n^2) per iteration is still
  # cheap at this scale (a few seconds at most even at n~150), so trade a
  # bit of runtime for actually converging.
  iter <- if (n > 80) 300 else if (n > 40) 450 else iter
  w <- nchar(labels) * 0.017 * cex + 0.015
  h <- 0.030 * cex + 0.012
  set.seed(1)
  px <- x + runif(n, -2e-3, 2e-3); py <- y + runif(n, -2e-3, 2e-3)
  bound <- 3 * max(abs(c(x, y)), 1)
  for (it in seq_len(iter)) {
    moved_any <- FALSE
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i >= j) next
      dx <- px[i] - px[j]; dy <- py[i] - py[j]
      ox <- (w[i] + w[j]) / 2 - abs(dx); oy <- (h[i] + h[j]) / 2 - abs(dy)
      if (is.na(ox) || is.na(oy)) next
      if (ox > 0 && oy > 0) {
        moved_any <- TRUE
        if (ox < oy) { s <- ox * 0.5 + 1e-4; sg <- if (dx >= 0) 1 else -1; px[i] <- px[i] + sg*s; px[j] <- px[j] - sg*s }
        else         { s <- oy * 0.5 + 1e-4; sg <- if (dy >= 0) 1 else -1; py[i] <- py[i] + sg*s; py[j] <- py[j] - sg*s }
      }
    }
    px <- pmin(pmax(px, -bound), bound); py <- pmin(pmax(py, -bound), bound)
    if (!moved_any) break
  }
  px[!is.finite(px)] <- x[!is.finite(px)]; py[!is.finite(py)] <- y[!is.finite(py)]
  list(x = px, y = py, moved = sqrt((px-x)^2 + (py-y)^2) > 0.015)
}

# Post-layout node de-overlap pass -- see AR06/AR11b history for why FR +
# this collision pass beats the untuneable alternative layouts.
declutter_nodes <- function(coords, radius, iter = 1500) {
  n <- nrow(coords); px <- coords[, 1]; py <- coords[, 2]
  for (it in seq_len(iter)) {
    moved_any <- FALSE
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i >= j) next
      dx <- px[i] - px[j]; dy <- py[i] - py[j]
      d <- sqrt(dx^2 + dy^2)
      min_d <- radius[i] + radius[j]
      if (d < min_d) {
        moved_any <- TRUE
        if (d > 1e-6) { ux <- dx / d; uy <- dy / d } else { a <- runif(1, 0, 2*pi); ux <- cos(a); uy <- sin(a) }
        s <- (min_d - d) * 0.5 + 1e-4
        px[i] <- px[i] + ux*s; py[i] <- py[i] + uy*s
        px[j] <- px[j] - ux*s; py[j] <- py[j] - uy*s
      }
    }
    if (!moved_any) break
  }
  cbind(px, py)
}

# vsize (igraph vertex.size units) for a gene node so its circle's diameter
# comfortably contains its own label at the given cex. Converts the needed
# radius (in plot-coordinate units) to vsize via *200 -- NOT the caller's
# `divisor` (used elsewhere for declutter_nodes' collision radius, a
# separate, looser "personal space" knob). Bug fixed 2026-08-18: this used
# to multiply by `divisor` (70-80 across callers) instead of 200, the fixed
# constant igraph actually uses to turn vertex.size into a rendered radius
# (see draw_edge_direction_markers()'s comment below) -- so every label-fit
# circle came out ~2.3x too small for its own text regardless of cex, most
# visible on long gene names (e.g. "FPGT-TNNI3K") at larger cex. pad=1.6,
# not ~1: declutter_nodes only guarantees circles don't OVERLAP (can end up
# exactly tangent, zero gap), so a label sized to just barely fit its own
# circle will visually touch a neighbour's label the moment two circles
# meet. The extra headroom here, plus the explicit extra gap
# declutter_radius_pad() adds below, is what actually keeps adjacent
# labels apart.
gene_label_vsize <- function(labels, cex, pad = 1.6, min_vsize = 4.5) {
  w <- nchar(labels) * 0.017 * cex + 0.02
  radius_needed <- (w / 2) * pad
  pmax(radius_needed * 200, min_vsize)
}

# Same idea, for TF nodes: floor the degree-based size so a TF with few
# targets (small circle) still fits its own (bold) label -- `degree_vsize`
# is the caller's existing out-degree-scaled size, used as-is whenever it's
# already bigger than what the label needs.
tf_label_vsize <- function(labels, cex, degree_vsize, pad = 1.6, min_vsize = 4.5) {
  pmax(degree_vsize, gene_label_vsize(labels, cex * 1.08, pad, min_vsize))
}

# Collision radius for declutter_nodes(), in the same plot-coordinate units
# as the layout -- vsize/200 (the true rendered radius, see
# draw_edge_direction_markers()'s comment below) plus a multiplicative
# margin so tangent circles keep a visible gap instead of touching at zero
# distance. Bug fixed 2026-08-18 alongside gene_label_vsize()/
# tf_label_vsize(): this used to divide by the caller's `divisor` (70-80)
# instead of 200, which happened to roughly cancel against those two
# functions' matching bug (both used the same wrong divisor) as long as
# vsize was ALSO wrong -- fixing vsize's conversion alone, without fixing
# this one too, made declutter reserve ~2.5x too much space per node
# (dividing an now-correctly-bigger vsize by 80 instead of 200), inflating
# the whole layout/canvas far more than the bigger circles actually need.
declutter_radius_pad <- function(vsize, extra = 1.18) (vsize / 200) * extra

# Thicker, clearer log2FC gradient legend. `pal` is the caller's LFC_PAL.
draw_gradient_legend <- function(x0, x1, y0, y1, lo_label, hi_label, title, pal,
                                  cex_title = 0.85, cex_lab = 0.78, border_lwd = 1.1) {
  n <- length(pal)
  xs <- seq(x0, x1, length.out = n + 1)
  rect(xs[-length(xs)], y0, xs[-1], y1, col = pal, border = NA)
  rect(x0, y0, x1, y1, border = "black", lwd = border_lwd)
  xmid <- (x0 + x1) / 2
  segments(xmid, y0, xmid, y1, col = "gray15", lwd = 0.9)
  text(x0, y1 + 0.045, title, cex = cex_title, adj = c(0, 0), col = "black", font = 2)
  text(x0, y0 - 0.035, lo_label, cex = cex_lab, adj = c(0, 1), col = "black")
  text(x1, y0 - 0.035, hi_label, cex = cex_lab, adj = c(1, 1), col = "black")
  text(xmid, y0 - 0.035, "0", cex = cex_lab, adj = c(0.5, 1), col = "black")
}

# Manual edge-endpoint direction markers: pointed arrow = Activating, flat
# T-bar = Repressive -- drawn at the target node's ACTUAL circle border so
# direction is readable from shape alone, not just edge colour. igraph's
# plot() should be called with edge.arrow.mode = 0 alongside this (no native
# arrowheads, so nothing is drawn twice).
#
# Node radius in plot-coordinate units is vertex.size/200 -- this is
# igraph's own internal constant (confirmed directly from
# igraph:::shapes("circle")$clip and $plot, both of which compute
# `vertex.size <- 1/200 * params("vertex","size")` before doing anything
# geometric with it), NOT the vsize/70-vsize/80 divisors used elsewhere in
# these scripts for declutter_nodes() -- those are an unrelated, looser
# "how much breathing room to reserve between nodes" tuning knob, not the
# true rendered radius, and using them here was the original bug (markers
# floating away from the node instead of sitting on its border). Marker
# size scales with each edge's own target-node radius (arrow_len_f etc.),
# not a fixed absolute size, so a marker on a big hub-TF node and one on a
# small gene node both look proportionate. Edge curvature (edge.curved) is
# ignored for placement -- these networks use small curvature (<=0.1), a
# minor, acceptable offset right at the endpoint where the curve is nearly
# tangent to the straight line anyway.
draw_edge_direction_markers <- function(g, lay, vsize, is_activating, col, lwd,
                                         arrow_len_f = 1.0, arrow_width_f = 0.55,
                                         bar_half_f = 0.65, gap_f = 0.12,
                                         marker_vsize = NULL) {
  node_r <- vsize / 200
  # marker_vsize (optional): a FIXED vsize used only to size the arrow/T-bar
  # glyph itself, decoupled from `node_r`, which must stay based on the
  # REAL vsize since it's also what places the marker at the node's actual
  # border. Conflating the two (as an earlier version of this call did, by
  # passing a fake uniform `vsize` for both) put the marker's border_pt at
  # where a small node's edge would be, landing INSIDE any node whose true
  # circle is bigger than that -- e.g. every arrow into a big TF node.
  marker_r_fixed <- if (!is.null(marker_vsize)) marker_vsize / 200 else NULL
  em <- igraph::ends(g, igraph::E(g), names = FALSE)
  for (i in seq_len(nrow(em))) {
    a <- em[i, 1]; b <- em[i, 2]
    p_from <- lay[a, ]; p_to <- lay[b, ]
    d <- p_to - p_from; dist <- sqrt(sum(d^2))
    if (!is.finite(dist) || dist < 1e-8) next
    u <- d / dist
    perp <- c(-u[2], u[1])
    target_r <- node_r[b]
    marker_r <- if (is.null(marker_r_fixed)) target_r else marker_r_fixed
    border_pt <- p_to - u * (target_r * (1 + gap_f))
    if (isTRUE(is_activating[i])) {
      arrow_len <- marker_r * arrow_len_f
      arrow_width <- marker_r * arrow_width_f
      base_c <- border_pt - u * arrow_len
      b1 <- base_c + perp * arrow_width; b2 <- base_c - perp * arrow_width
      polygon(c(border_pt[1], b1[1], b2[1]), c(border_pt[2], b1[2], b2[2]),
              col = col[i], border = NA)
    } else {
      bar_half <- marker_r * bar_half_f
      bar_c <- border_pt - u * (bar_half * 0.1)
      e1 <- bar_c + perp * bar_half; e2 <- bar_c - perp * bar_half
      segments(e1[1], e1[2], e2[1], e2[2], col = col[i], lwd = lwd[i] * 1.5)
    }
  }
}
