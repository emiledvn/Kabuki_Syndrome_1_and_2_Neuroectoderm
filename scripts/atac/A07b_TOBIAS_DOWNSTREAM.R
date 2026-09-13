#!/usr/bin/env Rscript
# A07b_TOBIAS_DOWNSTREAM.R -- DAR/expression attribution for TOBIAS differential
# TF binding: which TFs change footprinting AND show concordant expression change?
# This is the "05d" post-processing step 05_tobias_consensus.sh's header describes
# as not-yet-written -- it did not exist anywhere before this script.
#
# Adapted from data_import_local/local_scripts/09_INTEGRATION_tobias_analysis.R.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/A07b_TOBIAS_DOWNSTREAM.R
#
# CRITICAL FIX vs the original: the original pointed at
# data/NE_ATAC/TOBIAS_Results_Strict -- the OLD circular strict-combined TOBIAS run
# (peaks pre-filtered by differential accessibility). It predates and was never
# updated for 05_tobias_consensus.sh's fix. This version reads from A06's
# consensus-peak BINDetect output instead (data/tobias_output/<contrast>/BINDetect/),
# which is what makes the whole exercise statistically valid in the first place.
#
# Self-checkpointing: skips entirely if results/tables/A07b_concordant_tfs.csv exists.
# Requires A06_TOBIAS.sh and R02_DESeq2.R to have completed first.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(ggVennDiagram)
  library(eulerr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(yaml)
})

## Shared helper for both family-grouped TF-binding heatmaps below.
## ComplexHeatmap (not pheatmap) specifically because row_split's row_title
## prints the family NAME directly next to its block of rows -- no color
## legend needed for family membership at all, per request (pheatmap's
## annotation_row can only show a color strip + a separate legend).
## cluster_row_slices = FALSE preserves the exact row order passed in
## (already family-block-then-alphabetical) rather than ComplexHeatmap
## re-sorting blocks itself.
plot_family_heatmap <- function(mat, star_mat, family_vec, title, out_file, height, width = 3.5) {
  family_factor <- factor(family_vec, levels = unique(family_vec))
  col_fun <- colorRamp2(c(-max(abs(mat), na.rm = TRUE), 0, max(abs(mat), na.rm = TRUE)),
                         c("#4575B4", "white", "#D73027"))
  ht <- Heatmap(mat, name = "Binding\nchange", col = col_fun,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_split = family_factor, cluster_row_slices = FALSE,
                row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
                row_gap = unit(2, "mm"), border = TRUE,
                row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 11), column_names_rot = 45,
                column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                cell_fun = function(j, i, x, y, w, h, col) {
                  if (star_mat[i, j] == "*") grid.text("*", x, y, gp = gpar(fontsize = 13, col = "black"))
                })
  pdf(out_file, width = width, height = height)
  draw(ht)
  dev.off()
}

## Real MEME-minimal-format PWM parser + base-graphics sequence-logo renderer,
## copied from AR12_TOBIAS_method_illustration.R (same real motif file this
## pipeline's BINDetect used -- data/reference/motifs.meme -- not invented).
read_meme_pwm <- function(meme_lines, motif_id) {
  idx <- grep(paste0("^MOTIF ", motif_id, " "), meme_lines)
  if (length(idx) == 0) stop("motif not found in motifs.meme: ", motif_id)
  header_line <- meme_lines[idx[1] + 1]
  w <- as.integer(sub(".*w= *([0-9]+).*", "\\1", header_line))
  rows <- meme_lines[(idx[1] + 2):(idx[1] + 1 + w)]
  mat <- do.call(rbind, lapply(rows, function(r) as.numeric(strsplit(trimws(r), "\\s+")[[1]])))
  colnames(mat) <- c("A", "C", "G", "T")
  mat
}
draw_seqlogo <- function(pwm, xpos, y0 = 0, max_height = 1, slot_width = 0.85) {
  base_col <- c(A = "#009E73", C = "#0072B2", G = "#E69F00", T = "#D55E00")
  ic <- 2 + rowSums(pwm * log2(pmax(pwm, 1e-9)))  # bits, DNA alphabet max = 2
  nominal_h <- strheight("A", cex = 1); nominal_w <- strwidth("A", cex = 1)
  for (i in seq_len(nrow(pwm))) {
    heights <- pwm[i, ] * ic[i] / 2 * max_height
    ord <- order(heights)  # smallest first -> tallest letter ends up on top
    ycum <- y0
    for (b in ord) {
      h <- heights[b]
      if (h < 1e-4) next
      cex_b <- min(h / nominal_h, slot_width / nominal_w)
      text(xpos[i], ycum + h / 2, names(heights)[b], col = base_col[names(heights)[b]], font = 2, cex = cex_b)
      ycum <- ycum + h
    }
  }
}
## Row order matches `mat`'s rownames exactly (required by anno_image).
## Representative motif_id per TF = highest total_tfbs across all contrasts'
## bindetect_results.txt (a TF can match >1 JASPAR matrix; most-bound-sites
## is the most representative model, independent of which contrast's row
## slice_max(abs(binding_change)) happened to keep upstream).
render_tf_motif_logos <- function(mat, tobias_dir, contrasts) {
  meme_lines <- readLines("data/reference/motifs.meme")
  bindetect_all <- bind_rows(lapply(contrasts, function(ct) {
    read_tsv(file.path(tobias_dir, ct$name, "BINDetect", "bindetect_results.txt"), show_col_types = FALSE) %>%
      dplyr::select(name, motif_id, total_tfbs)
  }))
  motif_map <- bindetect_all %>% filter(name %in% rownames(mat)) %>%
    group_by(name) %>% slice_max(total_tfbs, n = 1, with_ties = FALSE) %>% ungroup()
  missing_motif <- setdiff(rownames(mat), motif_map$name)
  if (length(missing_motif) > 0) stop("[A07b] no motif_id found for: ", paste(missing_motif, collapse = ", "))
  motif_id_for <- setNames(motif_map$motif_id, motif_map$name)

  # Each TF's logo is rendered on its own canvas, sized to its own motif
  # width and autoscaled to its own tallest position -- shorter/narrower
  # motifs end up with visually bigger letters than long ones, since each
  # gets to use the full annotation box on its own. height/width ratio
  # matched to the anno_image box's actual aspect ratio (measured from a
  # rendered preview) so each logo fills the box with no pillar/letterboxing.
  logo_dir <- tempfile("A07b_logos_"); dir.create(logo_dir)
  logo_paths <- character(nrow(mat))
  for (i in seq_len(nrow(mat))) {
    tf <- rownames(mat)[i]
    pwm <- read_meme_pwm(meme_lines, motif_id_for[tf])
    w <- nrow(pwm)
    ic <- 2 + rowSums(pwm * log2(pmax(pwm, 1e-9)))
    y_top <- max(ic) / 2  # tallest position's letter-stack height at max_height=1 -- crops dead headroom
    out_png <- file.path(logo_dir, paste0(tf, ".png"))
    png(out_png, width = 55 * w, height = round(55 * w / 3.75), res = 300, bg = "white")
    par(mar = c(0, 0, 0, 0))
    plot(NA, xlim = c(0.5, w + 0.5), ylim = c(0, y_top), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    draw_seqlogo(pwm, xpos = seq_len(w), y0 = 0, max_height = 1, slot_width = 0.98)
    dev.off()
    logo_paths[i] <- out_png
  }
  setNames(logo_paths, rownames(mat))
}

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A07b] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TOBIAS_DIR  <- "data/tobias_output"
TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/atac"
RDS_DIR     <- "results/RDS"
SESSION_DIR <- "results/Session_info"
CONCORDANT_OUT <- file.path(TABLES_DIR, "A07b_concordant_tfs.csv")

for (d in c(TABLES_DIR, PLOTS_DIR, RDS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(CONCORDANT_OUT)) {
  cat("[A07b] Already complete (", CONCORDANT_OUT, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS   <- cfg$contrasts
PADJ_CUTOFF <- cfg$thresholds$fdr
MIN_EXPR    <- cfg$thresholds$min_mean_expression
TOBIAS_PVAL <- cfg$thresholds$tobias_binding_pvalue
TOBIAS_CHANGE <- cfg$thresholds$tobias_binding_change

if (!dir.exists(TOBIAS_DIR)) stop("[A07b] ERROR: ", TOBIAS_DIR, " not found -- run A06_TOBIAS.sh first.")
norm_counts_path <- file.path(RDS_DIR, "R02_normalized_counts.rds")
if (!file.exists(norm_counts_path)) stop("[A07b] ERROR: ", norm_counts_path, " not found -- run R02_DESeq2.R first.")

theme_pub <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())

cat("================================================================================\n")
cat("A07b_TOBIAS_DOWNSTREAM -- DAR/expression attribution for differential TF binding\n")
cat(sprintf("THRESHOLDS: binding pvalue < %.2f AND |change| > %.2f | expression padj < %.2f (no LFC floor -- see note below)\n",
            TOBIAS_PVAL, TOBIAS_CHANGE, PADJ_CUTOFF))
cat("================================================================================\n\n")

########
## Expression filter universe: mean normalized count > MIN_EXPR across all samples
########

norm_counts <- readRDS(norm_counts_path)
gene_map    <- readRDS(file.path(RDS_DIR, "R02_gene_id_name_map.rds"))
mean_counts <- tibble(gene_id = rownames(norm_counts), mean_count = rowMeans(norm_counts)) %>%
  left_join(gene_map, by = "gene_id")
expressed_genes <- mean_counts %>% filter(!is.na(mean_count), mean_count > MIN_EXPR) %>% pull(gene_name)
cat(sprintf("[A07b] Expressed genes (mean count > %d): %d\n\n", MIN_EXPR, length(expressed_genes)))

########
## Load BINDetect results from A06's consensus-peak output
########

load_tobias_data <- function(contrast_name) {
  f <- file.path(TOBIAS_DIR, contrast_name, "BINDetect", "bindetect_results.txt")
  cat("[A07b] Loading:", f, "\n")
  if (!file.exists(f)) { cat("  WARNING: not found\n"); return(NULL) }
  tobias <- read_tsv(f, show_col_types = FALSE)
  change_col  <- grep("_change$", colnames(tobias), value = TRUE)[1]
  pvalue_col  <- grep("_pvalue$", colnames(tobias), value = TRUE)[1]
  if (is.na(change_col) || is.na(pvalue_col)) { cat("  WARNING: cannot find change/pvalue columns\n"); return(NULL) }
  # TOBIAS names this column "{reference}_{treatment}_change" (A06 passes
  # --cond_names REF TREAT) but the VALUE is reference minus treatment, not
  # treatment minus reference -- confirmed numerically, not just from the
  # name: RFX1 has WT_mean_score=0.179 < KDM6A_ko_mean_score=0.200 (higher in
  # the mutant) yet WT_KDM6A_ko_change=-0.219 (negative); SOX4 has
  # WT=0.138 > KDM6A_ko=0.126 (higher in WT) yet change=+0.223 (positive).
  # Negating here makes the working `binding_change` treatment-relative
  # (positive = higher in the mutant), consistent with every other
  # fold-change-style quantity in this pipeline (RNA/ATAC log2FoldChange) --
  # without this, "Increased"/"Decreased" below would be exactly backwards.
  tobias %>%
    dplyr::select(TF_name = name, binding_change = !!sym(change_col), binding_pvalue = !!sym(pvalue_col)) %>%
    mutate(binding_change = -binding_change,
           binding_sig = binding_pvalue < TOBIAS_PVAL & abs(binding_change) > TOBIAS_CHANGE,
           binding_direction = case_when(binding_change > 0 ~ "Increased", binding_change < 0 ~ "Decreased", TRUE ~ "No change"),
           contrast = contrast_name)
}

filter_and_dedup <- function(df) {
  if (is.null(df)) return(NULL)
  df %>% filter(TF_name %in% expressed_genes) %>%
    group_by(TF_name, contrast) %>% slice_max(abs(binding_change), n = 1, with_ties = FALSE) %>% ungroup()
}

tobias_by_contrast <- list()
for (ct in CONTRASTS) tobias_by_contrast[[ct$name]] <- filter_and_dedup(load_tobias_data(ct$name))

########
## Integrate with RNA-seq: is this TF's binding change accompanied by a
## concordant expression change?
########

integrate_tf <- function(tobias_data, rna_data) {
  if (is.null(tobias_data)) return(NULL)
  tobias_data %>%
    left_join(rna_data %>% dplyr::select(gene_name, log2FoldChange, padj) %>% rename(TF_name = gene_name), by = "TF_name") %>%
    # No LFC floor here, unlike the rest of the pipeline's lfc_threshold
    # convention -- deliberate, not an oversight. "Concordant TF" already
    # requires two independent significant signals to agree (TOBIAS binding
    # pvalue/magnitude AND RNA padj, matching direction), which is a much
    # stronger evidentiary bar than a single-axis LFC cutoff; adding a
    # magnitude floor on top would be redundant and would cost real hits,
    # since TFs are typically low-abundance/low-fold-change genes where a
    # small expression shift can still be functionally consequential (same
    # point already documented in R02_DESeq2.R's no-lfcShrink rationale).
    mutate(expression_sig = !is.na(padj) & padj < PADJ_CUTOFF,
           expression_direction = case_when(log2FoldChange > 0 ~ "Increased", log2FoldChange < 0 ~ "Decreased", TRUE ~ "No change"),
           tf_class = case_when(
             binding_sig & expression_sig & binding_direction == expression_direction ~ "Concordant TF",
             binding_sig ~ "Differential binding",
             TRUE ~ "NS"))
}

tobias_integrated <- bind_rows(lapply(CONTRASTS, function(ct) {
  rna_data <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)), show_col_types = FALSE)
  integrate_tf(tobias_by_contrast[[ct$name]], rna_data)
}))

cat("\n[A07b] TF classification summary:\n")
print(tobias_integrated %>% group_by(contrast, tf_class) %>% summarise(count = n(), .groups = "drop"))

for (ct in CONTRASTS) {
  subset_df <- tobias_integrated %>% filter(contrast == ct$name)
  if (nrow(subset_df) > 0) write_csv(subset_df, file.path(TABLES_DIR, sprintf("A07b_tf_binding_%s.csv", ct$name)))
}

concordant_tfs <- tobias_integrated %>% filter(tf_class == "Concordant TF") %>% arrange(contrast, desc(abs(binding_change)))
write_csv(concordant_tfs, CONCORDANT_OUT)
all_diff_tfs <- tobias_integrated %>% filter(binding_sig) %>% arrange(contrast, desc(abs(binding_change)))
write_csv(all_diff_tfs, file.path(TABLES_DIR, "A07b_all_differential_tfs.csv"))

########
## Plots
########

COLORS <- setNames(c("#E69F00", "#56B4E9"), c(CONTRASTS[[1]]$name, CONTRASTS[[2]]$name))

## Top differential TFs bar plot
top_tfs <- tobias_integrated %>% filter(binding_sig) %>% group_by(contrast) %>%
  arrange(desc(abs(binding_change))) %>% slice_head(n = 20) %>% ungroup()
if (nrow(top_tfs) > 0) {
  p1 <- ggplot(top_tfs, aes(x = reorder(TF_name, binding_change), y = binding_change, fill = contrast)) +
    geom_col(color = "black", linewidth = 0.2) + coord_flip() +
    facet_wrap(~contrast, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = COLORS) +
    labs(title = "Top 20 Differential TF Binding Changes", x = "Transcription Factor", y = "Binding Change Score") +
    theme_pub + theme(axis.text.y = element_text(size = 8), legend.position = "none")
  ggsave(file.path(PLOTS_DIR, "A07b_top_tfs_barplot.pdf"), p1, width = 10, height = 8)
}

## TF overlap Venn -- area-proportional (eulerr), not ggVennDiagram's default
## of equal-size circles regardless of how lopsided the two set sizes are
## (KDM6A_ko_vs_WT: 130 sig. TFs; KMT2D_Het_vs_WT: 177 -- a real size
## difference an equal-circle Venn was hiding). Circle labels use the short
## treatment name (matches every other plot's convention); COLORS itself
## stays keyed by full contrast name since that's what tobias_integrated$contrast holds.
## unique() is required, not cosmetic -- bindetect_results.txt has multiple
## motif entries per TF (e.g. several JASPAR variants of the same symbol),
## so pull(TF_name) contains repeats; euler() (unlike ggVennDiagram) errors
## on a set with duplicate elements ("vectors in `combinations` cannot
## contain duplicates"), since it treats input as strict sets, not multisets.
venn_data <- setNames(lapply(CONTRASTS, function(ct) tobias_integrated %>% filter(contrast == ct$name, binding_sig) %>% pull(TF_name) %>% unique()),
                       sapply(CONTRASTS, function(ct) ct$name))
treatment_labels <- sapply(CONTRASTS, function(ct) ct$treatment)

fit <- euler(venn_data)
totals <- fit$original.values
pct_labels <- sprintf("%d\n(%.0f%%)", totals, 100 * totals / sum(totals))

pdf(file.path(PLOTS_DIR, "A07b_tf_overlap_venn.pdf"), width = 7, height = 6)
print(plot(fit,
           fills = list(fill = unname(COLORS[names(venn_data)]), alpha = 0.65),
           edges = list(col = "black", lwd = 1.5),
           labels = list(labels = treatment_labels, col = "black", font = 2, cex = 1.15),
           quantities = list(labels = pct_labels, col = "black", font = 2, cex = 1),
           main = "Differential TF Overlap Between Conditions (TOBIAS, area-proportional)"))
dev.off()

## Binding vs expression scatter + volcano, per contrast
for (ct in CONTRASTS) {
  plot_data <- tobias_integrated %>% filter(contrast == ct$name, !is.na(log2FoldChange)) %>%
    mutate(tf_class = factor(tf_class, levels = c("NS", "Differential binding", "Concordant TF")))
  if (nrow(plot_data) == 0) next
  label_data <- plot_data %>% filter(tf_class == "Concordant TF")

  p_scatter <- ggplot(plot_data, aes(x = binding_change, y = log2FoldChange, color = tf_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_vline(xintercept = c(-TOBIAS_CHANGE, TOBIAS_CHANGE), linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(alpha = 0.6) +
    geom_text_repel(data = label_data, aes(label = TF_name), size = 3, max.overlaps = 25, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c(NS = "#CCCCCC", `Differential binding` = "#0072B2", `Concordant TF` = "#D55E00"), name = NULL) +
    labs(title = paste(ct$name, "-- TF Binding vs Expression"), x = "TOBIAS Binding Change", y = "RNA-seq log2 Fold Change") +
    theme_pub + theme(legend.position = "bottom")
  ggsave(file.path(PLOTS_DIR, sprintf("A07b_binding_vs_expression_%s.pdf", ct$name)), p_scatter, width = 8, height = 6)

  plot_data_v <- plot_data %>% mutate(neg_log10_p = -log10(binding_pvalue + 1e-300))
  label_data_v <- plot_data_v %>% filter(tf_class %in% c("Concordant TF", "Differential binding")) %>%
    distinct(TF_name, .keep_all = TRUE)
  p_volcano <- ggplot(plot_data_v, aes(x = binding_change, y = neg_log10_p, color = tf_class)) +
    geom_hline(yintercept = -log10(TOBIAS_PVAL), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-TOBIAS_CHANGE, TOBIAS_CHANGE), linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7) +
    geom_text_repel(data = label_data_v, aes(label = TF_name), size = 3, max.overlaps = Inf, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c(NS = "#CCCCCC", `Differential binding` = "#0072B2", `Concordant TF` = "#D55E00"), name = NULL) +
    labs(title = paste("Volcano:", ct$name), x = "Binding Change Score", y = expression(-log[10](p - value))) +
    theme_pub + theme(legend.position = "bottom")
  ggsave(file.path(PLOTS_DIR, sprintf("A07b_volcano_%s.pdf", ct$name)), p_volcano, width = 7, height = 7)
}

## Concordant TF count summary
if (nrow(concordant_tfs) > 0) {
  conc_summary <- concordant_tfs %>% group_by(contrast) %>% summarise(count = n(), .groups = "drop")
  p5 <- ggplot(conc_summary, aes(x = contrast, y = count, fill = contrast)) +
    geom_col(color = "black", width = 0.6) + geom_text(aes(label = count), vjust = -0.5, fontface = "bold") +
    scale_fill_manual(values = COLORS) +
    labs(title = "Concordant TFs by Condition", y = "Number of Concordant TFs", x = NULL) +
    theme_pub + theme(legend.position = "none")
  ggsave(file.path(PLOTS_DIR, "A07b_concordant_tfs_summary.pdf"), p5, width = 5, height = 5)
}

## Heatmap of top-variable TF binding patterns across contrasts.
## chromVAR (A10) was tried as cross-validation here and dropped: its
## differentialDeviations() t-test has wildly different power between the
## two contrasts (2/612 motifs genome-wide reach padj<0.05 in KDM6A_ko_vs_WT
## vs 407/612 in KMT2D_Het_vs_WT, a >100x difference in baseline positive
## rate) -- not a strong, even-handed line of support, whichever threshold
## is picked. Reverted to marking each cell with TOBIAS's OWN binding_sig
## flag instead (this same contrast's own pvalue<TOBIAS_PVAL AND
## |change|>TOBIAS_CHANGE call) -- self-consistent across both contrasts by
## construction, unlike chromVAR's cross-contrast power asymmetry.
##
## Selection of which TFs enter this pool at all is still based on
## binding_sig in at least one contrast, but the VALUE shown for every
## included TF is its real binding_change in every contrast, not just the
## ones where it happened to be significant. Filtering the value itself to
## binding_sig (an earlier approach) silently substituted 0 for any
## non-significant cell -- which reads as "no change" on the heatmap when
## the true value might be a real, just-not-significant change in the same
## direction. values_fill stays NA (not 0) for the genuine edge case of a TF
## entirely absent from one contrast's BINDetect output, so that's visually
## distinct from "no change".
sig_tfs <- tobias_integrated %>% filter(binding_sig) %>% pull(TF_name) %>% unique()
heatmap_data <- tobias_integrated %>% filter(TF_name %in% sig_tfs) %>% dplyr::select(TF_name, contrast, binding_change) %>%
  group_by(TF_name, contrast) %>% summarise(binding_change = mean(binding_change, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = contrast, values_from = binding_change, values_fill = NA)
if (nrow(heatmap_data) > 0 && ncol(heatmap_data) > 2) {
  mat <- heatmap_data %>% column_to_rownames("TF_name") %>% as.matrix()
  tf_var <- apply(mat, 1, var, na.rm = TRUE)
  top_var_tfs <- names(sort(tf_var, decreasing = TRUE))[1:min(50, length(tf_var))]
  mat <- mat[top_var_tfs, , drop = FALSE]

  ## Star = TOBIAS's own binding_sig for THIS TF in THIS contrast (not a
  ## cross-method check) -- re-derived from tobias_integrated by exact
  ## (TF_name, contrast) match, since that's the authoritative per-cell
  ## significance call already computed above.
  sig_lookup <- tobias_integrated %>% distinct(TF_name, contrast, binding_sig)
  star_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      tf <- rownames(mat)[i]; ct_name <- colnames(mat)[j]
      if (is.na(mat[i, j])) next
      hit <- sig_lookup %>% filter(TF_name == tf, contrast == ct_name)
      if (nrow(hit) == 1 && isTRUE(hit$binding_sig)) star_mat[i, j] <- "*"
    }
  }

  ## TF family grouping (curatorial, documented here -- not an automated
  ## classifier): standard TF structural-family assignments (Lambert et al.
  ## 2018 human TF classification / JASPAR family conventions) for exactly
  ## the 50 TFs in this pool, built by hand once and checked against every
  ## name actually present (see the full top50 list this was built from,
  ## printed by [A07b] above) -- NOT a generic prefix rule (would wrongly
  ## split e.g. nuclear receptors, which don't share a common prefix, and
  ## wrongly lump e.g. ZBED/ZBTB/generic-ZNF as one "ZNF" bucket when they're
  ## structurally distinct zinc-finger subclasses). Extend this table by
  ## hand if a future rerun's top-50 set includes a TF not listed here
  ## (falls back to "Other" rather than erroring, logged below).
  tf_family <- c(
    KLF10="SP/KLF", KLF11="SP/KLF", KLF16="SP/KLF", KLF4="SP/KLF", KLF5="SP/KLF", KLF7="SP/KLF",
    SP3="SP/KLF", SP4="SP/KLF", SP8="SP/KLF", SP9="SP/KLF",
    RFX1="RFX", RFX2="RFX", RFX3="RFX", RFX5="RFX", RFX7="RFX",
    MAZ="C2H2-ZNF (misc)", ZNF135="C2H2-ZNF (misc)", ZNF148="C2H2-ZNF (misc)", ZNF211="C2H2-ZNF (misc)",
    ZNF281="C2H2-ZNF (misc)", ZNF610="C2H2-ZNF (misc)", ZNF8="C2H2-ZNF (misc)", ZNF93="C2H2-ZNF (misc)",
    ZBED1="ZBED", ZBED4="ZBED",
    ZBTB14="ZBTB", ZBTB17="ZBTB",
    NR1D1="Nuclear receptor", NR2F1="Nuclear receptor", RXRB="Nuclear receptor",
    CUX1="CUT homeodomain", CUX2="CUT homeodomain",
    ONECUT2="ONECUT homeodomain", ONECUT3="ONECUT homeodomain",
    OTX2="Other homeodomain", PBX2="Other homeodomain", POU4F1="Other homeodomain",
    SNAI1="Snail/ZEB (EMT)", SNAI2="Snail/ZEB (EMT)", ZEB1="Snail/ZEB (EMT)",
    TCF3="bHLH (E-protein)", TCF4="bHLH (E-protein)", TCF12="bHLH (E-protein)", TCFL5="bHLH (E-protein)",
    TEAD3="TEAD", TEAD4="TEAD",
    DBP="PAR-bZIP", TEF="PAR-bZIP",
    SATB1="SATB/AT-hook", TBP="Basal (TBP)"
  )
  family_vec <- tf_family[rownames(mat)]
  n_unclassified <- sum(is.na(family_vec))
  if (n_unclassified > 0) {
    cat(sprintf("[A07b] WARNING: %d TF(s) in the heatmap have no family assignment in tf_family -- shown as 'Other': %s\n",
                n_unclassified, paste(rownames(mat)[is.na(family_vec)], collapse = ", ")))
  }
  family_vec[is.na(family_vec)] <- "Other"

  ## Order rows by family block (alphabetical family, alphabetical TF within
  ## family), not hierarchical clustering -- clustering by binding_change
  ## alone would scatter same-family TFs across the heatmap by chance
  ## similarity in only 2 columns, defeating the point of grouping by category.
  row_order <- order(family_vec, rownames(mat))
  mat <- mat[row_order, , drop = FALSE]
  star_mat <- star_mat[row_order, , drop = FALSE]
  family_vec <- family_vec[row_order]
  gaps_row <- which(diff(as.integer(factor(family_vec, levels = unique(family_vec)))) != 0)

  name_to_treatment <- setNames(sapply(CONTRASTS, `[[`, "treatment"), sapply(CONTRASTS, `[[`, "name"))
  colnames(mat) <- name_to_treatment[colnames(mat)]
  colnames(star_mat) <- colnames(mat)

  plot_family_heatmap(mat, star_mat, family_vec,
                      "TOBIAS TF binding change, top 50 by variance\n(* = significant in this contrast)",
                      file.path(PLOTS_DIR, "A07b_heatmap_top50.pdf"), height = 13)
}

########
## Extended heatmap: ALL TOBIAS-significant TFs (binding_sig in >=1
## contrast, not just the top 50 by variance), user-curated family
## annotation (coarser than the top-50 heatmap's -- all homeodomain
## subclasses collapsed to one "Homeodomain" bucket, all zinc-finger
## subclasses (C2H2-misc/ZBED/ZBTB/PATZ1) collapsed to one "Zinc Finger"
## bucket, TBP/SATB1/DBP/TEF/NFYC collapsed to "Basal / Chromatin" -- a
## deliberately different, coarser grouping than the top-50 plot's, kept as
## a separate file rather than replacing it). 112 TFs total, exact mapping
## supplied and checked against the full binding_sig TF list -- covers all
## of it, nothing falls to "Other".
########
extended_family <- c(
  FOS="AP-1 (Fos/Jun)", FOSL1="AP-1 (Fos/Jun)", FOSL2="AP-1 (Fos/Jun)", JDP2="AP-1 (Fos/Jun)",
  JUNB="AP-1 (Fos/Jun)", JUND="AP-1 (Fos/Jun)",
  BACH1="bZIP (BACH/Maf)", BACH2="bZIP (BACH/Maf)", MAFK="bZIP (BACH/Maf)", NFE2="bZIP (BACH/Maf)",
  TEAD1="TEAD", TEAD2="TEAD", TEAD3="TEAD", TEAD4="TEAD",
  SNAI1="Snail/ZEB (EMT)", SNAI2="Snail/ZEB (EMT)", ZEB1="Snail/ZEB (EMT)",
  SOX10="SOX (HMG-box)", SOX13="SOX (HMG-box)", SOX15="SOX (HMG-box)", SOX2="SOX (HMG-box)",
  SOX4="SOX (HMG-box)", SOX8="SOX (HMG-box)", SOX9="SOX (HMG-box)", SRY="SOX (HMG-box)",
  FOXA3="FOX (Forkhead)", FOXD1="FOX (Forkhead)", FOXK1="FOX (Forkhead)", FOXK2="FOX (Forkhead)",
  FOXL1="FOX (Forkhead)", FOXN3="FOX (Forkhead)", FOXO4="FOX (Forkhead)", FOXO6="FOX (Forkhead)",
  FOXP1="FOX (Forkhead)", FOXP4="FOX (Forkhead)",
  KLF10="SP/KLF", KLF11="SP/KLF", KLF12="SP/KLF", KLF15="SP/KLF", KLF16="SP/KLF", KLF4="SP/KLF",
  KLF5="SP/KLF", KLF7="SP/KLF", SP1="SP/KLF", SP2="SP/KLF", SP3="SP/KLF", SP4="SP/KLF", SP8="SP/KLF", SP9="SP/KLF",
  TCF12="bHLH (E-protein)", TCF3="bHLH (E-protein)", TCF4="bHLH (E-protein)", TCFL5="bHLH (E-protein)",
  ASCL1="bHLH (Proneural)", NEUROG1="bHLH (Proneural)", NHLH1="bHLH (Proneural)", NHLH2="bHLH (Proneural)",
  PAX3="Homeodomain", PAX8="Homeodomain", POU2F1="Homeodomain", POU4F1="Homeodomain", POU5F1B="Homeodomain",
  CUX1="Homeodomain", CUX2="Homeodomain", ONECUT2="Homeodomain", ONECUT3="Homeodomain",
  MEIS1="Homeodomain", MEIS2="Homeodomain", PBX2="Homeodomain", OTX2="Homeodomain",
  RFX1="RFX", RFX2="RFX", RFX3="RFX", RFX4="RFX", RFX5="RFX", RFX7="RFX",
  ESRRA="Nuclear receptor", ESRRB="Nuclear receptor", NR1D1="Nuclear receptor", NR2C1="Nuclear receptor",
  NR2C2="Nuclear receptor", NR2F1="Nuclear receptor", NR2F2="Nuclear receptor", NR5A1="Nuclear receptor",
  NR6A1="Nuclear receptor", PPARD="Nuclear receptor", RARB="Nuclear receptor", RARG="Nuclear receptor",
  RXRB="Nuclear receptor", RXRG="Nuclear receptor", THRB="Nuclear receptor",
  BNC2="Zinc Finger", HINFP="Zinc Finger", MAZ="Zinc Finger", ZNF135="Zinc Finger", ZNF148="Zinc Finger",
  ZNF211="Zinc Finger", ZNF281="Zinc Finger", ZNF417="Zinc Finger", ZNF610="Zinc Finger", ZNF8="Zinc Finger",
  ZNF93="Zinc Finger", ZBED1="Zinc Finger", ZBED4="Zinc Finger", PATZ1="Zinc Finger",
  ZBTB14="Zinc Finger", ZBTB17="Zinc Finger",
  TBP="Basal / Chromatin", SATB1="Basal / Chromatin", DBP="Basal / Chromatin", TEF="Basal / Chromatin",
  NFYC="Basal / Chromatin"
)

ext_sig_tfs <- tobias_integrated %>% filter(binding_sig) %>% pull(TF_name) %>% unique()
ext_data <- tobias_integrated %>% filter(TF_name %in% ext_sig_tfs) %>% dplyr::select(TF_name, contrast, binding_change) %>%
  group_by(TF_name, contrast) %>% summarise(binding_change = mean(binding_change, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = contrast, values_from = binding_change, values_fill = NA)
if (nrow(ext_data) > 0 && ncol(ext_data) > 2) {
  ext_mat <- ext_data %>% column_to_rownames("TF_name") %>% as.matrix()

  ext_missing <- setdiff(rownames(ext_mat), names(extended_family))
  if (length(ext_missing) > 0) {
    cat(sprintf("[A07b] WARNING: %d TF(s) not in extended_family mapping -- shown as 'Other': %s\n",
                length(ext_missing), paste(ext_missing, collapse = ", ")))
  }
  ext_fam_vec <- extended_family[rownames(ext_mat)]
  ext_fam_vec[is.na(ext_fam_vec)] <- "Other"

  ext_sig_lookup <- tobias_integrated %>% distinct(TF_name, contrast, binding_sig)
  ext_star <- matrix("", nrow = nrow(ext_mat), ncol = ncol(ext_mat), dimnames = dimnames(ext_mat))
  for (i in seq_len(nrow(ext_mat))) {
    for (j in seq_len(ncol(ext_mat))) {
      tf <- rownames(ext_mat)[i]; ct_name <- colnames(ext_mat)[j]
      if (is.na(ext_mat[i, j])) next
      hit <- ext_sig_lookup %>% filter(TF_name == tf, contrast == ct_name)
      if (nrow(hit) == 1 && isTRUE(hit$binding_sig)) ext_star[i, j] <- "*"
    }
  }

  ext_row_order <- order(ext_fam_vec, rownames(ext_mat))
  ext_mat <- ext_mat[ext_row_order, , drop = FALSE]
  ext_star <- ext_star[ext_row_order, , drop = FALSE]
  ext_fam_vec <- ext_fam_vec[ext_row_order]
  ext_gaps_row <- which(diff(as.integer(factor(ext_fam_vec, levels = unique(ext_fam_vec)))) != 0)

  colnames(ext_mat) <- name_to_treatment[colnames(ext_mat)]
  colnames(ext_star) <- colnames(ext_mat)

  plot_family_heatmap(ext_mat, ext_star, ext_fam_vec,
                      sprintf("TOBIAS TF binding change, all %d significant TFs\n(* = significant in this contrast)", nrow(ext_mat)),
                      file.path(PLOTS_DIR, "A07b_heatmap_extended.pdf"), height = 24)
  cat(sprintf("[A07b] Extended heatmap: %d TFs across %d families -> A07b_heatmap_extended.pdf\n", nrow(ext_mat), length(unique(ext_fam_vec))))

  ########
  ## Strictest tier: TFs significant in BOTH contrasts (binding_sig TRUE for
  ## every column, not just >=1) -- every cell in this heatmap is starred by
  ## construction, so the star here is kept only for visual consistency with
  ## the other two heatmaps, not because it's doing any real filtering work.
  ########
  both_sig_tfs <- tobias_integrated %>% distinct(TF_name, contrast, binding_sig) %>%
    group_by(TF_name) %>% filter(all(binding_sig), n() == length(CONTRASTS)) %>% ungroup() %>% pull(TF_name) %>% unique()
  both_mat <- ext_mat[rownames(ext_mat) %in% both_sig_tfs, , drop = FALSE]
  both_star <- ext_star[rownames(ext_mat) %in% both_sig_tfs, , drop = FALSE]

  ## Exact user-supplied mapping for this specific 36-TF set (NOT derived
  ## from extended_family + a generic singleton-collapse rule -- that
  ## approach dumped PAX3 and KLF15 into a generic "Other" bucket; here
  ## KLF15 is reassigned to "Zinc Finger" (KLF proteins are structurally
  ## C2H2 zinc fingers, a defensible reclassification into an existing,
  ## non-singleton bucket rather than a catch-all) and PAX3 is kept as its
  ## own "Homeodomain" singleton rather than merged away. Checked to cover
  ## all 36 TFs exactly -- no fallback needed.
  both_sig_family <- c(
    FOS="AP-1 (Fos/Jun)", FOSL1="AP-1 (Fos/Jun)", FOSL2="AP-1 (Fos/Jun)", JDP2="AP-1 (Fos/Jun)",
    JUNB="AP-1 (Fos/Jun)", JUND="AP-1 (Fos/Jun)",
    BACH1="bZIP (BACH/Maf)", BACH2="bZIP (BACH/Maf)", NFE2="bZIP (BACH/Maf)",
    FOXA3="FOX (Forkhead)", FOXD1="FOX (Forkhead)", FOXK1="FOX (Forkhead)", FOXK2="FOX (Forkhead)",
    FOXL1="FOX (Forkhead)", FOXO6="FOX (Forkhead)", FOXP4="FOX (Forkhead)",
    ESRRA="Nuclear receptor", ESRRB="Nuclear receptor", NR2F2="Nuclear receptor", NR6A1="Nuclear receptor",
    RFX1="RFX", RFX2="RFX", RFX3="RFX", RFX5="RFX",
    SOX10="SOX (HMG-box)", SOX13="SOX (HMG-box)", SOX15="SOX (HMG-box)", SOX2="SOX (HMG-box)",
    SOX4="SOX (HMG-box)", SOX8="SOX (HMG-box)", SOX9="SOX (HMG-box)",
    BNC2="Zinc Finger", KLF15="Zinc Finger", ZBED4="Zinc Finger", ZNF610="Zinc Finger",
    PAX3="Homeodomain"
  )
  both_missing <- setdiff(rownames(both_mat), names(both_sig_family))
  if (length(both_missing) > 0) {
    cat(sprintf("[A07b] WARNING: %d TF(s) not in both_sig_family mapping -- shown as 'Other': %s\n",
                length(both_missing), paste(both_missing, collapse = ", ")))
  }
  both_fam_vec <- both_sig_family[rownames(both_mat)]
  both_fam_vec[is.na(both_fam_vec)] <- "Other"

  both_row_order <- order(both_fam_vec, rownames(both_mat))
  both_mat <- both_mat[both_row_order, , drop = FALSE]
  both_star <- both_star[both_row_order, , drop = FALSE]
  both_fam_vec <- both_fam_vec[both_row_order]
  if (nrow(both_mat) > 0) {
    plot_family_heatmap(both_mat, both_star, both_fam_vec,
                        sprintf("TOBIAS TF binding change, %d TFs significant in BOTH contrasts", nrow(both_mat)),
                        file.path(PLOTS_DIR, "A07b_heatmap_both_sig.pdf"), height = 9)
    cat(sprintf("[A07b] Both-contrasts-significant heatmap: %d TFs across %d families -> A07b_heatmap_both_sig.pdf\n",
                nrow(both_mat), length(unique(both_fam_vec))))

    ## Same heatmap, PLUS each TF's real JASPAR motif logo as a row
    ## annotation -- separate file, does not replace A07b_heatmap_both_sig.pdf.
    logo_paths <- render_tf_motif_logos(both_mat, TOBIAS_DIR, CONTRASTS)
    both_family_factor <- factor(both_fam_vec, levels = unique(both_fam_vec))
    both_col_fun <- colorRamp2(c(-max(abs(both_mat), na.rm = TRUE), 0, max(abs(both_mat), na.rm = TRUE)),
                                c("#4575B4", "white", "#D73027"))
    logo_anno <- rowAnnotation(
      Motif = anno_image(logo_paths[rownames(both_mat)], width = unit(1.6, "cm"), space = unit(1, "mm"), border = TRUE),
      annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 0)
    ht_motif <- Heatmap(both_mat, name = "Binding\nchange", col = both_col_fun,
                width = unit(3, "cm"),
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_split = both_family_factor, cluster_row_slices = FALSE,
                row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
                row_gap = unit(2, "mm"), border = TRUE,
                row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 11), column_names_rot = 45,
                column_title = sprintf("TF binding change, %d TFs significant in BOTH contrasts -- with JASPAR motif", nrow(both_mat)),
                column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                right_annotation = logo_anno)
    pdf(file.path(PLOTS_DIR, "A07b_heatmap_both_sig_with_motifs.pdf"), width = 6.3, height = 9)
    draw(ht_motif)
    dev.off()
    cat("[A07b] Both-contrasts-significant heatmap + motif logos -> A07b_heatmap_both_sig_with_motifs.pdf\n")
  }
}

saveRDS(tobias_integrated, file.path(RDS_DIR, "A07b_tobias_integrated.rds"))
writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "A07b_TOBIAS_downstream_session_info.txt"))

cat(sprintf("\n[DONE] A07b_TOBIAS_DOWNSTREAM complete. Concordant TFs: %d\n", nrow(concordant_tfs)))
