#!/usr/bin/env Rscript
# A05_DESeq2.R -- standalone reporting layer for the ATAC differential-accessibility
# results (Loess-normalized arm, A04b) -- a PCA of the peak-level accessibility
# landscape, a sample-distance heatmap, volcano plots at both LFC tiers, tiered
# sig tables, and an MA plot -- mirroring what R02_DESeq2.R already gives the
# RNA arm, so the two assay arms are reported at the same depth.
#
# Loess (A04b), not TMM (A04), per the normalization-methodology decision in
# docs/A04b_normalization_methodology.md: Loess is the arm actually reported
# throughout the pipeline (A08, AR01 downstream); TMM (A04) is kept solely as
# the comparison baseline that justifies preferring Loess, not fed downstream.
# A04b already runs DESeq2 internally (estimateDispersions + nbinomWaldTest) --
# A05 is not a second, independent fit, it loads A04b's own dds_loess object
# directly, same precedent A05 always followed for A04's DEdata (see the
# count-source bug note below, still relevant: never rebuild counts from
# dba.peakset()).
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/A05_DESeq2.R
#
# Self-checkpointing: skips entirely if results/RDS/A05_dds_object.rds exists.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A05] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

METADATA    <- "data/sample_metadata.csv"
LOESS_RDS   <- "data/diffbind_output/CsawLoess_Norm/dds_loess.rds"
PLOTS_DIR   <- "results/atac"
TABLES_DIR  <- "results/tables"
RDS_DIR     <- "results/RDS"
SESSION_DIR <- "results/Session_info"
DDS_RDS     <- file.path(RDS_DIR, "A05_dds_object.rds")

for (d in c(PLOTS_DIR, TABLES_DIR, RDS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(DDS_RDS)) {
  cat("[A05] Already complete (", DDS_RDS, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!file.exists(METADATA)) stop("[A05] ERROR: ", METADATA, " not found.")
if (!file.exists(LOESS_RDS)) stop("[A05] ERROR: ", LOESS_RDS, " not found -- run A04b_Csaw_Loess_Norm.R first.")

PADJ_CUTOFF <- cfg$thresholds$fdr
LFC_THR     <- cfg$thresholds$lfc_threshold
CONTRASTS   <- cfg$contrasts

COLORS <- setNames(c("#56B4E9", "#E69F00"),
                    c(CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment))

theme_pub <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "right")

cat("================================================================================\n")
cat("A05_DESeq2 -- ATAC standalone reporting (Loess-normalized arm, A04b)\n")
cat("================================================================================\n\n")

########
## Reuse A04b's own fitted DESeqDataSet directly (dds_loess) -- already
## dispersion-fit and Wald-tested there (estimateDispersions(fitType="local")
## + nbinomWaldTest, matching the fitType DiffBind itself uses for the TMM arm,
## so the two remain comparable on identical fitting machinery). No second fit
## needed or wanted. Peak IDs are already "chr:start-end" as rownames (A04b
## sets these directly from dba.peakset()'s GRanges at construction time), so
## unlike the old TMM-sourced version of this script, no separate peaks_gr
## re-derivation/alignment step is needed here.
########

cat("[A05] Loading Loess-normalized DESeqDataSet:", LOESS_RDS, "\n")
dds <- readRDS(LOESS_RDS)
sample_ids <- colnames(dds)
dds$genotype <- factor(dds$Condition, levels = c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment))
cat(sprintf("[A05] DESeqDataSet (from A04b): %d peaks x %d samples\n\n", nrow(dds), ncol(dds)))

########
## Samplesheet directly from sample_metadata.csv, for the printed sample table only
## -- dds itself is already ordered/labeled from A04b, not rebuilt from this.
########

meta      <- read_csv(METADATA, show_col_types = FALSE)
atac_meta <- meta %>%
  filter(assay == "ATAC") %>%
  mutate(genotype = factor(genotype, levels = c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment))) %>%
  arrange(match(sample_id, sample_ids))

if (!identical(atac_meta$sample_id, sample_ids)) {
  stop("[A05] ERROR: sample_metadata.csv ATAC rows do not match the sample set in ", LOESS_RDS)
}
if (!identical(colnames(dds), atac_meta$sample_id)) stop("[A05] ERROR: dds_loess sample order does not match sample_metadata.csv.")

cat("Samples in analysis:\n")
print(atac_meta %>% dplyr::select(sample_id, genotype, replicate))
cat("\n")

norm_counts <- counts(dds, normalized = TRUE)
saveRDS(norm_counts, file.path(RDS_DIR, "A05_normalized_counts.rds"))

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(RDS_DIR, "A05_vst_counts.rds"))

########
## PCA + sample correlation
########

pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3.5, box.padding = 0.5, point.padding = 0.3) +
  scale_color_manual(values = c("WT" = "grey50", COLORS), name = "Genotype") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_pub + labs(title = "PCA - ATAC-seq peak accessibility")
ggsave(file.path(PLOTS_DIR, "A05_pca.pdf"), p_pca, width = 7, height = 5)

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- atac_meta$sample_id
colnames(sample_dist_matrix) <- atac_meta$sample_id
pdf(file.path(PLOTS_DIR, "A05_sample_correlation.pdf"), width = 8, height = 7)
pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distances")
dev.off()

########
## Differential accessibility, per contrast (config-driven) -- volcano + MA plots
########

extract_results <- function(dds, treatment, reference) {
  cat(sprintf("Analyzing: %s vs %s\n", treatment, reference))
  # "Condition", not "genotype" -- resultsNames(dds) reflects DiffBind's original
  # ~Condition design baked in at fit time; the dds$genotype alias column (added
  # above for PCA/plotting readability) doesn't change that.
  res <- results(dds, contrast = c("Condition", treatment, reference))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("peak_id") %>%
    left_join(as.data.frame(norm_counts) %>% rownames_to_column("peak_id"), by = "peak_id") %>%
    arrange(padj)

  sig <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_THR)
  cat(sprintf("  DARs (|log2FC| > %.2f): %d (Up: %d, Down: %d)\n\n", LFC_THR, nrow(sig),
              sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
  list(full = res_df, sig = sig)
}

make_volcano <- function(res_df, name, lfc_thresh, suffix) {
  res_plot <- res_df %>%
    mutate(sig_category = case_when(
             is.na(padj) | padj >= PADJ_CUTOFF | abs(log2FoldChange) < lfc_thresh ~ "NS",
             log2FoldChange > lfc_thresh ~ "Up",
             log2FoldChange < -lfc_thresh ~ "Down",
             TRUE ~ "NS"),
           neg_log10_padj = -log10(padj))
  n_up <- sum(res_plot$sig_category == "Up", na.rm = TRUE)
  n_down <- sum(res_plot$sig_category == "Down", na.rm = TRUE)

  p <- ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = sig_category), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#e31a1c", "Down" = "#56B4E9", "NS" = "grey70"),
                        name = "Regulation", breaks = c("Up", "Down", "NS")) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey30") +
    theme_pub +
    labs(title = name, subtitle = sprintf("padj < %.2f, |log2FC| > %.2f | Up: %d, Down: %d", PADJ_CUTOFF, lfc_thresh, n_up, n_down),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  ggsave(file.path(PLOTS_DIR, paste0("A05_volcano_", suffix, ".pdf")), p, width = 8, height = 7)
}

theme_pub_rawp <- theme_minimal(base_size = 15) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "right")

# y only is capped (x is left uncapped, but both contrasts share the same x_limit
# so the two plots are directly comparable); shape (not text) marks capped points
# -- the wording for what's capped/shared belongs in the figure caption, not on
# the plot itself.
CAP_Y_RAWP <- 25

make_volcano_rawp <- function(res_df, name, lfc_thresh, suffix, x_limit) {
  # y-axis is the raw (unadjusted) p-value, but significance coloring still comes
  # from padj -- BH-adjusted significance has no fixed raw-p cutoff line (the BH
  # threshold is rank-dependent, i.e. i/m*alpha), so no geom_hline here, unlike
  # make_volcano()'s padj-axis plot where the padj cutoff is a single fixed line.
  res_plot <- res_df %>%
    mutate(sig_category = case_when(
             is.na(padj) | padj >= PADJ_CUTOFF | abs(log2FoldChange) < lfc_thresh ~ "NS",
             log2FoldChange > lfc_thresh ~ "Up",
             log2FoldChange < -lfc_thresh ~ "Down",
             TRUE ~ "NS"),
           neg_log10_pvalue = -log10(pvalue),
           y_capped = neg_log10_pvalue > CAP_Y_RAWP,
           neg_log10_pvalue = pmin(neg_log10_pvalue, CAP_Y_RAWP))
  n_up <- sum(res_plot$sig_category == "Up", na.rm = TRUE)
  n_down <- sum(res_plot$sig_category == "Down", na.rm = TRUE)

  p <- ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
    geom_point(aes(color = sig_category, shape = y_capped), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#e31a1c", "Down" = "#56B4E9", "NS" = "grey70"), name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey30") +
    theme_pub_rawp +
    labs(title = name,
         subtitle = sprintf("padj < %.2f, |log2FC| > %.2f | Up: %d, Down: %d", PADJ_CUTOFF, lfc_thresh, n_up, n_down),
         x = "log2 Fold Change", y = "-log10(p-value)")
  ggsave(file.path(PLOTS_DIR, paste0("A05_volcano_rawp_", suffix, ".pdf")), p, width = 5, height = 3.5)
  ggsave(file.path(PLOTS_DIR, paste0("A05_volcano_rawp_", suffix, ".png")), p, width = 5, height = 3.5, dpi = 300)
}

make_ma_plot <- function(res_df, name, suffix) {
  res_plot <- res_df %>%
    mutate(sig = !is.na(padj) & padj < PADJ_CUTOFF,
           mean_count = rowMeans(dplyr::select(., all_of(atac_meta$sample_id))))
  p <- ggplot(res_plot, aes(x = mean_count, y = log2FoldChange)) +
    geom_point(aes(color = sig), alpha = 0.5, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c("TRUE" = "#e31a1c", "FALSE" = "grey70"), name = sprintf("padj < %.2f", PADJ_CUTOFF)) +
    geom_hline(yintercept = 0, color = "grey30") +
    theme_pub +
    labs(title = name, subtitle = "MA plot", x = "Mean normalized counts", y = "log2 Fold Change")
  ggsave(file.path(PLOTS_DIR, paste0("A05_MA_", suffix, ".pdf")), p, width = 7, height = 5)
}

results_by_contrast <- list()
for (ct in CONTRASTS) {
  res <- extract_results(dds, ct$treatment, ct$reference)
  results_by_contrast[[ct$name]] <- res

  write_csv(res$full, file.path(TABLES_DIR, sprintf("A05_%s_full.csv", ct$name)))
  write_csv(res$sig, file.path(TABLES_DIR, sprintf("A05_%s_sig.csv", ct$name)))

  make_volcano(res$full, ct$name, LFC_THR, ct$name)
  make_ma_plot(res$full, ct$name, ct$name)
}

# Shared x-axis limit across both contrasts so the rawp volcano plots are directly
# comparable (x is left uncapped, unlike y, so this is the only thing keeping the
# two plots' scales aligned).
x_limit_rawp <- max(abs(unlist(lapply(results_by_contrast, function(r) r$full$log2FoldChange))), na.rm = TRUE)
for (ct in CONTRASTS) {
  make_volcano_rawp(results_by_contrast[[ct$name]]$full, ct$name, LFC_THR, ct$name, x_limit_rawp)
}

summary_df <- bind_rows(lapply(names(results_by_contrast), function(nm) {
  r <- results_by_contrast[[nm]]
  tibble(Comparison = nm,
         Threshold = sprintf("|log2FC| > %.2f", LFC_THR),
         Total_DARs = nrow(r$sig),
         Upregulated = sum(r$sig$log2FoldChange > 0),
         Downregulated = sum(r$sig$log2FoldChange < 0))
}))
write_csv(summary_df, file.path(TABLES_DIR, "A05_DAR_summary.csv"))
cat("\n[A05] DAR summary:\n"); print(summary_df)

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "A05_DESeq2_session_info.txt"))

# Checkpoint marker written last, once every output above has succeeded -- writing
# it earlier (as an earlier version of this script did) meant a crash partway
# through left DDS_RDS present but tables/plots missing, and the next run would
# silently skip everything as "already complete".
saveRDS(dds, DDS_RDS)
cat("\n[DONE] A05_DESeq2 complete\n")
