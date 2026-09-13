#!/usr/bin/env Rscript
# AR03_ATAC_Overlap.R -- compare the two genotypes' Loess-normalized ATAC
# differential accessibility: Venn diagrams, convergence scatter, genome-wide
# LFC correlation (+ permutation-null overlap-significance test, Fisher OR as
# effect size only), and a final text summary. Direct ATAC-side analogue of
# R04_RNA_Overlap.R -- same statistical design, including the same fix: a
# baseMean-binned permutation null replaces Fisher's own p-value for
# significance, since Fisher assumes region-independent Bernoulli trials,
# violated by co-accessible/co-regulated regions (same reasoning as R04's
# gene-independence violation).
#
# Both contrasts are tested on the IDENTICAL 146,769-peak consensus set
# (verified: A04b's two saved full-results tables have identical peak_id
# sets), so overlap is a direct peak_id intersection -- no genomic-range
# merge needed, the same simplification R04 relies on for genes.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR03_ATAC_Overlap.R
# Requires: A04b_Csaw_Loess_Norm.R and A05_DESeq2.R already run. Loess, not
# TMM: per the normalization-methodology decision in
# docs/A04b_normalization_methodology.md, Loess is the arm reported
# throughout the pipeline.
#
# Self-checkpointing: skips entirely if results/ATAC_ANALYSIS_SUMMARY.txt exists.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2)
  library(VennDiagram); library(grid)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR03] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/atac"
SESSION_DIR <- "results/Session_info"
SUMMARY     <- file.path("results", "ATAC_ANALYSIS_SUMMARY.txt")

if (file.exists(SUMMARY)) {
  cat("[AR03] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS <- cfg$contrasts
if (length(CONTRASTS) != 2) stop("[AR03] ERROR: AR03 assumes exactly 2 contrasts (pairwise Venn/scatter/correlation); config has ", length(CONTRASTS), ".")
NAME_A <- CONTRASTS[[1]]$name; NAME_B <- CONTRASTS[[2]]$name
LABEL_A <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment
PADJ_CUTOFF <- cfg$thresholds$fdr
LFC_THR     <- cfg$thresholds$lfc_threshold

full_a_path <- file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", NAME_A))
full_b_path <- file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", NAME_B))
if (!file.exists(full_a_path) || !file.exists(full_b_path)) stop("[AR03] ERROR: missing A04b output -- run A04b_Csaw_Loess_Norm.R first.")

COLORS <- list(A = "#E69F00", B = "#56B4E9", shared_up = "#1f78b4", shared_down = "#e31a1c")

res_a <- read_csv(full_a_path, show_col_types = FALSE)
res_b <- read_csv(full_b_path, show_col_types = FALSE)
stopifnot(setequal(res_a$peak_id, res_b$peak_id))
cat(sprintf("[AR03] %s: %d peaks | %s: %d peaks | shared consensus set: %s\n",
            NAME_A, nrow(res_a), NAME_B, nrow(res_b), setequal(res_a$peak_id, res_b$peak_id)))

########
## 1. Venn diagrams
########

make_venn <- function(list1, list2, name1, name2, title, filename) {
  if (length(list1) == 0 || length(list2) == 0) { warning("Empty list for Venn: ", title); return(invisible(NULL)) }
  overlap <- length(intersect(list1, list2))
  pdf(file.path(PLOTS_DIR, paste0(filename, ".pdf")), width = 6, height = 6)
  grid.newpage()
  VennDiagram::draw.pairwise.venn(
    area1 = length(list1), area2 = length(list2), cross.area = overlap,
    category = c(name1, name2), fill = c(COLORS$A, COLORS$B), alpha = 0.5,
    cat.pos = c(0, 0), cat.dist = c(0.05, 0.05), scaled = TRUE, lwd = 2, cat.cex = 1.2, cex = 1.5)
  grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
  intersect(list1, list2)
}

dar_ids <- function(res, lfc_thresh, direction) {
  res %>% filter(!is.na(padj), padj < PADJ_CUTOFF,
                 if (direction == "up") log2FoldChange > lfc_thresh else log2FoldChange < -lfc_thresh) %>%
    pull(peak_id)
}

shared_up   <- make_venn(dar_ids(res_a, LFC_THR, "up"),   dar_ids(res_b, LFC_THR, "up"),   LABEL_A, LABEL_B,
                          sprintf("Shared GAINED (|LFC| > %.2f)", LFC_THR), "AR03_venn_up")
shared_down <- make_venn(dar_ids(res_a, LFC_THR, "down"), dar_ids(res_b, LFC_THR, "down"), LABEL_A, LABEL_B,
                          sprintf("Shared LOST (|LFC| > %.2f)", LFC_THR), "AR03_venn_down")

########
## 2. Genome-wide LFC correlation + permutation-null overlap-significance test
##
## Same design/rationale as R04_RNA_Overlap.R Section 3 -- see that script's
## header comment for the full justification. Fisher's OR is kept as the
## effect size (a ratio of observed counts, much less sensitive to the
## region-independence violation than its own p-value); the permutation null
## is the significance test.
########

lfc_merged <- inner_join(
  res_a %>% transmute(peak_id, A_lfc = log2FoldChange, A_padj = padj, baseMean),
  res_b %>% transmute(peak_id, B_lfc = log2FoldChange, B_padj = padj),
  by = "peak_id") %>%
  mutate(sig_a = !is.na(A_padj) & A_padj < PADJ_CUTOFF & abs(A_lfc) > LFC_THR,
         sig_b = !is.na(B_padj) & B_padj < PADJ_CUTOFF & abs(B_lfc) > LFC_THR,
         category = case_when(
           sig_a & sig_b & sign(A_lfc) == sign(B_lfc) ~ "Convergent",
           sig_a & sig_b & sign(A_lfc) != sign(B_lfc) ~ "Divergent",
           sig_a & !sig_b ~ paste(LABEL_A, "only"),
           !sig_a & sig_b ~ paste(LABEL_B, "only"),
           TRUE ~ "NS"))

n_total <- nrow(lfc_merged); n_a <- sum(lfc_merged$sig_a); n_b <- sum(lfc_merged$sig_b); n_both <- sum(lfc_merged$sig_a & lfc_merged$sig_b)
fisher_res <- fisher.test(matrix(c(n_both, n_a - n_both, n_b - n_both, n_total - n_a - n_b + n_both), nrow = 2))
cor_all <- cor.test(lfc_merged$A_lfc, lfc_merged$B_lfc, method = "spearman")

N_PERM <- 10000
N_BASEMEAN_BINS <- 20
set.seed(1)  # fixed seed: permutation p-value must be reproducible run-to-run
bin <- as.integer(cut(rank(lfc_merged$baseMean, ties.method = "first"), breaks = N_BASEMEAN_BINS, labels = FALSE))
bin_idx <- split(seq_len(n_total), bin)
sig_a_vec <- lfc_merged$sig_a; sig_b_vec <- lfc_merged$sig_b
perm_overlaps <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  sig_b_perm <- sig_b_vec
  for (idx in bin_idx) sig_b_perm[idx] <- sample(sig_b_vec[idx])
  perm_overlaps[i] <- sum(sig_a_vec & sig_b_perm)
}
perm_p <- (1 + sum(perm_overlaps >= n_both)) / (N_PERM + 1)

cat(sprintf("[AR03] Genome-wide correlation (%d peaks): rho = %.4f, p = %s\n", n_total, cor_all$estimate, format.pval(cor_all$p.value)))
cat(sprintf("[AR03] Overlap: %s DARs=%d, %s DARs=%d, both=%d, OR=%.2f [%.2f, %.2f] (Fisher; effect size, not the significance test below)\n",
            LABEL_A, n_a, LABEL_B, n_b, n_both, fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2]))
cat(sprintf("[AR03] Overlap significance: permutation null (%d reps, %d baseMean bins), mean null overlap = %.1f, observed = %d, empirical p = %s\n",
            N_PERM, N_BASEMEAN_BINS, mean(perm_overlaps), n_both, format.pval(perm_p)))

p_corr <- ggplot(lfc_merged, aes(x = B_lfc, y = A_lfc, color = category)) +
  geom_point(data = filter(lfc_merged, category == "NS"), size = 0.3, alpha = 0.2) +
  geom_point(data = filter(lfc_merged, category != "NS"), size = 1.2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") + geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  scale_color_manual(values = setNames(c("grey85", COLORS$B, COLORS$A, "#2C3E50", "#E74C3C"),
                                        c("NS", paste(LABEL_B, "only"), paste(LABEL_A, "only"), "Convergent", "Divergent")), name = NULL) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("rho = %.3f\nn = %s", cor_all$estimate, format(n_total, big.mark = ",")),
           hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic") +
  labs(x = sprintf("%s log2FC vs WT (ATAC)", LABEL_B), y = sprintf("%s log2FC vs WT (ATAC)", LABEL_A),
       title = sprintf("Genome-wide accessibility correlation: %s vs %s", LABEL_A, LABEL_B)) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.15, 0.85), panel.grid.minor = element_blank())
ggsave(file.path(PLOTS_DIR, "AR03_genome_wide_lfc_correlation.pdf"), p_corr, width = 8, height = 8)

########
## 3. Convergence scatter, matching R04's Section 2 style
########

merged_data <- res_b %>%
  transmute(peak_id, LFC_B = log2FoldChange, P_B = padj) %>%
  inner_join(res_a %>% transmute(peak_id, LFC_A = log2FoldChange, P_A = padj), by = "peak_id") %>%
  mutate(Sig_Score = P_A * P_B,
         Category = case_when(
           (P_A < PADJ_CUTOFF & P_B < PADJ_CUTOFF) & (LFC_A > LFC_THR & LFC_B > LFC_THR) ~ "Shared Gained",
           (P_A < PADJ_CUTOFF & P_B < PADJ_CUTOFF) & (LFC_A < -LFC_THR & LFC_B < -LFC_THR) ~ "Shared Lost",
           TRUE ~ "Other"))

p_scatter <- ggplot(merged_data, aes(x = LFC_B, y = LFC_A)) +
  geom_point(data = subset(merged_data, Category == "Other"), color = "grey90", size = 0.5, alpha = 0.5) +
  geom_point(data = subset(merged_data, Category == "Shared Lost"), color = COLORS$shared_down, size = 1.5, alpha = 0.5) +
  geom_point(data = subset(merged_data, Category == "Shared Gained"), color = COLORS$shared_up, size = 1.5, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "grey50") + geom_hline(yintercept = 0, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(title = sprintf("Chromatin Accessibility Convergence: %s vs %s", LABEL_A, LABEL_B),
       subtitle = sprintf("Shared DARs: %d gained, %d lost (padj < %.2f, |LFC| > %.2f)",
                           sum(merged_data$Category == "Shared Gained"), sum(merged_data$Category == "Shared Lost"), PADJ_CUTOFF, LFC_THR),
       x = sprintf("log2 Fold Change (%s vs WT)", LABEL_B), y = sprintf("log2 Fold Change (%s vs WT)", LABEL_A))
ggsave(file.path(PLOTS_DIR, "AR03_scatter_convergence.pdf"), p_scatter, width = 8, height = 8)

########
## 4. Save overlap peak lists
########

save_overlap_list <- function(peak_list, name) {
  if (length(peak_list) == 0) return(invisible(NULL))
  merged_data %>% filter(peak_id %in% peak_list) %>% arrange(Sig_Score) %>%
    dplyr::select(peak_id, LFC_A, P_A, LFC_B, P_B, Sig_Score) %>%
    write_csv(file.path(TABLES_DIR, sprintf("AR03_overlap_%s_peaks.csv", name)))
}
save_overlap_list(shared_up, "gained")
save_overlap_list(shared_down, "lost")

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR03_ATAC_overlap_session_info.txt"))

########
## 5. Final text summary
########

sink(SUMMARY)
cat("================================================================================\n")
cat("ATAC-SEQ DIFFERENTIAL ACCESSIBILITY OVERLAP SUMMARY (Loess-normalized)\n")
cat("Project: KS1_2_Neuro-ectoderm\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

dar_summary_path <- file.path(TABLES_DIR, "A05_DAR_summary.csv")
if (file.exists(dar_summary_path)) {
  dar_summary <- read_csv(dar_summary_path, show_col_types = FALSE)
  cat("DAR counts by comparison and threshold (Loess):\n\n"); print(dar_summary); cat("\n")
}

cat("Overlap between genotypes:\n")
cat(sprintf("  Shared GAINED (|LFC| > %.2f): %d peaks\n", LFC_THR, length(shared_up)))
cat(sprintf("  Shared LOST (|LFC| > %.2f): %d peaks\n\n", LFC_THR, length(shared_down)))
cat(sprintf("Genome-wide correlation: rho = %.4f, p = %s\n", cor_all$estimate, format.pval(cor_all$p.value)))
cat(sprintf("DAR overlap effect size: OR = %.2f [%.2f, %.2f] (Fisher's exact test)\n", fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2]))
cat(sprintf("DAR overlap significance: permutation null, %d reps, baseMean-binned; observed = %d, empirical p = %s\n", N_PERM, n_both, format.pval(perm_p)))
cat("  (Fisher's own p-value is not reported as the significance test: it assumes\n")
cat("  region-independent Bernoulli trials, which co-accessible regions violate, making\n")
cat("  it anti-conservative. The OR is a ratio of observed counts and far less sensitive\n")
cat("  to that violation, so it is kept as the effect size; the permutation p above is\n")
cat("  the significance test.)\n\n")

cat("Analysis parameters:\n")
cat(sprintf("  Adjusted p-value: < %.2f\n", PADJ_CUTOFF))
cat(sprintf("  Log2 fold-change: |log2FC| > %.2f\n", LFC_THR))
cat("  Normalization: Loess (csaw normOffsets); see docs/A04b_normalization_methodology.md\n\n")
cat("================================================================================\n")
cat("END OF SUMMARY\n")
sink()

cat("\n[DONE] AR03_ATAC_Overlap complete -- summary: ", SUMMARY, "\n", sep = "")
