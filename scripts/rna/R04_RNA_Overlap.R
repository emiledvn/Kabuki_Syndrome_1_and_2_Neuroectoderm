#!/usr/bin/env Rscript
# R04_RNA_Overlap.R -- compare the two genotypes' RNA differential expression:
# Venn diagrams, convergence scatter, genome-wide LFC correlation (+ Fisher's
# exact test for overlap significance), concordant-DEG heatmap, GO enrichment
# on shared genes, and a final text summary report.
#
# Adapted from data_import_local/local_scripts/05_RNA_overlap_analysis.R, the
# overlap-heatmap part of 04_RNA_downstream_analysis.R, and
# 06_RNA_generate_summary.R (folded in at the end -- there's no separate "R05"
# in the pipeline structure).
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R04_RNA_Overlap.R
#
# Fix applied vs the original: GO enrichment on shared/overlap genes now passes
# an explicit `universe` (same background used in R02), not the clusterProfiler
# default of the whole genome.
#
# Self-checkpointing: skips entirely if results/tables/RNA_ANALYSIS_SUMMARY.txt exists.
# Requires R02_DESeq2.R to have completed first (only 2 contrasts supported --
# Venn/scatter/correlation are inherently pairwise).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(VennDiagram); library(grid)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(DESeq2)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R04] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/rna"
RDS_DIR     <- "results/RDS"
SESSION_DIR <- "results/Session_info"
SUMMARY     <- file.path("results", "RNA_ANALYSIS_SUMMARY.txt")

if (file.exists(SUMMARY)) {
  cat("[R04] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS <- cfg$contrasts
if (length(CONTRASTS) != 2) stop("[R04] ERROR: R04 assumes exactly 2 contrasts (pairwise Venn/scatter/correlation); config has ", length(CONTRASTS), ".")
NAME_A <- CONTRASTS[[1]]$name; NAME_B <- CONTRASTS[[2]]$name
LABEL_A <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment
PADJ_CUTOFF <- cfg$thresholds$fdr
LFC_THR     <- cfg$thresholds$lfc_threshold

full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
if (!file.exists(full_a_path) || !file.exists(full_b_path)) stop("[R04] ERROR: missing R02 output -- run R02_DESeq2.R first.")

# shared_up/shared_down match the Up/Down/NS convention used everywhere DEG
# significance is colour-coded in this pipeline (R02/A05 volcano plots, AR07's
# small volcano panels): Up = red, Down = blue. (COLORS$A/B, for the Venn
# diagrams, are a separate genotype-identity palette, not a significance one.)
COLORS <- list(A = "#E69F00", B = "#56B4E9", shared_up = "#e31a1c", shared_down = "#56B4E9")

res_a <- read_csv(full_a_path, show_col_types = FALSE)
res_b <- read_csv(full_b_path, show_col_types = FALSE)

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

deg_ids <- function(res, lfc_thresh, direction) {
  res %>% filter(!is.na(padj), padj < PADJ_CUTOFF,
                 if (direction == "up") log2FoldChange > lfc_thresh else log2FoldChange < -lfc_thresh) %>%
    pull(gene_id)
}

shared_up   <- make_venn(deg_ids(res_a, LFC_THR, "up"),   deg_ids(res_b, LFC_THR, "up"),   LABEL_A, LABEL_B,
                          sprintf("Shared Upregulated (|LFC| > %.2f)", LFC_THR), "R04_venn_up")
shared_down <- make_venn(deg_ids(res_a, LFC_THR, "down"), deg_ids(res_b, LFC_THR, "down"), LABEL_A, LABEL_B,
                          sprintf("Shared Downregulated (|LFC| > %.2f)", LFC_THR), "R04_venn_down")

########
## 2. Convergence scatter
########

merged_data <- res_b %>%
  dplyr::select(gene_id, gene_name, LFC_B = log2FoldChange, P_B = padj) %>%
  inner_join(res_a %>% dplyr::select(gene_id, LFC_A = log2FoldChange, P_A = padj), by = "gene_id") %>%
  mutate(Unique_Name = coalesce(gene_name, gene_id),
         Sig_Score = P_A * P_B,
         Category = case_when(
           (P_A < PADJ_CUTOFF & P_B < PADJ_CUTOFF) & (LFC_A > LFC_THR & LFC_B > LFC_THR) ~ "Shared Up",
           (P_A < PADJ_CUTOFF & P_B < PADJ_CUTOFF) & (LFC_A < -LFC_THR & LFC_B < -LFC_THR) ~ "Shared Down",
           TRUE ~ "Other"))

# Labels: only genes with a real symbol. Two distinct ways a gene can lack
# one, both excluded here: gene_name is NA (Unique_Name's coalesce() above
# falls back to the bare gene_id for these -- fine for the CSV export via
# save_overlap_list below, not for an axis label); and GENCODE v47 itself
# sets gene_name equal to the bare Ensembl ID (no version suffix) for ~36k
# under-characterised loci (mostly lncRNAs) -- confirmed directly in the GTF,
# e.g. ENSG00000254277 carries gene_name "ENSG00000254277" verbatim, so
# !is.na() alone doesn't catch it. Filtered before ranking, not after, so a
# real top-50 by Sig_Score isn't displaced by an unlabelable gene.
top_labels <- bind_rows(
  merged_data %>% filter(Category == "Shared Down", !is.na(gene_name), !str_detect(gene_name, "^ENSG\\d")) %>% arrange(Sig_Score) %>% head(50),
  merged_data %>% filter(Category == "Shared Up", !is.na(gene_name), !str_detect(gene_name, "^ENSG\\d")) %>% arrange(Sig_Score) %>% head(50))

# Correlation: Spearman (not Pearson) on the full gene set behind THIS plot
# (merged_data, gene_id inner join -- not lfc_merged from Section 3 below,
# which joins on gene_name+gene_id too and so covers a very slightly
# different gene set) -- rank-based, so it isn't distorted by the handful of
# extreme-LFC low-count DESeq2 outliers this dataset has (see AR07 Task 17's
# +/-22 outliers), and it's computed on the FULL data, not the axis-clipped
# window below, since coord_cartesian only changes what's drawn, not what
# the reported statistic should describe. Matches the method already used
# for the genome-wide correlation in Section 3.
scatter_cor <- cor.test(merged_data$LFC_A, merged_data$LFC_B, method = "spearman")

# Axis fixed at -10/7: a compromise between the full data range (+/-~24, a
# handful of near-zero-count DESeq2 outliers) and clipping so tight the
# labeled genes (near the origin) get squeezed.
axis_lo <- -10
axis_hi <- 7
p_scatter <- ggplot(merged_data, aes(x = LFC_B, y = LFC_A)) +
  geom_point(data = subset(merged_data, Category == "Other"), color = "grey70", size = 0.5, alpha = 0.5) +
  geom_point(data = subset(merged_data, Category == "Shared Down"), color = COLORS$shared_down, size = 1.5, alpha = 0.5) +
  geom_point(data = subset(merged_data, Category == "Shared Up"), color = COLORS$shared_up, size = 1.5, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "grey50") + geom_hline(yintercept = 0, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_labels, aes(label = gene_name), box.padding = 0.4, max.overlaps = 100,
                   min.segment.length = 0, size = 2.6, fontface = "italic", color = "grey20",
                   segment.size = 0.3, segment.color = "grey60") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 3.8, fontface = "italic",
           label = sprintf("Spearman rho = %.3f\np = %s, n = %s genes",
                            scatter_cor$estimate, format.pval(scatter_cor$p.value, digits = 2), format(nrow(merged_data), big.mark = ","))) +
  coord_cartesian(xlim = c(axis_lo, axis_hi), ylim = c(axis_lo, axis_hi)) +
  theme_minimal(base_size = 14) +
  labs(title = sprintf("Transcriptional Convergence: %s vs %s", LABEL_A, LABEL_B),
       subtitle = sprintf("Shared DEGs: %d up, %d down (padj < %.2f, |LFC| > %.2f)",
                           sum(merged_data$Category == "Shared Up"), sum(merged_data$Category == "Shared Down"), PADJ_CUTOFF, LFC_THR),
       x = sprintf("log2 Fold Change (%s vs WT)", LABEL_B), y = sprintf("log2 Fold Change (%s vs WT)", LABEL_A))
ggsave(file.path(PLOTS_DIR, "R04_scatter_convergence.pdf"), p_scatter, width = 8, height = 8)

########
## 3. Genome-wide LFC correlation + overlap-significance test
##
## Fisher's exact test (kept below for its odds ratio) assumes genes are
## independent Bernoulli trials, which is false -- genes co-vary in regulatory
## modules -- so with ~n_total genes tested its p-value is badly anti-conservative
## (comes out ~1e-266, not a number to put in a methods section as-is). The OR is
## a ratio of observed counts and much less sensitive to that independence
## violation than the p-value is, so it stays as the reported effect size. For
## significance, a permutation null is used instead: shuffle which genes are
## "sig_b" within baseMean-quantile bins (preserving the real relationship between
## expression level and power to reach significance, which an unbinned shuffle
## would ignore) and rebuild the overlap-count null empirically.
########

lfc_merged <- inner_join(
  res_a %>% transmute(gene_name, gene_id, A_lfc = log2FoldChange, A_padj = padj, baseMean),
  res_b %>% transmute(gene_name, gene_id, B_lfc = log2FoldChange, B_padj = padj),
  by = c("gene_name", "gene_id")) %>%
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

cat(sprintf("[R04] Genome-wide correlation (%d genes): rho = %.4f, p = %s\n", n_total, cor_all$estimate, format.pval(cor_all$p.value)))
cat(sprintf("[R04] Overlap: %s DEGs=%d, %s DEGs=%d, both=%d, OR=%.2f [%.2f, %.2f] (Fisher; effect size, not the significance test below)\n",
            LABEL_A, n_a, LABEL_B, n_b, n_both, fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2]))
cat(sprintf("[R04] Overlap significance: permutation null (%d reps, %d baseMean bins), mean null overlap = %.1f, observed = %d, empirical p = %s\n",
            N_PERM, N_BASEMEAN_BINS, mean(perm_overlaps), n_both, format.pval(perm_p)))

# Directional breakdown (for AR07's directional-concordance-summary barplot):
# same sig_a/sig_b/category classification as above, just split by the sign of
# each gene's own LFC instead of collapsed into Convergent/Divergent/only/NS.
concordance_counts <- bind_rows(
  lfc_merged %>% filter(category == "Convergent") %>%
    transmute(category = "Shared", direction = if_else(A_lfc > 0, "Up", "Down")),
  lfc_merged %>% filter(category == "Divergent") %>%
    transmute(category = "Discordant", direction = NA_character_),
  lfc_merged %>% filter(category == paste(LABEL_A, "only")) %>%
    transmute(category = paste(LABEL_A, "only"), direction = if_else(A_lfc > 0, "Up", "Down")),
  lfc_merged %>% filter(category == paste(LABEL_B, "only")) %>%
    transmute(category = paste(LABEL_B, "only"), direction = if_else(B_lfc > 0, "Up", "Down"))) %>%
  count(category, direction, name = "n")
write_csv(concordance_counts, file.path(TABLES_DIR, "R04_directional_concordance_counts.csv"))

write_csv(tibble(OR = fisher_res$estimate, OR_CI_low = fisher_res$conf.int[1], OR_CI_high = fisher_res$conf.int[2],
                  perm_p = perm_p, n_perm = N_PERM, n_basemean_bins = N_BASEMEAN_BINS,
                  n_total = n_total, n_a = n_a, n_b = n_b, n_both = n_both),
          file.path(TABLES_DIR, "R04_overlap_effect_size.csv"))

p_corr <- ggplot(lfc_merged, aes(x = B_lfc, y = A_lfc, color = category)) +
  geom_point(data = filter(lfc_merged, category == "NS"), size = 0.3, alpha = 0.2) +
  geom_point(data = filter(lfc_merged, category != "NS"), size = 1.2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") + geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  scale_color_manual(values = setNames(c("grey85", COLORS$B, COLORS$A, "#2C3E50", "#E74C3C"),
                                        c("NS", paste(LABEL_B, "only"), paste(LABEL_A, "only"), "Convergent", "Divergent")), name = NULL) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("rho = %.3f\nn = %s", cor_all$estimate, format(n_total, big.mark = ",")),
           hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic") +
  labs(x = sprintf("%s log2FC vs WT", LABEL_B), y = sprintf("%s log2FC vs WT", LABEL_A),
       title = sprintf("Genome-wide transcriptional correlation: %s vs %s", LABEL_A, LABEL_B)) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.15, 0.85), panel.grid.minor = element_blank())
ggsave(file.path(PLOTS_DIR, "R04_genome_wide_lfc_correlation.pdf"), p_corr, width = 8, height = 8)

########
## 4. Concordant-DEG heatmap (top 200 up/down per genotype, intersected)
########

vsd <- readRDS(file.path(RDS_DIR, "R02_vst_counts.rds"))

top_ids <- function(res, direction) {
  res %>% filter(!is.na(padj), padj < PADJ_CUTOFF, if (direction == "up") log2FoldChange > LFC_THR else log2FoldChange < -LFC_THR) %>%
    arrange(padj) %>% slice_head(n = 200) %>% pull(gene_id)
}
overlap_blocks <- list(UP_UP = intersect(top_ids(res_a, "up"), top_ids(res_b, "up")),
                        DOWN_DOWN = intersect(top_ids(res_a, "down"), top_ids(res_b, "down")))
overlap_blocks <- Filter(function(x) length(x) >= 2, overlap_blocks)

if (length(overlap_blocks) > 0) {
  gene_order <- unlist(overlap_blocks)
  block_labels <- setNames(rep(names(overlap_blocks), lengths(overlap_blocks)), gene_order)
  id_to_name <- bind_rows(res_a, res_b) %>% distinct(gene_id, gene_name) %>% deframe()

  genotype_order <- c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment)
  col_order <- order(match(colData(vsd)$genotype, genotype_order))
  mat <- assay(vsd)[gene_order, col_order]
  rownames(mat) <- id_to_name[rownames(mat)]
  mat <- mat - rowMeans(mat)

  row_ann <- data.frame(Block = block_labels, row.names = id_to_name[gene_order])
  col_ann <- data.frame(Genotype = colData(vsd)$genotype[col_order], row.names = colnames(mat))
  ann_colors <- list(Genotype = setNames(c("grey50", COLORS$A, COLORS$B), genotype_order),
                      Block = c(UP_UP = "#B40426", DOWN_DOWN = "#3B4CC0"))
  gaps_row <- cumsum(lengths(overlap_blocks))[-length(overlap_blocks)]

  pdf(file.path(PLOTS_DIR, "R04_heatmap_concordant_overlap.pdf"), width = 8, height = 12)
  pheatmap::pheatmap(mat, annotation_row = row_ann, annotation_col = col_ann, annotation_colors = ann_colors,
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                      cluster_rows = FALSE, cluster_cols = TRUE, gaps_row = gaps_row, border_color = "grey",
                      show_rownames = TRUE, fontsize_row = 8,
                      main = sprintf("Concordant DEGs - %s & %s vs WT (|LFC| > %.2f)", LABEL_A, LABEL_B, LFC_THR))
  dev.off()
}

########
## 5. GO enrichment on shared/overlap genes -- explicit universe, same as R02
########

gene_map <- readRDS(file.path(RDS_DIR, "R02_gene_id_name_map.rds"))
vst_gene_ids <- rownames(vsd)
universe_entrez <- tryCatch(
  bitr(gene_map$gene_name[gene_map$gene_id %in% vst_gene_ids], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  error = function(e) character(0))

gene_names <- res_a %>% dplyr::select(gene_id, gene_name) %>% deframe()

run_overlap_go <- function(gene_list, direction) {
  if (length(gene_list) < 10) return(NULL)
  symbols <- na.omit(gene_names[gene_list])
  entrez <- tryCatch(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
  if (length(entrez) < 10) return(NULL)
  ego <- enrichGO(gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF, readable = TRUE)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>% filter(p.adjust < PADJ_CUTOFF) %>% mutate(Direction = direction)
}

go_overlap <- bind_rows(run_overlap_go(shared_up, "Upregulated"), run_overlap_go(shared_down, "Downregulated"))
if (nrow(go_overlap) > 0) write_csv(go_overlap, file.path(TABLES_DIR, "R04_GO_overlap.csv"))

########
## 6. Save overlap gene lists
########

save_overlap_list <- function(gene_list, name) {
  if (length(gene_list) == 0) return(invisible(NULL))
  merged_data %>% filter(gene_id %in% gene_list) %>% arrange(Sig_Score) %>%
    dplyr::select(gene_id, gene_name, LFC_A, P_A, LFC_B, P_B, Sig_Score) %>%
    write_csv(file.path(TABLES_DIR, sprintf("R04_overlap_%s_genes.csv", name)))
}
save_overlap_list(shared_up, "up")
save_overlap_list(shared_down, "down")

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "R04_RNA_overlap_session_info.txt"))

########
## 7. Final text summary (from 06_RNA_generate_summary.R)
########

sink(SUMMARY)
cat("================================================================================\n")
cat("RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY\n")
cat("Project: KS1_2_Neuro-ectoderm\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

cat("Sample breakdown by genotype:\n")
meta <- read_csv("data/sample_metadata.csv", show_col_types = FALSE) %>% filter(assay == "RNA")
print(meta %>% dplyr::count(genotype))
cat("\n")

deg_summary <- read_csv(file.path(TABLES_DIR, "R02_DEG_summary.csv"), show_col_types = FALSE)
cat("DEG counts by comparison and threshold:\n\n"); print(deg_summary); cat("\n")

cat("Overlap between genotypes:\n")
cat(sprintf("  Shared upregulated (|LFC| > %.2f): %d genes\n", LFC_THR, length(shared_up)))
cat(sprintf("  Shared downregulated (|LFC| > %.2f): %d genes\n\n", LFC_THR, length(shared_down)))
cat(sprintf("Genome-wide correlation: rho = %.4f, p = %s\n", cor_all$estimate, format.pval(cor_all$p.value)))
cat(sprintf("DEG overlap effect size: OR = %.2f [%.2f, %.2f] (Fisher's exact test)\n", fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2]))
cat(sprintf("DEG overlap significance: permutation null, %d reps, baseMean-binned; observed = %d, empirical p = %s\n", N_PERM, n_both, format.pval(perm_p)))
cat("  (Fisher's own p-value is not reported as the significance test: it assumes\n")
cat("  gene-independent Bernoulli trials, which co-regulated genes violate, making it\n")
cat("  anti-conservative. The OR is a ratio of observed counts and far less sensitive\n")
cat("  to that violation, so it is kept as the effect size; the permutation p above is\n")
cat("  the significance test.)\n\n")

cat("Analysis parameters:\n")
cat(sprintf("  Adjusted p-value: < %.2f\n", PADJ_CUTOFF))
cat(sprintf("  Log2 fold-change: |log2FC| > %.2f\n", LFC_THR))
cat(sprintf("  Pre-filtering: genes with >= %d total counts\n", cfg$thresholds$rna_min_total_counts))
cat("  GO enrichment: Biological Process, BH-adjusted, explicit background universe\n\n")
cat("================================================================================\n")
cat("END OF SUMMARY\n")
sink()

cat("\n[DONE] R04_RNA_Overlap complete -- summary: ", SUMMARY, "\n", sep = "")
