#!/usr/bin/env Rscript
# R02_DESeq2.R -- RNA differential expression (DESeq2), PCA/QC, volcano plots,
# NE-identity marker validation, per-genotype heatmaps, and GO enrichment.
# Adapted from data_import_local/local_scripts/03_RNA_differential_expression.R and
# the per-genotype (non-overlap) parts of 04_RNA_downstream_analysis.R.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R02_DESeq2.R
#
# Fixes applied vs the originals:
#   - GO enrichment now passes an explicit `universe` (all genes that passed the
#     count pre-filter, i.e. everything DESeq2 actually tested) instead of
#     defaulting to the whole genome as background.
#   - Gene counts + samplesheet are built directly from data/RNA_nfcore_output and
#     data/sample_metadata.csv (single source of truth) instead of separately
#     maintained samplesheet_NE.csv / RNA_gene_counts_filtered.tsv files.
#
# Self-checkpointing: skips entirely if results/RDS/R02_dds_object.rds exists.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R02] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

METADATA   <- "data/sample_metadata.csv"
COUNTS     <- "data/RNA_nfcore_output/star_salmon/salmon.merged.gene_counts.tsv"
GTF        <- cfg$reference$gtf
PLOTS_DIR  <- "results/rna"
TABLES_DIR <- "results/tables"
RDS_DIR    <- "results/RDS"
SESSION_DIR <- "results/Session_info"
DDS_RDS    <- file.path(RDS_DIR, "R02_dds_object.rds")

for (d in c(PLOTS_DIR, TABLES_DIR, RDS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(DDS_RDS)) {
  cat("[R02] Already complete (", DDS_RDS, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!file.exists(METADATA)) stop("[R02] ERROR: ", METADATA, " not found.")
if (!file.exists(COUNTS))   stop("[R02] ERROR: ", COUNTS, " not found -- run R01_RUN_nfcore_RNA.sh first.")
if (!file.exists(GTF))      stop("[R02] ERROR: ", GTF, " not found.")

PADJ_CUTOFF   <- cfg$thresholds$fdr
LFC_THR       <- cfg$thresholds$lfc_threshold
MIN_TOTAL     <- cfg$thresholds$rna_min_total_counts
CONTRASTS     <- cfg$contrasts

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
cat("R02_DESeq2 -- RNA differential expression\n")
cat("================================================================================\n\n")

########
## Gene-ID -> gene-name map from our pinned GTF (single source of truth for the
## reference; strip Ensembl version suffixes since GTF gene_id and salmon/tximport
## gene_id both carry them but must match exactly for the join below).
########

cat("[R02] Building gene_id -> gene_name map from GTF...\n")
gtf_lines <- system(sprintf("zcat %s | awk -F'\t' '$3==\"gene\"{print $9}'", GTF), intern = TRUE)
gene_map <- tibble(raw = gtf_lines) %>%
  mutate(
    gene_id   = str_remove(str_match(raw, 'gene_id "([^"]+)"')[, 2], "\\.\\d+$"),
    gene_name = str_match(raw, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()
cat(sprintf("  %d gene_id -> gene_name mappings\n\n", nrow(gene_map)))
saveRDS(gene_map, file.path(RDS_DIR, "R02_gene_id_name_map.rds"))

########
## Samplesheet directly from sample_metadata.csv
########

meta     <- read_csv(METADATA, show_col_types = FALSE)
rna_meta <- meta %>%
  filter(assay == "RNA") %>%
  mutate(genotype = factor(genotype, levels = c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment)))

cat("Samples in analysis:\n")
print(rna_meta %>% dplyr::select(sample_id, genotype, rin))
cat("\n")

########
## Load nf-core/rnaseq gene counts, filter, build DESeq2 object
########

counts_raw <- read_tsv(COUNTS, show_col_types = FALSE)
id_col <- if ("gene_id" %in% names(counts_raw)) "gene_id" else names(counts_raw)[1]
counts_raw <- counts_raw %>% rename(gene_id = !!id_col) %>% mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

count_matrix <- counts_raw %>%
  dplyr::select(gene_id, all_of(rna_meta$sample_id)) %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round()
count_matrix <- count_matrix[, rna_meta$sample_id]  # enforce metadata order

n_before <- nrow(count_matrix)
count_matrix <- count_matrix[rowSums(count_matrix) >= MIN_TOTAL, ]
cat(sprintf("[R02] Pre-filter: %d / %d genes kept (>= %d total counts)\n\n", nrow(count_matrix), n_before, MIN_TOTAL))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = rna_meta, design = ~ genotype)

cat("[R02] Running DESeq2...\n")
dds <- DESeq(dds)
saveRDS(dds, DDS_RDS)

norm_counts <- counts(dds, normalized = TRUE)
saveRDS(norm_counts, file.path(RDS_DIR, "R02_normalized_counts.rds"))

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(RDS_DIR, "R02_vst_counts.rds"))

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
  theme_pub + labs(title = "PCA - RNA-seq samples")
ggsave(file.path(PLOTS_DIR, "R02_pca.pdf"), p_pca, width = 7, height = 5)

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- rna_meta$sample_id
colnames(sample_dist_matrix) <- rna_meta$sample_id
pdf(file.path(PLOTS_DIR, "R02_sample_correlation.pdf"), width = 8, height = 7)
pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distances")
dev.off()

########
## Differential expression, per contrast (config-driven)
########

# Deliberately no lfcShrink() here, despite apeglm/ashr being pinned in ks_1_2_r.yml.
# Tested: apeglm shrinkage was evaluated against this exact dds object and rejected
# for DEG calling. It over-suppresses low-baseMean neuroectoderm patterning/
# pluripotency TFs -- e.g. EMX1, LHX2, SIX3, FOXG1, POU5F1, HOXB1 in KDM6A_ko_vs_WT
# had raw log2FC -2.1 to -4.4 (padj as low as 0.001) collapsed to -0.05 to -0.24
# post-shrinkage, because TFs are typically low-abundance and apeglm's prior
# dominates the MLE at n=3 per genotype. That is the opposite of a noise correction
# for this dataset: TFs driving neuroectoderm identity are exactly the genes where a
# small absolute expression change is expected to be functionally consequential, so
# collapsing their LFC toward zero would suppress real signal, not artifacts. padj
# is unaffected by shrinkage either way (DESeq2 computes it from the MLE before any
# shrinkage step), so raw log2FoldChange + padj is used throughout for DEG calling.
# Any future rank-based analysis (GSEA/RRHO) should rank by the Wald statistic
# (results()$stat), not by LFC of either kind, for the same reason.
extract_results <- function(dds, treatment, reference) {
  cat(sprintf("Analyzing: %s vs %s\n", treatment, reference))
  res <- results(dds, contrast = c("genotype", treatment, reference))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_map, by = "gene_id") %>%
    left_join(as.data.frame(norm_counts) %>% rownames_to_column("gene_id"), by = "gene_id") %>%
    arrange(padj)

  sig <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_THR)
  cat(sprintf("  DEGs (|log2FC| > %.2f): %d (Up: %d, Down: %d)\n\n", LFC_THR, nrow(sig),
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
  top_genes <- res_plot %>% filter(sig_category %in% c("Up", "Down")) %>% arrange(padj) %>% head(20)

  p <- ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = sig_category), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#e31a1c", "Down" = "#56B4E9", "NS" = "grey70"),
                        name = "Regulation", breaks = c("Up", "Down", "NS")) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey30") +
    geom_text_repel(data = top_genes, aes(label = gene_name), size = 3, max.overlaps = 20) +
    theme_pub +
    labs(title = name, subtitle = sprintf("padj < %.2f, |log2FC| > %.2f | Up: %d, Down: %d", PADJ_CUTOFF, lfc_thresh, n_up, n_down),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  ggsave(file.path(PLOTS_DIR, paste0("R02_volcano_", suffix, ".pdf")), p, width = 8, height = 7)
}

results_by_contrast <- list()
for (ct in CONTRASTS) {
  res <- extract_results(dds, ct$treatment, ct$reference)
  results_by_contrast[[ct$name]] <- res

  write_csv(res$full, file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)))
  write_csv(res$sig, file.path(TABLES_DIR, sprintf("R02_%s_sig.csv", ct$name)))

  make_volcano(res$full, ct$name, LFC_THR, ct$name)
}

summary_df <- bind_rows(lapply(names(results_by_contrast), function(nm) {
  r <- results_by_contrast[[nm]]
  tibble(Comparison = nm,
         Threshold = sprintf("|log2FC| > %.2f", LFC_THR),
         Total_DEGs = nrow(r$sig),
         Upregulated = sum(r$sig$log2FoldChange > 0),
         Downregulated = sum(r$sig$log2FoldChange < 0))
}))
write_csv(summary_df, file.path(TABLES_DIR, "R02_DEG_summary.csv"))
cat("\n[R02] DEG summary:\n"); print(summary_df)

########
## Per-genotype heatmaps of top DEGs
########

make_deg_heatmap <- function(res_df, name, n_genes = 50, lfc_thresh = LFC_THR) {
  top_degs <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > lfc_thresh) %>%
    arrange(padj) %>% head(n_genes) %>% pull(gene_id)
  if (length(top_degs) == 0) { warning("No DEGs for ", name); return(invisible(NULL)) }

  mat <- assay(vsd)[top_degs, , drop = FALSE]
  gene_names <- res_df %>% filter(gene_id %in% top_degs) %>% dplyr::select(gene_id, gene_name) %>% deframe()
  rownames(mat) <- gene_names[rownames(mat)]
  mat_scaled <- t(scale(t(mat)))

  annotation_col <- data.frame(Genotype = colData(vsd)$genotype, row.names = colnames(mat))
  ann_colors <- list(Genotype = c("WT" = "grey50", COLORS))

  pdf(file.path(PLOTS_DIR, sprintf("R02_heatmap_%s_top%d.pdf", name, n_genes)), width = 8, height = 10)
  pheatmap(mat_scaled, annotation_col = annotation_col, annotation_colors = ann_colors,
           color = colorRampPalette(c("#4393c3", "white", "#d6604d"))(100), breaks = seq(-2, 2, length.out = 100),
           show_rownames = TRUE, fontsize_row = 8,
           main = sprintf("Top %d DEGs - %s (|LFC| > %.2f)", n_genes, name, lfc_thresh))
  dev.off()
}

for (ct in CONTRASTS) {
  make_deg_heatmap(results_by_contrast[[ct$name]]$full, ct$name, 50)
}

########
## GO enrichment, per genotype and direction -- with an explicit background
## universe (all genes that passed the count pre-filter, i.e. everything DESeq2
## actually tested), not the clusterProfiler default of the whole genome.
########

universe_entrez <- tryCatch(
  bitr(gene_map$gene_name[gene_map$gene_id %in% rownames(count_matrix)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  error = function(e) character(0)
)
cat(sprintf("\n[R02] GO background universe: %d Entrez IDs (all genes passing the count pre-filter)\n", length(universe_entrez)))

run_go <- function(gene_symbols, name, direction) {
  if (length(gene_symbols) < 10) { warning(sprintf("Too few genes (%d) for GO: %s %s", length(gene_symbols), name, direction)); return(NULL) }
  entrez_ids <- tryCatch(bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
  if (length(entrez_ids) < 10) return(NULL)
  ego <- enrichGO(gene = entrez_ids, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF, readable = TRUE)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>% filter(p.adjust < PADJ_CUTOFF) %>% mutate(Comparison = name, Direction = direction)
}

go_all <- bind_rows(lapply(CONTRASTS, function(ct) {
  res <- results_by_contrast[[ct$name]]$full
  up   <- res %>% filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange >  LFC_THR) %>% pull(gene_name) %>% na.omit()
  down <- res %>% filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange < -LFC_THR) %>% pull(gene_name) %>% na.omit()
  bind_rows(run_go(up, ct$name, "Up"), run_go(down, ct$name, "Down"))
}))

if (nrow(go_all) > 0) {
  write_csv(go_all, file.path(TABLES_DIR, "R02_GO_enrichment_per_genotype.csv"))
  cat(sprintf("[R02] GO enrichment: %d significant terms across genotypes/directions\n", nrow(go_all)))
}

########
## NE differentiation marker barplot
########

# Two conceptually distinct groups, not one undifferentiated list:
#   - NE identity (expected uniformly HIGH, regardless of genotype): PAX6,
#     NES (Nestin) are neural-progenitor-restricted. SOX2 is included here,
#     not with the pluripotency-exit group, despite being one of the core
#     Yamanaka pluripotency factors -- it is pluripotency-ASSOCIATED but not
#     pluripotency-RESTRICTED, since it is also required for and maintained
#     in neural stem/progenitor identity. Confirmed empirically: SOX2 is the
#     single most genotype-uniform gene of any candidate checked (WT/KDM6A_ko/
#     KMT2D_Het mean-count fold-range 1.29x, not a DEG in either contrast) --
#     consistent with it marking the shared NE lineage rather than tracking
#     pluripotency exit.
#   - Pluripotency exit (expected uniformly LOW/near-off, regardless of
#     genotype): POU5F1 (OCT4) and NANOG are pluripotency-restricted and
#     genuinely near-zero in all three genotypes here, confirming exit
#     occurred everywhere -- their between-genotype fold-change looks large
#     only because DESeq2 fold-change is noisy near the count floor, not
#     because of a real differentiation difference. ZFP42 (REX1) is added as
#     a better-behaved pluripotency-restricted marker away from that floor
#     (fold-range 1.25x, not a DEG).
marker_genes_ne     <- c("PAX6", "SOX2", "NES")
marker_genes_pluri  <- c("POU5F1", "NANOG", "ZFP42")
marker_genes <- c(marker_genes_ne, marker_genes_pluri)
marker_ids <- gene_map %>% filter(gene_name %in% marker_genes)

build_marker_data <- function(value_matrix) {
  as.data.frame(value_matrix[intersect(marker_ids$gene_id, rownames(value_matrix)), , drop = FALSE]) %>%
    rownames_to_column("gene_id") %>%
    left_join(marker_ids, by = "gene_id") %>%
    dplyr::select(gene_name, everything(), -gene_id) %>%
    pivot_longer(-gene_name, names_to = "sample_id", values_to = "value") %>%
    left_join(rna_meta %>% dplyr::select(sample_id, genotype), by = "sample_id") %>%
    mutate(gene_name = factor(gene_name, levels = marker_genes))
}

make_marker_plot <- function(marker_data, y_label, free_scale) {
  marker_summary <- marker_data %>% group_by(gene_name, genotype) %>%
    summarise(mean_value = mean(value), sd_value = sd(value), .groups = "drop")
  ggplot(marker_summary, aes(x = genotype, y = mean_value, fill = genotype)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = pmax(mean_value - sd_value, 0), ymax = mean_value + sd_value), width = 0.25, linewidth = 0.4) +
    geom_point(data = marker_data, aes(x = genotype, y = value), size = 1.5, shape = 21, fill = "white", color = "black",
               position = position_jitter(width = 0.1, seed = 42)) +
    facet_wrap(~ gene_name, nrow = 1, scales = if (free_scale) "free_y" else "fixed") +
    scale_fill_manual(values = c("WT" = "grey60", COLORS)) +
    labs(y = y_label, x = NULL,
         subtitle = sprintf("NE identity (uniformly high): %s  |  Pluripotency exit (uniformly low): %s",
                             paste(marker_genes_ne, collapse = ", "), paste(marker_genes_pluri, collapse = ", "))) +
    theme_pub + theme(legend.position = "none", strip.text = element_text(face = "italic", size = 11),
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                       plot.subtitle = element_text(size = 8, color = "grey30"))
}

# Original: DESeq2-normalized counts, each gene its own y-range -- best for
# spotting a within-gene genotype effect, not for comparing magnitude across
# genes (e.g. "is POU5F1 much lower than SOX2 overall" is invisible here).
marker_data_counts <- build_marker_data(norm_counts)
p_markers <- make_marker_plot(marker_data_counts, "Normalized counts", free_scale = TRUE)
ggsave(file.path(PLOTS_DIR, "R02_NE_markers_barplot.pdf"), p_markers, width = 10, height = 4)

# Requested alternates: a shared y-axis across genes (so absolute magnitude
# is directly comparable), in both DESeq2-normalized-count and TPM space,
# each raw and log10(x+1)-transformed. The log transform is applied to the
# underlying values before summarising -- not via a log ggplot axis on
# linear geom_col bars, since bars-from-zero don't render sensibly on a true
# log axis and mean +/- SD error bars can dip negative under a naive axis
# transform.
tpm_path <- "data/RNA_nfcore_output/star_salmon/salmon.merged.gene_tpm.tsv"
tpm_matrix <- read_tsv(tpm_path, show_col_types = FALSE) %>%
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
tpm_matrix <- tpm_matrix[, rna_meta$sample_id]
marker_data_tpm <- build_marker_data(tpm_matrix)

marker_variants <- list(
  normcounts_linear = list(data = marker_data_counts,
                            label = "Normalized counts"),
  normcounts_log    = list(data = marker_data_counts %>% mutate(value = log10(value + 1)),
                            label = expression(log[10]("normalized counts" + 1))),
  TPM_linear        = list(data = marker_data_tpm,
                            label = "TPM"),
  TPM_log           = list(data = marker_data_tpm %>% mutate(value = log10(value + 1)),
                            label = expression(log[10](TPM + 1)))
)
for (variant in names(marker_variants)) {
  v <- marker_variants[[variant]]
  p <- make_marker_plot(v$data, v$label, free_scale = FALSE)
  ggsave(file.path(PLOTS_DIR, sprintf("R02_NE_markers_barplot_%s.pdf", variant)), p, width = 10, height = 4)
}

########
## NE identity-acquisition heatmap -- a broader curated marker panel, distinct
## from the top-DEG heatmaps elsewhere in this script (those show what
## differs between genotype and WT; this shows whether the expected identity
## was acquired at all, regardless of genotype). Genes grouped into fixed
## biological blocks rather than clustered, so the narrative reads directly
## off the block order -- plus mesoderm/endoderm blocks as negative controls,
## to show the differentiation was specifically neuroectodermal, not just
## "some genes moved". TPM, not normalized counts: several of these genes
## (TBXT, SOX17, GATA4, ...) fall below R02's own count pre-filter and are
## entirely absent from norm_counts/dds, but are legitimately present
## (near-zero, as expected for a negative control) in TPM, which covers the
## full annotated gene set.
##
## Deliberately NOT row-scaled (z-score): row-scaling normalizes each gene to
## its own across-sample SD, so a gene that only wobbles by noise (e.g. a
## near-zero-count gene, or a borderline case like SOX1/OTX2) gets stretched
## to the same color intensity as a gene with a real, large fold-change --
## exactly backwards for a panel meant to show which genes genuinely changed.
## Plotted as log2(TPM+1) directly on one shared sequential color scale
## instead, so magnitude is preserved: a real change (e.g. POU5F1 dropping to
## near-zero) looks dramatic, a flat gene (e.g. SOX2, uniformly retained)
## looks flat, same principle as the log-vs-linear barplot fix above.
########

identity_panel <- list(
  Pluripotency = c("POU5F1", "NANOG", "ZFP42", "LIN28A", "DPPA4", "SALL4", "DNMT3B"),
  `NE/Neural`  = c("PAX6", "SOX1", "SOX2", "NES", "OTX2", "SOX3", "VIM"),
  Mesoderm     = c("TBXT", "MESP1", "MIXL1"),
  Endoderm     = c("SOX17", "FOXA2", "GATA4", "GATA6")
)
identity_genes <- unlist(identity_panel, use.names = FALSE)
identity_ids <- gene_map %>% filter(gene_name %in% identity_genes) %>% distinct(gene_name, .keep_all = TRUE)

identity_mat <- tpm_matrix[intersect(identity_ids$gene_id, rownames(tpm_matrix)), , drop = FALSE]
rownames(identity_mat) <- identity_ids$gene_name[match(rownames(identity_mat), identity_ids$gene_id)]
identity_mat <- identity_mat[identity_genes[identity_genes %in% rownames(identity_mat)], , drop = FALSE]

genotype_order <- c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment)
col_order <- order(match(rna_meta$genotype[match(colnames(identity_mat), rna_meta$sample_id)], genotype_order))
identity_mat <- identity_mat[, col_order]
identity_mat_log <- log2(identity_mat + 1)

block_labels <- setNames(rep(names(identity_panel), lengths(identity_panel)), identity_genes)
final_blocks <- block_labels[rownames(identity_mat_log)]
gaps_row <- cumsum(rle(final_blocks)$lengths); gaps_row <- gaps_row[-length(gaps_row)]

col_genotypes <- as.character(rna_meta$genotype[match(colnames(identity_mat_log), rna_meta$sample_id)])
# Columns were already ordered by genotype above (col_order) -- cluster_cols=TRUE
# would let pheatmap's own dendrogram override that with data-driven clustering,
# which on a 21-gene panel can easily interleave samples across genotype.
# cluster_cols=FALSE forces the genotype grouping to actually hold; gaps_col
# (from the same real column composition, not a hardcoded replicate count)
# reinforces it visually.
gaps_col <- cumsum(rle(col_genotypes)$lengths); gaps_col <- gaps_col[-length(gaps_col)]

row_ann <- data.frame(Category = final_blocks, row.names = rownames(identity_mat_log))
col_ann <- data.frame(Genotype = col_genotypes, row.names = colnames(identity_mat_log))
ann_colors <- list(Genotype = setNames(c("grey50", COLORS), genotype_order),
                    Category = setNames(c("#8E44AD", "#2980B9", "#D68910", "#229954"), names(identity_panel)))

pdf(file.path(PLOTS_DIR, "R02_identity_acquisition_heatmap.pdf"), width = 8, height = 9)
pheatmap(identity_mat_log, annotation_row = row_ann, annotation_col = col_ann, annotation_colors = ann_colors,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         breaks = seq(0, max(identity_mat_log), length.out = 100),
         cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps_row, gaps_col = gaps_col, border_color = "grey",
         show_rownames = TRUE, fontsize_row = 9,
         main = "NE identity acquisition: curated marker panel (log2(TPM+1))")
dev.off()

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "R02_DESeq2_session_info.txt"))
cat("\n[DONE] R02_DESeq2 complete\n")
