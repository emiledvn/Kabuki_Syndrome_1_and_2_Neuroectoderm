#!/usr/bin/env Rscript
# AR01_ATAC_RNA_Integration.R -- GREAT peak-to-gene annotation of ATAC DARs,
# ATAC-RNA concordance classification, and joint visualization.
# Adapted from data_import_local/local_scripts/07_INTEGRATION_peak_annotation.R
# and 08_INTEGRATION_atac_rna_correlation.R.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR01_ATAC_RNA_Integration.R
#
# Fixes applied vs the originals:
#   - Every GO enrichment call now passes an explicit `universe` background
#     (genes with a GREAT regulatory-domain assignment in this dataset, or the
#     R02 count-filter universe where appropriate) -- 07's ATAC-target GO call
#     and 08's second ("Custom GO Dotplot") call both lacked one originally.
#   - Config-driven contrasts throughout instead of hardcoded "KDM6A"/"KMT2D".
#   - The three originally-concatenated scripts in 07 (main + a
#     "07b_PLOT_separated_GO_enrichment.R" + an appended DAR-counts section) and
#     08's dead duplicate heatmap block are cleaned up rather than ported verbatim.
#
# Deliberately NOT ported: the original's hardcoded 40-gene "functional
# category" curated heatmap and its hardcoded 12-term GO dotplot. Both were
# curated by hand after inspecting results from a DIFFERENT prior dataset --
# porting those exact gene/term names here would be meaningless on new data.
# Once AR01_GO_shared_concordant.csv exists for real, curate a new version by hand.
#
# Self-checkpointing: skips entirely if results/tables/AR01_shared_targets_only.csv exists.
# Requires A04b (Csaw Loess), A04f (Loess DAR BEDs), and R02 (DESeq2) to have
# completed first. Loess, not TMM (A04): per the normalization-methodology
# decision in docs/A04b_normalization_methodology.md, Loess is the arm
# reported throughout the pipeline; TMM is kept only as the comparison
# baseline that justifies preferring Loess, not fed downstream.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2)
  library(GenomicRanges)
  library(rGREAT)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggVennDiagram); library(UpSetR)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR01] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

DAR_DIR     <- "data/diffbind_output/CsawLoess_Norm"
TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/integration"
RDS_DIR     <- "results/RDS"
SESSION_DIR <- "results/Session_info"
DONE_MARKER <- file.path(TABLES_DIR, "AR01_shared_targets_only.csv")

for (d in c(TABLES_DIR, PLOTS_DIR, RDS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(DONE_MARKER)) {
  cat("[AR01] Already complete (", DONE_MARKER, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS   <- cfg$contrasts
PADJ_CUTOFF <- cfg$thresholds$fdr
LFC_THR     <- cfg$thresholds$lfc_threshold

if (!dir.exists(DAR_DIR)) stop("[AR01] ERROR: ", DAR_DIR, " not found -- run A04b_Csaw_Loess_Norm.R and A04f_Loess_DAR_BED.R first.")
r02_files_exist <- all(sapply(CONTRASTS, function(ct) file.exists(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)))))
if (!r02_files_exist) stop("[AR01] ERROR: R02 output missing -- run R02_DESeq2.R first.")

theme_pub <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())

cat("================================================================================\n")
cat("AR01_ATAC_RNA_Integration\n")
cat("ANNOTATION STRATEGY (rGREAT): basal plus extension (GREAT standard),\n")
cat("  basal domain -5kb/+1kb from TSS, extension up to 1000kb, truncated at\n")
cat("  neighboring basal domains.\n")
cat("================================================================================\n\n")

########
## 1. GREAT peak-to-gene annotation of each contrast's DARs
########

annotate_peaks_great <- function(bed_file, name, direction) {
  cat("[AR01] Annotating:", bed_file, "\n")
  if (!file.exists(bed_file) || file.info(bed_file)$size == 0) { cat("  (empty or missing, skipping)\n"); return(NULL) }
  peaks <- read.table(bed_file, header = FALSE, sep = "\t")
  names(peaks)[1:3] <- c("chr", "start", "end")
  peaks_gr <- makeGRangesFromDataFrame(peaks, seqnames.field = "chr", start.field = "start", end.field = "end")

  res <- great(peaks_gr, "GO:BP", "hg38")
  assoc_gr <- getRegionGeneAssociations(res)

  as.data.frame(assoc_gr) %>%
    mutate(peak_id = paste0(name, "_", direction, "_", row_number()), contrast = name, direction = direction,
           width = end - start + 1) %>%
    tidyr::unnest(c(annotated_genes, dist_to_TSS)) %>%
    dplyr::rename(SYMBOL = annotated_genes, distanceToTSS = dist_to_TSS) %>%
    mutate(annotation_class = case_when(
      abs(distanceToTSS) <= 2000 ~ "Promoter", abs(distanceToTSS) <= 20000 ~ "Proximal",
      abs(distanceToTSS) <= 100000 ~ "Distal", TRUE ~ "Intergenic"))
}

annotated_peaks <- list()
for (ct in CONTRASTS) {
  for (direction in c("GAINED", "LOST")) {
    bed <- file.path(DAR_DIR, sprintf("%s_%s.bed", ct$name, direction))
    key <- paste(ct$name, direction, sep = "_")
    annotated_peaks[[key]] <- annotate_peaks_great(bed, ct$name, direction)
  }
}
all_annotations <- bind_rows(annotated_peaks)
cat(sprintf("\n[AR01] Total peak-to-gene links: %d, unique genes: %d\n\n", nrow(all_annotations), length(unique(all_annotations$SYMBOL))))

write_csv(all_annotations %>% dplyr::select(peak_id, seqnames, start, end, width, SYMBOL, annotation_class, distanceToTSS, contrast, direction),
          file.path(TABLES_DIR, "AR01_all_great_peaks.csv"))
saveRDS(all_annotations, file.path(RDS_DIR, "AR01_all_annotations.rds"))

## Distance-to-TSS and peaks-per-gene plots
distance_data <- all_annotations %>% filter(abs(distanceToTSS) <= 100000) %>% mutate(distance_kb = distanceToTSS / 1000)
p_dist <- ggplot(distance_data, aes(x = distance_kb, fill = paste(contrast, direction))) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-20, 20), linetype = "dotted", alpha = 0.5) +
  labs(title = "Distance to Nearest TSS (GREAT Assignment)",
       subtitle = "Dashed: +/-2kb (Promoter) | Dotted: +/-20kb (Proximal/Distal boundary)",
       x = "Distance to TSS (kb)", y = "Number of Peak-Gene Links", fill = NULL) +
  theme_pub + theme(legend.position = "bottom")
ggsave(file.path(PLOTS_DIR, "AR01_great_distance_to_TSS.pdf"), p_dist, width = 10, height = 6)

peaks_per_gene <- all_annotations %>% filter(!is.na(SYMBOL)) %>% count(contrast, direction, SYMBOL, name = "peak_count") %>%
  mutate(peak_count_capped = pmin(peak_count, 10))
p_ppg <- ggplot(peaks_per_gene, aes(x = factor(peak_count_capped), fill = paste(contrast, direction))) +
  geom_bar(position = "dodge", color = "black", linewidth = 0.3) +
  labs(title = "Differential Peaks per Gene (GREAT Model)", x = "Peaks per Gene (10 = >=10)", y = "Number of Genes", fill = NULL) +
  theme_pub + theme(legend.position = "bottom")
ggsave(file.path(PLOTS_DIR, "AR01_great_peaks_per_gene.pdf"), p_ppg, width = 10, height = 6)

## DAR counts + genomic annotation distribution
dar_counts <- all_annotations %>% distinct(peak_id, contrast, direction) %>% count(contrast, direction)
p_counts <- ggplot(dar_counts, aes(x = contrast, y = n, fill = contrast)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  facet_wrap(~ direction) +
  labs(x = NULL, y = "Number of DARs") + theme_pub + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOTS_DIR, "AR01_DAR_counts.pdf"), p_counts, width = 6, height = 5)

anno_dist <- all_annotations %>% distinct(peak_id, contrast, direction, annotation_class) %>%
  mutate(annotation_class = factor(annotation_class, levels = c("Promoter", "Proximal", "Distal", "Intergenic"))) %>%
  count(contrast, direction, annotation_class) %>% group_by(contrast, direction) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
# facet by direction (same layout as AR01_DAR_counts.pdf above) instead of
# cramming both contrast and direction into one interaction() x-label --
# that produced two unrotated 2-line labels per bar that collided with their
# neighbours in a 6in-wide plot.
p_anno <- ggplot(anno_dist, aes(x = contrast, y = pct, fill = annotation_class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  facet_wrap(~ direction) +
  scale_fill_manual(values = c(Promoter = "#E41A1C", Proximal = "#FF7F00", Distal = "#4DAF4A", Intergenic = "#377EB8"), name = "Annotation") +
  labs(x = NULL, y = "Percentage of DARs (%)") + theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOTS_DIR, "AR01_DAR_genomic_annotation.pdf"), p_anno, width = 7, height = 5)

########
## 2. GO enrichment on ATAC target genes -- explicit universe (all GREAT-assigned
## genes in this dataset), fixing the original's whole-genome-default bug.
## Also: uses rGREAT only for peak->gene mapping (as the original did), so this
## is still a plain gene-set enrichment on the assigned genes, not rGREAT's own
## region-based test -- documented here rather than silently assumed correct.
########

atac_genes_for_go <- all_annotations %>% mutate(group = paste(contrast, direction, sep = "_")) %>%
  distinct(group, SYMBOL) %>% filter(!is.na(SYMBOL))
atac_background <- unique(all_annotations$SYMBOL); atac_background <- atac_background[!is.na(atac_background)]

atac_go_compare <- compareCluster(SYMBOL ~ group, data = atac_genes_for_go, fun = "enrichGO",
                                   universe = atac_background, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                   pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF, qvalueCutoff = PADJ_CUTOFF, readable = TRUE)
if (!is.null(atac_go_compare) && nrow(as.data.frame(atac_go_compare)) > 0) {
  write_csv(as.data.frame(atac_go_compare), file.path(TABLES_DIR, "AR01_great_GO_results.csv"))
  p_go <- dotplot(atac_go_compare, showCategory = 8, font.size = 10) +
    labs(title = "GO Enrichment: Genes with Differential Chromatin Accessibility", x = NULL) +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(PLOTS_DIR, "AR01_great_GO_dotplot.pdf"), p_go, width = 10, height = 8)
}

########
## 3. ATAC-RNA integration: report every individual peak-gene link (all_peaks
## table below), and separately collapse to one row per gene by taking its
## single STRONGEST assigned peak (max |Fold|) -- not sum(Fold)/min(FDR) as
## originally written here. Two independent problems with that original
## aggregation, found during review:
##   - sum(log2FoldChange) across peaks models the peaks' effects as
##     MULTIPLYING together (sum of logs = log of a product), not adding --
##     the wrong operator even if a real cumulative-enhancer-dose effect
##     exists, and it makes gene-level magnitude incomparable across genes
##     with different peak counts (a gene assigned 3 modest peaks would
##     dwarf a gene assigned 1 genuinely larger peak). GREAT's basal+extension
##     model (up to 1Mb) also assigns the same peak to multiple genes and
##     assigns more peaks to genes simply sitting in gene-dense regions, so
##     peak count isn't a clean count of independent regulatory elements to
##     begin with.
##   - min(FDR) looked like an uncorrected multiple-comparisons selection, but
##     is not one in practice: every peak reaching this function already
##     individually passed A04f's DAR call (padj<FDR & |LFC|>threshold) --
##     annotate_peaks_great() is only ever run on A04f's GAINED/LOST BEDs, not
##     the full peak set -- so there is no non-significant candidate in the
##     pool to select against; min() here is a no-op over already-passing
##     values, not a source of inflated significance.
## Picking the single strongest peak (Fold and FDR from the SAME peak, unlike
## the original which could draw Fold and FDR from two different peaks)
## mirrors A07b's existing precedent for the same kind of many-to-one
## collapse (slice_max(abs(binding_change), n=1) per TF). Nothing is
## discarded: AR01_peak_gene_links_{contrast}_{direction}.csv below retains
## every peak assigned to every gene, for anyone who wants multi-peak/dose
## analysis at full resolution (e.g. a peak-resolved network) rather than
## the one-row-per-gene summary the rest of this script's plots/heatmaps use.
##
## diffbind_data is now A04b's Loess results table (peak_id="chr:start-end" +
## DESeq2's own log2FoldChange/padj columns), not DiffBind's dba.report()
## output (which had native seqnames/start/end + Fold/FDR columns) -- so the
## GRanges has to be built by parsing peak_id, and the field names read from
## it are DESeq2's, not DiffBind's. Kept the working column names below
## (Fold/FDR) unchanged past the read-in point to minimize the diff -- they
## still mean the same thing (fold-change, adjusted p-value), just sourced
## from A04b instead of A04 now.
########

integrate_atac_rna <- function(annotations, rna_data, diffbind_data, contrast_name, direction) {
  cat("[AR01] Integrating:", contrast_name, direction, "\n")
  peaks <- annotations %>% filter(contrast == !!contrast_name, direction == !!direction, !is.na(SYMBOL))
  if (nrow(peaks) == 0) return(NULL)

  peaks_gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  diff_m   <- regmatches(diffbind_data$peak_id, regexec("^(.+):(\\d+)-(\\d+)$", diffbind_data$peak_id))
  diff_gr  <- GRanges(seqnames = sapply(diff_m, `[`, 2),
                       ranges = IRanges(start = as.numeric(sapply(diff_m, `[`, 3)), end = as.numeric(sapply(diff_m, `[`, 4))),
                       Fold = diffbind_data$log2FoldChange, FDR = diffbind_data$padj)
  overlaps <- findOverlaps(peaks_gr, diff_gr)
  peaks$Fold <- NA; peaks$FDR <- NA
  peaks$Fold[queryHits(overlaps)] <- diff_gr$Fold[subjectHits(overlaps)]
  peaks$FDR[queryHits(overlaps)]  <- diff_gr$FDR[subjectHits(overlaps)]

  peaks_with_evidence <- peaks %>% filter(!is.na(Fold), !is.na(FDR))

  # Full peak-level detail, one row per peak-gene link -- not collapsed.
  write_csv(peaks_with_evidence %>%
              dplyr::select(peak_id, SYMBOL, seqnames, start, end, annotation_class, distanceToTSS, Fold, FDR),
            file.path(TABLES_DIR, sprintf("AR01_peak_gene_links_%s_%s.csv", contrast_name, direction)))

  # Gene-level collapse: single strongest peak per gene (max |Fold|); Fold and
  # FDR both come from that same peak. See the section-3 header comment above
  # for why this replaces the original sum(Fold)/min(FDR).
  strongest_peak <- peaks_with_evidence %>% group_by(SYMBOL) %>%
    slice_max(abs(Fold), n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(SYMBOL, ATAC_Fold = Fold, ATAC_FDR = FDR)

  annotation_mode <- peaks_with_evidence %>% count(SYMBOL, annotation_class) %>%
    group_by(SYMBOL) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(SYMBOL, primary_annotation = annotation_class)

  consolidated <- strongest_peak %>%
    left_join(peaks_with_evidence %>% count(SYMBOL, name = "n_peaks"), by = "SYMBOL") %>%
    left_join(annotation_mode, by = "SYMBOL")

  consolidated %>%
    left_join(rna_data %>% dplyr::select(gene_name, log2FoldChange, padj), by = c("SYMBOL" = "gene_name")) %>%
    # abs(ATAC_Fold) > LFC_THR: the selected representative peak already
    # passed |log2FC| > LFC_THR upstream (A04f's GAINED/LOST BEDs), so this is
    # normally satisfied automatically -- kept explicit so atac_sig/rna_sig
    # apply the identical magnitude+significance rule rather than relying on
    # that upstream invariant silently.
    mutate(atac_sig = !is.na(ATAC_FDR) & ATAC_FDR < PADJ_CUTOFF & abs(ATAC_Fold) > LFC_THR,
           rna_sig = !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_THR,
           is_super_target = n_peaks >= 3,
           concordance = case_when(
             atac_sig & rna_sig & sign(ATAC_Fold) == sign(log2FoldChange) ~ "Concordant",
             atac_sig & rna_sig & sign(ATAC_Fold) != sign(log2FoldChange) ~ "Discordant",
             atac_sig & !rna_sig ~ "ATAC-only",
             !atac_sig & rna_sig ~ "RNA-only",
             TRUE ~ "Not significant"))
}

integration_results <- list()
for (ct in CONTRASTS) {
  rna_data <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)), show_col_types = FALSE)
  diffbind_data <- read_csv(file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", ct$name)), show_col_types = FALSE)
  for (direction in c("GAINED", "LOST")) {
    key <- paste(ct$name, direction, sep = "_")
    integration_results[[key]] <- integrate_atac_rna(all_annotations, rna_data, diffbind_data, ct$name, direction)
  }
}

for (key in names(integration_results)) {
  if (!is.null(integration_results[[key]]))
    write_csv(integration_results[[key]], file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s.csv", key)))
}

## Combined scatter plots per contrast
##
## Scope note: `gained`/`lost` (and therefore `combined`) only contain genes
## with an ATAC-significant DAR already (integrate_atac_rna() above filters
## `all_annotations` to contrast/direction-matched peaks before this point) --
## there is no non-DAR/background gene in this correlation. So this rho
## answers "given a DAR exists near this gene, does its RNA change track the
## ATAC change," not "is ATAC-RNA coupling stronger than chance genome-wide."
## Do not report/interpret it as the latter -- it is a within-DAR concordance
## statistic, not a null-calibrated genome-wide one.
for (ct in CONTRASTS) {
  gained <- integration_results[[paste(ct$name, "GAINED", sep = "_")]]
  lost   <- integration_results[[paste(ct$name, "LOST", sep = "_")]]
  combined <- bind_rows(gained, lost) %>% filter(!is.na(ATAC_Fold), !is.na(log2FoldChange))
  if (nrow(combined) < 4) next
  cor_result <- cor.test(combined$ATAC_Fold, combined$log2FoldChange, method = "spearman", exact = FALSE)
  p <- ggplot(combined, aes(x = ATAC_Fold, y = log2FoldChange, color = concordance)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(alpha = 0.6, size = 2) +
    labs(title = paste(ct$name, "- ATAC vs RNA Correlation, within ATAC DARs (GREAT Model)"),
         subtitle = sprintf("Spearman rho = %.3f, p = %s | genes with an ATAC DAR only, not genome-wide",
                             cor_result$estimate, format.pval(cor_result$p.value, digits = 2)),
         x = "Aggregated ATAC Fold Change (log2)", y = "RNA log2 Fold Change", color = "Category") +
    theme_pub + theme(legend.position = "bottom")
  ggsave(file.path(PLOTS_DIR, sprintf("AR01_great_scatter_%s_combined.pdf", ct$name)), p, width = 8, height = 8)
}

########
## 4. GO enrichment on strict concordant genes -- explicit universe (this one
## already had it right in the original; kept as-is).
########

go_targets <- bind_rows(lapply(names(integration_results), function(key) {
  df <- integration_results[[key]]; if (is.null(df)) return(NULL)
  df %>% filter(concordance == "Concordant") %>% mutate(group = key) %>% dplyr::select(SYMBOL, group)
})) %>% filter(!is.na(SYMBOL))

if (nrow(go_targets) > 0) {
  go_compare <- compareCluster(SYMBOL ~ group, data = go_targets, fun = "enrichGO", universe = atac_background,
                                OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = PADJ_CUTOFF, qvalueCutoff = PADJ_CUTOFF, readable = TRUE)
  if (!is.null(go_compare) && nrow(as.data.frame(go_compare)) > 0) {
    write_csv(as.data.frame(go_compare), file.path(TABLES_DIR, "AR01_great_GO_concordant.csv"))
    p_go2 <- dotplot(go_compare, showCategory = 8, font.size = 10) +
      labs(title = "Biological Pathways Driven by ATAC-RNA Concordance", x = NULL) +
      theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(PLOTS_DIR, "AR01_great_GO_concordant_dotplot.pdf"), p_go2, width = 11, height = 8)
  }
}

########
## 5. Genotype overlap (Venn + UpSet) of concordant target genes
########

gene_lists <- setNames(lapply(names(integration_results), function(key) {
  df <- integration_results[[key]]; if (is.null(df)) return(character(0))
  df %>% filter(concordance == "Concordant") %>% pull(SYMBOL)
}), names(integration_results))

p_venn <- ggVennDiagram(gene_lists, label_alpha = 0, edge_size = 0.5, set_color = "black") +
  scale_fill_gradient(low = "white", high = "#56B4E9") +
  labs(title = "Overlap of Concordant Target Genes", fill = "Gene Count")
ggsave(file.path(PLOTS_DIR, "AR01_overlap_venn.pdf"), p_venn, width = 8, height = 6)

pdf(file.path(PLOTS_DIR, "AR01_overlap_upset.pdf"), width = 10, height = 6, onefile = FALSE)
upset(fromList(gene_lists), nsets = length(gene_lists), nintersects = 15, order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Overlapping Genes", sets.x.label = "Total Genes per Condition")
invisible(dev.off())

all_unique_genes <- unique(unlist(gene_lists))
overlap_matrix <- data.frame(SYMBOL = all_unique_genes)
for (key in names(gene_lists)) overlap_matrix[[paste0("In_", key)]] <- overlap_matrix$SYMBOL %in% gene_lists[[key]]
overlap_matrix$Total_Conditions <- rowSums(overlap_matrix[, paste0("In_", names(gene_lists)), drop = FALSE])

## "Shared" = concordant in at least one GAINED/LOST group per contrast, for both contrasts
in_any <- function(name_prefix) rowSums(overlap_matrix[, grep(paste0("^In_", name_prefix), names(overlap_matrix)), drop = FALSE]) > 0
shared_targets <- overlap_matrix[in_any(CONTRASTS[[1]]$name) & in_any(CONTRASTS[[2]]$name), ]

write_csv(overlap_matrix, file.path(TABLES_DIR, "AR01_overlap_matrix_full.csv"))
write_csv(shared_targets, DONE_MARKER)

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR01_ATAC_RNA_Integration_session_info.txt"))

cat(sprintf("\n[DONE] AR01_ATAC_RNA_Integration complete. Shared concordant targets: %d\n", nrow(shared_targets)))
cat("[NOTE] The original analysis's hand-curated functional-category heatmap and\n")
cat("       hardcoded-GO-term dotplot were not ported (see script header). Curate\n")
cat("       fresh versions from AR01_great_GO_concordant.csv once real results exist.\n")
