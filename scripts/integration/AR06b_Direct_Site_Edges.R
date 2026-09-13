#!/usr/bin/env Rscript
# AR06b_Direct_Site_Edges.R -- TF -> peak -> gene "direct-site" edge table,
# extracted from the source working repo's AR06_TF_Regulatory_Network.R.
#
# That script built TWO things: (1) a peak-resolved TF-target edge table from
# TOBIAS's own per-site footprint calls, and (2) an igraph-based network
# visualization/hub-TF/module analysis on top of it. Only (2) was superseded
# (by the CollecTRI-based network actually used in the paper); (1) is a
# genuinely separate computation that AR11d_CollecTRI_dense_with_candidates.R
# still depends on (its `AR06_direct_site_edges_<contrast>.csv` input). This
# script keeps only the "Direct-site" edge-construction pathway (TF's own
# binding site must itself be meaningfully differential; candidate genes
# restricted to DEGs; peak->gene assignment via rGREAT) -- the "DAR-anchored"
# pathway (site must also fall in a separately-called DAR peak) that the
# original script built alongside it is dropped here, since nothing
# downstream in this repository reads it.
#
# Requires (all produced by scripts already in this repository):
#   scripts/atac/A06_TOBIAS.sh              -- data/tobias_output/, consensus_peaks.bed
#   scripts/integration/AR01_ATAC_RNA_Integration.R  -- (not read directly; DEG/consensus-peak
#                                              re-annotation is redone here via rGREAT, matching
#                                              AR01's own basal+extension model)
#   scripts/rna/R02_DESeq2.R                 -- results/RDS/R02_normalized_counts.rds,
#                                              results/RDS/R02_gene_id_name_map.rds,
#                                              results/tables/R02_<contrast>_full.csv
#   scripts/atac/A07b_TOBIAS_DOWNSTREAM.R    -- results/tables/A07b_concordant_tfs.csv
#                                              (optional, cross-check only)
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/integration/AR06b_Direct_Site_Edges.R
# Self-checkpointing: skips if every contrast's output CSV already exists.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(data.table)
  library(GenomicRanges)
  library(parallel)
  library(yaml)
  library(rGREAT)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR06b] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TOBIAS_DIR <- "data/tobias_output"
TABLES_DIR <- "results/tables"
RDS_DIR    <- "results/RDS"
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

CONTRASTS     <- cfg$contrasts
PADJ_CUTOFF   <- cfg$thresholds$fdr
LFC_THR       <- cfg$thresholds$lfc_threshold
MIN_MEAN_EXPR <- cfg$thresholds$min_mean_expression
TOBIAS_PVAL   <- cfg$thresholds$tobias_binding_pvalue
TOBIAS_CHANGE <- cfg$thresholds$tobias_binding_change
N_CORES <- if (!is.null(cfg$resources$network_build_parallel)) cfg$resources$network_build_parallel else 4

DONE_MARKERS <- file.path(TABLES_DIR, sprintf("AR06_direct_site_edges_%s.csv", sapply(CONTRASTS, function(ct) ct$name)))
if (all(file.exists(DONE_MARKERS))) {
  cat("[AR06b] Already complete (all AR06_direct_site_edges_*.csv exist). Skipping. Delete them to force a rerun.\n")
  quit(save = "no", status = 0)
}

if (!dir.exists(TOBIAS_DIR)) stop("[AR06b] ERROR: ", TOBIAS_DIR, " not found -- run scripts/atac/A06_TOBIAS.sh first.")
norm_counts_path <- file.path(RDS_DIR, "R02_normalized_counts.rds")
if (!file.exists(norm_counts_path)) stop("[AR06b] ERROR: ", norm_counts_path, " not found -- run scripts/rna/R02_DESeq2.R first.")

cat("================================================================================\n")
cat("AR06b_Direct_Site_Edges -- TF -> peak -> gene edge table (Direct-site pathway only)\n")
cat(sprintf("THRESHOLDS: FDR < %.2f, |log2FC| > %.2f | TOBIAS binding pvalue < %.2f, |change| > %.2f | min expr > %d\n",
            PADJ_CUTOFF, LFC_THR, TOBIAS_PVAL, TOBIAS_CHANGE, MIN_MEAN_EXPR))
cat("================================================================================\n\n")

########
## Step 1: expressed-TF universe + canonical symbol map.
########

norm_counts <- readRDS(norm_counts_path)
gene_map    <- readRDS(file.path(RDS_DIR, "R02_gene_id_name_map.rds"))
mean_expr_by_symbol <- tibble(gene_id = rownames(norm_counts), mean_expr = rowMeans(norm_counts)) %>%
  left_join(gene_map, by = "gene_id") %>%
  filter(!is.na(gene_name)) %>%
  group_by(gene_name) %>% summarise(mean_expr = max(mean_expr), .groups = "drop")
expressed_symbols_upper <- toupper(mean_expr_by_symbol$gene_name[mean_expr_by_symbol$mean_expr > MIN_MEAN_EXPR])

symbol_lookup <- setNames(mean_expr_by_symbol$gene_name, toupper(mean_expr_by_symbol$gene_name))

tf_is_expressed <- function(tf_name) {
  parts <- toupper(strsplit(tf_name, "::", fixed = TRUE)[[1]])
  all(parts %in% expressed_symbols_upper)
}
canonical_tf_symbol <- function(tf_name) {
  parts <- strsplit(tf_name, "::", fixed = TRUE)[[1]]
  resolved <- symbol_lookup[toupper(parts)]
  if (any(is.na(resolved))) return(NA_character_)
  paste(resolved, collapse = "::")
}

bindetect_ref_path <- file.path(TOBIAS_DIR, CONTRASTS[[1]]$name, "BINDetect", "bindetect_results.txt")
if (!file.exists(bindetect_ref_path)) stop("[AR06b] ERROR: ", bindetect_ref_path, " not found -- run scripts/atac/A06_TOBIAS.sh first.")
bindetect_ref <- read_tsv(bindetect_ref_path, show_col_types = FALSE)

for (ct in CONTRASTS[-1]) {
  other_path <- file.path(TOBIAS_DIR, ct$name, "BINDetect", "bindetect_results.txt")
  other <- read_tsv(other_path, show_col_types = FALSE)
  if (!setequal(other$output_prefix, bindetect_ref$output_prefix)) {
    stop("[AR06b] ERROR: ", ct$name, "'s BINDetect TF panel differs from ", CONTRASTS[[1]]$name, "'s -- motif DB mismatch.")
  }
}

tf_universe <- bindetect_ref %>%
  transmute(TF_name = name, motif_id, output_prefix, motif_cluster = cluster,
            is_composite = str_detect(TF_name, "::"),
            expressed = vapply(TF_name, tf_is_expressed, logical(1)),
            TF_symbol_canonical = vapply(TF_name, canonical_tf_symbol, character(1)))

n_expressed <- sum(tf_universe$expressed, na.rm = TRUE)
cat(sprintf("[AR06b] TF-expression filter (mean normalized count > %d): %d/%d motifs kept.\n",
            MIN_MEAN_EXPR, n_expressed, nrow(tf_universe)))
write_csv(tf_universe, file.path(TABLES_DIR, "AR06_tf_expressed_universe.csv"))

expressed_tfs <- tf_universe %>% filter(expressed, !is.na(TF_symbol_canonical))

########
## Step 2: TF-level binding stats (sign-corrected reference-relative ->
## treatment-relative, same flip A07b applies at the TF level).
########

build_tf_binding_summary <- function(ct) {
  bindetect_path <- file.path(TOBIAS_DIR, ct$name, "BINDetect", "bindetect_results.txt")
  bindetect <- read_tsv(bindetect_path, show_col_types = FALSE)
  change_col <- grep("_change$", names(bindetect), value = TRUE)[1]
  pvalue_col <- grep("_pvalue$", names(bindetect), value = TRUE)[1]

  rna_data <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)), show_col_types = FALSE)

  tb <- bindetect %>%
    transmute(TF_name = name, motif_id, output_prefix,
              binding_change = -.data[[change_col]],
              binding_pvalue = .data[[pvalue_col]]) %>%
    left_join(tf_universe %>% dplyr::select(output_prefix, motif_cluster, is_composite, expressed, TF_symbol_canonical), by = "output_prefix")

  tb <- tb %>% mutate(binding_padj_all = p.adjust(binding_pvalue, method = "BH"))
  cluster_rep <- tb %>% filter(!is.na(motif_cluster)) %>% group_by(motif_cluster) %>%
    slice_min(binding_pvalue, n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(motif_cluster, cluster_rep_pvalue = binding_pvalue)
  cluster_rep$binding_padj_cluster <- p.adjust(cluster_rep$cluster_rep_pvalue, method = "BH")
  tb <- tb %>% left_join(cluster_rep %>% dplyr::select(motif_cluster, binding_padj_cluster), by = "motif_cluster")

  tb %>%
    mutate(binding_sig = binding_padj_all < TOBIAS_PVAL & abs(binding_change) > TOBIAS_CHANGE,
           binding_direction = case_when(binding_change > 0 ~ "Increased", binding_change < 0 ~ "Decreased", TRUE ~ "No change")) %>%
    left_join(rna_data %>% dplyr::select(gene_name, log2FoldChange, padj), by = c("TF_symbol_canonical" = "gene_name")) %>%
    mutate(expression_sig = !is.na(padj) & padj < PADJ_CUTOFF,
           expression_direction = case_when(log2FoldChange > 0 ~ "Increased", log2FoldChange < 0 ~ "Decreased", TRUE ~ "No change"),
           tf_class = case_when(
             binding_sig & expression_sig & binding_direction == expression_direction ~ "Concordant TF",
             binding_sig ~ "Differential binding",
             TRUE ~ "NS"),
           contrast = ct$name)
}

a07b_concordant_path <- file.path(TABLES_DIR, "A07b_concordant_tfs.csv")
a07b_concordant_names <- if (file.exists(a07b_concordant_path)) {
  read_csv(a07b_concordant_path, show_col_types = FALSE)$TF_name
} else character(0)

tf_binding_by_contrast <- list()
for (ct in CONTRASTS) {
  tb <- build_tf_binding_summary(ct) %>%
    mutate(a07b_flagged_concordant = TF_name %in% a07b_concordant_names)
  tf_binding_by_contrast[[ct$name]] <- tb
  write_csv(tb, file.path(TABLES_DIR, sprintf("AR06_tf_binding_summary_%s.csv", ct$name)))
  cat(sprintf("[AR06b] %s: %d Concordant TF, %d Differential binding, %d NS (of %d total TFs)\n",
              ct$name, sum(tb$tf_class == "Concordant TF"), sum(tb$tf_class == "Differential binding"),
              sum(tb$tf_class == "NS"), nrow(tb)))
}

########
## Step 3: consensus-peak GREAT annotation + "direct-site" whitelist -- a TF
## site must itself be meaningfully differential (Step 4); the candidate gene
## must be a DEG; peak->gene assignment via rGREAT basal+extension (same
## model AR01_ATAC_RNA_Integration.R uses for the DAR-anchored side),
## computed once on the full consensus peak set (shared across contrasts).
########

CONSENSUS_PEAKS_PATH <- file.path(TOBIAS_DIR, "consensus_peaks.bed")
if (!file.exists(CONSENSUS_PEAKS_PATH)) stop("[AR06b] ERROR: ", CONSENSUS_PEAKS_PATH, " not found -- run scripts/atac/A06_TOBIAS.sh first.")
consensus_peaks <- read_tsv(CONSENSUS_PEAKS_PATH, col_names = c("chr", "start", "end"), show_col_types = FALSE)
consensus_gr <- GRanges(consensus_peaks$chr, IRanges(consensus_peaks$start, consensus_peaks$end))
cat(sprintf("[AR06b] Consensus peak set (A06, shared across both contrasts): %d peaks\n", length(consensus_gr)))

great_cache_path <- file.path(RDS_DIR, "AR06_great_consensus_peak_gene_assoc.rds")
if (file.exists(great_cache_path)) {
  cat("[AR06b] loading cached consensus peak->gene associations from", great_cache_path, "\n")
  consensus_great_assoc <- readRDS(great_cache_path)
} else {
  cat(sprintf("[AR06b] running rGREAT basal+extension on the full %d-peak consensus set (one-time, cached to RDS)...\n", length(consensus_gr)))
  great_res <- great(consensus_gr, "GO:BP", "hg38")
  great_assoc_gr <- getRegionGeneAssociations(great_res)
  consensus_great_assoc <- as.data.frame(great_assoc_gr) %>%
    tidyr::unnest(c(annotated_genes, dist_to_TSS)) %>%
    dplyr::rename(SYMBOL = annotated_genes, distanceToTSS = dist_to_TSS) %>%
    mutate(annotation_class = case_when(
      abs(distanceToTSS) <= 2000 ~ "Promoter", abs(distanceToTSS) <= 20000 ~ "Proximal",
      abs(distanceToTSS) <= 100000 ~ "Distal", TRUE ~ "Intergenic"))
  saveRDS(consensus_great_assoc, great_cache_path)
}
cat(sprintf("[AR06b] %d peak-gene associations across %d consensus peaks, %d unique genes\n",
            nrow(consensus_great_assoc), length(consensus_gr), length(unique(consensus_great_assoc$SYMBOL))))

build_direct_site_whitelist <- function(ct) {
  rna_data <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)), show_col_types = FALSE)
  degs <- rna_data %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_THR) %>%
    dplyr::select(gene_name) %>% distinct()
  links <- consensus_great_assoc %>% filter(SYMBOL %in% degs$gene_name) %>%
    transmute(seqnames = as.character(seqnames), start, end, SYMBOL, distanceToTSS, annotation_class,
              Fold = NA_real_, FDR = NA_real_, dar_direction = NA_character_)
  if (nrow(links) == 0) return(list(gr = GRanges(), links = tibble(), n_degs = nrow(degs)))
  gr <- GRanges(links$seqnames, IRanges(links$start, links$end))
  list(gr = gr, links = links, n_degs = nrow(degs))
}

direct_site_whitelist <- setNames(lapply(CONTRASTS, build_direct_site_whitelist), sapply(CONTRASTS, function(ct) ct$name))
for (ct in CONTRASTS) {
  w <- direct_site_whitelist[[ct$name]]
  cat(sprintf("[AR06b] %s: %d DEGs -> %d candidate GREAT-assigned peak-gene links (direct-site whitelist)\n",
              ct$name, w$n_degs, nrow(w$links)))
}

########
## Step 4: per-TF site extraction from TOBIAS's per-site overview files,
## restricted to expressed TFs, filtered to bound-in-either-condition AND
## meaningfully differential at the site itself (no DAR-peak requirement).
########

extract_tf_sites <- function(output_prefix, ct, whitelist_gr, min_abs_log2fc) {
  f <- file.path(TOBIAS_DIR, ct$name, "BINDetect", output_prefix, paste0(output_prefix, "_overview.txt"))
  if (!file.exists(f)) return(NULL)
  treat <- ct$treatment
  select_cols <- c("TFBS_start", "TFBS_end", "TFBS_score", "peak_chr", "peak_start", "peak_end",
                    "WT_score", paste0(treat, "_score"),
                    "WT_bound", paste0(treat, "_bound"),
                    paste0("WT_", treat, "_log2fc"))
  sites <- tryCatch(data.table::fread(f, sep = "\t", select = select_cols, showProgress = FALSE),
                     error = function(e) NULL)
  if (is.null(sites) || nrow(sites) == 0) return(NULL)
  data.table::setnames(sites,
    old = c(paste0(treat, "_score"), paste0(treat, "_bound"), paste0("WT_", treat, "_log2fc")),
    new = c("treatment_score", "treatment_bound", "raw_log2fc"))

  sites <- sites[(WT_bound == 1 | treatment_bound == 1) & abs(raw_log2fc) > min_abs_log2fc]
  if (nrow(sites) == 0) return(NULL)

  sites_gr <- GRanges(sites$peak_chr, IRanges(sites$peak_start, sites$peak_end))
  ov <- findOverlaps(sites_gr, whitelist_gr)
  if (length(ov) == 0) return(NULL)

  out <- sites[queryHits(ov)]
  out[, `:=`(output_prefix = output_prefix, dar_row = subjectHits(ov))]
  out
}

extract_all_sites_for_contrast <- function(ct, whitelist_gr, min_abs_log2fc, label) {
  tf_list <- expressed_tfs$output_prefix
  cat(sprintf("[AR06b] %s: extracting per-site TOBIAS data for %d expressed TFs vs %s whitelist (%d cores)...\n",
              ct$name, length(tf_list), label, N_CORES))
  results <- parallel::mclapply(tf_list, extract_tf_sites, ct = ct, whitelist_gr = whitelist_gr,
                                  min_abs_log2fc = min_abs_log2fc, mc.cores = N_CORES)
  n_with_hits <- sum(!vapply(results, is.null, logical(1)))
  combined <- data.table::rbindlist(Filter(Negate(is.null), results))
  cat(sprintf("[AR06b] %s: %d/%d expressed TFs contributed >=1 site overlapping a %s; %d raw filtered site-rows total\n",
              ct$name, n_with_hits, length(tf_list), label, nrow(combined)))
  combined
}

sites_by_contrast_direct <- setNames(
  lapply(CONTRASTS, function(ct) extract_all_sites_for_contrast(ct, direct_site_whitelist[[ct$name]]$gr, min_abs_log2fc = TOBIAS_CHANGE, label = "DEG-proximal peak")),
  sapply(CONTRASTS, function(ct) ct$name))

########
## Step 5: sign convention flip (reference-relative -> treatment-relative),
## verified numerically before applying.
########

apply_sign_check_and_flip <- function(sites_list, source_label) {
  for (ct_name in names(sites_list)) {
    sites <- sites_list[[ct_name]]
    if (is.null(sites) || nrow(sites) == 0) next
    diff <- sites$treatment_score - sites$WT_score
    big <- abs(diff) > 0.05
    if (any(big)) {
      ok <- all(sign(diff[big]) == sign(-sites$raw_log2fc[big]))
      if (!ok) stop("[AR06b] ERROR: ", ct_name, " (", source_label, ") -- per-site TOBIAS log2fc sign-flip check FAILED. ",
                    "Raw WT_{treatment}_log2fc no longer appears to be reference-relative -- ",
                    "do not blindly apply the -1 flip; re-derive the convention before proceeding.")
      cat(sprintf("[AR06b] %s (%s): sign-flip check PASSED (%d/%d sites with |score diff| > 0.05 agree)\n",
                  ct_name, source_label, sum(big), length(big)))
    }
    sites[, site_log2fc_corrected := -raw_log2fc]
    sites_list[[ct_name]] <- sites
  }
  sites_list
}

sites_by_contrast_direct <- apply_sign_check_and_flip(sites_by_contrast_direct, "Direct-site")

########
## Step 5b: per-site significance -- robust z-score (median/MAD) against that
## TF's own site_log2fc_corrected distribution, converted to a two-sided
## normal p-value.
########

add_site_significance <- function(sites_list) {
  for (ct_name in names(sites_list)) {
    sites <- sites_list[[ct_name]]
    if (is.null(sites) || nrow(sites) == 0) next
    sites[, site_zscore := {
      m <- stats::median(site_log2fc_corrected)
      s <- stats::mad(site_log2fc_corrected)
      if (is.na(s) || s == 0) rep(NA_real_, .N) else (site_log2fc_corrected - m) / s
    }, by = output_prefix]
    sites[, site_pvalue := 2 * stats::pnorm(-abs(site_zscore))]
    sites[, site_sig := !is.na(site_pvalue) & site_pvalue < TOBIAS_PVAL]
    cat(sprintf("[AR06b] %s: per-site z-score/p-value assigned to %d sites (%d site_sig at p<%.2f)\n",
                ct_name, nrow(sites), sum(sites$site_sig, na.rm = TRUE), TOBIAS_PVAL))
    sites_list[[ct_name]] <- sites
  }
  sites_list
}

sites_by_contrast_direct <- add_site_significance(sites_by_contrast_direct)

########
## Step 6: TF -> peak -> gene edge assembly (one-to-many: one site's peak can
## be GREAT-assigned to multiple genes).
########

build_edges_for_contrast <- function(ct, sites, links, edge_source) {
  if (is.null(sites) || nrow(sites) == 0) return(NULL)
  rna_data <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", ct$name)), show_col_types = FALSE)
  tf_binding <- tf_binding_by_contrast[[ct$name]]

  gene_side <- links[sites$dar_row, c("SYMBOL", "annotation_class", "distanceToTSS", "Fold", "FDR", "dar_direction")]
  names(gene_side) <- c("SYMBOL", "annotation_class", "distanceToTSS", "ATAC_Fold", "ATAC_FDR", "dar_direction")

  dplyr::bind_cols(as.data.frame(sites), gene_side) %>%
    filter(!is.na(SYMBOL)) %>%
    left_join(tf_universe %>% dplyr::select(output_prefix, TF_symbol_canonical, motif_cluster), by = "output_prefix") %>%
    filter(!is.na(TF_symbol_canonical)) %>%
    left_join(rna_data %>% dplyr::select(gene_name, gene_log2FoldChange = log2FoldChange, gene_padj = padj),
              by = c("SYMBOL" = "gene_name")) %>%
    left_join(tf_binding %>% dplyr::select(output_prefix, binding_change, binding_sig, binding_direction, tf_class, expression_sig),
              by = "output_prefix") %>%
    mutate(contrast = ct$name, edge_source = edge_source)
}

edges_raw_direct <- setNames(
  lapply(CONTRASTS, function(ct) build_edges_for_contrast(ct, sites_by_contrast_direct[[ct$name]], direct_site_whitelist[[ct$name]]$links, "Direct-site")),
  sapply(CONTRASTS, function(ct) ct$name))

########
## Step 7: collapse to one edge per (TF, gene) pair -- single strongest site
## by |site_log2fc_corrected|.
########

collapse_edges <- function(joined) {
  if (is.null(joined) || nrow(joined) == 0) return(NULL)
  joined %>%
    group_by(contrast, TF_symbol_canonical, SYMBOL) %>%
    mutate(n_sites_collapsed = n()) %>%
    slice_max(abs(site_log2fc_corrected), n = 1, with_ties = FALSE) %>%
    ungroup()
}

edges_by_contrast_direct <- lapply(edges_raw_direct, collapse_edges)

########
## Step 8: confidence tiering (existing config thresholds only) + per-TF/
## per-gene target ranking, then write.
########

tier_edges <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0) return(NULL)
  edges <- edges %>% mutate(
    gene_rna_sig = !is.na(gene_padj) & gene_padj < PADJ_CUTOFF & abs(gene_log2FoldChange) > LFC_THR,
    site_direction_matches_gene = sign(site_log2fc_corrected) == sign(gene_log2FoldChange),
    regulatory_mode = case_when(
      is.na(gene_log2FoldChange) | site_log2fc_corrected == 0 ~ "Ambiguous",
      site_direction_matches_gene ~ "Activating",
      TRUE ~ "Repressive"),
    confidence_tier = case_when(
      tf_class == "Concordant TF" & gene_rna_sig & site_direction_matches_gene & site_sig ~ "High-Activating",
      tf_class == "Concordant TF" & gene_rna_sig & !site_direction_matches_gene & site_sig ~ "High-Repressive",
      binding_sig & gene_rna_sig ~ "Medium",
      TRUE ~ "Low"),
    high_confidence_proximal = confidence_tier %in% c("High-Activating", "High-Repressive") &
                                annotation_class %in% c("Promoter", "Proximal"))

  edges <- edges %>%
    group_by(contrast, TF_symbol_canonical) %>%
    arrange(desc(abs(site_log2fc_corrected)), desc(abs(site_zscore)), .by_group = TRUE) %>%
    mutate(tf_target_rank = row_number(), n_targets_for_tf = n()) %>%
    ungroup() %>%
    group_by(contrast, SYMBOL) %>%
    arrange(desc(abs(site_log2fc_corrected)), desc(abs(site_zscore)), .by_group = TRUE) %>%
    mutate(gene_regulator_rank = row_number(), n_regulators_for_gene = n()) %>%
    ungroup()

  edges
}

edges_by_contrast_direct <- lapply(edges_by_contrast_direct, tier_edges)

for (ct_name in names(edges_by_contrast_direct)) {
  ed <- edges_by_contrast_direct[[ct_name]]
  if (is.null(ed)) { cat(sprintf("[AR06b] %s (Direct-site): no edges survived construction.\n", ct_name)); next }
  cat(sprintf("[AR06b] %s (Direct-site): %d edges (High-Activating: %d, High-Repressive: %d, Medium: %d, Low: %d) across %d TFs -> %d DEGs\n",
              ct_name, nrow(ed), sum(ed$confidence_tier == "High-Activating"), sum(ed$confidence_tier == "High-Repressive"),
              sum(ed$confidence_tier == "Medium"), sum(ed$confidence_tier == "Low"),
              length(unique(ed$TF_symbol_canonical)), length(unique(ed$SYMBOL))))
  ed_out <- ed %>% dplyr::select(TF = TF_symbol_canonical, motif_cluster, SYMBOL, contrast, annotation_class, distanceToTSS,
                                    site_log2fc_corrected, site_zscore, site_pvalue, site_sig, TFBS_score, n_sites_collapsed,
                                    gene_log2FoldChange, gene_padj, tf_binding_change = binding_change, tf_class,
                                    regulatory_mode, confidence_tier,
                                    tf_target_rank, n_targets_for_tf, gene_regulator_rank, n_regulators_for_gene)
  saveRDS(ed_out, file.path(RDS_DIR, sprintf("AR06_direct_site_edges_full_%s.rds", ct_name)))
  write_csv(ed_out %>% filter(confidence_tier != "Low"),
            file.path(TABLES_DIR, sprintf("AR06_direct_site_edges_%s.csv", ct_name)))
}

cat("\n[AR06b] Done.\n")
