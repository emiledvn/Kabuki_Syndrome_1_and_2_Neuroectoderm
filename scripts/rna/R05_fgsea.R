#!/usr/bin/env Rscript
# R05_fgsea.R -- rank-based GSEA against three gene-set databases (GO
# Biological Process via clusterProfiler::gseGO(); MSigDB Hallmark and
# Reactome via msigdbr + fgsea directly), per contrast, ranked by the DESeq2
# Wald statistic (results()$stat) -- NOT by LFC of any kind. See
# R02_DESeq2.R's no-shrinkage decision comment (extract_results()): apeglm
# shrinkage was tested and rejected because it over-suppresses low-baseMean
# neuroectoderm TFs; ranking by raw or shrunk LFC would carry the same
# TF-suppression risk into GSEA's rank ordering, since GSEA's result depends
# on relative rank across the whole gene list, not just a hard threshold. The
# Wald stat incorporates effect size and precision without that low-count
# blind spot, and is what the DESeq2/fgsea documentation itself recommends.
#
# Hallmark/Reactome added alongside GO BP (2026-07-20) specifically because
# GO BP is highly redundant at the top of any ranked list -- many near-
# duplicate terms describing the same underlying gene module from different
# GO sub-branches (e.g. "electron transport chain" / "aerobic electron
# transport chain" / "respiratory electron transport chain" all describe
# essentially the same OXPHOS genes). Hallmark (50 sets) and Reactome (1839
# curated, far less redundant than GO BP) collapse that into one line and
# surface a cleaner top-N. All three databases share the same ranked vector
# and the same curated KS-relevant term filter below.
#
# Complements R02/R04's threshold-based enrichGO() ORA: GSEA doesn't discard
# effect-size/rank information the way ORA's hard LFC cutoff does, and isn't
# sensitive to the specific |LFC| threshold chosen.
#
# Duplicate gene symbols (multiple Ensembl IDs mapping to one SYMBOL) are
# resolved by keeping the entry with the largest |stat| per symbol -- fgsea
# requires a uniquely-named ranked vector, and this keeps the most informative
# (least ambiguous) value per symbol rather than an arbitrary first-match.
#
# Plotting is NOT done here -- AR02_GO_Plots.R renders all of R02/R04/R05/AR01's
# GO-shaped tables in one consistent style; an earlier version of this script
# had its own raw enrichplot::dotplot() call, removed to avoid two different-
# looking GSEA plots for the same result floating around.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R05_fgsea.R
# Requires: R02_DESeq2.R already run.
#
# Self-checkpointing: skips entirely if results/tables/R05_fgsea_summary.csv exists.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(stringr)
  library(DESeq2)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(msigdbr); library(fgsea)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R05] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

RDS_DIR     <- "results/RDS"
TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/rna"
SESSION_DIR <- "results/Session_info"
SUMMARY     <- file.path(TABLES_DIR, "R05_fgsea_summary.csv")
for (d in c(TABLES_DIR, PLOTS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(SUMMARY)) {
  cat("[R05] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

DDS_RDS <- file.path(RDS_DIR, "R02_dds_object.rds")
GENEMAP_RDS <- file.path(RDS_DIR, "R02_gene_id_name_map.rds")
if (!file.exists(DDS_RDS)) stop("[R05] ERROR: ", DDS_RDS, " not found -- run R02_DESeq2.R first.")

PADJ_CUTOFF <- cfg$thresholds$fdr
CONTRASTS   <- cfg$contrasts

cat("================================================================================\n")
cat("R05_fgsea -- GO BP GSEA (clusterProfiler::gseGO), ranked by Wald statistic\n")
cat("================================================================================\n\n")

dds <- readRDS(DDS_RDS)
gene_map <- readRDS(GENEMAP_RDS)

build_ranked_vector <- function(dds, treatment, reference) {
  res <- results(dds, contrast = c("genotype", treatment, reference))
  df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(stat)) %>%
    left_join(gene_map, by = "gene_id") %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    group_by(gene_name) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))
  setNames(df$stat, df$gene_name)
}

# Curated view: same GSEA results, filtered to terms/pathways plausibly
# relevant to this project's biology (neuroectoderm/forebrain patterning,
# chromatin regulation, the Wnt/BMP/FGF/Notch developmental signalling axes,
# EMT) -- reported ALONGSIDE the full table, not instead of it. Standard,
# legitimate practice for surfacing biologically-prioritized hits out of a
# long ranked list; matches GO's lowercase "forebrain regionalization" style
# descriptions and Hallmark/Reactome's ALL_CAPS_WITH_UNDERSCORES style
# equally since matching is case-insensitive.
KS_RELEVANT_PATTERN <- paste(
  "neuro", "neural", "forebrain", "ectoderm", "chromatin", "histone",
  "methylat", "acetyl", "polycomb", "wnt", "\\bbmp\\b", "fgf", "notch",
  "epithelial.*mesenchymal", "mesenchymal.*epithelial", "\\bemt\\b",
  "patterning", "neural.crest", "pluripotenc", "stem.cell", "differentiation",
  sep = "|"
)
curate_ks_terms <- function(df) df %>% filter(str_detect(Description, regex(KS_RELEVANT_PATTERN, ignore_case = TRUE)))

summary_rows <- list()

run_and_save <- function(res_df, contrast_name, label, tested_universe_note) {
  res_df <- res_df %>% arrange(p.adjust)
  write_csv(res_df, file.path(TABLES_DIR, sprintf("R05_%s_%s.csv", label, contrast_name)))
  curated <- curate_ks_terms(res_df)
  write_csv(curated, file.path(TABLES_DIR, sprintf("R05_%s_%s_curated.csv", label, contrast_name)))

  n_sig <- sum(res_df$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
  n_sig_curated <- sum(curated$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
  cat(sprintf("[R05] %s %s: %d/%d significant (FDR < %.2f)%s | curated subset: %d/%d significant\n",
              contrast_name, label, n_sig, nrow(res_df), PADJ_CUTOFF, tested_universe_note, n_sig_curated, nrow(curated)))
  data.frame(Comparison = contrast_name, GeneSetDB = label, N_tested = nrow(res_df), N_significant = n_sig,
             N_curated = nrow(curated), N_curated_significant = n_sig_curated)
}

for (ct in CONTRASTS) {
  cat(sprintf("\n[R05] === %s ===\n", ct$name))
  cat("[R05] Building Wald-stat-ranked gene list...\n")
  ranked <- build_ranked_vector(dds, ct$treatment, ct$reference)
  cat(sprintf("[R05] %d uniquely-named genes.\n", length(ranked)))

  ## GO Biological Process, via clusterProfiler::gseGO()
  gsea_go <- tryCatch(
    gseGO(geneList = ranked, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 1, seed = TRUE, verbose = FALSE),
    error = function(e) { warning("[R05] gseGO failed for ", ct$name, ": ", conditionMessage(e)); NULL }
  )
  if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
    summary_rows[[length(summary_rows) + 1]] <- run_and_save(
      gsea_go@result %>% mutate(Comparison = ct$name), ct$name, "fgsea", "")
  } else cat(sprintf("[R05] %s GO BP -- no GSEA result.\n", ct$name))

  ## Hallmark and Reactome, via msigdbr + fgsea directly (gseGO() is GO-specific)
  for (msig in list(list(collection = "H", subcollection = NULL, label = "hallmark"),
                     list(collection = "C2", subcollection = "CP:REACTOME", label = "reactome"))) {
    gene_sets <- msigdbr(species = "Homo sapiens", collection = msig$collection, subcollection = msig$subcollection)
    pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
    set.seed(1)  # fgsea's permutation null is stochastic -- fixed seed for reproducibility, same convention as gseGO(seed=TRUE) above
    raw <- suppressWarnings(fgsea(pathways = pathways, stats = ranked, eps = 0))
    if (nrow(raw) == 0) { cat(sprintf("[R05] %s %s -- no GSEA result.\n", ct$name, msig$label)); next }
    res_df <- as.data.frame(raw) %>%
      mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
      transmute(ID = pathway, Description = pathway, setSize = size, NES = NES,
                pvalue = pval, p.adjust = padj, qvalue = padj, leading_edge = leadingEdge,
                Comparison = ct$name)
    summary_rows[[length(summary_rows) + 1]] <- run_and_save(
      res_df, ct$name, msig$label, sprintf(" of %d %s gene sets", nrow(res_df), msig$label))
  }
}

if (length(summary_rows) > 0) {
  write_csv(bind_rows(summary_rows), SUMMARY)
  cat("\n[R05] Summary:\n"); print(bind_rows(summary_rows))
}

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "R05_fgsea_session_info.txt"))
cat("\n[DONE] R05_fgsea complete\n")
