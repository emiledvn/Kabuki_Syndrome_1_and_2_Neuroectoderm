#!/usr/bin/env Rscript
# AR07b_save_103gene_GO_table.R -- one-off: reruns the exact GO BP enrichment
# used for AR07_GO_103gene_concordant.pdf (AR07_Publication_Figures.R TASK 3,
# 103-gene ATAC-RNA concordant set from AR05) and writes the full result table
# to disk, since AR07_Publication_Figures.R plots it but never saves it.
# Same universe, same enrichGO call, same redundancy-collapse step -- read-only
# over already-computed tables, no thresholds changed.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR07b_save_103gene_GO_table.R

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(yaml)
})

cfg <- yaml::read_yaml("config/pipeline_config.yaml")
FDR <- cfg$thresholds$fdr

TABLES_DIR <- "results/tables"

norm_counts <- readRDS("results/RDS/R02_normalized_counts.rds")
gene_map    <- readRDS("results/RDS/R02_gene_id_name_map.rds")
universe_symbols <- gene_map$gene_name[gene_map$gene_id %in% rownames(norm_counts)]
universe_entrez  <- bitr(universe_symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
cat(sprintf("GO background universe: %d expressed genes\n", length(universe_symbols)))

run_go_bp <- function(genes, label) {
  entrez <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = FDR, readable = TRUE)
  ego@result %>% filter(p.adjust < FDR) %>% mutate(Group = label)
}

collapse_redundant_go_terms <- function(go_df, jaccard_thresh = 0.7) {
  jaccard <- function(g1, g2) length(intersect(g1, g2)) / length(union(g1, g2))
  go_df %>%
    mutate(genes = str_split(geneID, "/")) %>%
    group_by(Group) %>%
    group_modify(function(d, key) {
      n <- nrow(d)
      edges <- integer(0)
      if (n > 1) for (i in seq_len(n - 1)) for (j in seq((i + 1), n)) {
        if (jaccard(d$genes[[i]], d$genes[[j]]) >= jaccard_thresh) edges <- c(edges, i, j)
      }
      g <- igraph::make_empty_graph(n = n, directed = FALSE)
      if (length(edges) > 0) g <- igraph::add_edges(g, edges)
      d$cluster_id <- igraph::components(g)$membership
      d %>% group_by(cluster_id) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup() %>%
        dplyr::select(-cluster_id)
    }) %>% ungroup() %>% dplyr::select(-genes)
}

ar05_data <- read_csv(file.path(TABLES_DIR, "AR05_concordant_heatmap_data.csv"), show_col_types = FALSE)
gained_genes <- ar05_data %>% filter(Direction == "Up")   %>% pull(SYMBOL)
lost_genes   <- ar05_data %>% filter(Direction == "Down") %>% pull(SYMBOL)
cat(sprintf("gained (n=%d), lost (n=%d)\n", length(gained_genes), length(lost_genes)))

go_gained <- run_go_bp(gained_genes, "Gained (concordant UP)")
go_lost   <- run_go_bp(lost_genes,   "Lost (concordant DOWN)")
go_all_full <- bind_rows(go_gained, go_lost)
cat(sprintf("GO BP terms (uncollapsed): %d\n", nrow(go_all_full)))

go_all_collapsed <- collapse_redundant_go_terms(go_all_full)
cat(sprintf("GO BP terms (redundancy-collapsed, Jaccard>=0.7): %d\n", nrow(go_all_collapsed)))

out_full <- file.path(TABLES_DIR, "AR07_GO_103gene_concordant_full.csv")
out_collapsed <- file.path(TABLES_DIR, "AR07_GO_103gene_concordant_collapsed.csv")
write_csv(go_all_full, out_full)
write_csv(go_all_collapsed, out_collapsed)
cat(sprintf("[OK] wrote %s (%d rows)\n", out_full, nrow(go_all_full)))
cat(sprintf("[OK] wrote %s (%d rows)\n", out_collapsed, nrow(go_all_collapsed)))
