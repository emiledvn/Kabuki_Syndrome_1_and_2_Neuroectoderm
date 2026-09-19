#!/usr/bin/env Rscript
# AR11c_CollecTRI_candidate_uncurated_links.R -- high-confidence AR06 edges
# (chromatin-accessibility-inferred regulatory links, High-Activating/
# High-Repressive confidence tier) where the TF IS well-curated in CollecTRI
# generally (has >=1 curated target somewhere) but this SPECIFIC (TF, gene)
# pair is not curated. TF-level coverage (not just this one pair) is the
# relevance filter: it argues the gap is more likely a genuine CollecTRI
# curation gap for this particular gene than the TF simply being obscure/
# understudied.
#
# These are candidates, not validated relationships -- used downstream
# (AR11d) as a visually distinct, clearly-flagged "not yet curated" layer,
# never treated as equivalent to an orthogonally supported (AR11b) edge.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/integration/AR11c_CollecTRI_candidate_uncurated_links.R
# Requires: AR06b_Direct_Site_Edges.R (AR06_graph_direct_site_<contrast>.rds), AR11_CollecTRI_fetch_regulons.R
# Self-checkpointing: skips entirely if results/tables/AR11_candidate_uncurated_links.csv exists.

suppressPackageStartupMessages({ library(igraph); library(dplyr); library(readr); library(yaml) })

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR11c] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)
CONTRASTS <- sapply(cfg$contrasts, function(ct) ct$name)

RDS_DIR    <- "results/RDS"
TABLES_DIR <- "results/tables"
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
DONE_MARKER <- file.path(TABLES_DIR, "AR11_candidate_uncurated_links.csv")

if (file.exists(DONE_MARKER)) {
  cat("[AR11c] Already complete (", DONE_MARKER, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

collectri_path <- file.path(RDS_DIR, "AR11_collectri_signed.rds")
if (!file.exists(collectri_path)) stop("[AR11c] ERROR: ", collectri_path, " not found -- run AR11_CollecTRI_fetch_regulons.R first.")
graph_paths_exist <- all(sapply(CONTRASTS, function(ct) file.exists(file.path(RDS_DIR, sprintf("AR06_graph_direct_site_%s.rds", ct)))))
if (!graph_paths_exist) stop("[AR11c] ERROR: AR06_graph_direct_site_*.rds missing -- run AR06b_Direct_Site_Edges.R first.")

collectri <- readRDS(collectri_path)
collectri_tfs <- unique(collectri$source)
collectri_pairs <- paste(collectri$source, collectri$target)

out <- list()
for (ct in CONTRASTS) {
  g <- readRDS(file.path(RDS_DIR, sprintf("AR06_graph_direct_site_%s.rds", ct)))
  el <- igraph::as_data_frame(g, what = "edges")
  cand <- el %>%
    filter(from %in% collectri_tfs,                        # TF is well-curated in general
           !(paste(from, to) %in% collectri_pairs),          # this specific pair is absent
           confidence_tier %in% c("High-Activating", "High-Repressive")) %>%
    arrange(desc(abs(site_log2fc_corrected)))
  cat(sprintf("[AR11c] %s: %d high-confidence chromatin-accessibility-inferred edges where the TF is CollecTRI-covered but this pair is not curated\n", ct, nrow(cand)))
  out[[ct]] <- cand
}
combined <- bind_rows(out) %>%
  select(from, to, contrast, regulatory_mode, confidence_tier, site_log2fc_corrected, tf_binding_change, gene_log2FoldChange, gene_padj)
write_csv(combined, DONE_MARKER)

cat("\n[DONE] AR11c_CollecTRI_candidate_uncurated_links complete --", nrow(combined), "candidate edges saved to", DONE_MARKER, "\n")
