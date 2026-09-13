#!/usr/bin/env Rscript
# AR04_FourWay_Venn.R -- 4-way Venn (RNA DEG x2 genotypes, ATAC DAR-gene x2
# genotypes) showing where transcriptional and chromatin changes converge
# across BOTH genotypes simultaneously, plus GO enrichment on each Venn's
# full 4-way intersection.
#
# ATAC sets are DAR-ASSOCIATED GENES (via AR01's GREAT peak-to-gene
# annotation), not peaks directly: RNA DEGs and ATAC DARs live in different
# ID spaces (gene_id vs peak_id), so the only way to intersect them at all is
# through the same peak-to-gene mapping AR01 already uses -- these sets are
# read directly from AR01_great_correlation_{contrast}_{GAINED,LOST}.csv, not
# recomputed.
#
# Three Venns: (1) any direction (union of up+down per assay/genotype), (2)
# UP-UP-UP-UP (RNA up + ATAC GAINED, both genotypes), (3) DOWN-DOWN-DOWN-DOWN
# (RNA down + ATAC LOST, both genotypes).
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR04_FourWay_Venn.R
# Requires: R02_DESeq2.R and AR01_ATAC_RNA_Integration.R already run.
#
# Self-checkpointing: skips entirely if results/tables/AR04_venn4_any_intersection.csv exists.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(ggVennDiagram)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR04] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/integration"
SESSION_DIR <- "results/Session_info"
DONE <- file.path(TABLES_DIR, "AR04_venn4_any_intersection.csv")
for (d in c(TABLES_DIR, PLOTS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(DONE)) {
  cat("[AR04] Already complete (", DONE, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS <- cfg$contrasts
NAME_A <- CONTRASTS[[1]]$name; NAME_B <- CONTRASTS[[2]]$name
LABEL_A <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment
FDR <- cfg$thresholds$fdr; LFC_THR <- cfg$thresholds$lfc_threshold

cat("================================================================================\n")
cat("AR04_FourWay_Venn -- RNA DEG x ATAC DAR-gene overlap, both genotypes\n")
cat("================================================================================\n\n")

read_rna <- function(name) read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", name)), show_col_types = FALSE)
rna_a <- read_rna(NAME_A); rna_b <- read_rna(NAME_B)

rna_set <- function(df, dir) {
  if (dir == "up") df %>% filter(!is.na(padj), padj < FDR, log2FoldChange > LFC_THR) %>% pull(gene_name)
  else df %>% filter(!is.na(padj), padj < FDR, log2FoldChange < -LFC_THR) %>% pull(gene_name)
}
RNA_A_up <- rna_set(rna_a, "up"); RNA_A_down <- rna_set(rna_a, "down")
RNA_B_up <- rna_set(rna_b, "up"); RNA_B_down <- rna_set(rna_b, "down")

atac_set <- function(name, dir) {
  f <- file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_%s.csv", name, ifelse(dir == "up", "GAINED", "LOST")))
  if (!file.exists(f)) { cat("[AR04] WARNING: ", f, " not found -- empty set.\n"); return(character(0)) }
  read_csv(f, show_col_types = FALSE) %>% pull(SYMBOL) %>% unique()
}
ATAC_A_up <- atac_set(NAME_A, "up"); ATAC_A_down <- atac_set(NAME_A, "down")
ATAC_B_up <- atac_set(NAME_B, "up"); ATAC_B_down <- atac_set(NAME_B, "down")

# Background universe: genes DESeq2 actually tested (R02 count pre-filter)
# AND with a GREAT regulatory-domain assignment in this dataset -- explicit,
# matching R02/AR01's own convention, not clusterProfiler's whole-genome
# default.
norm_counts <- readRDS("results/RDS/R02_normalized_counts.rds")
gene_map    <- readRDS("results/RDS/R02_gene_id_name_map.rds")
rna_universe <- gene_map$gene_name[gene_map$gene_id %in% rownames(norm_counts)]
all_annotations <- readRDS("results/RDS/AR01_all_annotations.rds")
atac_universe <- unique(all_annotations$SYMBOL)
universe_symbols <- intersect(rna_universe, atac_universe)
universe_entrez <- tryCatch(bitr(universe_symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
cat(sprintf("[AR04] GO background universe: %d genes (RNA-tested AND GREAT-assigned)\n\n", length(universe_symbols)))

# ColorBrewer Set1 (red/blue/green/purple) -- the standard, widely-recognized
# categorical Venn palette. Regions are colored by alpha-blending the colors
# of every set they belong to (classic translucent-circle look), so the
# 4-way intersection is automatically the most saturated region on the plot
# with no separate outline needed.
SET_COLORS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
BLEND_ALPHA <- 0.42

blend_colors <- function(hexes, alpha = BLEND_ALPHA, bg = "#FFFFFF") {
  bg_rgb <- grDevices::col2rgb(bg) / 255
  for (h in hexes) {
    c_rgb <- grDevices::col2rgb(h) / 255
    bg_rgb <- alpha * c_rgb + (1 - alpha) * bg_rgb
  }
  grDevices::rgb(bg_rgb[1], bg_rgb[2], bg_rgb[3])
}

run_venn <- function(sets, title, filename) {
  vdata <- process_data(Venn(sets))
  region_edge  <- vdata$regionEdge
  region_label <- vdata$regionLabel
  set_label    <- vdata$setLabel

  ids <- unique(region_edge$id)
  fill_lookup <- vapply(ids, function(id) blend_colors(SET_COLORS[as.integer(strsplit(id, "/")[[1]])]),
                         character(1))

  full_id <- paste(seq_along(sets), collapse = "/")
  is_full <- region_label$id == full_id

  p <- ggplot() +
    geom_polygon(data = region_edge, aes(X, Y, group = id, fill = id), color = NA) +
    scale_fill_manual(values = fill_lookup, guide = "none") +
    geom_text(data = region_label[!is_full, ], aes(X, Y, label = count), size = 5, color = "grey15") +
    geom_text(data = region_label[is_full, ], aes(X, Y, label = count),
              size = 6.5, fontface = "bold", color = "white") +
    geom_text(data = set_label, aes(X, Y, label = name),
              size = 5, fontface = "bold", color = SET_COLORS[as.integer(set_label$id)]) +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(10, 70, 10, 40),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 12)))
  ggsave(file.path(PLOTS_DIR, paste0(filename, ".pdf")), p, width = 8.5, height = 7)
  Reduce(intersect, sets)
}

run_go <- function(genes, label) {
  if (length(genes) < 5) { cat(sprintf("[AR04] %s: too few genes (%d) for GO, skipping.\n", label, length(genes))); return(NULL) }
  entrez <- tryCatch(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
  if (length(entrez) < 5) return(NULL)
  ego <- enrichGO(gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = FDR, readable = TRUE)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>% filter(p.adjust < FDR) %>% mutate(Set = label)
}

########
## 1. Any-direction 4-way Venn
########
sets_any <- list(union(RNA_A_up, RNA_A_down), union(RNA_B_up, RNA_B_down),
                  union(ATAC_A_up, ATAC_A_down), union(ATAC_B_up, ATAC_B_down))
names(sets_any) <- c(sprintf("%s RNA", LABEL_A), sprintf("%s RNA", LABEL_B),
                      sprintf("%s ATAC", LABEL_A), sprintf("%s ATAC", LABEL_B))
inter_any <- run_venn(sets_any, sprintf("RNA DEG x ATAC DAR-gene overlap, %s & %s (any direction)", LABEL_A, LABEL_B), "AR04_venn4_any")
write_csv(tibble(SYMBOL = inter_any), DONE)
cat(sprintf("[AR04] Any-direction 4-way intersection: %d genes\n", length(inter_any)))

########
## 2. UP-UP-UP-UP
########
sets_up <- list(RNA_A_up, RNA_B_up, ATAC_A_up, ATAC_B_up)
names(sets_up) <- c(sprintf("%s RNA up", LABEL_A), sprintf("%s RNA up", LABEL_B),
                     sprintf("%s ATAC gained", LABEL_A), sprintf("%s ATAC gained", LABEL_B))
inter_up <- run_venn(sets_up, sprintf("Concordant UP: RNA up + ATAC gained, %s & %s", LABEL_A, LABEL_B), "AR04_venn4_up")
write_csv(tibble(SYMBOL = inter_up), file.path(TABLES_DIR, "AR04_venn4_up_intersection.csv"))
cat(sprintf("[AR04] UP-UP-UP-UP intersection: %d genes: %s\n", length(inter_up), paste(inter_up, collapse = ", ")))

########
## 3. DOWN-DOWN-DOWN-DOWN
########
sets_down <- list(RNA_A_down, RNA_B_down, ATAC_A_down, ATAC_B_down)
names(sets_down) <- c(sprintf("%s RNA down", LABEL_A), sprintf("%s RNA down", LABEL_B),
                       sprintf("%s ATAC lost", LABEL_A), sprintf("%s ATAC lost", LABEL_B))
inter_down <- run_venn(sets_down, sprintf("Concordant DOWN: RNA down + ATAC lost, %s & %s", LABEL_A, LABEL_B), "AR04_venn4_down")
write_csv(tibble(SYMBOL = inter_down), file.path(TABLES_DIR, "AR04_venn4_down_intersection.csv"))
cat(sprintf("[AR04] DOWN-DOWN-DOWN-DOWN intersection: %d genes: %s\n", length(inter_down), paste(inter_down, collapse = ", ")))

########
## 4. GO enrichment on each Venn's full 4-way intersection
########
go_results <- bind_rows(
  run_go(inter_any, "Any-direction 4-way intersection"),
  run_go(inter_up, "UP-UP-UP-UP intersection"),
  run_go(inter_down, "DOWN-DOWN-DOWN-DOWN intersection")
)
if (!is.null(go_results) && nrow(go_results) > 0) {
  write_csv(go_results, file.path(TABLES_DIR, "AR04_GO_venn4_intersections.csv"))
  cat(sprintf("\n[AR04] GO enrichment: %d significant terms across the 3 intersections\n", nrow(go_results)))
} else cat("\n[AR04] No significant GO terms for any intersection (likely too few genes).\n")

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR04_FourWay_Venn_session_info.txt"))
cat("\n[DONE] AR04_FourWay_Venn complete\n")
