#!/usr/bin/env Rscript
# AR02_GO_Plots.R -- publication-style GO Biological Process dotplots for the
# GO enrichment tables already produced by R02, R04, and AR01: fixed point
# size, colour = -log10(FDR), gene count printed inside each point,
# enrichment ratio (GeneRatio / BgRatio) on the x-axis, terms ordered by that
# ratio. Style adapted from a prior (unrelated) project's GO figures, per
# actual file paths/column names in this repo -- not copied verbatim.
#
# Read-only: does not rerun any upstream enrichGO()/compareCluster() call,
# does not touch ATAC BAMs/peaks/TOBIAS. Self-checkpointing per input: each
# plot is skipped individually (with a message) if its source table doesn't
# exist yet -- e.g. AR01's tables, since AR01 hasn't run in this repo yet.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR02_GO_Plots.R

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(readr)
  library(ggplot2); library(scales)
})

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/integration"
SESSION_DIR <- "results/Session_info"
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SESSION_DIR, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        axis.text.y = element_text(size = 9))

# enrichGO()/compareCluster() report GeneRatio/BgRatio as "k/n" strings --
# parsed directly rather than relying on a FoldEnrichment column, since not
# every GO table here carries one (R02/R04 get it for free from
# clusterProfiler's own as.data.frame(); AR01's compareCluster() output does
# not), so parsing keeps all three tables on the same code path.
parse_ratio <- function(x) {
  vapply(x, function(r) {
    v <- as.numeric(strsplit(r, "/")[[1]])
    v[1] / v[2]
  }, numeric(1))
}

GO_COL_TYPES <- c(ID = "character", Description = "character",
                   GeneRatio = "character", BgRatio = "character",
                   pvalue = "numeric", p.adjust = "numeric",
                   qvalue = "numeric", geneID = "character", Count = "integer")

read_go_csv <- function(path) {
  # Explicit colClasses: a GO table can legitimately be header-only with zero
  # rows (nothing significant under the current universe/threshold).
  # read.csv() can't infer types from an empty column and defaults it to
  # logical, which then breaks bind_rows() against a non-empty table of the
  # same shape -- so pin the enrichGO()/compareCluster() output schema
  # explicitly. Named colClasses matches by column name (per ?read.table),
  # so this is robust to the extra columns (Comparison/Direction/Cluster/...)
  # that differ across R02/R04/AR01's tables.
  all_cols <- names(read.csv(path, nrows = 0))
  extra <- setdiff(all_cols, names(GO_COL_TYPES))
  col_types <- c(GO_COL_TYPES, setNames(rep("character", length(extra)), extra))
  read.csv(path, colClasses = col_types)
}

go_dotplot <- function(go_df, title, subtitle, n = 15, facet_var = NULL) {
  df <- go_df %>% filter(!is.na(p.adjust)) %>%
    mutate(enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

  if (!is.null(facet_var)) {
    df <- df %>% group_by(.data[[facet_var]]) %>% arrange(p.adjust) %>% slice_head(n = n) %>% ungroup()
  } else {
    df <- df %>% arrange(p.adjust) %>% slice_head(n = n)
  }
  df <- df %>% mutate(label = str_wrap(Description, 40))

  p <- ggplot(df, aes(x = enrichment_ratio, y = reorder(label, enrichment_ratio), colour = -log10(p.adjust))) +
    geom_point(size = 8, alpha = 0.9) +
    geom_text(aes(label = Count), colour = "white", size = 2.8, fontface = "bold") +
    scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
                            name = expression(-log[10](FDR)),
                            labels = scales::label_number(accuracy = 0.1)) +
    labs(title = title, subtitle = subtitle, x = "Enrichment ratio", y = NULL) +
    theme_pub

  if (!is.null(facet_var)) p <- p + facet_wrap(vars(.data[[facet_var]]), scales = "free")
  p
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(PLOTS_DIR, paste0(name, ".pdf")), p, width = width, height = height)
  cat(sprintf("[AR02] Saved %s.pdf\n", name))
}

# GSEA (R05) is a different result shape from ORA (enrichGO/compareCluster) --
# no GeneRatio/BgRatio, NES instead (signed, so "top" means |NES| not just
# large positive), setSize instead of Count. Same visual language as
# go_dotplot() (fixed dot size, same colour gradient/legend, gene-set-size as
# white text inside the point, wrapped labels) for consistency across every
# plot this script produces, but its own function since the x-axis and
# top-N/ordering logic genuinely differ.
gsea_dotplot <- function(gsea_df, title, subtitle, n = 15, fdr = 0.05) {
  # R05 deliberately runs gseGO/fgsea with pvalueCutoff=1 (keeps every tested
  # term, not just significant ones, so the curated view can still show a
  # relevant-but-subthreshold term) -- so unlike the ORA tables above, whose
  # source enrichGO() call already restricts rows to significant ones, this
  # table can contain non-significant rows. Filter on fdr explicitly before
  # ranking by |NES|, or a big-but-non-significant effect could end up in the
  # "top 15" despite the subtitle claiming FDR-significance.
  df <- gsea_df %>% filter(!is.na(p.adjust), p.adjust < fdr) %>%
    mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
    group_by(direction) %>% arrange(desc(abs(NES))) %>% slice_head(n = n) %>% ungroup() %>%
    mutate(label = str_wrap(Description, 40))
  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = NES, y = reorder(label, abs(NES)), colour = -log10(p.adjust))) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(size = 8, alpha = 0.9) +
    geom_text(aes(label = setSize), colour = "white", size = 2.8, fontface = "bold") +
    scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
                            name = expression(-log[10](FDR)),
                            labels = scales::label_number(accuracy = 0.1)) +
    labs(title = title, subtitle = subtitle, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_pub +
    facet_wrap(vars(direction), scales = "free")  # free x too -- Down/Up NES ranges don't overlap, a shared x-axis wastes half of each panel on empty space
}

made_any <- FALSE

########
## R02 -- per-genotype, per-direction RNA GO (one figure per contrast, Up/Down faceted)
########
r02_path <- file.path(TABLES_DIR, "R02_GO_enrichment_per_genotype.csv")
if (file.exists(r02_path)) {
  go_rna <- read_go_csv(r02_path)
  for (ct in unique(go_rna$Comparison)) {
    sub <- go_rna %>% filter(Comparison == ct)
    if (nrow(sub) == 0) next
    p <- go_dotplot(sub, sprintf("GO Biological Process enrichment: %s vs WT", ct),
                     "Top 15 terms per direction | FDR < 0.05, |LFC| > 1", facet_var = "Direction")
    save_plot(p, sprintf("AR02_GO_RNA_%s", ct), width = 13, height = 7)
    made_any <- TRUE
  }
} else cat("[AR02] Skipping R02 GO plot -- ", r02_path, " not found.\n")

########
## R04 -- shared/overlap RNA GO (Up vs Down faceted)
########
r04_path <- file.path(TABLES_DIR, "R04_GO_overlap.csv")
if (file.exists(r04_path)) {
  go_overlap <- read_go_csv(r04_path)
  if (nrow(go_overlap) > 0) {
    p <- go_dotplot(go_overlap, "GO Biological Process enrichment: shared RNA DEGs",
                     "Top 15 terms per direction | FDR < 0.05, |LFC| > 1", facet_var = "Direction")
    save_plot(p, "AR02_GO_RNA_overlap", width = 13, height = 7)
    made_any <- TRUE
  }
} else cat("[AR02] Skipping R04 GO plot -- ", r04_path, " not found.\n")

########
## AR01 -- ATAC target-gene GO (compareCluster, Cluster = contrast_direction)
########
ar01_path <- file.path(TABLES_DIR, "AR01_great_GO_results.csv")
if (file.exists(ar01_path)) {
  go_atac <- read_go_csv(ar01_path)
  if (nrow(go_atac) > 0) {
    p <- go_dotplot(go_atac, "GO Biological Process enrichment: ATAC target genes",
                     "Top 15 terms per group | FDR < 0.05, |LFC| > 1", facet_var = "Cluster")
    save_plot(p, "AR02_GO_ATAC_targets", width = 15, height = 9)
    made_any <- TRUE
  }
} else cat("[AR02] Skipping AR01 target-gene GO plot -- ", ar01_path, " not found (run AR01 first).\n")

########
## AR01 -- ATAC-RNA concordant-gene GO
########
ar01_conc_path <- file.path(TABLES_DIR, "AR01_great_GO_concordant.csv")
if (file.exists(ar01_conc_path)) {
  go_conc <- read_go_csv(ar01_conc_path)
  if (nrow(go_conc) > 0) {
    p <- go_dotplot(go_conc, "GO Biological Process enrichment: ATAC-RNA concordant genes",
                     "Top 15 terms per group | FDR < 0.05, |LFC| > 1", facet_var = "Cluster")
    save_plot(p, "AR02_GO_ATAC_RNA_concordant", width = 15, height = 9)
    made_any <- TRUE
  }
} else cat("[AR02] Skipping AR01 concordant GO plot -- ", ar01_conc_path, " not found (run AR01 first).\n")

########
## AR04 -- 4-way RNA DEG x ATAC DAR-gene Venn intersections (any-direction,
## UP-UP-UP-UP, DOWN-DOWN-DOWN-DOWN)
########
ar04_path <- file.path(TABLES_DIR, "AR04_GO_venn4_intersections.csv")
if (file.exists(ar04_path)) {
  go_venn4 <- read_go_csv(ar04_path)
  if (nrow(go_venn4) > 0) {
    p <- go_dotplot(go_venn4, "GO Biological Process enrichment: RNA DEG x ATAC DAR-gene 4-way intersections",
                     "Top 15 terms per intersection | FDR < 0.05, |LFC| > 1", facet_var = "Set")
    save_plot(p, "AR02_GO_venn4_intersections", width = 16, height = 9)
    made_any <- TRUE
  }
} else cat("[AR02] Skipping AR04 4-way-Venn GO plot -- ", ar04_path, " not found (run AR04 first).\n")

########
## R05 -- GSEA against three gene-set databases (GO BP, Hallmark, Reactome),
## ranked by Wald statistic, plus each one's KS-relevant curated subset. The
## curated subset CAN legitimately be empty (0 rows, e.g. no Hallmark term
## matches the KS-relevant regex for a given contrast) -- handled by the
## nrow==0 skip below, so a plain readr::read_csv() is fine here (no
## empty-table colClasses gotcha like the ORA tables above, which can be
## empty AND still need specific column types for bind_rows()). Excludes
## R05_fgsea_summary.csv, which matches the "fgsea" glob pattern too but
## isn't a results table. gsea_dotplot() returns NULL if nothing survives its
## own FDR filter -- skipped rather than saving an empty/broken plot.
########
R05_DBS <- c(fgsea = "GO BP", hallmark = "Hallmark", reactome = "Reactome")
for (db in names(R05_DBS)) {
  db_label <- R05_DBS[[db]]
  files <- Sys.glob(file.path(TABLES_DIR, sprintf("R05_%s_*.csv", db)))
  files <- files[!grepl("_summary\\.csv$", files)]  # R05_fgsea_summary.csv matches the fgsea glob too -- not a results table (no p.adjust column)
  for (f in files) {
    gsea_df <- read_csv(f, show_col_types = FALSE)
    if (nrow(gsea_df) == 0) next
    ct <- unique(gsea_df$Comparison)[1]
    is_curated <- grepl("_curated\\.csv$", f)
    p <- gsea_dotplot(gsea_df, sprintf("%s GSEA%s: %s vs WT", db_label, if (is_curated) " (KS-relevant curated terms)" else "", ct),
                       sprintf("Top %s by |NES| | FDR < 0.05, ranked by Wald statistic", if (is_curated) "curated terms" else "15 terms per direction"),
                       n = if (is_curated) 50 else 15)
    if (is.null(p)) { cat(sprintf("[AR02] Skipping %s -- no FDR-significant rows.\n", basename(f))); next }
    save_plot(p, sprintf("AR02_GSEA_%s_%s%s", db, ct, if (is_curated) "_curated" else ""), width = 13, height = 7)
    made_any <- TRUE
  }
  if (length(files) == 0) cat(sprintf("[AR02] Skipping R05 %s GSEA plots -- no R05_%s_*.csv found.\n", db_label, db))
}

if (!made_any) cat("[AR02] No GO tables found yet -- nothing to plot. Run R02/R04/R05/AR01 first.\n")

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR02_GO_Plots_session_info.txt"))
cat("\n[DONE] AR02_GO_Plots complete\n")
