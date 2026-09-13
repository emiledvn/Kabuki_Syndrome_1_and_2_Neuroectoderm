#!/usr/bin/env Rscript
# R06_Generate_Report.R -- assembles a single comprehensive reference document
# from the pipeline's own already-computed output tables: every number here is
# read directly from a results/tables/*.csv|tsv file, never hand-transcribed,
# so the report can be regenerated any time thresholds/methods change instead
# of manually re-verified against the underlying data (the exact problem that
# came up updating docs/RESULTS_DRAFT.md after the LFC-threshold unification).
#
# Sections (per request): 1 figures index + associated metrics, 2 top 100
# DEGs up/down per contrast + intersection, 3 GO/GSEA terms, 4 ATAC numbers,
# 5 ATAC motif enrichment, 6 ATAC-RNA intersection, 7 TOBIAS/chromVAR
# differential motif exploration.
#
# Not self-checkpointing (unlike the data-generating steps): this only reads
# already-computed tables, cheap to rerun, and should always reflect current
# output -- always overwrites docs/COMPLETE_REPORT.md.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R06_Generate_Report.R
# Requires: the full pipeline (R02/R04/A04b/A04f/A05/A07b/A09/A10/AR01/AR02/AR03)
# already run.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R06] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR <- "results/tables"
PLOTS_DIR  <- "results/rna"
DOCS_DIR   <- "docs"
OUT_MD     <- file.path(DOCS_DIR, "COMPLETE_REPORT.md")
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

CONTRASTS   <- cfg$contrasts
NAME_A      <- CONTRASTS[[1]]$name;  NAME_B <- CONTRASTS[[2]]$name
LABEL_A     <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment
FDR         <- cfg$thresholds$fdr
LFC_THR     <- cfg$thresholds$lfc_threshold

read_if <- function(path, ...) if (file.exists(path)) read_csv(path, show_col_types = FALSE, ...) else NULL
read_tsv_if <- function(path, ...) if (file.exists(path)) read_tsv(path, show_col_types = FALSE, ...) else NULL
fmt <- function(x, digits = 3) if (is.numeric(x)) formatC(x, digits = digits, format = "g") else as.character(x)
md_table <- function(df, max_rows = NULL) {
  if (is.null(df) || nrow(df) == 0) return("*(no rows)*")
  if (!is.null(max_rows)) df <- head(df, max_rows)
  df <- df %>% mutate(across(where(is.numeric), ~ fmt(.x)))
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

lines <- c(
  "# KS1_2_Neuro-ectoderm -- Complete Results Reference",
  "",
  sprintf("Generated %s by `scripts/R06_Generate_Report.R` -- every number below is read directly from its source table (paths given), not transcribed by hand.", format(Sys.Date(), "%Y-%m-%d")),
  "",
  sprintf("Shared thresholds: FDR (padj) < %.2f | |log2FC| > %.2f (single threshold, applied identically to RNA DEGs, ATAC DARs, and all downstream gene/peak lists).", FDR, LFC_THR),
  ""
)

########################################################################
## 1. Figures index + associated metrics
########################################################################
lines <- c(lines, "## 1. Figures index and associated metrics", "")

deg_summary <- read_if(file.path(TABLES_DIR, "R02_DEG_summary.csv"))
dar_summary <- read_if(file.path(TABLES_DIR, "A05_DAR_summary.csv"))
rna_summary_txt <- if (file.exists("results/RNA_ANALYSIS_SUMMARY.txt")) readLines("results/RNA_ANALYSIS_SUMMARY.txt") else character(0)
atac_summary_txt <- if (file.exists("results/ATAC_ANALYSIS_SUMMARY.txt")) readLines("results/ATAC_ANALYSIS_SUMMARY.txt") else character(0)
grab <- function(txt, pattern) { m <- grep(pattern, txt, value = TRUE); if (length(m) > 0) trimws(m[1]) else NA_character_ }

lines <- c(lines,
  "### R02 -- RNA differential expression",
  md_table(deg_summary), "",
  "Figures: `R02_pca.pdf`, `R02_sample_correlation.pdf`, `R02_volcano_{contrast}.pdf`, `R02_heatmap_{contrast}_top50.pdf`, `R02_NE_markers_barplot*.pdf` (5 variants), `R02_identity_acquisition_heatmap.pdf`.", "",
  "### R04 -- RNA cross-genotype overlap", "",
  grab(rna_summary_txt, "Genome-wide correlation"), "  ",
  grab(rna_summary_txt, "DEG overlap effect size"), "  ",
  grab(rna_summary_txt, "DEG overlap significance"), "",
  "Figures: `R04_venn_{up,down}.pdf`, `R04_scatter_convergence.pdf`, `R04_genome_wide_lfc_correlation.pdf`, `R04_heatmap_concordant_overlap.pdf`.", "",
  "### A04b/A04f/A05 -- ATAC differential accessibility (Loess-normalized)",
  md_table(dar_summary), "",
  "Figures: `A04b_MA_{contrast}_{Loess,TMM}.pdf`, `A04d_*_TMM_vs_Loess.pdf` (normalization QC), `A04e_MA_WTvsWT_Null_*.pdf` (null check), `A05_pca.pdf`, `A05_volcano_{contrast}.pdf`, `A05_MA_{contrast}.pdf`.", "",
  "### A08 -- ATAC signal visualization",
  "Figures: `A08_{contrast}_{persample,merged}.pdf`, `A08_shared_DARs_merged.pdf`, `A08_union_DARs_merged.pdf`, `A08_union_DARs_by_distance.pdf`.", "",
  "### AR03 -- ATAC cross-genotype overlap",
  grab(atac_summary_txt, "Genome-wide correlation"), "  ",
  grab(atac_summary_txt, "DAR overlap effect size"), "  ",
  grab(atac_summary_txt, "DAR overlap significance"), "",
  "Figures: `AR03_venn_{up,down}.pdf`, `AR03_scatter_convergence.pdf`, `AR03_genome_wide_lfc_correlation.pdf`.", "",
  "### A09 -- motif enrichment (monaLisa)", "",
  "Figures: `A09_{contrast}_{GAINED,LOST}_dotplot.pdf`.", "",
  "### A10 -- chromVAR differential motif accessibility",
  "Figures: `A10_{contrast}_dotplot.pdf`.", "",
  "### A07b -- TOBIAS differential TF binding",
  "Figures: `A07b_top_tfs_barplot.pdf`, `A07b_tf_overlap_venn.pdf`, `A07b_binding_vs_expression_{contrast}.pdf`, `A07b_volcano_{contrast}.pdf`, `A07b_heatmap_top50.pdf`, `A07b_heatmap_extended.pdf`, `A07b_heatmap_both_sig.pdf`, `A07b_concordant_tfs_summary.pdf`.", "",
  "### AR01 -- ATAC-RNA integration (GREAT)",
  "Figures: `AR01_great_distance_to_TSS.pdf`, `AR01_great_peaks_per_gene.pdf`, `AR01_DAR_counts.pdf`, `AR01_DAR_genomic_annotation.pdf`, `AR01_great_scatter_{contrast}_combined.pdf`, `AR01_overlap_{venn,upset}.pdf`.", "",
  "### AR02 -- GO/GSEA plotting layer",
  "Figures: `AR02_GO_RNA_{contrast}.pdf`, `AR02_GO_RNA_overlap.pdf`, `AR02_GO_ATAC_targets.pdf`, `AR02_GO_ATAC_RNA_concordant.pdf`, `AR02_GO_venn4_intersections.pdf`, `AR02_GSEA_{fgsea,hallmark,reactome}_{contrast}[_curated].pdf`.", "",
  "### AR04/AR05 -- 4-way RNA-ATAC Venn and concordant-gene heatmap",
  "Figures: `AR04_venn4_{any,up,down}.pdf`, `AR05_concordant_heatmap.pdf` (full detail in Section 6).", ""
)

########################################################################
## 2. Top 100 DEGs up/down per contrast + intersection
########################################################################
lines <- c(lines, "## 2. Top DEGs (up to 100 per direction), by contrast and shared", "")

top_deg_block <- function(name) {
  full <- read_if(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", name)))
  if (is.null(full)) return(c(sprintf("*(%s: no data)*", name), ""))
  up   <- full %>% filter(!is.na(padj), padj < FDR, log2FoldChange >  LFC_THR) %>% arrange(padj) %>% head(100) %>%
    dplyr::select(gene_name, log2FoldChange, padj, baseMean)
  down <- full %>% filter(!is.na(padj), padj < FDR, log2FoldChange < -LFC_THR) %>% arrange(padj) %>% head(100) %>%
    dplyr::select(gene_name, log2FoldChange, padj, baseMean)
  n_up_total   <- sum(!is.na(full$padj) & full$padj < FDR & full$log2FoldChange >  LFC_THR)
  n_down_total <- sum(!is.na(full$padj) & full$padj < FDR & full$log2FoldChange < -LFC_THR)
  c(sprintf("### %s (%d total up, %d total down; showing up to top 100 each by padj)", name, n_up_total, n_down_total), "",
    "**Top UP:**", "", md_table(up), "",
    "**Top DOWN:**", "", md_table(down), "")
}
lines <- c(lines, top_deg_block(NAME_A), top_deg_block(NAME_B))

overlap_block <- function(direction, file) {
  df <- read_if(file.path(TABLES_DIR, file))
  if (is.null(df)) return(c(sprintf("*(shared %s: no data)*", direction), ""))
  df <- df %>% arrange(Sig_Score) %>% head(100) %>%
    dplyr::select(gene_name, LFC_A, P_A, LFC_B, P_B)
  names(df)[names(df) == "LFC_A"] <- sprintf("LFC_%s", NAME_A); names(df)[names(df) == "P_A"] <- sprintf("padj_%s", NAME_A)
  names(df)[names(df) == "LFC_B"] <- sprintf("LFC_%s", NAME_B); names(df)[names(df) == "P_B"] <- sprintf("padj_%s", NAME_B)
  c(sprintf("### Shared %s DEGs (both contrasts; up to top 100 by combined significance)", direction), "", md_table(df), "")
}
lines <- c(lines,
  overlap_block("UP", "R04_overlap_up_genes.csv"),
  overlap_block("DOWN", "R04_overlap_down_genes.csv"))

########################################################################
## 3. GO / GSEA terms
########################################################################
lines <- c(lines, "## 3. GO enrichment and GSEA terms", "")
go_block <- function(title, path, n = 15, comparison_val = NULL, group_col = "Direction") {
  df <- read_if(file.path(TABLES_DIR, path))
  if (is.null(df) || nrow(df) == 0) return(c(sprintf("### %s", title), "", "*(no data)*", ""))
  if (!is.null(comparison_val) && "Comparison" %in% names(df)) df <- df %>% filter(Comparison == comparison_val)
  if (nrow(df) == 0) return(c(sprintf("### %s", title), "", "*(no significant terms)*", ""))
  out <- c(sprintf("### %s", title), "")
  for (g in sort(unique(df[[group_col]]))) {
    sub <- df %>% filter(.data[[group_col]] == g) %>% arrange(p.adjust) %>% head(n) %>%
      dplyr::select(Description, GeneRatio, p.adjust, Count)
    out <- c(out, sprintf("**%s** (top %d by FDR):", g, n), "", md_table(sub), "")
  }
  out
}

lines <- c(lines, go_block(sprintf("R02 -- RNA GO BP, %s", NAME_A), "R02_GO_enrichment_per_genotype.csv", comparison_val = NAME_A))
lines <- c(lines, go_block(sprintf("R02 -- RNA GO BP, %s", NAME_B), "R02_GO_enrichment_per_genotype.csv", comparison_val = NAME_B))
lines <- c(lines, go_block("R04 -- RNA GO BP, shared/overlap DEGs", "R04_GO_overlap.csv"))
lines <- c(lines, go_block("AR01 -- ATAC target genes GO BP (GREAT)", "AR01_great_GO_results.csv", group_col = "Cluster"))
lines <- c(lines, go_block("AR01 -- ATAC-RNA concordant genes GO BP", "AR01_great_GO_concordant.csv", group_col = "Cluster"))

gsea_block <- function(title, path, n = 15) {
  df <- read_if(file.path(TABLES_DIR, path))
  if (is.null(df) || nrow(df) == 0) return(c(sprintf("### %s", title), "", "*(no significant terms)*", ""))
  df <- df %>% filter(!is.na(p.adjust), p.adjust < FDR) %>% mutate(Direction = ifelse(NES > 0, "Up (NES>0)", "Down (NES<0)"))
  out <- c(sprintf("### %s", title), "")
  for (g in c("Up (NES>0)", "Down (NES<0)")) {
    sub <- df %>% filter(Direction == g) %>% arrange(desc(abs(NES))) %>% head(n) %>%
      dplyr::select(Description, NES, p.adjust, setSize)
    out <- c(out, sprintf("**%s** (top %d by |NES| among FDR-significant):", g, n), "", md_table(sub), "")
  }
  out
}
for (ct in list(list(n = NAME_A), list(n = NAME_B))) {
  nm <- ct$n
  lines <- c(lines, gsea_block(sprintf("R05 -- GSEA GO BP (curated, KS-relevant), %s", nm), sprintf("R05_fgsea_%s_curated.csv", nm)))
  lines <- c(lines, gsea_block(sprintf("R05 -- GSEA Hallmark (curated), %s", nm), sprintf("R05_hallmark_%s_curated.csv", nm)))
  lines <- c(lines, gsea_block(sprintf("R05 -- GSEA Reactome (curated), %s", nm), sprintf("R05_reactome_%s_curated.csv", nm)))
}

########################################################################
## 4. ATAC numbers
########################################################################
lines <- c(lines, "## 4. ATAC accessibility numbers", "")

a04b <- read_tsv_if(file.path(TABLES_DIR, "A04b_summary.tsv"))
lines <- c(lines,
  "### Normalization comparison (A04b, FDR-only, both norm methods -- context for the Loess-vs-TMM methodology decision, NOT the final DAR definition)",
  "", md_table(a04b), "",
  sprintf("### Final DAR calling (Loess-normalized, padj < %.2f AND |log2FC| > %.2f) -- A05", FDR, LFC_THR),
  "", md_table(dar_summary), "")

for (ct in CONTRASTS) {
  gbed <- sprintf("data/diffbind_output/CsawLoess_Norm/%s_GAINED.bed", ct$name)
  lbed <- sprintf("data/diffbind_output/CsawLoess_Norm/%s_LOST.bed", ct$name)
  n_g <- if (file.exists(gbed)) length(readLines(gbed)) else NA
  n_l <- if (file.exists(lbed)) length(readLines(lbed)) else NA
  lines <- c(lines, sprintf("A04f GAINED/LOST BED files, %s: %d GAINED, %d LOST (matches A05 Up/Down exactly by construction -- same padj+LFC filter).", ct$name, n_g, n_l))
}
lines <- c(lines, "",
  "### Cross-genotype ATAC overlap (AR03)", "",
  grab(atac_summary_txt, "Genome-wide correlation"), "  ",
  grab(atac_summary_txt, "DAR overlap effect size"), "  ",
  grab(atac_summary_txt, "DAR overlap significance"), "  ",
  grab(atac_summary_txt, "Shared GAINED"), "  ",
  grab(atac_summary_txt, "Shared LOST"), "")

########################################################################
## 5. ATAC motif enrichment (A09, monaLisa)
########################################################################
lines <- c(lines, "", "## 5. ATAC motif enrichment (A09, monaLisa: sequence enrichment in GAINED/LOST DARs vs background)", "")
a09_summary <- read_if(file.path(TABLES_DIR, "A09_summary.tsv"))  # write_csv() despite .tsv extension in A09 itself -- reading as CSV to match actual content
lines <- c(lines, md_table(a09_summary), "")
for (ct in CONTRASTS) for (dir in c("GAINED", "LOST")) {
  df <- read_if(file.path(TABLES_DIR, sprintf("A09_%s_%s_vs_Background.csv", ct$name, dir)))
  if (is.null(df)) next
  sig <- df %>% filter(negLog10Padj > -log10(FDR)) %>% arrange(desc(abs(log2enr))) %>% head(15) %>%
    dplyr::select(motif_name, log2enr, negLog10Padj, n_fg_with_hit, n_bg_with_hit)
  lines <- c(lines, sprintf("### %s, %s (top 15 by |log2 enrichment| among FDR-significant motifs)", ct$name, dir), "", md_table(sig), "")
}

########################################################################
## 6. ATAC-RNA intersection (AR01, GREAT)
########################################################################
lines <- c(lines, "## 6. ATAC-RNA intersection (AR01, GREAT peak-to-gene annotation)", "")
all_peaks <- read_if(file.path(TABLES_DIR, "AR01_all_great_peaks.csv"))
if (!is.null(all_peaks)) {
  lines <- c(lines, sprintf("Total peak-to-gene links: %d | Unique genes with a GREAT assignment: %d", nrow(all_peaks), length(unique(all_peaks$SYMBOL))), "")
}
for (ct in CONTRASTS) for (dir in c("GAINED", "LOST")) {
  df <- read_if(file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_%s.csv", ct$name, dir)))
  if (is.null(df)) next
  tab <- df %>% count(concordance, name = "N") %>% arrange(desc(N))
  lines <- c(lines, sprintf("### %s %s -- ATAC-RNA concordance classification", ct$name, dir), "", md_table(tab), "")
}
shared <- read_if(file.path(TABLES_DIR, "AR01_shared_targets_only.csv"))
if (!is.null(shared)) {
  lines <- c(lines, sprintf("### Shared concordant target genes across BOTH contrasts: %d", nrow(shared)), "",
             md_table(shared %>% dplyr::select(SYMBOL, Total_Conditions) %>% arrange(desc(Total_Conditions)), max_rows = 200), "")
}

lines <- c(lines, "### AR04 -- 4-way RNA DEG x ATAC DAR-gene Venn (both genotypes simultaneously)", "",
  "Figures: `AR04_venn4_any.pdf` (any direction), `AR04_venn4_up.pdf` (UP-UP-UP-UP), `AR04_venn4_down.pdf` (DOWN-DOWN-DOWN-DOWN). Figure `AR05_concordant_heatmap.pdf`: two-panel (RNA log2FC | ATAC Fold) heatmap of the shared concordant target genes, row-grouped by category.", "")
for (f in c("AR04_venn4_any_intersection.csv", "AR04_venn4_up_intersection.csv", "AR04_venn4_down_intersection.csv")) {
  df <- read_if(file.path(TABLES_DIR, f))
  if (is.null(df)) next
  lines <- c(lines, sprintf("**%s**: %d genes -- %s", f, nrow(df), paste(df$SYMBOL, collapse = ", ")), "")
}
go_venn4 <- read_if(file.path(TABLES_DIR, "AR04_GO_venn4_intersections.csv"))
if (!is.null(go_venn4)) {
  for (s in unique(go_venn4$Set)) {
    sub <- go_venn4 %>% filter(Set == s) %>% arrange(p.adjust) %>% head(15) %>% dplyr::select(Description, GeneRatio, p.adjust, Count)
    lines <- c(lines, sprintf("**GO BP, %s** (top 15 by FDR):", s), "", md_table(sub), "")
  }
}

ar05_data <- read_if(file.path(TABLES_DIR, "AR05_concordant_heatmap_data.csv"))
if (!is.null(ar05_data)) {
  cat_counts <- ar05_data %>% count(category, Direction, name = "N") %>% arrange(Direction, desc(N))
  lines <- c(lines, "### AR05 -- concordant heatmap category breakdown", "",
             sprintf("%d shared concordant genes (both genotypes), %d/%d matched a reference-curated functional category, %d assigned \"Other\":",
                     nrow(ar05_data), sum(ar05_data$category != paste0("Other (", ar05_data$Direction, ")")), nrow(ar05_data),
                     sum(ar05_data$category == paste0("Other (", ar05_data$Direction, ")"))),
             "", md_table(cat_counts), "",
             "Full gene-level table (category, RNA log2FC and ATAC Fold per genotype): `AR05_concordant_heatmap_data.csv`.", "")
}

########################################################################
## 7. TOBIAS + chromVAR differential motif exploration
########################################################################
lines <- c(lines, "## 7. TOBIAS and chromVAR differential TF/motif exploration", "")

for (ct in CONTRASTS) {
  full_tf <- read_if(file.path(TABLES_DIR, sprintf("A07b_tf_binding_%s.csv", ct$name)))
  if (is.null(full_tf)) next
  tab <- full_tf %>% count(tf_class, name = "N") %>% arrange(desc(N))
  lines <- c(lines, sprintf("### A07b TF classification, %s", ct$name), "", md_table(tab), "")
}
concordant <- read_if(file.path(TABLES_DIR, "A07b_concordant_tfs.csv"))
if (!is.null(concordant)) {
  lines <- c(lines, "### A07b Concordant TFs (all, both contrasts -- binding_sig AND expression_sig, matching direction)", "",
             md_table(concordant %>% dplyr::select(TF_name, contrast, binding_change, binding_direction, log2FoldChange, padj) %>% arrange(contrast, desc(abs(binding_change)))), "")
}
all_diff <- read_if(file.path(TABLES_DIR, "A07b_all_differential_tfs.csv"))
if (!is.null(all_diff)) {
  for (ct in CONTRASTS) {
    sub <- all_diff %>% filter(contrast == ct$name) %>% arrange(desc(abs(binding_change))) %>% head(15) %>%
      dplyr::select(TF_name, binding_change, binding_direction, binding_pvalue, tf_class)
    lines <- c(lines, sprintf("### A07b top 15 differential-binding TFs by |change|, %s", ct$name), "", md_table(sub), "")
  }
}

a10_summary <- read_if(file.path(TABLES_DIR, "A10_summary.tsv"))  # same write_csv()/.tsv-extension mismatch as A09
lines <- c(lines, "### A10 chromVAR differential motif accessibility summary (real per-sample n=3 replicate structure)", "", md_table(a10_summary), "")
for (ct in CONTRASTS) {
  df <- read_if(file.path(TABLES_DIR, sprintf("A10_%s.csv", ct$name)))
  if (is.null(df)) next
  sig <- df %>% filter(p_value_adjusted < FDR) %>% arrange(p_value_adjusted) %>% head(15) %>%
    dplyr::select(motif_name, motif_id, mean_deviation_diff, Direction, p_value_adjusted)
  lines <- c(lines, sprintf("### A10 top 15 significant motifs by FDR, %s", ct$name), "", md_table(sig), "")
}

# Cross-validation: does chromVAR agree with TOBIAS wherever both have power?
# Same analysis as the ad-hoc check earlier in this session, folded into the
# report so it's reproducible rather than a one-off scratch result.
lines <- c(lines, "### Cross-validation: TOBIAS vs chromVAR direction agreement", "",
  "Among TF/motif entries significant in BOTH A07b (TOBIAS differential binding) AND A10 (chromVAR, FDR < 0.05), case-insensitive gene-symbol match, dimer motifs expanded to their component TFs:", "")
if (!is.null(all_diff)) {
  for (ct in CONTRASTS) {
    a10ct <- read_if(file.path(TABLES_DIR, sprintf("A10_%s.csv", ct$name)))
    if (is.null(a10ct)) next
    tobias_ct <- all_diff %>% filter(contrast == ct$name) %>% mutate(TF_upper = toupper(TF_name))
    a10_long <- a10ct %>% tidyr::separate_rows(motif_name, sep = "::") %>% mutate(motif_name_upper = toupper(motif_name))
    merged <- tobias_ct %>% inner_join(a10_long, by = c("TF_upper" = "motif_name_upper"), relationship = "many-to-many")
    both_sig <- merged %>% filter(p_value_adjusted < FDR) %>% mutate(agree = binding_direction == Direction)
    n_agree <- sum(both_sig$agree); n_total <- nrow(both_sig)
    lines <- c(lines, sprintf("- **%s**: %d TF/motif pairs significant in both; direction agreement %d/%d (%s%%)",
                               ct$name, n_total, n_agree, n_total, if (n_total > 0) sprintf("%.0f", 100 * n_agree / n_total) else "NA"))
  }
}
lines <- c(lines, "")

writeLines(lines, OUT_MD)
cat(sprintf("[R06] Wrote %s (%d lines)\n", OUT_MD, length(lines)))
