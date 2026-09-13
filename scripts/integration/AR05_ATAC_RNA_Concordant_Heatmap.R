#!/usr/bin/env Rscript
# AR05_ATAC_RNA_Concordant_Heatmap.R -- two-panel (RNA log2FC | ATAC Fold)
# ComplexHeatmap of genes concordantly changed in ATAC accessibility AND RNA
# expression, in BOTH genotypes (AR01_shared_targets_only.csv, n=103), row-
# split by category.
#
# Style/structure adapted from a reference script the user supplied (from a
# prior, different project/dataset) -- adapted here to this project's actual
# paths, column names, and ATAC_Fold sign convention (already
# genotype-relative and correctly signed here: GAINED > 0, LOST < 0 by
# construction in AR01, so unlike the reference no direction-based sign flip
# is needed). The reference script's functional-category labels (BMP/WNT
# signaling, ECM/Mesenchymal, Forebrain/NE, etc.) are reused ONLY for genes
# that are actually present in this dataset's own shared-concordant set (many
# are -- independent cross-validation of the earlier hand curation, done on a
# different dataset). Genes not covered by that hand-curated list are labeled
# "Other (Up)"/"Other (Down)" rather than assigning an unverified category.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR05_ATAC_RNA_Concordant_Heatmap.R
# Requires: R02_DESeq2.R and AR01_ATAC_RNA_Integration.R already run.
#
# Self-checkpointing: skips entirely if results/integration/AR05_concordant_heatmap.pdf exists.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
  library(ComplexHeatmap); library(circlize)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR05] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/integration"
SESSION_DIR <- "results/Session_info"
OUT_PDF <- file.path(PLOTS_DIR, "AR05_concordant_heatmap.pdf")
for (d in c(PLOTS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(OUT_PDF)) {
  cat("[AR05] Already complete (", OUT_PDF, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

CONTRASTS <- cfg$contrasts
NAME_A <- CONTRASTS[[1]]$name; NAME_B <- CONTRASTS[[2]]$name
LABEL_A <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment

cat("================================================================================\n")
cat("AR05_ATAC_RNA_Concordant_Heatmap -- shared concordant targets, both genotypes\n")
cat("================================================================================\n\n")

shared <- read_csv(file.path(TABLES_DIR, "AR01_shared_targets_only.csv"), show_col_types = FALSE)
cat(sprintf("[AR05] %d shared concordant genes (both genotypes, matching direction)\n", nrow(shared)))

rna_a <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A)), show_col_types = FALSE)
rna_b <- read_csv(file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B)), show_col_types = FALSE)

atac_tabs <- list(
  A_GAINED = read_csv(file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_GAINED.csv", NAME_A)), show_col_types = FALSE),
  A_LOST   = read_csv(file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_LOST.csv",   NAME_A)), show_col_types = FALSE),
  B_GAINED = read_csv(file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_GAINED.csv", NAME_B)), show_col_types = FALSE),
  B_LOST   = read_csv(file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_LOST.csv",   NAME_B)), show_col_types = FALSE)
)

# Each shared gene is concordant in exactly one direction per genotype (the
# In_*_GAINED / In_*_LOST flags are mutually exclusive within a genotype by
# AR01's own construction) -- pick the matching ATAC_Fold from whichever
# direction table it's flagged TRUE in.
get_atac_fold <- function(symbol, name, gained_col, lost_col) {
  if (isTRUE(shared[[gained_col]][shared$SYMBOL == symbol])) {
    tab <- atac_tabs[[paste0(name, "_GAINED")]]
  } else if (isTRUE(shared[[lost_col]][shared$SYMBOL == symbol])) {
    tab <- atac_tabs[[paste0(name, "_LOST")]]
  } else return(NA_real_)
  v <- tab$ATAC_Fold[tab$SYMBOL == symbol]
  if (length(v) == 0) NA_real_ else v[1]
}

hm_data <- shared %>%
  rowwise() %>%
  mutate(
    KDM6A_ATAC_fold = get_atac_fold(SYMBOL, "A", "In_KDM6A_ko_vs_WT_GAINED", "In_KDM6A_ko_vs_WT_LOST"),
    KMT2D_ATAC_fold = get_atac_fold(SYMBOL, "B", "In_KMT2D_Het_vs_WT_GAINED", "In_KMT2D_Het_vs_WT_LOST"),
    Direction = if (In_KDM6A_ko_vs_WT_GAINED || In_KMT2D_Het_vs_WT_GAINED) "Up" else "Down"
  ) %>%
  ungroup() %>%
  left_join(rna_a %>% transmute(SYMBOL = gene_name, KDM6A_RNA_lfc = log2FoldChange), by = "SYMBOL") %>%
  left_join(rna_b %>% transmute(SYMBOL = gene_name, KMT2D_RNA_lfc = log2FoldChange), by = "SYMBOL") %>%
  filter(!is.na(KDM6A_RNA_lfc), !is.na(KMT2D_RNA_lfc)) %>%
  filter(!SYMBOL %in% c("LRRC61", "PEG3"))

cat(sprintf("[AR05] %d genes with complete RNA+ATAC data for both genotypes\n\n", nrow(hm_data)))

# Reference categories -- revised curation (user-supplied, 2026-07-29),
# replacing the earlier cross-dataset hand curation. "Other" in the source
# list is direction-agnostic by name, but suffixed (Up)/(Down) here from
# hm_data's own computed Direction (not the source list's Direction column)
# so it still row-splits into the Up-block / Down-block the same way every
# other, direction-specific category already does.
ref_categories <- tribble(
  ~SYMBOL,    ~category,
  "BAMBI",    "BMP / TGF-beta Signalling", "TGFB2",  "BMP / TGF-beta Signalling", "BMP4",   "BMP / TGF-beta Signalling",
  "BMPR1B",   "BMP / TGF-beta Signalling", "ID4",    "BMP / TGF-beta Signalling",
  "WNT3A",    "EMT / Neural Plate Border", "ALX4",   "EMT / Neural Plate Border", "DACH2",  "EMT / Neural Plate Border",
  "MEF2C",    "EMT / Neural Plate Border", "PCGF5",  "EMT / Neural Plate Border", "SEMA3D", "EMT / Neural Plate Border",
  "TSHZ2",    "EMT / Neural Plate Border",
  "ADAMTS5",  "ECM / Mesenchymal", "CDH6",   "ECM / Mesenchymal", "CTNNA2", "ECM / Mesenchymal",
  "FN1",      "ECM / Mesenchymal", "ITGA1",  "ECM / Mesenchymal", "LAMA4",  "ECM / Mesenchymal",
  "ALPK2",    "Mechanosensing", "CAV1",   "Mechanosensing", "PIEZO2", "Mechanosensing", "RND3", "Mechanosensing",
  "ATOH1",    "Sensory / Peripheral Neurogenesis", "CNTN3",  "Sensory / Peripheral Neurogenesis",
  "GRIN3A",   "Sensory / Peripheral Neurogenesis", "KCNJ6",  "Sensory / Peripheral Neurogenesis",
  "NCAM2",    "Sensory / Peripheral Neurogenesis", "NEFM",   "Sensory / Peripheral Neurogenesis",
  "NTRK2",    "Sensory / Peripheral Neurogenesis", "SCN3A",  "Sensory / Peripheral Neurogenesis",
  "STMN2",    "Sensory / Peripheral Neurogenesis",
  "ANGPT1",   "Paracrine & Microenvironment", "FAT4",   "Paracrine & Microenvironment",
  "IL1R1",    "Paracrine & Microenvironment", "KITLG",  "Paracrine & Microenvironment",
  "CYP26A1",  "Retinoic Acid Clearance", "CYP26C1", "Retinoic Acid Clearance",
  "AFF3",     "Other", "BICC1",  "Other", "CBLN2",  "Other", "FRMD3",  "Other", "GCK",    "Other",
  "GMNC",     "Other", "JPH2",   "Other", "LYPD1",  "Other", "MMRN1",  "Other", "PDE11A", "Other",
  "PDE3A",    "Other", "PLK2",   "Other", "S100A6", "Other", "SLC13A4","Other", "TMEM163","Other",
  "VWC2",     "Other", "ZPLD1",  "Other",
  "DMBX1",    "Anterior Neural / Placodal Identity", "EN1",    "Anterior Neural / Placodal Identity",
  "IGFBP5",   "Anterior Neural / Placodal Identity", "LHX2",   "Anterior Neural / Placodal Identity",
  "DLX3",     "Anterior Neural / Placodal Identity", "DLX6",   "Anterior Neural / Placodal Identity",
  "CDH3",     "Epithelial Integrity", "GRHL2",  "Epithelial Integrity", "PRSS8",  "Epithelial Integrity",
  "RIPK4",    "Epithelial Integrity", "SFN",    "Epithelial Integrity",
  "CSPG5",    "CNS Lineage & Guidance", "CXCL14", "CNS Lineage & Guidance", "SEMA3E", "CNS Lineage & Guidance",
  "SLITRK1",  "CNS Lineage & Guidance", "WNT8A",  "CNS Lineage & Guidance", "ST8SIA6","CNS Lineage & Guidance",
  "CHRM4",    "GPCR & Synaptic Excitability", "GNAO1",  "GPCR & Synaptic Excitability", "GRM3",   "GPCR & Synaptic Excitability",
  "HRH3",     "GPCR & Synaptic Excitability", "HTR1D",  "GPCR & Synaptic Excitability", "KCNB1",  "GPCR & Synaptic Excitability",
  "KCNH2",    "GPCR & Synaptic Excitability", "SCN5A",  "GPCR & Synaptic Excitability", "GRIK3",  "GPCR & Synaptic Excitability",
  "DMTN",     "Cytoskeletal / Trafficking", "FHDC1",  "Cytoskeletal / Trafficking", "MYO3B",  "Cytoskeletal / Trafficking",
  "MYO5B",    "Cytoskeletal / Trafficking", "PITPNM3","Cytoskeletal / Trafficking",
  "A4GALT",   "Other", "CAPN6",  "Other", "CCM2L",  "Other", "COL13A1","Other", "DHRS2",  "Other",
  "ISM1",     "Other", "KLHDC8A","Other", "L1TD1",  "Other", "MACC1",  "Other",
  "NANOS1",   "Other", "OSBP2",  "Other", "PEG3",   "Other", "PPP1R17","Other", "PPP2R2C","Other",
  "RAB25",    "Other", "LRRC61", "Other", "PIK3R5", "Other"
)

hm_data <- hm_data %>%
  left_join(ref_categories, by = "SYMBOL") %>%
  mutate(category = ifelse(is.na(category) | category == "Other", paste0("Other (", Direction, ")"), category))

n_curated <- sum(hm_data$SYMBOL %in% ref_categories$SYMBOL)
cat(sprintf("[AR05] %d/%d genes matched a reference-curated category; %d assigned to Other (Up/Down)\n\n",
            n_curated, nrow(hm_data), nrow(hm_data) - n_curated))

cat_levels <- c("BMP / TGF-beta Signalling", "EMT / Neural Plate Border", "ECM / Mesenchymal", "Mechanosensing",
                "Sensory / Peripheral Neurogenesis", "Paracrine & Microenvironment", "Retinoic Acid Clearance", "Other (Up)",
                "Anterior Neural / Placodal Identity", "Epithelial Integrity", "CNS Lineage & Guidance",
                "GPCR & Synaptic Excitability", "Cytoskeletal / Trafficking", "Other (Down)")
cat_levels <- cat_levels[cat_levels %in% unique(hm_data$category)]
hm_data <- hm_data %>% mutate(category = factor(category, levels = cat_levels)) %>%
  arrange(category, desc(abs(KMT2D_RNA_lfc)))

# Underlying data table, not just the PDF -- every other figure in this
# pipeline has a backing table; this one didn't.
write_csv(hm_data %>% dplyr::select(SYMBOL, category, Direction, KMT2D_RNA_lfc, KDM6A_RNA_lfc, KMT2D_ATAC_fold, KDM6A_ATAC_fold),
          file.path(TABLES_DIR, "AR05_concordant_heatmap_data.csv"))

mat_rna <- as.matrix(hm_data[, c("KMT2D_RNA_lfc", "KDM6A_RNA_lfc")])
rownames(mat_rna) <- hm_data$SYMBOL; colnames(mat_rna) <- c(LABEL_B, LABEL_A)

mat_atac <- as.matrix(hm_data[, c("KMT2D_ATAC_fold", "KDM6A_ATAC_fold")])
rownames(mat_atac) <- hm_data$SYMBOL; colnames(mat_atac) <- c(LABEL_B, LABEL_A)
mat_atac[is.na(mat_atac)] <- 0

# One shared colour scale (global max across RNA and ATAC) so RNA and ATAC
# can carry a single "Log2FC" legend that doesn't need to say which panel
# it's from -- a per-matrix scale would make one shared legend meaningless
# (same colour, two different value ranges).
shared_max <- max(abs(mat_rna), abs(mat_atac), na.rm = TRUE)
col_shared <- colorRamp2(c(-shared_max, 0, shared_max), c("#2166AC", "white", "#B2182B"))

row_split <- hm_data$category
# 14 categories now (was 11) -- extended qualitative pool, ordered so that
# only ADJACENT categories in the fixed row_split order need to read apart
# (this is a stacked, ordered layout, not a scatter/all-pairs case);
# row_title text labels each block directly, so color reinforces the
# grouping rather than being the sole identity carrier.
cat_colors_pool <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b15928", "#999999",
                      "#fb9a99", "#fdbf6f", "#cab2d6", "#01665e", "#b2df8a", "#666666")
cat_colors <- setNames(cat_colors_pool[seq_along(cat_levels)], cat_levels)

# Category legend dropped (show_legend = FALSE) -- row_title already labels
# each block directly, so a separate colour-key legend is redundant; the
# colour bar itself stays as a grouping cue.
row_anno <- rowAnnotation(Category = hm_data$category, col = list(Category = cat_colors),
                           show_annotation_name = FALSE, show_legend = FALSE, width = unit(5, "mm"))

# Same `name` + same shared colour scale on both heatmaps -> ComplexHeatmap
# merges them into ONE "Log2FC" legend, with no RNA/ATAC distinction (per
# request: the legend shouldn't say which panel it's from).
ht_rna <- Heatmap(mat_rna, name = "Log2FC", col = col_shared,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   row_split = row_split, row_gap = unit(2, "mm"),
                   row_title_rot = 0, row_title_gp = gpar(fontsize = 7, fontface = "bold"),
                   column_title = "RNA-seq", column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                   show_row_names = FALSE, width = unit(18, "mm"), border = TRUE,
                   cell_fun = function(j, i, x, y, w, h, fill) grid.rect(x, y, w, h, gp = gpar(col = "grey80", lwd = 0.5, fill = fill)))

ht_atac <- Heatmap(mat_atac, name = "Log2FC", col = col_shared,
                    show_heatmap_legend = FALSE,
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    row_split = row_split, row_gap = unit(2, "mm"),
                    row_title_rot = 0, row_title_gp = gpar(fontsize = 7, fontface = "bold"),
                    column_title = "ATAC-seq", column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                    width = unit(18, "mm"), border = TRUE,
                    cell_fun = function(j, i, x, y, w, h, fill) grid.rect(x, y, w, h, gp = gpar(col = "grey80", lwd = 0.5, fill = fill)))

# Legend on the LEFT, ahead of the Category colour bar / RNA / ATAC / gene
# names reading order -- and no annotation legend (Category legend removed
# above).
pdf(OUT_PDF, width = 8, height = max(10, nrow(hm_data) * 0.16))
draw(row_anno + ht_rna + ht_atac, row_title = NULL,
     column_title = "Concordant ATAC-RNA targets (both genotypes)",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()
cat(sprintf("[AR05] Saved: %s (%d genes)\n", OUT_PDF, nrow(hm_data)))

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR05_ATAC_RNA_Concordant_Heatmap_session_info.txt"))
cat("\n[DONE] AR05_ATAC_RNA_Concordant_Heatmap complete\n")
