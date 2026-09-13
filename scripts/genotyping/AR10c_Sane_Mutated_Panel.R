#!/usr/bin/env Rscript
# AR10c_Sane_Mutated_Panel.R -- combines the KMT2D and KDM6A publication-ready
# allele-quantification panels (from AR10b_Sane_Transcript_DESeq2.R --locus
# <GENE>, run for both genes first) into one side-by-side figure.
#
#   conda activate ks_1_2_r
#   Rscript scripts/AR10b_Sane_Transcript_DESeq2.R --locus KMT2D
#   Rscript scripts/AR10b_Sane_Transcript_DESeq2.R --locus KDM6A
#   Rscript scripts/AR10c_Sane_Mutated_Panel.R

suppressMessages({
  library(patchwork)
  library(ggplot2)
})

TABLES_DIR <- "results/tables"
PLOTS_DIR  <- "results/figures/genotype_validation/AR10_Combined_Mutation_Loci"
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

kmt2d_rds <- file.path(TABLES_DIR, "AR10_KMT2D_sane_and_mutated_pub_plotobj.rds")
kdm6a_rds <- file.path(TABLES_DIR, "AR10_KDM6A_sane_and_mutated_pub_plotobj.rds")
if (!file.exists(kmt2d_rds) || !file.exists(kdm6a_rds)) {
  stop("Missing plot object(s) -- run AR10b_Sane_Transcript_DESeq2.R --locus KMT2D and --locus KDM6A first.")
}

p_kmt2d <- readRDS(kmt2d_rds)
p_kdm6a <- readRDS(kdm6a_rds)

panel <- (p_kmt2d | p_kdm6a) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(file.path(PLOTS_DIR, "AR10_KMT2D_KDM6A_allele_quantification_panel.png"), panel, width = 10, height = 5.2, dpi = 300)
ggsave(file.path(PLOTS_DIR, "AR10_KMT2D_KDM6A_allele_quantification_panel.pdf"), panel, width = 10, height = 5.2)
cat("Wrote", file.path(PLOTS_DIR, "AR10_KMT2D_KDM6A_allele_quantification_panel.png/.pdf"), "\n")
