#!/usr/bin/env Rscript
# AR07_Publication_Figures.R -- publication-ready figure layer, separate from
# the exploratory pipeline scripts. Reads already-computed result tables only
# (R02, R04, AR01, AR05) -- never reruns DESeq2 or any upstream analysis,
# never modifies an existing results file, never changes a threshold (FDR/LFC
# are read from config/pipeline_config.yaml like everywhere else in this
# pipeline). Each task is independent and wrapped so a missing input stops
# only that task, not the whole script.
#
# Language rules enforced in every axis label/title/legend below: ATAC signal
# is "chromatin accessibility", never "binding" (TOBIAS's actual
# binding-depth measurements aren't used in any of these 4 tasks, but the
# rule is kept in force for any ATAC-derived quantity); GREAT-assigned genes
# are labeled "GREAT-inferred target gene".
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR07_Publication_Figures.R
# Requires: R02_DESeq2.R, R04_RNA_Overlap.R, AR01_ATAC_RNA_Integration.R,
# AR05_ATAC_RNA_Concordant_Heatmap.R already run.
#
# Not self-checkpointing: always regenerates all 6 figures (cheap -- reads
# tables, no recomputation), so a re-run always reflects current upstream data.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(stringr)
  library(ggplot2); library(ggrepel)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(yaml); library(scales)
  library(DESeq2); library(pheatmap)
  library(ggtext)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR07] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR <- "results/tables"
PLOTS_DIR  <- "results/figures/rna_transcriptional_phenotype"
SESSION_DIR <- "results/Session_info"
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SESSION_DIR, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- cfg$contrasts
NAME_A <- CONTRASTS[[1]]$name; NAME_B <- CONTRASTS[[2]]$name   # KDM6A_ko_vs_WT, KMT2D_Het_vs_WT
LABEL_A <- CONTRASTS[[1]]$treatment; LABEL_B <- CONTRASTS[[2]]$treatment  # KDM6A_ko, KMT2D_Het
FDR     <- cfg$thresholds$fdr
LFC_THR <- cfg$thresholds$lfc_threshold

GENOTYPE_COLORS <- setNames(c("#56B4E9", "#E69F00"), c("KDM6A_ko", "KMT2D_Het"))

# Facet strip labels for the GO dotplot family's facet_wrap(~Comparison)
# strips -- plain black text (user preference, overriding the earlier
# genotype-coloured-text design). Kept as a function (not just using
# Comparison's own factor labels directly) so every call site stays
# structurally identical if colour-coding is ever wanted back; strip.text
# still renders via ggtext::element_markdown() elsewhere in this script, but
# plain text through that renderer looks identical to element_text().
genotype_strip_label <- function(name) name

theme_pub <- theme_bw(base_size = 10) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())

cat("================================================================================\n")
cat("AR07_Publication_Figures -- publication-ready figure layer\n")
cat("================================================================================\n\n")

########
## Shared helpers
########

check_input <- function(path) {
  if (!file.exists(path)) { cat(sprintf("  [SKIP] input not found: %s\n", path)); return(FALSE) }
  TRUE
}

save_plot_both <- function(p, basename, width, height) {
  pdf_path <- file.path(PLOTS_DIR, paste0(basename, ".pdf"))
  png_path <- file.path(PLOTS_DIR, paste0(basename, ".png"))
  ok <- tryCatch({
    ggsave(pdf_path, p, width = width, height = height)
    ggsave(png_path, p, width = width, height = height, dpi = 300)
    TRUE
  }, error = function(e) { cat(sprintf("  [ERROR] could not save %s: %s\n", basename, conditionMessage(e))); FALSE })
  if (ok) cat(sprintf("  [OK] wrote %s (.pdf + .png)\n", basename))
  ok
}

# GO BP dotplot -- same visual language as AR02_GO_Plots.R::go_dotplot():
# fixed point size, gene count printed as white text inside each point,
# colour = -log10(FDR), enrichment ratio (GeneRatio / BgRatio, i.e.
# background-corrected) on the x-axis, terms ordered by that ratio, labels
# wrapped. Kept as a separate copy (not sourced from AR02) since AR07 is
# read-only over already-computed tables and deliberately has no cross-script
# dependency on another analysis script.
go_dotplot_pub <- function(go_df, title, subtitle, n = 10, facet_var) {
  parse_ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
  df <- go_df %>% mutate(enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  df <- df %>% group_by(.data[[facet_var]]) %>% arrange(p.adjust) %>% slice_head(n = n) %>% ungroup()
  df <- df %>% mutate(label = str_wrap(Description, 40))
  ggplot(df, aes(x = enrichment_ratio, y = reorder(label, enrichment_ratio), colour = -log10(p.adjust))) +
    geom_point(size = 8, alpha = 0.9) +
    geom_text(aes(label = Count), colour = "white", size = 2.8, fontface = "bold") +
    scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
                            name = expression(-log[10](FDR)),
                            labels = scales::label_number(accuracy = 0.1)) +
    facet_wrap(vars(.data[[facet_var]]), scales = "free", ncol = 2) +
    labs(title = title, subtitle = subtitle, x = "Enrichment ratio", y = NULL) +
    theme_pub
}

# Universe = same convention as R02/AR01: genes DESeq2 actually tested (count
# pre-filter), converted to Entrez.
universe_entrez <- NULL
norm_counts_path <- "results/RDS/R02_normalized_counts.rds"
gene_map_path <- "results/RDS/R02_gene_id_name_map.rds"
if (check_input(norm_counts_path) && check_input(gene_map_path)) {
  norm_counts <- readRDS(norm_counts_path)
  gene_map    <- readRDS(gene_map_path)
  universe_symbols <- gene_map$gene_name[gene_map$gene_id %in% rownames(norm_counts)]
  universe_entrez <- tryCatch(bitr(universe_symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
  cat(sprintf("[AR07] GO background universe: %d expressed genes\n\n", length(universe_symbols)))
}

# Collapse near-duplicate GO BP terms via gene-set Jaccard similarity -- same
# approach as build_collapsed_shared_axis_plot()'s TASK 12 redundancy graph
# below (connected components of a Jaccard>=0.7 graph, canonical label =
# lowest-p.adjust member of each cluster), factored out here so TASK 3 can
# reuse it without depending on TASK 12's per-genotype-specific code. Jaccard
# is computed on geneID (readable=TRUE gene symbols) within each Group only --
# overlap across Gained/Lost isn't meaningful since concordant UP and DOWN
# genes are disjoint sets.
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

run_go_bp <- function(genes, label) {
  if (length(genes) < 5 || is.null(universe_entrez)) return(NULL)
  entrez <- tryCatch(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID, error = function(e) character(0))
  if (length(entrez) < 5) return(NULL)
  ego <- enrichGO(gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = FDR, readable = TRUE)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>% filter(p.adjust < FDR) %>% mutate(Group = label)
}

########################################################################
## TASK 1 -- GO BP dotplot for shared (RNA-only) concordant DEGs, UP vs DOWN
########################################################################
cat("=== TASK 1: GO dotplot, shared concordant DEGs (R04 up/down) ===\n")
tryCatch({
  up_path   <- file.path(TABLES_DIR, "R04_overlap_up_genes.csv")
  down_path <- file.path(TABLES_DIR, "R04_overlap_down_genes.csv")
  if (check_input(up_path) && check_input(down_path)) {
    shared_up   <- read_csv(up_path,   show_col_types = FALSE)$gene_name
    shared_down <- read_csv(down_path, show_col_types = FALSE)$gene_name
    cat(sprintf("  Input: %s (%d genes), %s (%d genes)\n", up_path, length(shared_up), down_path, length(shared_down)))

    go_up   <- run_go_bp(shared_up,   "Shared UP")
    go_down <- run_go_bp(shared_down, "Shared DOWN")
    go_all  <- bind_rows(go_up, go_down)

    if (!is.null(go_all) && nrow(go_all) > 0) {
      cat(sprintf("  GO BP terms found: %d (Shared UP: %d, Shared DOWN: %d)\n",
                  nrow(go_all), sum(go_all$Group == "Shared UP"), sum(go_all$Group == "Shared DOWN")))
      p <- go_dotplot_pub(go_all,
             "GO Biological Process: shared concordant DEGs (both genotypes)",
             sprintf("Top 10 terms per direction | FDR < %.2f, |log2FC| > %.2f | background: all expressed genes", FDR, LFC_THR),
             facet_var = "Group")
      save_plot_both(p, "AR07_GO_shared_concordant_up_down", width = 13, height = 6)
    } else {
      cat("  [SKIP] no significant GO BP terms found for either direction.\n")
    }
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 1 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 2 -- NPB / anterior identity gene barplot
########################################################################
cat("=== TASK 2: NPB / anterior identity gene barplot ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  if (check_input(full_a_path) && check_input(full_b_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)
    cat(sprintf("  Input: %s, %s\n", full_a_path, full_b_path))

    gene_groups <- tribble(
      ~gene_name, ~gene_group,
      "MSX1", "NPB / dorsal", "MSX2", "NPB / dorsal",
      "DLX3", "Placodal", "DLX5", "Placodal",
      "SNAI2", "EMT driver",
      "DMBX1", "Anterior", "LHX2", "Anterior", "RAX", "Anterior", "SIX3", "Anterior",
      "HOXA2", "Posterior drift", "TLX3", "Posterior drift",
      "CYP26A1", "RA metabolism", "CYP26C1", "RA metabolism"
    ) %>% mutate(gene_name = factor(gene_name, levels = gene_name),
                 gene_group = factor(gene_group, levels = unique(gene_group)))

    extract <- function(df, label) {
      df %>% filter(gene_name %in% gene_groups$gene_name) %>%
        transmute(gene_name, genotype = label, log2FC = log2FoldChange, padj = padj)
    }
    plot_data <- bind_rows(extract(rna_a, LABEL_A), extract(rna_b, LABEL_B)) %>%
      right_join(gene_groups, by = "gene_name") %>%
      mutate(gene_name = factor(gene_name, levels = levels(gene_groups$gene_name)),
             significant = !is.na(padj) & padj < FDR & abs(log2FC) > LFC_THR)

    missing_genes <- setdiff(as.character(gene_groups$gene_name), plot_data$gene_name[!is.na(plot_data$log2FC)])
    cat(sprintf("  %d/%d requested genes found in both RNA tables\n", length(unique(gene_groups$gene_name)) - length(unique(missing_genes)), length(unique(gene_groups$gene_name))))
    if (length(missing_genes) > 0) cat(sprintf("  Not found in R02 tables (dropped, not fabricated): %s\n", paste(unique(missing_genes), collapse = ", ")))

    p <- ggplot(plot_data, aes(x = gene_name, y = log2FC, fill = genotype, alpha = significant)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      facet_wrap(~gene_group, scales = "free_x", nrow = 2) +
      scale_fill_manual(values = GENOTYPE_COLORS, name = "Genotype") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), name = sprintf("Significant\n(padj<%.2f, |log2FC|>%.2f)", FDR, LFC_THR),
                          labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
      labs(title = "NPB, placodal, EMT-driver, and anterior/posterior identity gene panel",
           subtitle = sprintf("RNA log2FC vs WT | %s vs %s", LABEL_A, LABEL_B),
           x = NULL, y = "RNA log2 fold change (vs WT)") +
      theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "bottom") +
      guides(fill = guide_legend(override.aes = list(alpha = 1)))
    ok <- save_plot_both(p, "AR07_patterning_genes_barplot", width = 13, height = 8)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 2 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 3 -- GO enrichment of the 103-gene ATAC-RNA concordant set (AR05)
########################################################################
cat("=== TASK 3: GO dotplot, 103-gene ATAC-RNA concordant set (AR05) ===\n")
tryCatch({
  ar05_path <- file.path(TABLES_DIR, "AR05_concordant_heatmap_data.csv")
  if (check_input(ar05_path)) {
    ar05_data <- read_csv(ar05_path, show_col_types = FALSE)
    gained_genes <- ar05_data %>% filter(Direction == "Up")   %>% pull(SYMBOL)
    lost_genes   <- ar05_data %>% filter(Direction == "Down") %>% pull(SYMBOL)
    cat(sprintf("  Input: %s -- gained (n=%d), lost (n=%d)\n", ar05_path, length(gained_genes), length(lost_genes)))

    go_gained <- run_go_bp(gained_genes, "Gained (concordant UP)")
    go_lost   <- run_go_bp(lost_genes,   "Lost (concordant DOWN)")
    go_all    <- bind_rows(go_gained, go_lost)

    if (!is.null(go_all) && nrow(go_all) > 0) {
      n_before <- nrow(go_all)
      cat(sprintf("  GO BP terms found: %d (Gained: %d, Lost: %d)\n",
                  n_before, sum(go_all$Group == "Gained (concordant UP)"), sum(go_all$Group == "Lost (concordant DOWN)")))
      go_all <- collapse_redundant_go_terms(go_all)
      cat(sprintf("  Redundancy collapse (gene-set Jaccard >= 0.7): %d -> %d terms\n", n_before, nrow(go_all)))
      p <- go_dotplot_pub(go_all,
             "GO Biological Process: 103-gene ATAC-RNA concordant target set",
             sprintf("Top 10 terms per direction, redundancy-collapsed (Jaccard >= 0.7) | FDR < %.2f | background: all expressed genes", FDR),
             facet_var = "Group")
      save_plot_both(p, "AR07_GO_103gene_concordant", width = 13, height = 6)
    } else {
      cat("  [SKIP] no significant GO BP terms found for either direction.\n")
    }
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 3 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 4 -- ATAC (chromatin accessibility) vs RNA scatter, anterior loci
########################################################################
cat("=== TASK 4: chromatin-accessibility vs RNA scatter, anterior identity loci ===\n")
tryCatch({
  # Top 30 Concordant-Up + top 30 Concordant-Down genes (same direction in
  # both assays, NOT Discordant/ATAC-only/RNA-only; split by sign of RNA_lfc
  # so both directions get representation instead of one direction crowding
  # out the other) by combined ATAC+RNA significance (ATAC_FDR * RNA_padj,
  # same composite-score convention as R04_scatter_convergence.pdf's
  # Sig_Score), not a fixed curated gene list.
  N_LABELS_ATAC_RNA_PER_DIRECTION <- 30

  make_scatter <- function(contrast_name, contrast_label) {
    rna_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", contrast_name))
    gained_path <- file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_GAINED.csv", contrast_name))
    lost_path   <- file.path(TABLES_DIR, sprintf("AR01_great_correlation_%s_LOST.csv",   contrast_name))
    if (!check_input(rna_path) || !(check_input(gained_path) || check_input(lost_path))) return(NULL)

    atac <- bind_rows(
      if (file.exists(gained_path)) read_csv(gained_path, show_col_types = FALSE) else NULL,
      if (file.exists(lost_path))   read_csv(lost_path,   show_col_types = FALSE) else NULL
    ) %>%
      group_by(SYMBOL) %>% slice_max(abs(ATAC_Fold), n = 1, with_ties = FALSE) %>% ungroup()

    rna <- read_csv(rna_path, show_col_types = FALSE)
    merged <- atac %>%
      dplyr::select(SYMBOL, ATAC_Fold, ATAC_FDR, concordance) %>%
      left_join(rna %>% transmute(SYMBOL = gene_name, RNA_lfc = log2FoldChange, RNA_padj = padj), by = "SYMBOL")

    cat(sprintf("  %s: %d genes with a GREAT-inferred target-gene assignment (input: %s, %s)\n",
                contrast_name, nrow(merged), gained_path, lost_path))

    concordant_scored <- merged %>% filter(concordance == "Concordant", !is.na(RNA_padj), !is.na(ATAC_FDR)) %>%
      mutate(sig_score = RNA_padj * ATAC_FDR)
    label_data <- bind_rows(
      concordant_scored %>% filter(RNA_lfc > 0) %>% arrange(sig_score) %>% head(N_LABELS_ATAC_RNA_PER_DIRECTION),
      concordant_scored %>% filter(RNA_lfc < 0) %>% arrange(sig_score) %>% head(N_LABELS_ATAC_RNA_PER_DIRECTION))
    cat(sprintf("  %s: labelling top %d up + %d down Concordant genes (of %d/%d requested) by combined ATAC+RNA significance\n",
                contrast_name, sum(label_data$RNA_lfc > 0), sum(label_data$RNA_lfc < 0),
                N_LABELS_ATAC_RNA_PER_DIRECTION, N_LABELS_ATAC_RNA_PER_DIRECTION))

    ggplot(merged, aes(x = ATAC_Fold, y = RNA_lfc, color = concordance)) +
      geom_point(alpha = 0.5, size = 1.3) +
      geom_point(data = label_data, shape = 21, size = 2.6, stroke = 1, fill = "white", color = "black") +
      geom_text_repel(data = label_data, aes(label = SYMBOL), color = "black", fontface = "italic",
                       size = 2.6, box.padding = 0.4, point.padding = 0.15, force = 3, force_pull = 0.5,
                       max.overlaps = Inf, min.segment.length = 0, max.time = 5, max.iter = 100000, seed = 42) +
      geom_vline(xintercept = 0, color = "grey50") + geom_hline(yintercept = 0, color = "grey50") +
      geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "grey60") +
      geom_hline(yintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "grey60") +
      # Fewer colours: only Concordant/Discordant (the categories that carry the
      # actual "do ATAC and RNA agree" signal, and the only ones ever labelled)
      # get an accent colour; ATAC-only/RNA-only/Not significant -- the
      # uninformative majority of points -- are all flattened to one grey so
      # they read as background, not as three more categories to parse.
      scale_color_manual(values = c(Concordant = "#B2182B", Discordant = "#2166AC",
                                     "ATAC-only" = "grey75", "RNA-only" = "grey75", "Not significant" = "grey75"),
                          name = "Category") +
      labs(title = contrast_label,
           x = "Chromatin accessibility log2FC (strongest GREAT-inferred DAR)",
           y = "RNA log2 fold change") +
      theme_pub + theme(legend.position = "bottom")
  }

  p_a <- make_scatter(NAME_A, sprintf("%s vs WT", LABEL_A))
  p_b <- make_scatter(NAME_B, sprintf("%s vs WT", LABEL_B))

  if (!is.null(p_a) && !is.null(p_b)) {
    combined <- p_a + p_b
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      cat("  patchwork not installed -- saving the two contrast panels as separate files instead of one combined figure.\n")
      save_plot_both(p_a, "AR07_ATAC_RNA_scatter_anterior_loci_KDM6A_ko", width = 8, height = 7)
      save_plot_both(p_b, "AR07_ATAC_RNA_scatter_anterior_loci_KMT2D_Het", width = 8, height = 7)
    } else {
      library(patchwork)
      p_final <- (p_a + p_b) + patchwork::plot_annotation(
        title = "Chromatin accessibility vs RNA expression change, anterior identity loci",
        subtitle = sprintf("GREAT-inferred target genes | quadrant lines at 0 | dashed lines at |log2FC| = %.2f", LFC_THR))
      save_plot_both(p_final, "AR07_ATAC_RNA_scatter_anterior_loci", width = 15, height = 7.5)
    }
  } else {
    cat("  [SKIP] one or both contrasts had no usable input; figure not generated.\n")
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 4 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 5 -- BMP / WNT signalling gene barplot (Figure 4)
########################################################################
cat("=== TASK 5: BMP / WNT signalling gene barplot ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  if (check_input(full_a_path) && check_input(full_b_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)
    cat(sprintf("  Input: %s, %s\n", full_a_path, full_b_path))

    # WNT panel audited against the full R02 tables (2026-09-11): every WNT
    # ligand significant (padj<0.05) in >=1 genotype is now included -- the
    # previous list (WNT3A/WNT7A/WNT8A only) missed WNT1/WNT2B/WNT3/WNT4/
    # WNT5B/WNT8B/WNT9A. CTNNB1 dropped: beta-catenin activity is set
    # post-translationally by the destruction complex, so its own mRNA level
    # (log2FC ~0.14, padj 0.36-0.41 here, flat) is not a meaningful pathway
    # readout. TCF7 dropped (also flat/ns, redundant with LEF1 as the TF
    # readout). Added significant transcriptional targets (SP5, CD44, RNF43)
    # as the RNA-level pathway-activity readout, alongside AXIN2 (kept as the
    # standard reporter for reference even though non-significant here) and
    # LEF1.
    gene_groups <- tribble(
      ~gene_name, ~gene_group,
      "BMP4", "BMP signalling", "BMP7", "BMP signalling", "GDF7", "BMP signalling", "GDF10", "BMP signalling",
      "BMPR1B", "BMP signalling", "SMAD9", "BMP signalling",
      "ID1", "BMP signalling", "ID2", "BMP signalling", "ID3", "BMP signalling", "ID4", "BMP signalling",
      "NOG", "BMP signalling", "BAMBI", "BMP signalling",
      "WNT1", "WNT ligands", "WNT2B", "WNT ligands", "WNT3", "WNT ligands", "WNT3A", "WNT ligands",
      "WNT4", "WNT ligands", "WNT5B", "WNT ligands", "WNT7A", "WNT ligands", "WNT8A", "WNT ligands",
      "WNT8B", "WNT ligands", "WNT9A", "WNT ligands",
      "AXIN2", "WNT target genes", "SP5", "WNT target genes", "CD44", "WNT target genes",
      "RNF43", "WNT target genes", "DKK1", "WNT target genes", "LEF1", "WNT target genes"
    ) %>% mutate(gene_name = factor(gene_name, levels = gene_name),
                 gene_group = factor(gene_group, levels = unique(gene_group)))

    extract <- function(df, label) {
      df %>% filter(gene_name %in% gene_groups$gene_name) %>%
        transmute(gene_name, genotype = label, log2FC = log2FoldChange, padj = padj)
    }
    plot_data <- bind_rows(extract(rna_a, LABEL_A), extract(rna_b, LABEL_B)) %>%
      right_join(gene_groups, by = "gene_name") %>%
      mutate(gene_name = factor(gene_name, levels = levels(gene_groups$gene_name)),
             genotype = factor(genotype, levels = c(LABEL_B, LABEL_A)),
             significant = !is.na(padj) & padj < FDR & abs(log2FC) > LFC_THR)

    missing_genes <- setdiff(as.character(gene_groups$gene_name), plot_data$gene_name[!is.na(plot_data$log2FC)])
    cat(sprintf("  %d/%d requested genes found in both RNA tables\n", length(unique(gene_groups$gene_name)) - length(unique(missing_genes)), length(unique(gene_groups$gene_name))))
    if (length(missing_genes) > 0) cat(sprintf("  Not found in R02 tables (dropped, not fabricated): %s\n", paste(unique(missing_genes), collapse = ", ")))

    p <- ggplot(plot_data, aes(x = gene_name, y = log2FC, fill = genotype, alpha = significant)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
      geom_hline(yintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "grey40", linewidth = 0.4) +
      facet_grid(~gene_group, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = GENOTYPE_COLORS, name = "Genotype", breaks = c(LABEL_A, LABEL_B)) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), name = NULL,
                          labels = c(`TRUE` = sprintf("padj < %.2f and |log2FC| > %.2f", FDR, LFC_THR),
                                     `FALSE` = "not significant / below threshold")) +
      labs(title = "BMP and WNT signalling gene panel",
           x = NULL, y = "log2 fold change (vs WT)") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            legend.position = "bottom", legend.box = "vertical",
            panel.spacing = unit(0.4, "lines")) +
      guides(fill = guide_legend(override.aes = list(alpha = 1), order = 1),
             alpha = guide_legend(override.aes = list(fill = "grey40"), order = 2))
    ok <- save_plot_both(p, "AR07_BMP_WNT_signalling_barplot", width = 13, height = 6)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 5 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 6 -- shared upregulated / mesenchymal program barplot (Figure 3 companion)
########################################################################
cat("=== TASK 6: shared upregulated mesenchymal-program barplot ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  if (check_input(full_a_path) && check_input(full_b_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)
    cat(sprintf("  Input: %s, %s\n", full_a_path, full_b_path))

    gene_order <- c("SNAI2", "TWIST1", "FN1", "COL6A3", "TGFB2", "COL3A1",
                     "NCAM2", "PIEZO2", "MEF2C", "ALX4", "NTRK2")

    extract <- function(df, label) {
      df %>% filter(gene_name %in% gene_order) %>%
        transmute(gene_name, genotype = label, log2FC = log2FoldChange, padj = padj)
    }
    plot_data <- bind_rows(extract(rna_a, LABEL_A), extract(rna_b, LABEL_B)) %>%
      mutate(gene_name = factor(gene_name, levels = gene_order),
             significant = !is.na(padj) & padj < FDR & abs(log2FC) > LFC_THR)

    missing_genes <- setdiff(gene_order, plot_data$gene_name[!is.na(plot_data$log2FC)])
    cat(sprintf("  %d/%d requested genes found in both RNA tables\n", length(gene_order) - length(unique(missing_genes)), length(gene_order)))
    if (length(missing_genes) > 0) cat(sprintf("  Not found in R02 tables (dropped, not fabricated): %s\n", paste(unique(missing_genes), collapse = ", ")))

    p <- ggplot(plot_data, aes(x = gene_name, y = log2FC, fill = genotype, alpha = significant)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = GENOTYPE_COLORS, name = "Genotype") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), name = sprintf("Significant\n(padj<%.2f, |log2FC|>%.2f)", FDR, LFC_THR),
                          labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
      labs(title = "Shared upregulated - mesenchymal program",
           subtitle = sprintf("RNA log2FC vs WT | %s vs %s", LABEL_A, LABEL_B),
           x = NULL, y = "RNA log2 fold change (vs WT)") +
      theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "bottom") +
      guides(fill = guide_legend(override.aes = list(alpha = 1)))
    ok <- save_plot_both(p, "AR07_EMT_mesenchymal_barplot", width = 11, height = 6)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 6 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 7 -- Directional concordance summary: shared vs genotype-specific DEGs
########################################################################
cat("=== TASK 7: directional concordance summary barplot ===\n")
tryCatch({
  counts_path <- file.path(TABLES_DIR, "R04_directional_concordance_counts.csv")
  effect_path <- file.path(TABLES_DIR, "R04_overlap_effect_size.csv")
  if (check_input(counts_path) && check_input(effect_path)) {
    counts <- read_csv(counts_path, show_col_types = FALSE)
    effect <- read_csv(effect_path, show_col_types = FALSE)
    cat(sprintf("  Input: %s, %s\n", counts_path, effect_path))

    cat_levels <- c("Shared", "Discordant", paste(LABEL_A, "only"), paste(LABEL_B, "only"))
    plot_data <- counts %>%
      mutate(category = factor(category, levels = cat_levels),
             direction = factor(coalesce(direction, "n/a"), levels = c("Up", "Down", "n/a")))
    cat_totals <- plot_data %>% group_by(category) %>% summarise(total = sum(n), .groups = "drop")

    annotation <- sprintf("OR = %.1f (95%% CI: %.1f-%.1f)\npermutation p %s",
                           effect$OR, effect$OR_CI_low, effect$OR_CI_high,
                           ifelse(effect$perm_p < 1e-4, "< 1e-4", sprintf("= %.4f", effect$perm_p)))

    p <- ggplot(plot_data, aes(x = category, y = n, fill = direction)) +
      geom_col(color = "black", linewidth = 0.2, width = 0.7) +
      geom_text(data = cat_totals, aes(x = category, y = total, label = total), inherit.aes = FALSE,
                vjust = -0.4, fontface = "bold", size = 5) +
      annotate("label", x = 1, y = max(cat_totals$total) * 1.15, label = annotation,
               hjust = 0, size = 5.5, fill = "white") +
      scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC", `n/a` = "grey60"), name = "Direction") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title = "Directional concordance summary",
           subtitle = sprintf("Shared vs genotype-specific DEGs | %s vs %s | padj < %.2f, |log2FC| > %.2f",
                               LABEL_A, LABEL_B, FDR, LFC_THR),
           x = NULL, y = "Number of DEGs") +
      theme_pub + theme(legend.position = "right",
                         axis.text.y = element_text(size = 13, color = "black"),
                         axis.text.x = element_text(size = 13, color = "black", angle = 20, hjust = 1),
                         axis.title.y = element_text(size = 15, face = "bold"),
                         plot.title = element_text(size = 17, face = "bold"),
                         plot.subtitle = element_text(size = 12),
                         legend.text = element_text(size = 12), legend.title = element_text(size = 13, face = "bold"),
                         panel.border = element_blank(), axis.line = element_line(color = "black"))
    ok <- save_plot_both(p, "AR07_directional_concordance_summary_barplot", width = 9, height = 6.5)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 7 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 8 -- Combined top-DEG heatmap: top genes from BOTH genotypes
## (R02's per-contrast top50 heatmaps each rank by that contrast's own padj
## alone, so e.g. the KMT2D_Het_vs_WT one only shows genes that are top DEGs
## for KMT2D_Het -- a KDM6A_ko-specific top gene that isn't also strong in
## KMT2D_Het won't appear there. This unions the top N of each contrast and
## marks, per gene, which contrast(s) it was a top hit in.)
########################################################################
cat("=== TASK 8: combined top-DEG heatmap (KDM6A_ko + KMT2D_Het) ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  vst_path    <- "results/RDS/R02_vst_counts.rds"
  if (check_input(full_a_path) && check_input(full_b_path) && check_input(vst_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)
    vsd   <- readRDS(vst_path)
    cat(sprintf("  Input: %s, %s, %s\n", full_a_path, full_b_path, vst_path))

    N_TOP <- 50
    top_ids <- function(res_df) {
      res_df %>% filter(!is.na(padj), padj < FDR, abs(log2FoldChange) > LFC_THR) %>%
        arrange(padj) %>% head(N_TOP) %>% pull(gene_id)
    }
    top_a <- top_ids(rna_a); top_b <- top_ids(rna_b)
    combined_ids <- union(top_a, top_b)
    cat(sprintf("  Top %d DEGs per contrast: %s only=%d, %s only=%d, both=%d, combined=%d\n",
                N_TOP, LABEL_A, length(setdiff(top_a, top_b)), LABEL_B, length(setdiff(top_b, top_a)),
                length(intersect(top_a, top_b)), length(combined_ids)))

    id_to_name <- bind_rows(rna_a, rna_b) %>% distinct(gene_id, gene_name) %>% deframe()
    source_lab <- case_when(combined_ids %in% top_a & combined_ids %in% top_b ~ "Both",
                             combined_ids %in% top_a ~ paste(LABEL_A, "top"),
                             TRUE ~ paste(LABEL_B, "top"))

    mat <- assay(vsd)[combined_ids, , drop = FALSE]
    rownames(mat) <- id_to_name[rownames(mat)]
    mat_scaled <- t(scale(t(mat)))

    genotype_order <- c("WT", LABEL_A, LABEL_B)
    col_order <- order(match(colData(vsd)$genotype, genotype_order))
    mat_scaled <- mat_scaled[, col_order, drop = FALSE]

    row_order <- order(factor(source_lab, levels = c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top"))))
    mat_scaled <- mat_scaled[row_order, , drop = FALSE]
    row_ann <- data.frame(`Top in` = source_lab[row_order], check.names = FALSE, row.names = rownames(mat_scaled))
    col_ann <- data.frame(Genotype = colData(vsd)$genotype[col_order], row.names = colnames(mat_scaled))
    ann_colors <- list(Genotype = setNames(c("grey50", GENOTYPE_COLORS[LABEL_A], GENOTYPE_COLORS[LABEL_B]), genotype_order),
                        `Top in` = setNames(c("#56B4E9", "#984EA3", "#E69F00"),
                                             c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top"))))
    gaps_row <- cumsum(table(factor(source_lab, levels = c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top")))))
    gaps_row <- gaps_row[gaps_row > 0 & gaps_row < nrow(mat_scaled)]

    out_pdf <- file.path(PLOTS_DIR, "AR07_combined_top_DEG_heatmap.pdf")
    pdf(out_pdf, width = 8, height = max(10, nrow(mat_scaled) * 0.16))
    pheatmap(mat_scaled, annotation_row = row_ann, annotation_col = col_ann, annotation_colors = ann_colors,
             color = colorRampPalette(c("#4393c3", "white", "#d6604d"))(100), breaks = seq(-2, 2, length.out = 100),
             cluster_rows = FALSE, cluster_cols = TRUE, gaps_row = gaps_row, border_color = "grey",
             show_rownames = TRUE, fontsize_row = 7,
             main = sprintf("Top %d DEGs per genotype (union) - %s vs %s (padj < %.2f, |log2FC| > %.2f)",
                             N_TOP, LABEL_A, LABEL_B, FDR, LFC_THR))
    dev.off()
    cat(sprintf("  [OK] wrote AR07_combined_top_DEG_heatmap (.pdf, %d genes)\n", nrow(mat_scaled)))
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 8 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 9 -- Combined panel: directional concordance summary (TASK 7) +
## combined top-DEG heatmap (TASK 8) as one figure. Self-contained (rebuilds
## both plot objects from the same inputs, doesn't reuse TASK 7/8's local
## variables) -- same independence convention as go_dotplot_pub above.
########################################################################
cat("=== TASK 9: combined directional-concordance panel (barplot + heatmap) ===\n")
tryCatch({
  counts_path <- file.path(TABLES_DIR, "R04_directional_concordance_counts.csv")
  effect_path <- file.path(TABLES_DIR, "R04_overlap_effect_size.csv")
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  vst_path    <- "results/RDS/R02_vst_counts.rds"
  have_inputs <- check_input(counts_path) && check_input(effect_path) &&
    check_input(full_a_path) && check_input(full_b_path) && check_input(vst_path)

  if (have_inputs && !requireNamespace("patchwork", quietly = TRUE)) {
    cat("  [SKIP] patchwork not installed -- combined panel needs it (TASK 7/8 outputs remain available as separate files).\n")
  } else if (have_inputs) {
    library(patchwork)

    # -- panel A: directional concordance barplot (same as TASK 7) --
    counts <- read_csv(counts_path, show_col_types = FALSE)
    effect <- read_csv(effect_path, show_col_types = FALSE)
    cat_levels <- c("Shared", "Discordant", paste(LABEL_A, "only"), paste(LABEL_B, "only"))
    bar_data <- counts %>%
      mutate(category = factor(category, levels = cat_levels),
             direction = factor(coalesce(direction, "n/a"), levels = c("Up", "Down", "n/a")))
    bar_totals <- bar_data %>% group_by(category) %>% summarise(total = sum(n), .groups = "drop")
    annotation <- sprintf("OR = %.1f (95%% CI: %.1f-%.1f)\npermutation p %s",
                           effect$OR, effect$OR_CI_low, effect$OR_CI_high,
                           ifelse(effect$perm_p < 1e-4, "< 1e-4", sprintf("= %.4f", effect$perm_p)))
    p_bar <- ggplot(bar_data, aes(x = category, y = n, fill = direction)) +
      geom_col(color = "black", linewidth = 0.2, width = 0.7) +
      geom_text(data = bar_totals, aes(x = category, y = total, label = total), inherit.aes = FALSE,
                vjust = -0.4, fontface = "bold", size = 5) +
      annotate("label", x = 1, y = max(bar_totals$total) * 1.15, label = annotation,
               hjust = 0, size = 5.5, fill = "white") +
      scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC", `n/a` = "grey60"), name = "Direction") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title = "Directional concordance summary",
           subtitle = sprintf("Shared vs genotype-specific DEGs | padj < %.2f, |log2FC| > %.2f", FDR, LFC_THR),
           x = NULL, y = "Number of DEGs") +
      theme_pub + theme(legend.position = "right",
                         axis.text.y = element_text(size = 13, color = "black"),
                         axis.text.x = element_text(size = 13, color = "black", angle = 20, hjust = 1),
                         axis.title.y = element_text(size = 15, face = "bold"),
                         plot.title = element_text(size = 17, face = "bold"),
                         plot.subtitle = element_text(size = 12),
                         legend.text = element_text(size = 12), legend.title = element_text(size = 13, face = "bold"),
                         panel.border = element_blank(), axis.line = element_line(color = "black"))

    # -- panel B: combined top-DEG heatmap (same as TASK 8), silent = TRUE so
    # pheatmap returns a grob instead of drawing straight to a device --
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)
    vsd   <- readRDS(vst_path)
    N_TOP <- 50
    top_ids <- function(res_df) {
      res_df %>% filter(!is.na(padj), padj < FDR, abs(log2FoldChange) > LFC_THR) %>%
        arrange(padj) %>% head(N_TOP) %>% pull(gene_id)
    }
    top_a <- top_ids(rna_a); top_b <- top_ids(rna_b)
    combined_ids <- union(top_a, top_b)
    id_to_name <- bind_rows(rna_a, rna_b) %>% distinct(gene_id, gene_name) %>% deframe()
    source_lab <- case_when(combined_ids %in% top_a & combined_ids %in% top_b ~ "Both",
                             combined_ids %in% top_a ~ paste(LABEL_A, "top"),
                             TRUE ~ paste(LABEL_B, "top"))

    mat <- assay(vsd)[combined_ids, , drop = FALSE]
    rownames(mat) <- id_to_name[rownames(mat)]
    mat_scaled <- t(scale(t(mat)))

    genotype_order <- c("WT", LABEL_A, LABEL_B)
    col_order <- order(match(colData(vsd)$genotype, genotype_order))
    mat_scaled <- mat_scaled[, col_order, drop = FALSE]

    row_order <- order(factor(source_lab, levels = c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top"))))
    mat_scaled <- mat_scaled[row_order, , drop = FALSE]
    row_ann <- data.frame(`Top in` = source_lab[row_order], check.names = FALSE, row.names = rownames(mat_scaled))
    col_ann <- data.frame(Genotype = colData(vsd)$genotype[col_order], row.names = colnames(mat_scaled))
    ann_colors <- list(Genotype = setNames(c("grey50", GENOTYPE_COLORS[LABEL_A], GENOTYPE_COLORS[LABEL_B]), genotype_order),
                        `Top in` = setNames(c("#56B4E9", "#984EA3", "#E69F00"),
                                             c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top"))))
    gaps_row <- cumsum(table(factor(source_lab, levels = c(paste(LABEL_A, "top"), "Both", paste(LABEL_B, "top")))))
    gaps_row <- gaps_row[gaps_row > 0 & gaps_row < nrow(mat_scaled)]

    ht <- pheatmap(mat_scaled, annotation_row = row_ann, annotation_col = col_ann, annotation_colors = ann_colors,
                    color = colorRampPalette(c("#4393c3", "white", "#d6604d"))(100), breaks = seq(-2, 2, length.out = 100),
                    cluster_rows = FALSE, cluster_cols = TRUE, gaps_row = gaps_row, border_color = "grey",
                    show_rownames = TRUE, fontsize_row = 7,
                    main = sprintf("Top %d DEGs per genotype (union)", N_TOP), silent = TRUE)

    # -- combine: barplot on top, heatmap below (very different aspect ratios --
    # stacking keeps each panel's own text/legend legible instead of squashing
    # the tall heatmap into a side-by-side column) --
    p_final <- p_bar / patchwork::wrap_elements(full = ht$gtable) + patchwork::plot_layout(heights = c(1, 2.6))
    height_total <- 6.5 + max(10, nrow(mat_scaled) * 0.16)
    ok <- save_plot_both(p_final, "AR07_directional_concordance_panel", width = 10, height = height_total)
  } else {
    cat("  [SKIP] one or more inputs missing; combined panel not generated.\n")
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 9 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 10 -- RNA-seq PCA, paper-panel version: same DESeq2 VST object as
## R02_pca.pdf, but with per-sample point labels dropped (color + shape only)
## so it drops straight into a multi-panel figure. See R02_pca.pdf for the
## labelled per-sample QC version. Smaller panel (5x4, was 7x5) with larger
## text (own theme override, not AR07's shared theme_pub) so it still reads
## clearly once shrunk into a multi-panel layout.
########################################################################
cat("=== TASK 10: RNA-seq PCA (no point labels) ===\n")
tryCatch({
  vst_path <- file.path("results/RDS", "R02_vst_counts.rds")
  if (check_input(vst_path)) {
    vsd <- readRDS(vst_path)
    pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))
    pca_shapes <- setNames(c(16, 17, 15), c("WT", LABEL_A, LABEL_B))

    p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype, shape = genotype)) +
      geom_point(size = 4) +
      scale_color_manual(values = c("WT" = "grey50", GENOTYPE_COLORS), name = "Genotype") +
      scale_shape_manual(values = pca_shapes, name = "Genotype") +
      xlab(paste0("PC1: ", percent_var[1], "% variance")) +
      ylab(paste0("PC2: ", percent_var[2], "% variance")) +
      theme_pub + labs(title = "PCA - RNA-seq samples") +
      theme(axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            plot.title = element_text(size = 15),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 12))
    save_plot_both(p_pca, "AR07_pca_no_labels", width = 5, height = 4)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 10 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 11 -- GO BP dotplot, ORA (enrichGO, R02) -- single shared y-axis
## across both genotypes, NOT curated beyond redundancy collapsing (unlike
## Fig3_panels.R's GSEA version). Built from enrichGO's over-representation
## results (R02_GO_enrichment_per_genotype.csv, separate Up/Down gene lists
## per genotype).
##
## Pipeline, in order (order matters -- see below):
##  1. Per (genotype, direction), take the top N_PER_DIRECTION terms by FDR,
##     independently for all 4 groups -- NOT a single global top-N pooled
##     across all 4 (that crowds out KMT2D_Het/Down: only 11 terms total,
##     all far weaker than KDM6A_ko's much larger/stronger sets). A term
##     significant in BOTH directions for one genotype is NOT collapsed to
##     its stronger direction -- both points are plotted (verified this
##     happens: e.g. a term can be driven by distinct up- and down-regulated
##     gene subsets, which is real biology, not noise).
##  2. Cutoff-artifact fix: a term is force-included into a genotype's set if
##     it made the OTHER genotype's top-N cut, even if it ranked below this
##     genotype's own top-N (verified this matters in practice: "skin
##     development" ranks 77th of KDM6A_ko's 109 significant Down terms --
##     still FDR < 0.05, just weaker than KDM6A_ko's own top hits -- but IS
##     top-15 in KMT2D_Het/Down, so without this step it would silently
##     read as "KMT2D_Het-only" when it's actually shared). Bounded to
##     candidates from the other genotype's own top-N list, not "every
##     significant term ever", so this stays a targeted fix.
##  3. Compute cross-genotype term overlap (exact Description match) BEFORE
##     any redundancy collapsing -- this order matters: collapsing first and
##     comparing after risks silently breaking a real overlap (see step 4).
##  4. Collapse near-duplicate GO terms within each genotype's own list
##     (gene-set Jaccard >= JACCARD_THRESH, e.g. "forebrain neuron
##     differentiation" / "forebrain generation of neurons" are the same 10
##     genes under two labels). Within a cluster, prefer keeping whichever
##     member matches the OTHER genotype's list (from step 3) over just
##     picking lowest FDR -- verified necessary in practice: KMT2D_Het's
##     "transforming growth factor beta receptor superfamily signaling
##     pathway" / "cell surface receptor protein serine/threonine kinase
##     signaling pathway" cluster is an exact FDR tie (6.25e-06 both), so a
##     blind lowest-FDR pick could have arbitrarily dropped the very term
##     that overlaps with KDM6A_ko, undercounting the real shared biology.
##     (Up-row and Down-row gene sets for the same genotype are always
##     disjoint by construction, so this never accidentally merges a term's
##     own Up and Down points into one.)
##  5. Union the two genotypes' final (collapsed) term lists into ONE y-axis,
##     grouped into three blocks -- Shared / KDM6A_ko-only / KMT2D_Het-only
##     -- rather than one undifferentiated sorted list, so the shared
##     TGF-beta/mesenchyme/skeletal program the two genotypes have in common
##     is immediately visually separated from what's genotype-specific.
##
## x axis: signed log2(enrichment ratio) -- log2 (not raw ratio) because raw
## ratio is bounded below by ~1, so a raw signed version would leave a dead
## gap between -1 and +1; log2 puts "no enrichment" exactly at 0 and is
## symmetric in both directions. Positive = enriched among UP genes,
## negative = enriched among DOWN genes. Two columns (one per genotype,
## facet_wrap), gene count written in white inside each dot, colour =
## -log10(FDR) (same convention as go_dotplot_pub() above and
## AR02_GO_Plots.R::go_dotplot()).
########################################################################
cat("=== TASK 11: GO BP dotplot (enrichGO ORA, single shared y-axis) ===\n")
tryCatch({
  go_ora_path <- file.path(TABLES_DIR, "R02_GO_enrichment_per_genotype.csv")
  if (check_input(go_ora_path)) {
    N_PER_DIRECTION <- 15

    parse_ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
    go_ora <- read_csv(go_ora_path, show_col_types = FALSE) %>%
      mutate(enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
             genes = str_split(geneID, "/"))

    top_by_direction <- function(name) {
      go_ora %>% filter(Comparison == name) %>% group_by(Direction) %>% arrange(p.adjust) %>%
        slice_head(n = N_PER_DIRECTION) %>% ungroup()
    }
    a_top <- top_by_direction(NAME_A)
    b_top <- top_by_direction(NAME_B)

    force_include <- function(name, other_terms) go_ora %>% filter(Comparison == name, Description %in% other_terms)
    a_raw <- bind_rows(a_top, force_include(NAME_A, unique(b_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
    b_raw <- bind_rows(b_top, force_include(NAME_B, unique(a_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
    cat(sprintf("  %s: %d top-%d terms + %d force-included from cross-genotype overlap = %d total\n",
                LABEL_A, nrow(a_top), N_PER_DIRECTION, nrow(a_raw) - nrow(a_top), nrow(a_raw)))
    cat(sprintf("  %s: %d top-%d terms + %d force-included from cross-genotype overlap = %d total\n",
                LABEL_B, nrow(b_top), N_PER_DIRECTION, nrow(b_raw) - nrow(b_top), nrow(b_raw)))

    common_terms <- intersect(unique(a_raw$Description), unique(b_raw$Description))
    cat(sprintf("  Cross-genotype overlap: %d terms shared\n", length(common_terms)))

    # NOTE: this base version does NOT collapse redundant near-duplicate GO
    # terms (e.g. "forebrain neuron differentiation" / "forebrain generation
    # of neurons" both appear as separate rows) -- that collapsing is TASK 12
    # (AR07_GO_BP_ORA_dotplot_shared_axis_collapsed.pdf), a deliberately
    # separate output so the two can be compared side by side rather than
    # one silently overwriting the other.
    a_final <- a_raw
    b_final <- b_raw
    common_final <- common_terms

    combined <- bind_rows(a_final, b_final) %>%
      mutate(signed_score = ifelse(Direction == "Down", -log2(enrichment_ratio), log2(enrichment_ratio)),
             Comparison = factor(Comparison, levels = c(NAME_A, NAME_B),
                                  labels = c(genotype_strip_label(LABEL_A), genotype_strip_label(LABEL_B))))

    # single shared y-axis, grouped into 3 blocks so the overlap is visually
    # obvious rather than just implied by matching labels down a sorted list
    term_category <- tibble(Description = union(a_final$Description, b_final$Description)) %>%
      mutate(category = case_when(Description %in% common_final ~ "Shared",
                                   Description %in% a_final$Description ~ paste0(LABEL_A, " only"),
                                   TRUE ~ paste0(LABEL_B, " only")))
    term_order <- combined %>% left_join(term_category, by = "Description") %>%
      group_by(Description, category) %>% summarise(sort_key = mean(signed_score), .groups = "drop") %>%
      mutate(category = factor(category, levels = c("Shared", paste0(LABEL_A, " only"), paste0(LABEL_B, " only")))) %>%
      arrange(category, desc(sort_key)) %>% pull(Description)
    term_label <- setNames(str_wrap(str_to_sentence(term_order), 45), term_order)

    plot_data <- combined %>% mutate(label = factor(term_label[Description], levels = rev(term_label[term_order])))

    p_ora <- ggplot(plot_data, aes(x = signed_score, y = label)) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_point(aes(colour = -log10(p.adjust)), size = 6, alpha = 0.9) +
      geom_text(aes(label = Count), colour = "white", size = 2.4, fontface = "bold") +
      scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
                              name = expression(-log[10](FDR)),
                              labels = scales::label_number(accuracy = 0.1)) +
      facet_wrap(~Comparison) +
      labs(title = "GO Biological Process (ORA, enrichGO)",
           subtitle = sprintf("Top %d terms/direction/genotype, NOT redundancy-collapsed (see the _collapsed version) | FDR < %.2f\n+ = enriched in UP genes, - = enriched in DOWN genes | number in dot = gene count | rows grouped: Shared / genotype-specific",
                               N_PER_DIRECTION, FDR),
           x = "Signed log2(enrichment ratio)", y = NULL) +
      theme_pub + theme(panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.3),
                          strip.text = ggtext::element_markdown(face = "bold", size = 9))
    save_plot_both(p_ora, "AR07_GO_BP_ORA_dotplot_shared_axis", width = 11, height = 13)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 11 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 12 -- same single-shared-y-axis GO BP dotplot as TASK 11, but WITH
## redundant near-duplicate GO terms collapsed to one canonical label. Kept
## as a fully separate output (deliberately does not reuse/overwrite TASK
## 11's file) so the two can be compared directly.
##
## Redundancy collapsing is done as ONE GLOBAL pass across both genotypes
## combined, NOT two independent per-genotype passes -- verified this
## matters: two genotypes can each have their own internally-consistent
## "lowest FDR" pick within a redundant pair, but DISAGREE with each other
## (e.g. KDM6A_ko's own numbers favour "osteoblast differentiation" while
## KMT2D_Het's own numbers favour "regulation of osteoblast
## differentiation" for the same redundant pair) -- collapsing independently
## per genotype would silently rename the same GO concept to two different
## labels, which then fails to match as "Shared" even though it clearly is.
## A single global collapse forces both genotypes onto the SAME canonical
## label (lowest FDR anywhere it occurs, across both genotypes/directions)
## for a given redundant cluster.
##
## Factored into a function of N_PER_DIRECTION (top N terms per direction
## per genotype) and an output basename/height, since TASK 12/14/15 below are
## the same pipeline at N=15/8/5 -- three separate output files (a "look at
## everything" version and two smaller, more presentable cuts), not three
## copies of this logic.
########################################################################
build_collapsed_shared_axis_plot <- function(n_per_direction, output_name, plot_height, plot_width = 11, short_subtitle = FALSE) {
  go_ora_path <- file.path(TABLES_DIR, "R02_GO_enrichment_per_genotype.csv")
  if (!check_input(go_ora_path)) return(invisible(NULL))
  JACCARD_THRESH <- 0.7

  parse_ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
  go_ora <- read_csv(go_ora_path, show_col_types = FALSE) %>%
    mutate(enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
           genes = str_split(geneID, "/"))

  top_by_direction <- function(name) {
    go_ora %>% filter(Comparison == name) %>% group_by(Direction) %>% arrange(p.adjust) %>%
      slice_head(n = n_per_direction) %>% ungroup()
  }
  a_top <- top_by_direction(NAME_A)
  b_top <- top_by_direction(NAME_B)

  force_include <- function(name, other_terms) go_ora %>% filter(Comparison == name, Description %in% other_terms)
  a_raw <- bind_rows(a_top, force_include(NAME_A, unique(b_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
  b_raw <- bind_rows(b_top, force_include(NAME_B, unique(a_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)

  jaccard <- function(g1, g2) length(intersect(g1, g2)) / length(union(g1, g2))
  find_redundant_edges <- function(data) {
    n <- nrow(data)
    edges <- list()
    for (i in seq_len(n - 1)) for (j in seq((i + 1), n)) {
      if (jaccard(data$genes[[i]], data$genes[[j]]) >= JACCARD_THRESH) edges[[length(edges) + 1]] <- c(data$Description[i], data$Description[j])
    }
    edges
  }
  # gene-set Jaccard is only meaningful WITHIN a (genotype, direction) group
  # (different runs have different gene universes) -- edges are found per
  # genotype's own raw set, then pooled into one global graph over term
  # labels: if two labels are near-duplicate in EITHER genotype's own
  # analysis, they're the same GO concept everywhere.
  all_edges <- c(find_redundant_edges(a_raw), find_redundant_edges(b_raw))
  combined_raw <- bind_rows(a_raw, b_raw)
  all_terms <- unique(combined_raw$Description)
  redundancy_graph <- igraph::make_empty_graph(n = length(all_terms), directed = FALSE)
  igraph::V(redundancy_graph)$name <- all_terms
  if (length(all_edges) > 0) {
    edge_vec <- unlist(lapply(all_edges, function(e) match(e, all_terms)))
    redundancy_graph <- igraph::add_edges(redundancy_graph, edge_vec)
  }
  cluster_of <- setNames(igraph::components(redundancy_graph)$membership, all_terms)

  canonical_label <- combined_raw %>% group_by(Description) %>% summarise(min_padj = min(p.adjust), .groups = "drop") %>%
    mutate(cluster_id = cluster_of[Description]) %>%
    group_by(cluster_id) %>% mutate(canonical = Description[which.min(min_padj)]) %>% ungroup() %>%
    dplyr::select(Description, canonical)
  n_renamed <- sum(canonical_label$Description != canonical_label$canonical)
  cat(sprintf("  Global redundancy collapse (Jaccard >= %.1f): %d/%d terms merged into a shared canonical label\n",
              JACCARD_THRESH, n_renamed, length(all_terms)))

  combined <- combined_raw %>%
    mutate(Description = canonical_label$canonical[match(Description, canonical_label$Description)]) %>%
    group_by(Comparison, Direction, Description) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(signed_score = ifelse(Direction == "Down", -log2(enrichment_ratio), log2(enrichment_ratio)),
           Comparison = factor(Comparison, levels = c(NAME_A, NAME_B), labels = c(LABEL_A, LABEL_B)))

  a_final <- combined %>% filter(Comparison == LABEL_A)
  b_final <- combined %>% filter(Comparison == LABEL_B)
  common_final <- intersect(a_final$Description, b_final$Description)
  cat(sprintf("  After collapsing: %s %d terms | %s %d terms | %d shared\n",
              LABEL_A, n_distinct(a_final$Description), LABEL_B, n_distinct(b_final$Description), length(common_final)))

  term_category <- tibble(Description = union(a_final$Description, b_final$Description)) %>%
    mutate(category = case_when(Description %in% common_final ~ "Shared",
                                 Description %in% a_final$Description ~ paste0(LABEL_A, " only"),
                                 TRUE ~ paste0(LABEL_B, " only")))
  term_order <- combined %>% left_join(term_category, by = "Description") %>%
    group_by(Description, category) %>% summarise(sort_key = mean(signed_score), .groups = "drop") %>%
    mutate(category = factor(category, levels = c("Shared", paste0(LABEL_A, " only"), paste0(LABEL_B, " only")))) %>%
    arrange(category, desc(sort_key)) %>% pull(Description)
  term_label <- setNames(str_wrap(str_to_sentence(term_order), 45), term_order)

  plot_data <- combined %>%
    mutate(label = factor(term_label[Description], levels = rev(term_label[term_order])),
           Comparison = factor(Comparison, levels = c(LABEL_A, LABEL_B),
                                labels = c(genotype_strip_label(LABEL_A), genotype_strip_label(LABEL_B))))

  p_ora <- ggplot(plot_data, aes(x = signed_score, y = label)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(aes(colour = -log10(p.adjust)), size = 6, alpha = 0.9) +
    geom_text(aes(label = Count), colour = "white", size = 2.4, fontface = "bold") +
    scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
                            name = expression(-log[10](FDR)),
                            labels = scales::label_number(accuracy = 0.1)) +
    facet_wrap(~Comparison) +
    labs(title = "GO Biological Process (ORA, enrichGO)",
         subtitle = if (short_subtitle) sprintf("Top %d terms/direction/genotype, redundancy-collapsed | FDR < %.2f", n_per_direction, FDR)
                    else sprintf("Top %d terms/direction/genotype, redundancy-collapsed (Jaccard >= %.1f, global canonical label) | FDR < %.2f\n+ = enriched in UP genes, - = enriched in DOWN genes | number in dot = gene count | rows grouped: Shared / genotype-specific",
                                  n_per_direction, JACCARD_THRESH, FDR),
         x = "Signed log2(enrichment ratio)", y = NULL) +
    theme_pub + theme(panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.3),
                       strip.text = ggtext::element_markdown(face = "bold", size = 9))
  save_plot_both(p_ora, output_name, width = plot_width, height = plot_height)
}

cat("=== TASK 12: GO BP dotplot (enrichGO ORA, single shared y-axis, redundancy-collapsed) ===\n")
tryCatch(build_collapsed_shared_axis_plot(15, "AR07_GO_BP_ORA_dotplot_shared_axis_collapsed", 12),
         error = function(e) cat(sprintf("  [ERROR] TASK 12 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 13 -- pie-glyph variant of the same collapsed single-shared-y-axis GO
## BP data (TASK 12's term set/collapsing, reused): NOT split by Up/Down --
## one pie per (genotype, term) instead of up to two signed points. Pie size
## = total gene count (Up + Down combined), pie slices = Up (red) vs Down
## (blue) share of that count. A term significant in only one direction
## renders as a solid single-colour disc; a term significant in BOTH
## directions (e.g. KDM6A_ko's "forebrain development") renders as a
## genuinely split two-colour pie -- visually distinct from the shared-axis
## plots' two separate points for the same case.
##
## Requires scatterpie (geom_scatterpie/geom_scatterpie_legend) --
## coord_equal() is mandatory for the pies to render as true circles, not
## ellipses; the two genotype columns are placed far apart on x (1 and 4,
## not 1 and 2) purely so the coord_equal-constrained panel gets enough
## width to fit axis tick labels and the size legend without collision --
## the visible gap between the two columns is a layout artifact of this
## fixed-aspect requirement, not meaningful spacing.
########################################################################
cat("=== TASK 13: GO BP pie-glyph (enrichGO ORA, single shared y-axis, Up/Down as pie slices) ===\n")
tryCatch({
  go_ora_path <- file.path(TABLES_DIR, "R02_GO_enrichment_per_genotype.csv")
  if (check_input(go_ora_path) && requireNamespace("scatterpie", quietly = TRUE)) {
    N_PER_DIRECTION <- 15
    JACCARD_THRESH  <- 0.7

    parse_ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
    go_ora <- read_csv(go_ora_path, show_col_types = FALSE) %>%
      mutate(enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
             genes = str_split(geneID, "/"))

    top_by_direction <- function(name) {
      go_ora %>% filter(Comparison == name) %>% group_by(Direction) %>% arrange(p.adjust) %>%
        slice_head(n = N_PER_DIRECTION) %>% ungroup()
    }
    a_top <- top_by_direction(NAME_A)
    b_top <- top_by_direction(NAME_B)

    force_include <- function(name, other_terms) go_ora %>% filter(Comparison == name, Description %in% other_terms)
    a_raw <- bind_rows(a_top, force_include(NAME_A, unique(b_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
    b_raw <- bind_rows(b_top, force_include(NAME_B, unique(a_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)

    jaccard <- function(g1, g2) length(intersect(g1, g2)) / length(union(g1, g2))
    find_redundant_edges <- function(data) {
      n <- nrow(data)
      edges <- list()
      for (i in seq_len(n - 1)) for (j in seq((i + 1), n)) {
        if (jaccard(data$genes[[i]], data$genes[[j]]) >= JACCARD_THRESH) edges[[length(edges) + 1]] <- c(data$Description[i], data$Description[j])
      }
      edges
    }
    all_edges <- c(find_redundant_edges(a_raw), find_redundant_edges(b_raw))
    combined_raw <- bind_rows(a_raw, b_raw)
    all_terms <- unique(combined_raw$Description)
    redundancy_graph <- igraph::make_empty_graph(n = length(all_terms), directed = FALSE)
    igraph::V(redundancy_graph)$name <- all_terms
    if (length(all_edges) > 0) {
      edge_vec <- unlist(lapply(all_edges, function(e) match(e, all_terms)))
      redundancy_graph <- igraph::add_edges(redundancy_graph, edge_vec)
    }
    cluster_of <- setNames(igraph::components(redundancy_graph)$membership, all_terms)

    canonical_label <- combined_raw %>% group_by(Description) %>% summarise(min_padj = min(p.adjust), .groups = "drop") %>%
      mutate(cluster_id = cluster_of[Description]) %>%
      group_by(cluster_id) %>% mutate(canonical = Description[which.min(min_padj)]) %>% ungroup() %>%
      dplyr::select(Description, canonical)

    combined <- combined_raw %>%
      mutate(Description = canonical_label$canonical[match(Description, canonical_label$Description)]) %>%
      group_by(Comparison, Direction, Description) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup() %>%
      mutate(Comparison = factor(Comparison, levels = c(NAME_A, NAME_B), labels = c(LABEL_A, LABEL_B)))

    pie_data <- combined %>% group_by(Comparison, Description) %>%
      summarise(Up = sum(Count[Direction == "Up"]), Down = sum(Count[Direction == "Down"]), .groups = "drop") %>%
      mutate(total_count = Up + Down)

    a_terms <- pie_data$Description[pie_data$Comparison == LABEL_A]
    b_terms <- pie_data$Description[pie_data$Comparison == LABEL_B]
    common_final <- intersect(a_terms, b_terms)
    cat(sprintf("  %s %d terms | %s %d terms | %d shared\n", LABEL_A, length(a_terms), LABEL_B, length(b_terms), length(common_final)))

    term_category <- tibble(Description = union(a_terms, b_terms)) %>%
      mutate(category = case_when(Description %in% common_final ~ "Shared",
                                   Description %in% a_terms ~ paste0(LABEL_A, " only"), TRUE ~ paste0(LABEL_B, " only")))
    term_order <- pie_data %>% left_join(term_category, by = "Description") %>%
      group_by(Description, category) %>% summarise(sort_key = max(total_count), .groups = "drop") %>%
      mutate(category = factor(category, levels = c("Shared", paste0(LABEL_A, " only"), paste0(LABEL_B, " only")))) %>%
      arrange(category, desc(sort_key)) %>% pull(Description)

    n_terms <- length(term_order)
    y_lookup <- setNames(rev(seq_len(n_terms)), term_order)
    x_lookup <- setNames(c(1, 4), c(LABEL_A, LABEL_B))
    term_label <- setNames(str_wrap(str_to_sentence(term_order), 40), term_order)

    R_MIN <- 0.2; R_MAX <- 0.55
    max_total <- max(pie_data$total_count)
    pie_data <- pie_data %>% mutate(x = x_lookup[as.character(Comparison)], y = y_lookup[Description],
                                      r = R_MIN + (R_MAX - R_MIN) * sqrt(total_count / max_total))

    legend_counts <- c(5, 15, max_total)
    legend_radii <- R_MIN + (R_MAX - R_MIN) * sqrt(legend_counts / max_total)

    p_pie <- ggplot() +
      scatterpie::geom_scatterpie(aes(x = x, y = y, r = r), data = pie_data, cols = c("Up", "Down"), color = "grey30", linewidth = 0.15) +
      scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC"), name = "Direction") +
      scatterpie::geom_scatterpie_legend(legend_radii, x = 6.2, y = 4, n = 3, breaks = legend_radii,
                                          labeller = function(r) legend_counts[match(round(r, 6), round(legend_radii, 6))]) +
      annotate("text", x = 6.2, y = 7, label = "Gene count", fontface = "bold", size = 3, hjust = 0.5) +
      scale_x_continuous(breaks = x_lookup, labels = vapply(names(x_lookup), genotype_strip_label, character(1)), limits = c(0, 8)) +
      scale_y_continuous(breaks = y_lookup, labels = term_label[names(y_lookup)], limits = c(-1, n_terms + 1)) +
      coord_equal() +
      labs(title = "GO Biological Process (ORA, enrichGO) -- pie-glyph", x = NULL, y = NULL,
           subtitle = sprintf("Top %d terms/direction/genotype, redundancy-collapsed | FDR < %.2f | pie size = gene count, slices = Up/Down share | rows: Shared / genotype-specific",
                               N_PER_DIRECTION, FDR)) +
      theme_pub + theme(axis.text.x = ggtext::element_markdown(face = "bold", size = 11))
    save_plot_both(p_pie, "AR07_GO_BP_ORA_piechart_shared_axis", width = 10, height = 16)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 13 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 14 -- same design/pipeline as TASK 12, top 8 terms per direction per
## genotype instead of 15: a curated, paper-panel-sized cut. Kept as its own
## file (not overwriting TASK 12) since TASK 12 is the "look at everything"
## version and this is the "smaller, more presentable" one.
########################################################################
cat("=== TASK 14: GO BP dotplot (enrichGO ORA, single shared y-axis, redundancy-collapsed, top 8/direction/genotype) ===\n")
tryCatch(build_collapsed_shared_axis_plot(8, "AR07_GO_BP_ORA_dotplot_shared_axis_collapsed_top8", 8),
         error = function(e) cat(sprintf("  [ERROR] TASK 14 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 15 -- same as TASK 12/14, top 5 terms per direction per genotype:
## the most curated cut, closest to a main-figure panel.
########################################################################
cat("=== TASK 15: GO BP dotplot (enrichGO ORA, single shared y-axis, redundancy-collapsed, top 5/direction/genotype) ===\n")
tryCatch(build_collapsed_shared_axis_plot(5, "AR07_GO_BP_ORA_dotplot_shared_axis_collapsed_top5", 6, plot_width = 9, short_subtitle = TRUE),
         error = function(e) cat(sprintf("  [ERROR] TASK 15 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 16 -- patterning/lineage gene barplot, ported from Fig3_panels.R's
## Panel 3B (the fuller 17-gene, 6-category version -- adds "NC specifiers"
## and several Anterior-identity genes that TASK 2's smaller panel doesn't
## have). Brought in with AR07's own conventions, not Fig3's: GENOTYPE_COLORS
## (KDM6A_ko blue / KMT2D_Het orange, matching the PCA and every other
## genotype-coloured plot in this script) instead of Fig3_panels.R's
## deliberately-reversed FIG3_COLORS, and theme_pub instead of theme_fig3.
## Only the genotype-encoding bars get GENOTYPE_COLORS -- the gene_group
## facet strip text ("NPB dorsal", "Placodal", ...) and the plot title stay
## plain black via theme_pub's defaults, same as everywhere else in AR07;
## these are gene-category labels, not genotype, so they were never a
## candidate for the genotype-colour treatment applied to the GO dotplots.
########################################################################
cat("=== TASK 16: patterning/lineage gene barplot (from Fig3_panels.R Panel 3B) ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  if (check_input(full_a_path) && check_input(full_b_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)

    gene_groups <- tribble(
      ~gene_name, ~gene_group,
      "MSX1", "NPB dorsal", "MSX2", "NPB dorsal", "ZIC1", "NPB dorsal",
      "DLX3", "Placodal", "DLX5", "Placodal", "DLX6", "Placodal",
      "FOXD3", "NC specifiers", "TFAP2A", "NC specifiers", "TFAP2C", "NC specifiers", "SNAI2", "NC specifiers",
      "TWIST1", "NC specifiers", "MITF", "NC specifiers", "MEF2C", "NC specifiers",
      "SIX3", "Anterior identity", "RAX", "Anterior identity", "FOXG1", "Anterior identity",
      "OTX2", "Anterior identity", "LHX2", "Anterior identity", "DMBX1", "Anterior identity",
      "FGF17", "Anterior identity", "CALB2", "Anterior identity",
      "HOXA2", "Posterior drift", "TLX3", "Posterior drift",
      "CYP26A1", "RA metabolism", "CYP26C1", "RA metabolism"
    ) %>% mutate(gene_name = factor(gene_name, levels = gene_name),
                 gene_group = factor(gene_group, levels = unique(gene_group)))

    extract <- function(df, label) {
      df %>% filter(gene_name %in% gene_groups$gene_name) %>%
        transmute(gene_name, genotype = label, log2FC = log2FoldChange, padj = padj)
    }
    plot_data <- bind_rows(extract(rna_a, LABEL_A), extract(rna_b, LABEL_B)) %>%
      right_join(gene_groups, by = "gene_name") %>%
      mutate(gene_name = factor(gene_name, levels = levels(gene_groups$gene_name)),
             genotype = factor(genotype, levels = c(LABEL_B, LABEL_A)),
             significant = !is.na(padj) & padj < FDR & abs(log2FC) > LFC_THR)

    missing_genes <- setdiff(as.character(gene_groups$gene_name), plot_data$gene_name[!is.na(plot_data$log2FC)])
    cat(sprintf("  %d/%d genes found in both R02 tables\n",
                length(unique(gene_groups$gene_name)) - length(unique(missing_genes)), length(unique(gene_groups$gene_name))))
    if (length(missing_genes) > 0) cat(sprintf("  Not found (dropped, not fabricated): %s\n", paste(unique(missing_genes), collapse = ", ")))

    p <- ggplot(plot_data, aes(x = gene_name, y = log2FC, fill = genotype, alpha = significant)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
      geom_hline(yintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "grey40", linewidth = 0.4) +
      facet_grid(~gene_group, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = GENOTYPE_COLORS, name = "Genotype", breaks = c(LABEL_A, LABEL_B)) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), name = NULL,
                          labels = c(`TRUE` = sprintf("padj < %.2f and |log2FC| > %.2f", FDR, LFC_THR),
                                     `FALSE` = "not significant / below threshold")) +
      labs(title = "Patterning and lineage gene expression",
           x = NULL, y = "log2 fold change (vs WT)") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            legend.position = "bottom", legend.box = "vertical",
            panel.spacing = unit(0.4, "lines")) +
      guides(fill = guide_legend(override.aes = list(alpha = 1), order = 1),
             alpha = guide_legend(override.aes = list(fill = "grey40"), order = 2))
    save_plot_both(p, "AR07_patterning_lineage_genes_barplot", width = 14, height = 6)
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 16 failed: %s\n", conditionMessage(e))))
cat("\n")

########################################################################
## TASK 17 -- RNA-seq volcano plots, small paper-panel version, one per
## contrast, in BOTH y-axis flavours: raw p-value (matches
## A05_DESeq2.R::make_volcano_rawp() for the ATAC pairing, no fixed
## significance line since BH's threshold is rank- not value-dependent) and
## adjusted p-value/padj (matches A05_DESeq2.R::make_volcano() /
## R02_volcano_*.pdf, DOES get a fixed geom_hline at -log10(FDR)). Y-axis
## automatic; x-axis fixed at +/-15 with clipped points shown as triangles --
## a fully automatic x-axis let a single untested (padj=NA), near-zero-count
## gene per contrast (baseMean ~15-21, nominal log2FC ~+/-22) blow the scale
## out to +/-20 even though every real DEG sits within +/-13. Point color
## (Up/Down/NS) is always from padj regardless of which value is plotted.
##
## Point labels: top 15 DEGs by padj per contrast (overall ranking across
## Up+Down, not split evenly between the two directions) -- same convention
## as R02_volcano_*.pdf's top_genes (N=20 there; 15 here for this panel's
## smaller size).
########################################################################
N_LABELS_RAWP <- 50
CAP_X_RAWP <- 15
cat("=== TASK 17: RNA-seq volcano plots, small paper-panel, p and padj y-axes ===\n")
tryCatch({
  full_a_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_A))
  full_b_path <- file.path(TABLES_DIR, sprintf("R02_%s_full.csv", NAME_B))
  if (check_input(full_a_path) && check_input(full_b_path)) {
    rna_a <- read_csv(full_a_path, show_col_types = FALSE)
    rna_b <- read_csv(full_b_path, show_col_types = FALSE)

    theme_pub_rawp <- theme_minimal(base_size = 15) +
      theme(axis.text = element_text(color = "black"),
            axis.title = element_text(face = "bold"),
            axis.line = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "right")

    make_volcano_small <- function(res_df, name, suffix, y_col, y_label, hline_at, file_tag) {
      res_plot <- res_df %>%
        mutate(sig_category = case_when(
                 is.na(padj) | padj >= FDR | abs(log2FoldChange) < LFC_THR ~ "NS",
                 log2FoldChange > LFC_THR ~ "Up",
                 log2FoldChange < -LFC_THR ~ "Down",
                 TRUE ~ "NS"),
               neg_log10_y = -log10(.data[[y_col]]),
               x_capped = abs(log2FoldChange) > CAP_X_RAWP,
               log2FoldChange = pmin(pmax(log2FoldChange, -CAP_X_RAWP), CAP_X_RAWP))
      n_up <- sum(res_plot$sig_category == "Up", na.rm = TRUE)
      n_down <- sum(res_plot$sig_category == "Down", na.rm = TRUE)

      labels <- res_plot %>% filter(sig_category %in% c("Up", "Down")) %>% arrange(padj) %>% head(N_LABELS_RAWP)

      p <- ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_y)) +
        geom_point(aes(color = sig_category, shape = x_capped), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("Up" = "#e31a1c", "Down" = "#56B4E9", "NS" = "grey70"), name = NULL) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
        scale_x_continuous(limits = c(-CAP_X_RAWP, CAP_X_RAWP)) +
        geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "grey30") +
        { if (!is.null(hline_at)) geom_hline(yintercept = hline_at, linetype = "dashed", color = "grey30") } +
        geom_text_repel(data = labels, aes(label = gene_name), size = 2.8, fontface = "italic",
                         max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                         box.padding = 0.4, point.padding = 0.15, force = 3, force_pull = 0.5,
                         max.time = 5, max.iter = 100000, seed = 42) +
        theme_pub_rawp +
        labs(title = name,
             subtitle = sprintf("padj < %.2f, |log2FC| > %.2f | Up: %d, Down: %d", FDR, LFC_THR, n_up, n_down),
             x = "log2 Fold Change", y = y_label)
      save_plot_both(p, paste0("AR07_volcano_", file_tag, "_", suffix), width = 6.5, height = 5)
    }

    for (ct in list(list(df = rna_a, name = NAME_A, suffix = NAME_A),
                     list(df = rna_b, name = NAME_B, suffix = NAME_B))) {
      make_volcano_small(ct$df, ct$name, ct$suffix, "pvalue", "-log10(p-value)", NULL, "rawp")
      make_volcano_small(ct$df, ct$name, ct$suffix, "padj", "-log10(adjusted p-value)", -log10(FDR), "padj")
    }
  }
}, error = function(e) cat(sprintf("  [ERROR] TASK 17 failed: %s\n", conditionMessage(e))))
cat("\n")

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "AR07_Publication_Figures_session_info.txt"))
cat("[DONE] AR07_Publication_Figures complete\n")
