#!/usr/bin/env Rscript
# AR07b_Figure2_main_composite.R -- main-text Figure 2 as ONE composite for A4 print (190 x 185 mm, portrait).
# A  Venn of shared Up / Down DEGs between genotypes + Fisher OR and permutation p (from R04)
# B  KDM6A_ko volcano (RNA)         C  KMT2D_Het volcano (RNA)
# D  transcriptional convergence scatter       E  GO BP ORA dotplot (top 5 / direction / genotype, redundancy-collapsed)
# Rebuilt natively (not pasted PDFs) so font, sizes and colours are identical across panels: Arial/Helvetica-metric sans,
# 7 pt base, 6 pt gene labels, 5.5 pt minimum. Same data and thresholds as R04 / AR07 (padj < 0.05, |log2FC| > 1).
# Direction colours (red up / blue down) are used ONLY for direction; GO dot colour is a separate sequential purple.
# Writes results/Plots/AR07b_Figure2_main_composite.{pdf,png}; touches no existing file.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/publication_figures/AR07b_Figure2_main_composite.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2); library(ggrepel); library(patchwork); library(yaml); library(igraph) })

cfg <- yaml::yaml.load(iconv(paste(readLines("config/pipeline_config.yaml", warn = FALSE, encoding = "UTF-8"), collapse = "\n"), "UTF-8", "ASCII", sub = ""))
TAB <- "results/Tables"; PLT <- "results/Plots"
FDR <- cfg$thresholds$fdr; LFC <- cfg$thresholds$lfc_threshold
NAME_A <- "KDM6A_ko_vs_WT"; NAME_B <- "KMT2D_Het_vs_WT"; LAB_A <- "KDM6A_ko"; LAB_B <- "KMT2D_Het"

# ---- global style ---------------------------------------------------------------------------------------------
PT <- 72.27 / 25.4                                  # ggplot size unit (mm) -> pt conversion helper: size_mm = pt / .pt
.pt_mm <- function(pt) pt / ggplot2::.pt
FS_BASE <- 7; FS_AXIS <- 6.5; FS_TITLE <- 7.5; FS_LAB <- 5.5; FS_SMALL <- 5.5
COL_A <- "#56B4E9"; COL_B <- "#E69F00"   # genotype code used across the paper: KDM6A_ko blue, KMT2D_Het orange
UP <- "#D73027"; DOWN <- "#4575B4"   # direction colours (volcanoes, scatter); genotype colours only in the Venn
NSCOL <- "grey82"; DARK <- "grey15"
theme_fig <- theme_minimal(base_size = FS_BASE, base_family = "sans") +
  theme(text = element_text(colour = DARK), axis.text = element_text(colour = "black", size = FS_AXIS),
        axis.title = element_text(size = FS_BASE), plot.title = element_text(face = "bold", size = FS_TITLE, hjust = 0),
        plot.subtitle = element_text(size = FS_SMALL, colour = "grey30", hjust = 0),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        legend.text = element_text(size = FS_SMALL), legend.title = element_text(size = FS_SMALL), legend.key.size = unit(3, "mm"),
        plot.tag = element_text(face = "bold", size = 11), plot.margin = margin(2, 3, 2, 2))

# ---- data -----------------------------------------------------------------------------------------------------
rna <- lapply(c(A = NAME_A, B = NAME_B), function(ct) read_csv(file.path(TAB, sprintf("R02_%s_full.csv", ct)), show_col_types = FALSE))
cls <- function(d) d %>% mutate(cat = case_when(is.na(padj) | padj >= FDR | abs(log2FoldChange) < LFC ~ "NS", log2FoldChange > 0 ~ "Up", TRUE ~ "Down"))
rna <- lapply(rna, cls)
merged <- inner_join(rna$A %>% select(gene_id, gene_name, lfcA = log2FoldChange, catA = cat), rna$B %>% select(gene_id, lfcB = log2FoldChange, catB = cat), by = "gene_id") %>%
  mutate(Category = case_when(catA == "Up" & catB == "Up" ~ "Shared Up", catA == "Down" & catB == "Down" ~ "Shared Down", TRUE ~ "Other"))
ov <- read_csv(file.path(TAB, "R04_overlap_effect_size.csv"), show_col_types = FALSE)
rho <- cor(merged$lfcA, merged$lfcB, method = "spearman")
cat(sprintf("[AR07b] shared up %d down %d | n_a %d n_b %d n_both %d | OR %.1f perm p %.2g | rho %.3f (n=%d)\n", sum(merged$Category == "Shared Up"), sum(merged$Category == "Shared Down"), ov$n_a, ov$n_b, ov$n_both, ov$OR, ov$perm_p, rho, nrow(merged)))

# ---- A: area-proportional Venns (Up above, Down below) ------------------------------------------------------------
lens <- function(r1, r2, d) { if (d >= r1 + r2) return(0); if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1)) + r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2)) - 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)) }
circ <- function(cx, cy, r, id) tibble(x = cx + r * cos(seq(0, 2 * pi, length.out = 200)), y = cy + r * sin(seq(0, 2 * pi, length.out = 200)), id = id)
venn_layers <- function(na, nb, nboth, cy, dir_lab) {
  ra <- sqrt(na / pi); rb <- sqrt(nb / pi); d <- uniroot(function(d) lens(ra, rb, d) - nboth, c(abs(ra - rb) + 1e-6, ra + rb - 1e-6))$root
  cxa <- -d / 2 - (ra - rb) / 4; cxb <- cxa + d
  lo <- cxb - rb; hi <- cxa + ra
  list(polys = bind_rows(circ(cxa, cy, ra, "A"), circ(cxb, cy, rb, "B")),
       txt = tibble(x = c((cxa - ra + lo) / 2, (lo + hi) / 2, (hi + cxb + rb) / 2), y = cy, lab = c(na - nboth, nboth, nb - nboth)),
       hdr = tibble(x = c(cxa + 0.5, cxb - 0.5), hj = c(1, 0), y = cy + max(ra, rb) + 1.6, lab = c(LAB_A, LAB_B)), ytop = cy + max(ra, rb), ybot = cy - max(ra, rb))
}
na_u <- sum(rna$A$cat == "Up"); nb_u <- sum(rna$B$cat == "Up"); sh_u <- sum(merged$Category == "Shared Up")
na_d <- sum(rna$A$cat == "Down"); nb_d <- sum(rna$B$cat == "Down"); sh_d <- sum(merged$Category == "Shared Down")
vu <- venn_layers(na_u, nb_u, sh_u, cy = 17, "u"); vd <- venn_layers(na_d, nb_d, sh_d, cy = -8, "d")
ann <- sprintf("DEG overlap: %d genes\nFisher OR = %.1f (95%% CI %.0f\u2013%.0f)\nPermutation p < 10\u207b\u2074 (%s perm.)",
               ov$n_both, ov$OR, ov$OR_CI_low, ov$OR_CI_high, format(ov$n_perm, big.mark = ","))
pA <- ggplot() +
  geom_polygon(data = vu$polys, aes(x, y, group = id, fill = id, colour = id), alpha = 0.55, linewidth = 0.4) +
  geom_polygon(data = vd$polys, aes(x, y, group = id, fill = id, colour = id), alpha = 0.55, linewidth = 0.4) +
  scale_fill_manual(values = c(A = COL_A, B = COL_B), guide = "none") + scale_colour_manual(values = c(A = COL_A, B = COL_B), guide = "none") +
  geom_text(data = vu$txt, aes(x, y, label = lab), size = .pt_mm(7.5), fontface = "bold", colour = DARK) + geom_text(data = vd$txt, aes(x, y, label = lab), size = .pt_mm(7.5), fontface = "bold", colour = DARK) +
  geom_text(data = vu$hdr, aes(x, y, label = lab, hjust = hj), size = .pt_mm(FS_LAB), fontface = "bold", colour = DARK) + geom_text(data = vd$hdr, aes(x, y, label = lab, hjust = hj), size = .pt_mm(FS_LAB), fontface = "bold", colour = DARK) +
  annotate("text", x = -26.5, y = 17, label = "Upregulated", angle = 90, hjust = 0.5, size = .pt_mm(FS_LAB), fontface = "bold", colour = DARK) +
  annotate("text", x = -26.5, y = -8, label = "Downregulated", angle = 90, hjust = 0.5, size = .pt_mm(FS_LAB), fontface = "bold", colour = DARK) +
  annotate("label", x = -21, y = -20.5, label = ann, hjust = 0, vjust = 1, size = .pt_mm(FS_SMALL), lineheight = 1.1, colour = DARK, fill = "white", label.size = 0.3, label.r = unit(0, "pt"), label.padding = unit(1.8, "mm")) +
  coord_equal(xlim = c(-27, 15), ylim = c(-34, 31), clip = "off") + labs(title = "Shared DEGs") + theme_void(base_size = FS_BASE, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = FS_TITLE, hjust = 0), plot.tag = element_text(face = "bold", size = 11), plot.margin = margin(2, 3, 2, 2))

# ---- B/C: volcanoes ----------------------------------------------------------------------------------------------
XL <- 10; YC <- 22
genes_A <- c("SIX3","RAX","FOXG1","OTX2","WNT4","WNT5B","WNT8B","DKK1","FZD5","SFRP1","PAX7","BAMBI","BMPR1B","CALB2","CYP26A1","DLX3","DLX6","DMBX1","FGF17","GDF10","MEF2C","MITF","TLX3","WNT9A")
genes_B <- c("MSX1","MSX2","ZIC1","BMP7","SMAD9","ID1","ID2","NOG","WNT2B","WNT7A","SP5","CD44","RNF43","SOX10","BAMBI","BMPR1B","CYP26A1","DLX3","DLX6","FGF17","GDF10","MEF2C","MITF","TLX3","WNT9A")
volc <- function(d, lab_genes, ttl, show_ylab, pcol) {
  d <- d %>% mutate(y = -log10(pmax(pvalue, 1e-300)), cap = abs(log2FoldChange) > XL | y > YC, x = pmin(pmax(log2FoldChange, -XL), XL), y = pmin(y, YC), cat = factor(cat, c("NS", "Down", "Up")))
  n_up <- sum(d$cat == "Up"); n_dn <- sum(d$cat == "Down")
  lab <- d %>% filter(toupper(gene_name) %in% lab_genes) %>% distinct(gene_name, .keep_all = TRUE)
  cat(sprintf("[AR07b] %s: %d/%d label genes found\n", ttl, nrow(lab), length(lab_genes)))
  rl <- function(dd, xl) geom_text_repel(data = dd, aes(label = gene_name), size = .pt_mm(FS_LAB), fontface = "italic", colour = "black", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey45",
                                         box.padding = 0.25, point.padding = 0.15, force = 6, force_pull = 0.2, xlim = xl, direction = "both", max.time = 5, max.iter = 100000, seed = 42, bg.color = "white", bg.r = 0.12)
  ggplot(d %>% arrange(cat), aes(x, y, colour = cat, shape = cap)) + geom_point(size = 0.55, alpha = 0.6, stroke = 0) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey35", linewidth = 0.25) +
    scale_colour_manual(values = c(NS = NSCOL, Down = DOWN, Up = UP), guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    rl(lab %>% filter(x >= 0), c(0, XL)) + rl(lab %>% filter(x < 0), c(-XL, 0)) +
    scale_x_continuous(limits = c(-XL, XL), breaks = seq(-10, 10, 5)) + scale_y_continuous(limits = c(0, YC), breaks = seq(0, 20, 5)) +
    labs(title = ttl, subtitle = sprintf("Up %d | Down %d (padj < %.2f, |log2FC| > %g)", n_up, n_dn, FDR, LFC), x = expression(log[2]~fold~change), y = if (show_ylab) expression(-log[10]~italic(p)) else NULL) + theme_fig
}
pB <- volc(rna$A, genes_A, paste0(LAB_A, " vs WT"), TRUE, COL_A)
pC <- volc(rna$B, genes_B, paste0(LAB_B, " vs WT"), TRUE, COL_B)

# ---- D: convergence scatter ---------------------------------------------------------------------------------------
genes_r04 <- c("COL3A1","FN1","SPHK1","SNCG","JUN","ATOH1","NTRK2","GRIN3A","SCN3A","KCNJ6","CNTN3","NCAM2","NEFM","STMN2","WNT3A","BMP4","ID4","DLX5","LHX2","CYP26C1","BAMBI","BMPR1B","CALB2","CYP26A1","DLX3","DLX6","DMBX1","FGF17","GDF10","GDF7","HOXA2","MEF2C","MITF","SNAI2","TLX3","TWIST1","WNT8A","WNT9A")
genes_r04_vlight <- c("COL3A1","FN1","SPHK1","SNCG","JUN","ATOH1","NTRK2","GRIN3A","SCN3A","KCNJ6","CNTN3","NCAM2","NEFM","STMN2","WNT3A","BMP4","ID4","DLX5","LHX2","CYP26C1","HOXA2","SNAI2","TWIST1","WNT8A","GDF7","DLX3","DLX6","DMBX1","FGF17","CALB2","GRHL2","SFN","RIPK4","CDH3","PRSS8","HRH3","GRIK3","CHRM4","HTR1D","SCN5A","GNAO1","KCNH2","KCNB1","GRM3")
ld <- merged %>% filter(toupper(gene_name) %in% genes_r04_vlight) %>% filter(Category != "Shared Down" | toupper(gene_name) %in% genes_r04) %>% distinct(gene_name, .keep_all = TRUE)
cat(sprintf("[AR07b] scatter labels: %d/%d found\n", nrow(ld), length(genes_r04_vlight)))
AR <- c(-8, 8)
pD <- ggplot(merged, aes(lfcB, lfcA)) +
  geom_point(data = filter(merged, Category == "Other"), colour = "grey82", size = 0.25, alpha = 0.5, stroke = 0) +
  geom_point(data = filter(merged, Category == "Shared Down"), colour = DOWN, size = 0.8, alpha = 0.7, stroke = 0) + geom_point(data = filter(merged, Category == "Shared Up"), colour = UP, size = 0.8, alpha = 0.7, stroke = 0) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_text_repel(data = ld, aes(label = gene_name), size = .pt_mm(FS_LAB), fontface = "italic", colour = "black", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey45", box.padding = 0.3, point.padding = 0.15, force = 6, force_pull = 0.2, max.time = 5, max.iter = 100000, seed = 42, bg.color = "white", bg.r = 0.12) +
  annotate("text", x = AR[2], y = AR[1], hjust = 1, vjust = 0, size = .pt_mm(FS_SMALL), lineheight = 1.1, label = sprintf("Spearman ρ = %.2f, n = %s genes\nshared DEGs: %d up, %d down", rho, format(nrow(merged), big.mark = ","), sum(merged$Category == "Shared Up"), sum(merged$Category == "Shared Down"))) +
  coord_equal(xlim = AR, ylim = AR) + labs(title = "Transcriptional convergence", subtitle = sprintf("Shared DEGs: padj < %.2f and |log2FC| > %g in both genotypes", FDR, LFC), x = bquote(log[2]~fold~change~"("*.(LAB_B)~vs~WT*")"), y = bquote(log[2]~fold~change~"("*.(LAB_A)~vs~WT*")")) + theme_fig

# ---- E: GO BP ORA dotplot (top 5 / direction / genotype, global redundancy collapse) ---------------------------------
N_PER <- 5; JAC <- 0.7
parse_ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
go <- read_csv(file.path(TAB, "R02_GO_enrichment_per_genotype.csv"), show_col_types = FALSE) %>% mutate(er = parse_ratio(GeneRatio) / parse_ratio(BgRatio), genes = str_split(geneID, "/"))
topd <- function(nm) go %>% filter(Comparison == nm) %>% group_by(Direction) %>% arrange(p.adjust) %>% slice_head(n = N_PER) %>% ungroup()
a_top <- topd(NAME_A); b_top <- topd(NAME_B)
fi <- function(nm, other) go %>% filter(Comparison == nm, Description %in% other)
a_raw <- bind_rows(a_top, fi(NAME_A, unique(b_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
b_raw <- bind_rows(b_top, fi(NAME_B, unique(a_top$Description))) %>% distinct(Description, Direction, .keep_all = TRUE)
jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
edges <- function(d) { e <- list(); if (nrow(d) > 1) for (i in seq_len(nrow(d) - 1)) for (j in (i + 1):nrow(d)) if (jac(d$genes[[i]], d$genes[[j]]) >= JAC) e[[length(e) + 1]] <- c(d$Description[i], d$Description[j]); e }
all_e <- c(edges(a_raw), edges(b_raw)); craw <- bind_rows(a_raw, b_raw); terms <- unique(craw$Description)
g <- make_empty_graph(n = length(terms), directed = FALSE); V(g)$name <- terms; if (length(all_e)) g <- add_edges(g, unlist(lapply(all_e, function(e) match(e, terms))))
clu <- setNames(components(g)$membership, terms)
canon <- craw %>% group_by(Description) %>% summarise(mp = min(p.adjust), .groups = "drop") %>% mutate(cl = clu[Description]) %>% group_by(cl) %>% mutate(canonical = Description[which.min(mp)]) %>% ungroup()
comb <- craw %>% mutate(Description = canon$canonical[match(Description, canon$Description)]) %>% group_by(Comparison, Direction, Description) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(score = ifelse(Direction == "Down", -log2(er), log2(er)), Comparison = factor(Comparison, c(NAME_A, NAME_B), c(LAB_A, LAB_B)))
common <- intersect(comb$Description[comb$Comparison == LAB_A], comb$Description[comb$Comparison == LAB_B])
cat_of <- function(d) ifelse(d %in% common, "Shared", ifelse(d %in% comb$Description[comb$Comparison == LAB_A], "A", "B"))
ord <- comb %>% mutate(cat = factor(cat_of(Description), c("Shared", "A", "B"))) %>% group_by(Description, cat) %>% summarise(k = mean(score), .groups = "drop") %>% arrange(cat, desc(k)) %>% pull(Description)
short <- c("Transforming growth factor beta receptor superfamily signaling pathway" = "TGF-beta receptor signaling pathway", "Cell-cell adhesion via plasma-membrane adhesion molecules" = "Cell-cell adhesion (plasma membrane)")
lbl <- setNames(str_wrap(ifelse(str_to_sentence(ord) %in% names(short), short[str_to_sentence(ord)], str_to_sentence(ord)), 30), ord)
pdat <- comb %>% mutate(label = factor(lbl[Description], levels = rev(lbl[ord])), lg = -log10(p.adjust))
cat(sprintf("[AR07b] GO rows: %d terms, %d shared\n", length(ord), length(common)))
pE <- ggplot(pdat, aes(score, label)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_point(aes(colour = lg), size = 3.4) + geom_text(aes(label = Count), colour = "white", size = .pt_mm(FS_SMALL), fontface = "bold") +
  scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"), name = expression(-log[10](FDR))) +
  scale_x_continuous(expand = expansion(mult = 0.14)) + facet_wrap(~Comparison) +
  labs(title = "GO biological process (ORA)", subtitle = sprintf("Top %d terms per direction and genotype, redundancy-collapsed", N_PER), x = expression("Signed"~log[2]~"enrichment ("*"+ up, − down)"), y = NULL) +
  theme_fig + theme(panel.grid.major.y = element_line(colour = "grey92", linetype = "dashed", linewidth = 0.25), panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3),
                    strip.text = element_text(face = "bold", size = FS_BASE), axis.text.y = element_text(size = FS_SMALL, lineheight = 0.85), legend.position = "right", legend.key.height = unit(8, "mm"))

# ---- assemble --------------------------------------------------------------------------------------------------------
pA <- pA + theme(plot.margin = margin(2, 2, 2, 2))
fig <- ((pB | pC | pA) + plot_layout(widths = c(1, 1, 1.1))) / ((pD | pE) + plot_layout(widths = c(1.12, 0.88))) + plot_layout(heights = c(0.95, 1.15)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
W <- 190 / 25.4; H <- 185 / 25.4
ggsave(file.path(PLT, paste0("AR07b_Figure2_main_composite", Sys.getenv("OUT_SUFFIX", ""), ".pdf")), fig, width = W, height = H, device = cairo_pdf)
ggsave(file.path(PLT, paste0("AR07b_Figure2_main_composite", Sys.getenv("OUT_SUFFIX", ""), ".png")), fig, width = W, height = H, dpi = 300)
cat("[DONE] AR07b_Figure2_main_composite\n")
