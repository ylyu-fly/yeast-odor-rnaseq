## Subject: DESeq and GO
## Project: Yeast Odor RNA-seq
## Part: 
## Author: Yang Lyu
## Date created: 08/01/2026
## Date modified: 09/30/2026
## May require R 4.6.1

rm(list = ls())

library(DESeq2)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(gridExtra)

setwd("~/Research/Projects/RNA-seq/YeastOdor/")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

acute <- read.table("./ProcessedData/reads_count.txt", skip = 1, header = TRUE, row.names = 1)
age20 <- read.table("./ProcessedData/reads_count_age20.txt", skip = 1, header = TRUE, row.names = 1)

acute_counts <- acute[, 6:ncol(acute)]
age20_counts <- age20[, 6:ncol(age20)]

# Remove library prefix and sample numbers
colnames(acute_counts) <- gsub("^.*_\\d+_\\d+_", "", colnames(acute_counts))
colnames(age20_counts) <- gsub("^.*_\\d+_\\d+_", "", colnames(age20_counts))

# Remove BAM suffix
colnames(acute_counts) <- gsub("_Aligned.*", "", colnames(acute_counts))
colnames(age20_counts) <- gsub("_Aligned.*", "", colnames(age20_counts))

cat("Acute dataset gene x sample:", nrow(acute_counts), "x", ncol(acute_counts), "\n")
cat("Chronic dataset gene x sample:", nrow(age20_counts), "x", ncol(age20_counts), "\n")

# ============================================================================
# 2. PARSE METADATA
# ============================================================================

parse_sample_info <- function(names) {
  info <- data.frame(sample = names, tissue = NA, genotype = NA, treatment = NA)
  for (i in seq_along(names)) {
    if (grepl("Head", names[i])) info$tissue[i] <- "head"
    else if (grepl("Ovary", names[i])) info$tissue[i] <- "ovary"
    if (grepl("WT", names[i])) info$genotype[i] <- "WT"
    else if (grepl("5HT2A", names[i])) info$genotype[i] <- "5HT2A"
    if (grepl("Yeast", names[i])) info$treatment[i] <- "Yeast"
    else if (grepl("NoOdor", names[i])) info$treatment[i] <- "NoOdor"
  }
  return(info)
}

acute_meta <- parse_sample_info(colnames(acute_counts))
age20_meta <- parse_sample_info(colnames(age20_counts))

# Sample distribution
print(table(acute_meta$tissue, acute_meta$genotype, acute_meta$treatment))
print(table(age20_meta$tissue, age20_meta$genotype, age20_meta$treatment))

# ============================================================================
# 3. ACUTE PCA SETUP
# ============================================================================

coldata_acute_all <- data.frame(
  tissue = factor(acute_meta$tissue, levels = c("ovary", "head")),
  genotype = factor(acute_meta$genotype, levels = c("WT", "5HT2A")),
  treatment = factor(acute_meta$treatment, levels = c("NoOdor", "Yeast")),
  row.names = acute_meta$sample
)

acute_counts_all <- acute_counts[, rownames(coldata_acute_all)]
dds_acute_pca <- DESeqDataSetFromMatrix(
  countData = as.matrix(acute_counts_all),
  colData = coldata_acute_all,
  design = ~ tissue + genotype + treatment
)

keep <- rowSums(counts(dds_acute_pca)) >= 10
dds_acute_pca <- dds_acute_pca[keep, ]
dds_acute_pca <- DESeq(dds_acute_pca, quiet = TRUE)

vsd_acute <- varianceStabilizingTransformation(dds_acute_pca, blind = TRUE)
pca_acute <- prcomp(t(assay(vsd_acute)))
var_exp_acute <- pca_acute$sdev^2 / sum(pca_acute$sdev^2) * 100

# PCA variance
cat(sprintf("  PC1: %.2f%%\n", var_exp_acute[1]))
cat(sprintf("  PC2: %.2f%%\n", var_exp_acute[2]))
cat(sprintf("  PC3: %.2f%%\n", var_exp_acute[3]))

# ============================================================================
# 4. ACUTE ANALYSIS BY TISSUE
# ============================================================================

acute_results <- list()

for (tissue in c("head", "ovary")) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat(sprintf("ACUTE EXPOSURE (4 hours) - %s\n", toupper(tissue)))
  cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
  
  # Filter to tissue
  idx <- acute_meta$tissue == tissue
  counts_t <- acute_counts[, idx]
  meta_t <- acute_meta[idx, ]
  
  coldata_t <- data.frame(
    genotype = factor(meta_t$genotype, levels = c("WT", "5HT2A")),
    treatment = factor(meta_t$treatment, levels = c("NoOdor", "Yeast")),
    row.names = meta_t$sample
  )
  
  counts_t <- counts_t[, rownames(coldata_t)]
  
  # DESeq2 with factorial design
  dds_t <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts_t),
    colData = coldata_t,
    design = ~ genotype + treatment + genotype:treatment
  )
  
  keep_t <- rowSums(counts(dds_t)) >= 10
  dds_t <- dds_t[keep_t, ]
  dds_t <- DESeq(dds_t, quiet = TRUE)
  
  cat("Genes analyzed:", nrow(dds_t), "\n\n")
  
  # Test contrasts
  res_geno <- results(dds_t, name = "genotype_5HT2A_vs_WT")
  sig_geno <- (abs(res_geno$log2FoldChange) >= 0.5) & (res_geno$padj <= 0.05) & (!is.na(res_geno$padj))
  
  res_treat <- results(dds_t, name = "treatment_Yeast_vs_NoOdor")
  sig_treat <- (abs(res_treat$log2FoldChange) >= 0.5) & (res_treat$padj <= 0.05) & (!is.na(res_treat$padj))
  
  res_inter <- results(dds_t, name = "genotype5HT2A.treatmentYeast")
  sig_inter <- (abs(res_inter$log2FoldChange) >= 0.5) & (res_inter$padj <= 0.05) & (!is.na(res_inter$padj))
  
  acute_results[[tissue]] <- list(
    dds = dds_t,
    vsd = varianceStabilizingTransformation(dds_t, blind = TRUE),
    pca = prcomp(t(assay(varianceStabilizingTransformation(dds_t, blind = TRUE)))),
    coldata = coldata_t,
    main_genotype = list(res = res_geno, sig = sig_geno,
                         n = sum(sig_geno, na.rm = TRUE),
                         n_up = sum((res_geno$log2FoldChange >= 0.5) & sig_geno, na.rm = TRUE),
                         n_down = sum((res_geno$log2FoldChange <= -0.5) & sig_geno, na.rm = TRUE)),
    main_treatment = list(res = res_treat, sig = sig_treat,
                          n = sum(sig_treat, na.rm = TRUE),
                          n_up = sum((res_treat$log2FoldChange >= 0.5) & sig_treat, na.rm = TRUE),
                          n_down = sum((res_treat$log2FoldChange <= -0.5) & sig_treat, na.rm = TRUE)),
    interaction = list(res = res_inter, sig = sig_inter,
                       n = sum(sig_inter, na.rm = TRUE),
                       n_up = sum((res_inter$log2FoldChange >= 0.5) & sig_inter, na.rm = TRUE),
                       n_down = sum((res_inter$log2FoldChange <= -0.5) & sig_inter, na.rm = TRUE))
  )
  
  cat("Contrasts tested (FDR ≤ 0.05, |log2FC| ≥ 0.5):\n")
  cat(sprintf("  Genotype (5HT2A vs WT):        %3d DEGs (%3d up, %3d down)\n",
              acute_results[[tissue]]$main_genotype$n,
              acute_results[[tissue]]$main_genotype$n_up,
              acute_results[[tissue]]$main_genotype$n_down))
  cat(sprintf("  Treatment (Yeast vs NoOdor):   %3d DEGs (%3d up, %3d down)\n",
              acute_results[[tissue]]$main_treatment$n,
              acute_results[[tissue]]$main_treatment$n_up,
              acute_results[[tissue]]$main_treatment$n_down))
  cat(sprintf("  Interaction (genotype×treat):  %3d DEGs (%3d up, %3d down)\n",
              acute_results[[tissue]]$interaction$n,
              acute_results[[tissue]]$interaction$n_up,
              acute_results[[tissue]]$interaction$n_down))
}

# ============================================================================
# 5. ACUTE FIGURES
# ============================================================================

par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1), family = "sans")

pca_coords <- data.frame(
  PC1 = pca_acute$x[, 1],
  PC2 = pca_acute$x[, 2],
  PC3 = pca_acute$x[, 3],
  tissue = coldata_acute_all$tissue,
  genotype = coldata_acute_all$genotype,
  treatment = coldata_acute_all$treatment
)

# Shape by genotype: circle (WT) vs triangle (5HT2A)
pchs <- ifelse(pca_coords$genotype == "WT", 21, 24)

# Fill color by genotype + treatment
cols <- rep(NA, nrow(pca_coords))
cols[pca_coords$genotype == "WT" & pca_coords$treatment == "NoOdor"] <- "grey89"       # light grey
cols[pca_coords$genotype == "WT" & pca_coords$treatment == "Yeast"] <- "#4A4A4A"        # dark grey
cols[pca_coords$genotype == "5HT2A" & pca_coords$treatment == "NoOdor"] <- "#22BDD2"    # light blue
cols[pca_coords$genotype == "5HT2A" & pca_coords$treatment == "Yeast"] <- "#333a8c"     # dark blue

border_cols <- "#1a1a1a"  # Dark border color for all
cex_size <- 2.2

plot(pca_coords$PC1, pca_coords$PC3, pch = pchs, cex = cex_size, col = border_cols, bg = cols,
     main = "Young flies",
     xlab = sprintf("PC1 (%.1f%%)", var_exp_acute[1]),
     ylab = sprintf("PC3 (%.1f%%)", var_exp_acute[3]),
     xlim = c(-250, 250), ylim = c(-50, 50), 
     las = 1, cex.lab = 1.5, cex.main = 1.5)
abline(h = 0, v = 0, lty = 2, col = "gray70")

# Add tissue region labels based on PC1 range
pc1_min <- min(pca_coords$PC1)
pc1_max <- max(pca_coords$PC1)
pc1_mid <- (pc1_min + pc1_max) / 2
pc3_top <- max(pca_coords$PC3) * 0.9

text(pc1_min + (pc1_mid - pc1_min) * 0.5, 48, "Ovary", cex = 1.6, adj = 0.5, col = "gray30", font = 2)
text(pc1_mid + (pc1_max - pc1_mid) * 0.5, 48, "Head", cex = 1.6, adj = 0.5, col = "gray30", font = 2)

legend("bottomleft", 
       legend = c("Wildtype / No Odor",
                  "Wildtype / Yeast Odor",
                  bquote(italic("5-HT2A") ~ "/ No Odor"),
                  bquote(italic("5-HT2A") ~ "/ Yeast Odor")),
       pch = c(21, 21, 24, 24),
       pt.bg = c("grey89", "#4A4A4A", "#22BDD2", "#333a8c"),
       col = "gray30",
       cex = 0.8, bty = "n", pt.cex = 1.5)

# Volcano plots by tissue

for (tissue in c("head", "ovary")) {
  
  res <- acute_results[[tissue]]$main_genotype$res
  sig_mask <- acute_results[[tissue]]$main_genotype$sig
  
  x_not_sig <- res[!sig_mask, "log2FoldChange"]
  y_not_sig <- -log10(res[!sig_mask, "padj"] + 1e-300)
  
  up_mask <- sig_mask & (res$log2FoldChange >= 0.5)
  down_mask <- sig_mask & (res$log2FoldChange <= -0.5)
  
  plot(x_not_sig, y_not_sig, pch = 20, col = "gray60", cex = 0.4,
       xlim = c(-8, 8), ylim = c(0, 80),
       main = sprintf("%s: 5-HT2A vs Wildtype", tissue),
       xlab = "log2(Fold Change)", ylab = "-log10(adjusted p-value)",
       las = 1, cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2)
  
  points(res[up_mask, "log2FoldChange"], -log10(res[up_mask, "padj"] + 1e-300),
         pch = 20, col = "#C93735", cex = 0.5)
  points(res[down_mask, "log2FoldChange"], -log10(res[down_mask, "padj"] + 1e-300),
         pch = 20, col = "#22BDD2", cex = 0.5)
  
  abline(h = -log10(0.05), col = "gray11", lty = 2, lwd = 1.5)
  abline(v = c(-0.5, 0.5), col = "gray11", lty = 2, lwd = 1.5)
  
  # Add cutoff text annotations
  text(-7, -log10(0.05) + 3, "FDR = 0.05", cex = 1.2, col = "black", font = 2)
  
  legend("topleft",
         legend = c(sprintf("Up (n=%d)", acute_results[[tissue]]$main_genotype$n_up),
                    sprintf("Down (n=%d)", acute_results[[tissue]]$main_genotype$n_down), "NS"),
         col = c("#C93735", "#22BDD2", "gray60"), pch = 20, cex = 1.5, bty = "n")
}

## For paper main figure, plot head only
tissue <- "head"
res <- acute_results[[tissue]]$main_genotype$res
sig_mask <- acute_results[[tissue]]$main_genotype$sig

x_not_sig <- res[!sig_mask, "log2FoldChange"]
y_not_sig <- -log10(res[!sig_mask, "padj"] + 1e-300)

up_mask <- sig_mask & (res$log2FoldChange >= 0.5)
down_mask <- sig_mask & (res$log2FoldChange <= -0.5)

plot(x_not_sig, y_not_sig, pch = 20, col = "gray60", cex = 0.4,
     xlim = c(-8, 8), ylim = c(0, 80),
     main = sprintf("Heads of young flies: 5-HT2A vs Wildtype"),
     xlab = "log2(Fold Change)", ylab = "-log10(adjusted p-value)",
     las = 1, cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2)

points(res[up_mask, "log2FoldChange"], -log10(res[up_mask, "padj"] + 1e-300),
       pch = 20, col = "#C93735", cex = 0.5)
points(res[down_mask, "log2FoldChange"], -log10(res[down_mask, "padj"] + 1e-300),
       pch = 20, col = "#22BDD2", cex = 0.5)

abline(h = -log10(0.05), col = "gray11", lty = 2, lwd = 1.5)
abline(v = c(-0.5, 0.5), col = "gray11", lty = 2, lwd = 1.5)

# Add cutoff text annotations
text(-5, -log10(0.05) + 3, "FDR = 0.05", cex = 1.2, col = "black", font = 2)

legend("topleft",
       legend = c(sprintf("Up (n=%d)", acute_results[[tissue]]$main_genotype$n_up),
                  sprintf("Down (n=%d)", acute_results[[tissue]]$main_genotype$n_down), "NS"),
       col = c("#C93735", "#22BDD2", "gray60"), pch = 20, cex = 1.2, bty = "n")

# ============================================================================
# 6. CHRONIC ANALYSIS - HEAD ONLY
# ============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n", sep="")
cat("CHRONIC EXPOSURE (Day 20) - HEAD\n")
cat(paste(rep("=", 70), collapse=""), "\n\n", sep="")

coldata_age20 <- data.frame(
  genotype = factor(age20_meta$genotype, levels = c("WT", "5HT2A")),
  treatment = factor(age20_meta$treatment, levels = c("NoOdor", "Yeast")),
  row.names = age20_meta$sample
)

age20_counts <- age20_counts[, rownames(coldata_age20)]

dds_age20 <- DESeqDataSetFromMatrix(
  countData = as.matrix(age20_counts),
  colData = coldata_age20,
  design = ~ genotype + treatment + genotype:treatment
)

keep <- rowSums(counts(dds_age20)) >= 10
dds_age20 <- dds_age20[keep, ]
dds_age20 <- DESeq(dds_age20, quiet = TRUE)

cat("Genes analyzed:", nrow(dds_age20), "\n\n")

# PCA
vsd_age20 <- varianceStabilizingTransformation(dds_age20, blind = TRUE)
pca_age20 <- prcomp(t(assay(vsd_age20)))
var_exp_age20 <- pca_age20$sdev^2 / sum(pca_age20$sdev^2) * 100

cat("PCA Variance:\n")
cat(sprintf("  PC1: %.2f%%\n", var_exp_age20[1]))
cat(sprintf("  PC2: %.2f%%\n", var_exp_age20[2]))
cat(sprintf("  PC3: %.2f%%\n", var_exp_age20[3]))

# Test contrasts
res_geno <- results(dds_age20, name = "genotype_5HT2A_vs_WT")
sig_geno <- (abs(res_geno$log2FoldChange) >= 0.5) & (res_geno$padj <= 0.05) & (!is.na(res_geno$padj))

res_treat <- results(dds_age20, name = "treatment_Yeast_vs_NoOdor")
sig_treat <- (abs(res_treat$log2FoldChange) >= 0.5) & (res_treat$padj <= 0.05) & (!is.na(res_treat$padj))

res_inter <- results(dds_age20, name = "genotype5HT2A.treatmentYeast")
sig_inter <- (abs(res_inter$log2FoldChange) >= 0.5) & (res_inter$padj <= 0.05) & (!is.na(res_inter$padj))

age20_results <- list(
  main_genotype = list(res = res_geno, sig = sig_geno,
                       n = sum(sig_geno, na.rm = TRUE),
                       n_up = sum((res_geno$log2FoldChange >= 0.5) & sig_geno, na.rm = TRUE),
                       n_down = sum((res_geno$log2FoldChange <= -0.5) & sig_geno, na.rm = TRUE)),
  main_treatment = list(res = res_treat, sig = sig_treat,
                        n = sum(sig_treat, na.rm = TRUE),
                        n_up = sum((res_treat$log2FoldChange >= 0.5) & sig_treat, na.rm = TRUE),
                        n_down = sum((res_treat$log2FoldChange <= -0.5) & sig_treat, na.rm = TRUE)),
  interaction = list(res = res_inter, sig = sig_inter,
                     n = sum(sig_inter, na.rm = TRUE),
                     n_up = sum((res_inter$log2FoldChange >= 0.5) & sig_inter, na.rm = TRUE),
                     n_down = sum((res_inter$log2FoldChange <= -0.5) & sig_inter, na.rm = TRUE))
)

cat("\nContrasts tested (FDR ≤ 0.05, |log2FC| ≥ 0.5):\n")
cat(sprintf("  Genotype (5HT2A vs WT):        %3d DEGs (%3d up, %3d down)\n",
            age20_results$main_genotype$n, age20_results$main_genotype$n_up, age20_results$main_genotype$n_down))
cat(sprintf("  Treatment (Yeast vs NoOdor):   %3d DEGs (%3d up, %3d down)\n",
            age20_results$main_treatment$n, age20_results$main_treatment$n_up, age20_results$main_treatment$n_down))
cat(sprintf("  Interaction (genotype×treat):  %3d DEGs (%3d up, %3d down)\n",
            age20_results$interaction$n, age20_results$interaction$n_up, age20_results$interaction$n_down))

# ============================================================================
# 7. CHRONIC FIGURES
# ============================================================================

pca_coords <- data.frame(
  PC1 = pca_age20$x[, 1],
  PC2 = pca_age20$x[, 2],
  PC3 = pca_age20$x[, 3],
  genotype = coldata_age20$genotype,
  treatment = coldata_age20$treatment
)

# Shape by genotype: circle (WT) vs triangle (5HT2A)
pchs <- ifelse(pca_coords$genotype == "WT", 21, 24)

# Fill color by genotype + treatment
# WT: light grey (No Odor) vs dark grey (Yeast)
# 5HT2A: light blue (No Odor) vs dark blue (Yeast)
cols <- rep(NA, nrow(pca_coords))
cols[pca_coords$genotype == "WT" & pca_coords$treatment == "NoOdor"] <- "grey89"       # light grey
cols[pca_coords$genotype == "WT" & pca_coords$treatment == "Yeast"] <- "#4A4A4A"        # dark grey
cols[pca_coords$genotype == "5HT2A" & pca_coords$treatment == "NoOdor"] <- "#22BDD2"    # light blue
cols[pca_coords$genotype == "5HT2A" & pca_coords$treatment == "Yeast"] <- "#333a8c"     # dark blue

border_cols <- "#1a1a1a"  # Dark border color for all
cex_size <- 2.2

plot(pca_coords$PC1, pca_coords$PC3, pch = pchs, cex = cex_size, col = border_cols, bg = cols,
     main = "Middle-aged flies (Head only)",
     xlab = sprintf("PC1 (%.1f%%)", var_exp_age20[1]),
     ylab = sprintf("PC2 (%.1f%%)", var_exp_age20[2]),
     ylim = c(-12, 12), xlim = c(-20, 20),
     las = 1, cex.lab = 1.5, cex.main = 1.5)
abline(h = 0, v = 0, lty = 2, col = "gray70")

legend("bottomleft", 
       legend = c("Wildtype / No Odor",
                  "Wildtype / Yeast Odor",
                  bquote(italic("5-HT2A") ~ "/ No Odor"),
                  bquote(italic("5-HT2A") ~ "/ Yeast Odor")),
       pch = c(21, 21, 24, 24),
       pt.bg = c("grey89", "#4A4A4A", "#22BDD2", "#333a8c"),
       col = "gray30",
       cex = 0.8, bty = "n", pt.cex = 1.5)

# Volcano plot
res <- age20_results$main_genotype$res
sig_mask <- age20_results$main_genotype$sig

x_not_sig <- res[!sig_mask, "log2FoldChange"]
y_not_sig <- -log10(res[!sig_mask, "padj"] + 1e-300)

up_mask <- sig_mask & (res$log2FoldChange >= 0.5)
down_mask <- sig_mask & (res$log2FoldChange <= -0.5)

plot(x_not_sig, y_not_sig, pch = 20, col = "gray60", cex = 0.4,
     xlim = c(-8, 8), ylim = c(0, 80),
     main = sprintf("Head of middle-aged flies: 5-HT2A vs Wildtype"),
     xlab = "log2(fold change)", ylab = "-log10(adjusted p-value)",
     las = 1, cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2)

points(res[up_mask, "log2FoldChange"], -log10(res[up_mask, "padj"] + 1e-300),
       pch = 20, col = "#C93735", cex = 0.5)
points(res[down_mask, "log2FoldChange"], -log10(res[down_mask, "padj"] + 1e-300),
       pch = 20, col = "#22BDD2", cex = 0.5)

abline(h = -log10(0.05), col = "gray11", lty = 2, lwd = 1.5)
abline(v = c(-0.5, 0.5), col = "gray11", lty = 2, lwd = 1.5)

# Add cutoff text annotations
text(-5, -log10(0.05) + 3, "FDR = 0.05", cex = 1.2, col = "black", font = 2)

legend("topleft",
       legend = c(sprintf("Up (n=%d)", age20_results$main_genotype$n_up),
                  sprintf("Down (n=%d)", age20_results$main_genotype$n_down), "NS"),
       col = c("#C93735", "#22BDD2", "gray60"), pch = 20, cex = 1.2, bty = "n")

# ============================================================================
# 8. GO ENRICHMENT ANALYSIS
# ============================================================================

run_directional_go <- function(res, tissue_name, timepoint = "Acute") {
  
  cat("\n", paste(rep("-", 70), collapse=""), "\n", sep="")
  cat(sprintf("%s - %s\n", timepoint, tissue_name))
  cat(paste(rep("-", 70), collapse=""), "\n\n", sep="")
  
  # Extract UP and DOWN regulated genes (FDR <= 0.05, |log2FC| >= 0.5)
  up_genes <- rownames(res)[res$log2FoldChange >= 0.5 & res$padj <= 0.05 & !is.na(res$padj)]
  down_genes <- rownames(res)[res$log2FoldChange <= -0.5 & res$padj <= 0.05 & !is.na(res$padj)]
  
  cat(sprintf("UP-regulated genes: %d\n", length(up_genes)))
  cat(sprintf("DOWN-regulated genes: %d\n\n", length(down_genes)))
  
  results <- list(up = NULL, down = NULL)
  
  # ========== UP-REGULATED GENES ==========
  if (length(up_genes) > 0) {
    cat("Converting UP genes to Entrez IDs...\n")
    up_entrez <- bitr(up_genes, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
    cat(sprintf("  Mapped %d/%d genes successfully\n\n", nrow(up_entrez), length(up_genes)))
    
    # Biological Process (BP) enrichment for UP genes
    cat("Running GO enrichment (Biological Process) for UP genes:\n")
    go_up_bp <- enrichGO(gene = up_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    
    if (!is.null(go_up_bp) && nrow(go_up_bp) > 0) {
      cat(sprintf("  ✓ %d BP terms enriched\n", nrow(go_up_bp)))
      cat("  Top 5 terms:\n")
      for (i in 1:min(5, nrow(go_up_bp))) {
        genes_in_term <- as.numeric(strsplit(go_up_bp$GeneRatio[i], "/")[[1]][1])
        cat(sprintf("    %d. %s\n", i, go_up_bp$Description[i]))
        cat(sprintf("       p.adjust = %.2e, genes = %d/%s\n", 
                    go_up_bp$p.adjust[i], genes_in_term, go_up_bp$GeneRatio[i]))
      }
    } else {
      cat("  ✗ No significant BP terms enriched\n")
    }
    
    # Molecular Function (MF) enrichment for UP genes
    cat("\nRunning GO enrichment (Molecular Function) for UP genes:\n")
    go_up_mf <- enrichGO(gene = up_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                         ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    if (!is.null(go_up_mf) && nrow(go_up_mf) > 0) {
      cat(sprintf("  ✓ %d MF terms enriched\n", nrow(go_up_mf)))
    } else {
      cat("  ✗ No significant MF terms enriched\n")
    }
    
    # Cellular Component (CC) enrichment for UP genes
    cat("\nRunning GO enrichment (Cellular Component) for UP genes:\n")
    go_up_cc <- enrichGO(gene = up_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                         ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    if (!is.null(go_up_cc) && nrow(go_up_cc) > 0) {
      cat(sprintf("  ✓ %d CC terms enriched\n", nrow(go_up_cc)))
    } else {
      cat("  ✗ No significant CC terms enriched\n")
    }
    
    results$up <- list(bp = go_up_bp, mf = go_up_mf, cc = go_up_cc)
  }
  
  # ========== DOWN-REGULATED GENES ==========
  if (length(down_genes) > 0) {
    cat("\n")
    cat("Converting DOWN genes to Entrez IDs...\n")
    down_entrez <- bitr(down_genes, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
    cat(sprintf("  Mapped %d/%d genes successfully\n\n", nrow(down_entrez), length(down_genes)))
    
    # Biological Process (BP) enrichment for DOWN genes
    cat("Running GO enrichment (Biological Process) for DOWN genes:\n")
    go_down_bp <- enrichGO(gene = down_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                           ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    
    if (!is.null(go_down_bp) && nrow(go_down_bp) > 0) {
      cat(sprintf("  ✓ %d BP terms enriched\n", nrow(go_down_bp)))
      cat("  Top 5 terms:\n")
      for (i in 1:min(5, nrow(go_down_bp))) {
        genes_in_term <- as.numeric(strsplit(go_down_bp$GeneRatio[i], "/")[[1]][1])
        cat(sprintf("    %d. %s\n", i, go_down_bp$Description[i]))
        cat(sprintf("       p.adjust = %.2e, genes = %d/%s\n", 
                    go_down_bp$p.adjust[i], genes_in_term, go_down_bp$GeneRatio[i]))
      }
    } else {
      cat("  ✗ No significant BP terms enriched\n")
    }
    
    # Molecular Function (MF) enrichment for DOWN genes
    cat("\nRunning GO enrichment (Molecular Function) for DOWN genes:\n")
    go_down_mf <- enrichGO(gene = down_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                           ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    if (!is.null(go_down_mf) && nrow(go_down_mf) > 0) {
      cat(sprintf("  ✓ %d MF terms enriched\n", nrow(go_down_mf)))
    } else {
      cat("  ✗ No significant MF terms enriched\n")
    }
    
    # Cellular Component (CC) enrichment for DOWN genes
    cat("\nRunning GO enrichment (Cellular Component) for DOWN genes:\n")
    go_down_cc <- enrichGO(gene = down_entrez$ENTREZID, OrgDb = org.Dm.eg.db, 
                           ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1)
    if (!is.null(go_down_cc) && nrow(go_down_cc) > 0) {
      cat(sprintf("  ✓ %d CC terms enriched\n", nrow(go_down_cc)))
    } else {
      cat("  ✗ No significant CC terms enriched\n")
    }
    
    results$down <- list(bp = go_down_bp, mf = go_down_mf, cc = go_down_cc)
  }
  
  return(results)
}

# ============================================================================
# RUN DIRECTIONAL GO FOR ALL TISSUE/TIMEPOINT COMBINATIONS
# ============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n", sep="")
cat("DIRECTIONAL GO ENRICHMENT ANALYSIS\n")
cat("(Comparing UP-regulated vs DOWN-regulated genes)\n")
cat(paste(rep("=", 70), collapse=""), "\n", sep="")

all_directional_go <- list()

# ACUTE HEAD (4 hours)
all_directional_go[["acute_head"]] <- 
  run_directional_go(acute_results[["head"]]$main_genotype$res, "HEAD", "Acute (4 hours)")

# ACUTE OVARIES (4 hours)
all_directional_go[["acute_ovary"]] <- 
  run_directional_go(acute_results[["ovary"]]$main_genotype$res, "OVARIES", "Acute (4 hours)")

# CHRONIC HEAD (Day 20)
all_directional_go[["chronic_head"]] <- 
  run_directional_go(age20_results$main_genotype$res, "HEAD", "Chronic (Day 20)")


p1 <- barplot(all_directional_go[["acute_head"]]$up$bp, showCategory = 10, 
              title = "Young flies / Up-regulated genes")

p2 <- barplot(all_directional_go[["acute_head"]]$down$bp, showCategory = 10, 
              title = "Young flies / Down-regulated genes")

p3 <- barplot(all_directional_go[["chronic_head"]]$up$bp, showCategory = 10, 
              title = "Middle-aged flies / Up-regulated genes")

p4 <- barplot(all_directional_go[["chronic_head"]]$down$bp, showCategory = 10, 
              title = "Middle-aged flies / Down-regulated genes")

grid.arrange(grobs = list(p1, p2, p3, p4), ncol = 2, nrow = 2)


# ============================================================================
# CSV OUTPUT
# ============================================================================

# ACUTE EXPOSURE (Hour 4) =================================================

# UPREGULATED genes
write.csv(as.data.frame(all_directional_go[["acute_head"]]$up$bp), 
          "acute_head_UP_BP.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_head"]]$up$mf), 
          "acute_head_UP_MF.csv", row.names=TRUE)


write.csv(as.data.frame(all_directional_go[["acute_head"]]$up$cc), 
          "acute_head_UP_CC.csv", row.names=TRUE)


write.csv(as.data.frame(all_directional_go[["acute_head"]]$down$bp), 
          "acute_head_DOWN_BP.csv", row.names=TRUE)


write.csv(as.data.frame(all_directional_go[["acute_head"]]$down$mf), 
          "acute_head_DOWN_MF.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_head"]]$down$cc), 
          "acute_head_DOWN_CC.csv", row.names=TRUE)


write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$up$bp), 
          "acute_ovary_UP_BP.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$up$mf), 
          "acute_ovary_UP_MF.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$up$cc), 
          "acute_ovary_UP_CC.csv", row.names=TRUE)


# DOWNREGULATED genes
write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$down$bp), 
          "acute_ovary_DOWN_BP.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$down$mf), 
          "acute_ovary_DOWN_MF.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["acute_ovary"]]$down$cc), 
          "acute_ovary_DOWN_CC.csv", row.names=TRUE)


# CHRONIC HEAD (Day 20) =================================================

# UPREGULATED genes
write.csv(as.data.frame(all_directional_go[["chronic_head"]]$up$bp), 
          "chronic_head_UP_BP.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["chronic_head"]]$up$mf), 
          "chronic_head_UP_MF.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["chronic_head"]]$up$cc), 
          "chronic_head_UP_CC.csv", row.names=TRUE)


# DOWNREGULATED genes
write.csv(as.data.frame(all_directional_go[["chronic_head"]]$down$bp), 
          "chronic_head_DOWN_BP.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["chronic_head"]]$down$mf), 
          "chronic_head_DOWN_MF.csv", row.names=TRUE)

write.csv(as.data.frame(all_directional_go[["chronic_head"]]$down$cc), 
          "chronic_head_DOWN_CC.csv", row.names=TRUE)