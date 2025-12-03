if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


bioc_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "biomaRt",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "survival",
  "survminer",
  "ComplexHeatmap"
)

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    BiocManager::install(pkg, update = FALSE)
  }
}

# CRAN packages
cran_packages <- c(
  "glmnet",
  "dplyr",
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "corrplot",
  "RColorBrewer",
  "reshape2"
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg)
  }
}

cat("\n✓ All packages installed successfully\n")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/data", showWarnings = FALSE)

cat("✓ Output directories created:\n")
cat("  - results/figures/\n")
cat("  - results/tables/\n")
cat("  - results/data/\n")

library(TCGAbiolinks)
library(SummarizedExperiment)


download_tcga <- function(cancer_type) {
  
  cat(paste("\n>>> Downloading", cancer_type, "...\n"))
  
  query <- GDCquery(
    project = paste0("TCGA-", cancer_type),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  
  
  filename <- paste0("results/data/tcga_", tolower(cancer_type), "_expression.rds")
  saveRDS(data, filename)
  
  
  n_tumor <- sum(colData(data)$sample_type == "Primary Tumor")
  n_normal <- sum(colData(data)$sample_type == "Solid Tissue Normal")
  
  cat(paste("  ✓", cancer_type, "downloaded:\n"))
  cat(paste("    - Tumor samples:", n_tumor, "\n"))
  cat(paste("    - Normal samples:", n_normal, "\n"))
  cat(paste("    - Total:", n_tumor + n_normal, "\n"))
  cat(paste("    - Saved to:", filename, "\n"))
  
  return(data)
}


data_brca <- download_tcga("BRCA")
data_prad <- download_tcga("PRAD")


cat("\nSummary:\n")
cat(paste("  BRCA samples:", ncol(data_brca), "\n"))
cat(paste("  PRAD samples:", ncol(data_prad), "\n"))
cat(paste("  Total samples:", ncol(data_brca) + ncol(data_prad), "\n"))
cat("\nNext: Run 02_qc_and_filtering.R\n\n")

library(SummarizedExperiment)
library(corrplot)
library(pheatmap)
library(RColorBrewer)


data_brca <- readRDS("results/data/tcga_brca_expression.rds")
data_prad <- readRDS("results/data/tcga_prad_expression.rds")



cat("\n>>> Step 1: RMA Quality Control\n")

generate_qc_plots <- function(data, cancer_name) {
  
  cat(paste("\nGenerating QC for", cancer_name, "...\n"))
  
  
  tpm_mat <- assay(data, "tpm_unstrand")
  
  
  set.seed(42)
  subset_idx <- sample(1:nrow(tpm_mat), min(5000, nrow(tpm_mat)))
  tpm_subset <- tpm_mat[subset_idx, ]
  
  sample_corr <- cor(t(log2(tpm_subset + 1)))
  
  
  png(paste0("results/figures/Figure1_", cancer_name, "_correlation_heatmap.png"),
      width = 1000, height = 900, res = 100)
  
  corrplot(sample_corr, 
           method = "color",
           tl.cex = 0.3,
           title = paste("Array-Array Correlation after RMA -", cancer_name),
           mar = c(0, 0, 2, 0))
  
  dev.off()
  
  
  corr_values <- sample_corr[upper.tri(sample_corr)]
  
  png(paste0("results/figures/Figure1_", cancer_name, "_correlation_boxplot.png"),
      width = 800, height = 600, res = 100)
  
  boxplot(corr_values,
          main = paste("Sample Correlation Distribution -", cancer_name),
          ylab = "Correlation Coefficient",
          col = "lightblue",
          ylim = c(0.5, 1))
  
  stripchart(corr_values, method = "jitter", add = TRUE, 
             pch = 16, col = rgb(0, 0, 0, 0.3))
  
  dev.off()
  
  cat(paste("  ✓ Correlation stats for", cancer_name, ":\n"))
  cat(paste("    Mean:", round(mean(corr_values), 3), "\n"))
  cat(paste("    Median:", round(median(corr_values), 3), "\n"))
  
  return(data.frame(
    cancer = cancer_name,
    mean_corr = mean(corr_values),
    median_corr = median(corr_values),
    n_samples = ncol(tpm_mat)
  ))
}


qc_brca <- generate_qc_plots(data_brca, "BRCA")
qc_prad <- generate_qc_plots(data_prad, "PRAD")

qc_summary <- rbind(qc_brca, qc_prad)
write.csv(qc_summary, "results/tables/Table1_QC_summary.csv", row.names = FALSE)

cat("\n✓ QC plots generated\n")



cat("\n>>> Step 2: TPM Filtering\n")

filter_by_tpm <- function(data, cancer_name, tpm_min = 5, tpm_max = 1500) {
  
  cat(paste("\nFiltering", cancer_name, "...\n"))
  
  tpm_mat <- assay(data, "tpm_unstrand")
  mean_tpm <- rowMeans(tpm_mat)
  
  genes_keep <- mean_tpm > tpm_min & mean_tpm < tpm_max
  
  cat(paste("  Before filter:", nrow(tpm_mat), "genes\n"))
  cat(paste("  After filter (5 < TPM < 1500):", sum(genes_keep), "genes\n"))
  
  tpm_filtered <- tpm_mat[genes_keep, ]
  
  return(list(
    tpm = tpm_filtered,
    colData = colData(data)
  ))
}


filtered_brca <- filter_by_tpm(data_brca, "BRCA")
filtered_prad <- filter_by_tpm(data_prad, "PRAD")


common_genes <- intersect(rownames(filtered_brca$tpm), 
                          rownames(filtered_prad$tpm))

cat(paste("\n✓ Common genes between BRCA and PRAD:", length(common_genes), "\n"))


final_tpm_brca <- filtered_brca$tpm[common_genes, ]
final_tpm_prad <- filtered_prad$tpm[common_genes, ]


combined_tpm <- cbind(final_tpm_brca, final_tpm_prad)

cat(paste("✓ Combined matrix:", nrow(combined_tpm), "genes ×", 
          ncol(combined_tpm), "samples\n"))

cat("\n>>> Step 3: Metadata Preparation\n")

meta_brca <- as.data.frame(filtered_brca$colData[, c("barcode", "sample_type")])
meta_brca$cancer_type <- "BRCA"

meta_prad <- as.data.frame(filtered_prad$colData[, c("barcode", "sample_type")])
meta_prad$cancer_type <- "PRAD"

combined_meta <- rbind(meta_brca, meta_prad)
rownames(combined_meta) <- colnames(combined_tpm)

combined_meta$label <- ifelse(combined_meta$sample_type == "Primary Tumor", 1, 0)
combined_meta$cancer_binary <- ifelse(combined_meta$cancer_type == "BRCA", 1, 0)

cat("\nSample distribution:\n")
print(table(combined_meta$cancer_type, combined_meta$sample_type))


saveRDS(list(
  tpm = combined_tpm,
  meta = combined_meta,
  common_genes = common_genes
), "results/data/processed_brca_prad.rds")

cat("\n✓ Processed data saved\n")



library(biomaRt)
library(dplyr)

data_list <- readRDS("results/data/processed_brca_prad.rds")
combined_tpm <- data_list$tpm
combined_meta <- data_list$meta

cat("\n>>> Step 1: Gene Symbol Conversion\n")

ensembl_ids <- gsub("\\.\\d+$", "", rownames(combined_tpm))

cat("Connecting to biomaRt...\n")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

cat("Retrieving gene symbols...\n")
conv <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

conv <- conv[conv$hgnc_symbol != "", ]
conv <- conv[!duplicated(conv$hgnc_symbol), ]

# Map
idx <- match(ensembl_ids, conv$ensembl_gene_id)
genes_keep <- !is.na(idx)

expr_symbol <- combined_tpm[genes_keep, ]
rownames(expr_symbol) <- conv$hgnc_symbol[idx[genes_keep]]

cat(paste("✓ Genes after symbol conversion:", nrow(expr_symbol), "\n"))

cat("\n>>> Step 2: Log2 Transformation\n")

X_log <- log2(expr_symbol + 1)
X <- t(X_log)  # Transpose: samples × genes

cat(paste("✓ Matrix dimensions:", nrow(X), "samples ×", ncol(X), "genes\n"))


cat("\n>>> Step 3: Variance Filtering\n")

gene_vars <- apply(X, 2, var)
nzv <- gene_vars > 0

X_filtered <- X[, nzv]

cat(paste("  Before:", ncol(X), "genes\n"))
cat(paste("  After:", ncol(X_filtered), "genes\n"))


cat("\n>>> Step 4: PCA Analysis\n")

pca_result <- prcomp(X_filtered, scale. = TRUE)

var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
cum_var <- cumsum(var_explained)

cat(paste("  PC1 explains:", round(var_explained[1], 2), "%\n"))
cat(paste("  PC2 explains:", round(var_explained[2], 2), "%\n"))
cat(paste("  PCs for 90% variance:", which(cum_var >= 90)[1], "\n"))

# Plot 1: Variance explained
png("results/figures/Figure2_PCA_variance_explained.png", 
    width = 1200, height = 600, res = 100)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

plot(var_explained[1:50], type = "o", pch = 16,
     xlab = "Principal Component", ylab = "Variance Explained (%)",
     main = "Scree Plot", cex.main = 1.5)

plot(cum_var[1:50], type = "o", pch = 16,
     xlab = "Principal Component", ylab = "Cumulative Variance (%)",
     main = "Cumulative Variance", cex.main = 1.5)
abline(h = 90, col = "red", lty = 2)

dev.off()

# Plot 2: PCA biplot
png("results/figures/Figure3_PCA_biplot.png", 
    width = 1000, height = 800, res = 100)

colors <- c(BRCA = "red", PRAD = "blue")
col_vec <- colors[combined_meta$cancer_type]
pch_vec <- ifelse(combined_meta$label == 1, 16, 17)

plot(pca_result$x[, 1], pca_result$x[, 2],
     main = "PCA: BRCA vs PRAD",
     xlab = paste("PC1 (", round(var_explained[1], 1), "%)", sep = ""),
     ylab = paste("PC2 (", round(var_explained[2], 1), "%)", sep = ""),
     col = col_vec,
     pch = pch_vec,
     cex = 1.5)

legend("topright",
       legend = c("BRCA", "PRAD", "Tumor", "Normal"),
       col = c("red", "blue", "black", "black"),
       pch = c(16, 16, 16, 17),
       cex = 1.2)

dev.off()

cat("\n✓ PCA plots generated\n")


saveRDS(list(
  X = X_filtered,
  y_tumor = combined_meta$label,
  y_cancer = combined_meta$cancer_binary,
  meta = combined_meta,
  pca = pca_result,
  gene_symbols = colnames(X_filtered)
), "results/data/preprocessed_for_lasso.rds")

cat("\n✓ Preprocessed data saved\n")


cat("\nFiles generated:\n")
cat("  - results/figures/Figure2_PCA_variance_explained.png\n")
cat("  - results/figures/Figure3_PCA_biplot.png\n")
cat("  - results/data/preprocessed_for_lasso.rds\n")
cat("\nNext: Run 04_lasso_feature_selection.R\n\n")

library(glmnet)
library(dplyr)


data_list <- readRDS("results/data/preprocessed_for_lasso.rds")
X <- data_list$X
y_tumor <- data_list$y_tumor
y_cancer <- data_list$y_cancer
meta <- data_list$meta



cat("\n>>> Analysis 1: Tumor vs Normal Classification\n")

set.seed(42)

cvfit_tumor <- cv.glmnet(
  x = X,
  y = y_tumor,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  standardize = TRUE
)


png("results/figures/Figure4_LASSO_CV_tumor_vs_normal.png",
    width = 800, height = 600, res = 100)

plot(cvfit_tumor, main = "LASSO CV: Tumor vs Normal")

dev.off()


coef_tumor <- as.matrix(coef(cvfit_tumor, s = "lambda.1se"))
sel_idx <- which(coef_tumor[, 1] != 0)
sel_genes_tumor <- rownames(coef_tumor)[sel_idx]
sel_genes_tumor <- sel_genes_tumor[sel_genes_tumor != "(Intercept)"]

cat(paste("\n✓ Tumor vs Normal LASSO:\n"))
cat(paste("  Lambda (1se):", round(cvfit_tumor$lambda.1se, 6), "\n"))
cat(paste("  Selected genes:", length(sel_genes_tumor), "\n"))


top_genes_tumor <- sel_genes_tumor[order(abs(coef_tumor[sel_genes_tumor, 1]), decreasing = TRUE)][1:20]

cat("\nTop 20 genes for Tumor vs Normal:\n")
print(data.frame(
  gene = top_genes_tumor,
  coefficient = round(coef_tumor[top_genes_tumor, 1], 4)
))


cat("\n>>> Analysis 2: BRCA vs PRAD Classification\n")

set.seed(42)

cvfit_cancer <- cv.glmnet(
  x = X,
  y = y_cancer,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  standardize = TRUE
)


png("results/figures/Figure5_LASSO_CV_brca_vs_prad.png",
    width = 800, height = 600, res = 100)

plot(cvfit_cancer, main = "LASSO CV: BRCA vs PRAD")

dev.off()


coef_cancer <- as.matrix(coef(cvfit_cancer, s = "lambda.1se"))
sel_idx2 <- which(coef_cancer[, 1] != 0)
sel_genes_cancer <- rownames(coef_cancer)[sel_idx2]
sel_genes_cancer <- sel_genes_cancer[sel_genes_cancer != "(Intercept)"]

cat(paste("\n✓ BRCA vs PRAD LASSO:\n"))
cat(paste("  Lambda (1se):", round(cvfit_cancer$lambda.1se, 6), "\n"))
cat(paste("  Selected genes:", length(sel_genes_cancer), "\n"))

top_genes_cancer <- sel_genes_cancer[order(abs(coef_cancer[sel_genes_cancer, 1]), decreasing = TRUE)][1:20]

cat("\nTop 20 genes for BRCA vs PRAD:\n")
print(data.frame(
  gene = top_genes_cancer,
  coefficient = round(coef_cancer[top_genes_cancer, 1], 4)
))



cat("\n>>> Exporting for Python\n")

X_tumor_selected <- X[, sel_genes_tumor]
export_tumor <- as.data.frame(X_tumor_selected)
export_tumor$label <- y_tumor
export_tumor$cancerType <- meta$cancer_type
export_tumor$sampleID <- rownames(meta)

write.csv(export_tumor, "results/data/expression_tumor_vs_normal.csv", row.names = FALSE)

X_cancer_selected <- X[, sel_genes_cancer]
export_cancer <- as.data.frame(X_cancer_selected)
export_cancer$label <- y_cancer
export_cancer$tumorStatus <- ifelse(meta$label == 1, "Tumor", "Normal")
export_cancer$sampleID <- rownames(meta)

write.csv(export_cancer, "results/data/expression_brca_vs_prad.csv", row.names = FALSE)

lasso_tumor_genes <- data.frame(
  gene_name = sel_genes_tumor,
  coefficient = coef_tumor[sel_genes_tumor, 1]
) %>% arrange(desc(abs(coefficient)))

write.csv(lasso_tumor_genes, "results/tables/Table2_LASSO_genes_tumor_vs_normal.csv", row.names = FALSE)

lasso_cancer_genes <- data.frame(
  gene_name = sel_genes_cancer,
  coefficient = coef_cancer[sel_genes_cancer, 1]
) %>% arrange(desc(abs(coefficient)))

write.csv(lasso_cancer_genes, "results/tables/Table3_LASSO_genes_brca_vs_prad.csv", row.names = FALSE)

saveRDS(list(
  cvfit_tumor = cvfit_tumor,
  cvfit_cancer = cvfit_cancer,
  sel_genes_tumor = sel_genes_tumor,
  sel_genes_cancer = sel_genes_cancer,
  X = X,
  meta = meta
), "results/data/lasso_results.rds")

cat("\nFiles generated:\n")
cat("  - results/figures/Figure4_LASSO_CV_tumor_vs_normal.png\n")
cat("  - results/figures/Figure5_LASSO_CV_brca_vs_prad.png\n")
cat("  - results/data/expression_tumor_vs_normal.csv\n")
cat("  - results/data/expression_brca_vs_prad.csv\n")
cat("  - results/tables/Table2_LASSO_genes_tumor_vs_normal.csv\n")
cat("  - results/tables/Table3_LASSO_genes_brca_vs_prad.csv\n")
cat("\nNext: Run Python script 05_xgboost_shap.py\n\n")
