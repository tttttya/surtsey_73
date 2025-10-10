
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximport", "readr", "apeglm", "pheatmap"))

library(DESeq2)
library(tximport)
library(readr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(matrixStats)

setwd("E:/adrian2taisiia_2025.09.26/adrian2taisiia_2025.09.26/kallisto.100")

files <- c(
  "073_32/abundance.tsv",
  "073_33/abundance.tsv",
  "073_34/abundance.tsv",
  "073_35/abundance.tsv",
  "073_36/abundance.tsv",
  "073_37/abundance.tsv",
  "073_38/abundance.tsv",
  "073_39/abundance.tsv",
  "073_40/abundance.tsv",
  "073_41/abundance.tsv",
  "073_73/abundance.tsv",
  "073_74/abundance.tsv",
  "073_75/abundance.tsv",
  "073_76/abundance.tsv",
  "073_77/abundance.tsv",
  "073_78/abundance.tsv",
  "073_79/abundance.tsv",
  "073_80/abundance.tsv",
  "073_81/abundance.tsv",
  "073_82/abundance.tsv",
  "073_114/abundance.tsv",
  "073_115/abundance.tsv",
  "073_116/abundance.tsv",
  "073_117/abundance.tsv",
  "073_118/abundance.tsv",
  "073_119/abundance.tsv",
  "073_120/abundance.tsv",
  "073_121/abundance.tsv",
  "073_122/abundance.tsv",
  "073_123/abundance.tsv"
)

names(files) <- gsub("/abundance.tsv", "", files)

txi <- tximport(files, type = "kallisto", txOut = TRUE)
count_matrix <- txi$counts

meta_full <- read.table("E:/metadata.tsv", sep = "\t", header = TRUE, row.names = 2, check.names = FALSE)

meta <- meta_full[meta_full$run == 73, ]

common_samples <- intersect(colnames(count_matrix), rownames(meta))
count_matrix <- count_matrix[, common_samples]
meta <- meta[common_samples, ]

stopifnot(all(colnames(count_matrix) == rownames(meta)))

count_matrix_int <- round(count_matrix)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData   = meta,
  design    = ~ time + treatment
)

dds <- dds[rowSums(counts(dds)) > 5, ]

seta_indexes <- which(meta$treatment == "zero")
setb_indexes <- which(meta$treatment == "mix")

count_threshold <- 20
tpm_threshold <- 2

a <- counts(dds)[, seta_indexes]
b <- counts(dds)[, setb_indexes]
diff_median <- rowMedians(a) - rowMedians(b)
keep_counts <- abs(diff_median) >= count_threshold
dds <- dds[keep_counts, ]
cat(nrow(dds))

subset <- txi$abundance[names(dds), ]
a_tpm <- rowMedians(subset[, seta_indexes])
b_tpm <- rowMedians(subset[, setb_indexes])
keep_tpm <- pmax(a_tpm, b_tpm) >= tpm_threshold
dds <- dds[keep_tpm, ]
cat(nrow(dds))

dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", "mix", "zero"))
res <- res[order(res$padj), ]

effect_size_threshold <- log2(1.5)
deg_res <- res[!is.na(res$padj) & res$padj < 0.1 & abs(res$log2FoldChange) > effect_size_threshold, ]
cat(nrow(deg_res))

write.csv(as.data.frame(deg_res), "DESeq2_filtered_DEGs.csv")

plotMA(res, ylim = c(-5, 5))
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "treatment")

ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  geom_vline(xintercept = c(-1, 1), col = "blue") +
  theme_minimal()

deg_genes <- rownames(deg_res)
mat <- assay(vsd)[deg_genes, ]
mat_z <- t(scale(t(mat)))
mat_z <- mat_z[complete.cases(mat_z), ]
mat_z[!is.finite(mat_z)] <- 0
mat_z <- mat_z[rowSds(mat_z, na.rm = TRUE) > 0, ]
pheatmap(mat_z, annotation_col = meta, cluster_rows = TRUE, cluster_cols = TRUE)

