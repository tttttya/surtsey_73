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

txi <- tximport(files, type="kallisto", txOut=TRUE)
count_matrix <- txi$counts

meta_full <- read.table("E:/metadata.tsv", sep="\t", header=TRUE, row.names=2, check.names=FALSE)
meta <- meta_full[meta_full$run == 73, ]

stopifnot(all(colnames(count_matrix) == rownames(meta)))

count_matrix_int <- round(count_matrix)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_int,
                              colData   = meta,
                              design    = ~ time + treatment)

dds <- dds[rowSums(counts(dds)) > 10, ]
cat(nrow(dds))

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
res_mix_zero <- results(dds, contrast=c("treatment", "mix", "zero"))
res_mix_zero <- res_mix_zero[order(res_mix_zero$padj), ]

res_ilo_zero <- results(dds, contrast=c("treatment", "ilo_only", "zero"))
res_ilo_zero <- res_ilo_zero[order(res_ilo_zero$padj), ]

res_mix_mixplusilo <- results(dds, contrast=c("treatment", "mix", "mix_plus_ilo"))
res_mix_mixplusilo <- res_mix_mixplusilo[order(res_mix_mixplusilo$padj), ]

deg_mix_zero <- rownames(res_mix_zero[!is.na(res_mix_zero$padj) &
                                        res_mix_zero$padj < 0.05 &
                                        abs(res_mix_zero$log2FoldChange) > 1, ])

deg_ilo_zero <- rownames(res_ilo_zero[!is.na(res_ilo_zero$padj) &
                                        res_ilo_zero$padj < 0.05 &
                                        abs(res_ilo_zero$log2FoldChange) > 1, ])

deg_mix_mixplusilo <- rownames(res_mix_mixplusilo[!is.na(res_mix_mixplusilo$padj) &
                                                    res_mix_mixplusilo$padj < 0.05 &
                                                    abs(res_mix_mixplusilo$log2FoldChange) > 1, ])

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

deg_all <- unique(c(deg_mix_zero, deg_ilo_zero, deg_mix_mixplusilo))
cat("Total DEGs:", length(deg_all), "\n")

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[deg_all, ]
mat_z <- t(scale(t(mat)))
mat_z <- mat_z[complete.cases(mat_z), ]

annotation_col <- meta[, c("time", "treatment"), drop = FALSE]
annotation_col$time <- as.factor(annotation_col$time)
annotation_col$treatment <- as.factor(annotation_col$treatment)

gene_origin <- data.frame(
  gene = rownames(mat_z),
  origin = case_when(
    rownames(mat_z) %in% deg_mix_zero ~ "mix_zero",
    rownames(mat_z) %in% deg_ilo_zero ~ "ilo_zero",
    rownames(mat_z) %in% deg_mix_mixplusilo ~ "mix_mixplusilo",
    TRUE ~ "other"
  )
)
rownames(gene_origin) <- gene_origin$gene
gene_origin <- gene_origin["origin", drop = FALSE]

pheatmap(mat_z,
         annotation_col = annotation_col,
         annotation_row = gene_origin,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         fontsize_row = 6)

