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

meta <- read.table("E:/metadata.tsv", sep="\t", header=TRUE, row.names=2, check.names=FALSE)
meta <- meta[meta$run == 73, ]
stopifnot(all(colnames(count_matrix) == rownames(meta)))

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData   = meta,
                              design    = ~ time + treatment)
dds <- dds[rowSums(counts(dds)) > 10, ]
cat(nrow(dds))

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


# 4h filtration
meta_4h <- as.data.frame(colData(dds))
meta_4h <- meta_4h[meta_4h$time == 'four', ]
dds_4h <- dds[, rownames(meta_4h)]
dds_4h$time <- droplevels(dds_4h$time)
dds_4h$treatment <- droplevels(dds_4h$treatment)

design(dds_4h) <- ~ treatment
dds_4h <- DESeq(dds_4h)

# contrasts
res_mix_zero_4h <- results(dds_4h, contrast = c("treatment", "mix", "zero"))
res_ilo_zero_4h <- results(dds_4h, contrast = c("treatment", "ilo_only", "zero"))
res_mix_mixplusilo_4h <- results(dds_4h, contrast = c("treatment", "mix", "mix_plus_ilo"))

# DEGs
deg_mix_zero_4h <- rownames(res_mix_zero_4h[!is.na(res_mix_zero_4h$padj) &
                                              res_mix_zero_4h$padj < 0.05 &
                                              abs(res_mix_zero_4h$log2FoldChange) > 1, ])
deg_ilo_zero_4h <- rownames(res_ilo_zero_4h[!is.na(res_ilo_zero_4h$padj) &
                                              res_ilo_zero_4h$padj < 0.05 &
                                              abs(res_ilo_zero_4h$log2FoldChange) > 1, ])
deg_mix_mixplusilo_4h <- rownames(res_mix_mixplusilo_4h[!is.na(res_mix_mixplusilo_4h$padj) &
                                                          res_mix_mixplusilo_4h$padj < 0.05 &
                                                          abs(res_mix_mixplusilo_4h$log2FoldChange) > 1, ])
deg_all_4h <- unique(c(deg_mix_zero_4h, deg_ilo_zero_4h, deg_mix_mixplusilo_4h))

# heatmap matrix
vsd_4h <- vst(dds_4h, blind = FALSE)
mat_4h <- assay(vsd_4h)[deg_all_4h, ]
mat_z_4h <- t(scale(t(mat_4h)))
mat_z_4h <- mat_z_4h[complete.cases(mat_z_4h), ]

# annotations
annotation_col_4h <- as.data.frame(colData(dds_4h)[, "treatment", drop=FALSE])
gene_origin_4h <- data.frame(
  gene = rownames(mat_z_4h),
  origin = case_when(
    rownames(mat_z_4h) %in% deg_mix_zero_4h ~ "mix_zero",
    rownames(mat_z_4h) %in% deg_ilo_zero_4h ~ "ilo_zero",
    rownames(mat_z_4h) %in% deg_mix_mixplusilo_4h ~ "mix_mixplusilo",
    TRUE ~ "other"
  )
)
rownames(gene_origin_4h) <- gene_origin_4h$gene
gene_origin_4h <- gene_origin_4h["origin", drop=FALSE]

ann_colors <- list(
  treatment = c(zero="gray53", 
                mix="firebrick3", 
                ilo_only="green4", 
                mix_plus_ilo="#1F78B4"),
  
  origin = c(mix_zero="darkorange", 
             ilo_zero="palegreen3", 
             mix_mixplusilo="mediumorchid3", 
             other="#D3D3D3")
)

# heatmap
pheatmap(mat_z_4h,
         annotation_col = annotation_col_4h,
         annotation_row = gene_origin_4h,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         clustering_method = "ward.D2",
         fontsize_row = 6,
         border_color = NA)