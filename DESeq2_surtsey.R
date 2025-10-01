if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("readr")
BiocManager::install("apeglm")      
BiocManager::install("pheatmap")    

library(DESeq2)
library(tximport)
library(readr)
library(ggplot2)
library(pheatmap)

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

meta_full <- read.table("E:/metadata.tsv", 
                   sep="\t", header=TRUE, 
                   row.names=2,  
                   check.names=FALSE)
meta <- meta_full[meta_full$run == 73, ]
all(colnames(count_matrix) == rownames(meta))

count_matrix <- txi$counts
count_matrix_int <- round(count_matrix)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_int,
                              colData   = meta,
                              design    = ~ time + treatment)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

#results

res <- results(dds, contrast=c("treatment", "mix", "zero"))
res <- res[order(res$padj), ]  
write.csv(as.data.frame(res), file="DESeq2_results.csv")

plotMA(res, ylim=c(-5,5))

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

res$gene <- rownames(res)
ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), col="red") +
  geom_vline(xintercept=c(-1,1), col="blue") +
  theme_minimal()

topgenes <- head(order(res$padj), 20)
mat <- assay(vsd)[topgenes, ]
pheatmap(mat, annotation_col=meta)
