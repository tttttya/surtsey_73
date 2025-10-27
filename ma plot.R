
library(ggplot2)

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

# dataframe
res_df <- as.data.frame(res_mix_zero_4h)  # замените на ваш res
res_df$DEG <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

# highlighting DEGs
res_nonDEG <- res_df[!res_df$DEG, ]
res_DEG <- res_df[res_df$DEG, ]

# MA plot
ggplot() +
  geom_point(data = res_nonDEG, aes(x = baseMean, y = log2FoldChange),
             color = "grey70", alpha = 0.5, size = 1.2) +  
  geom_point(data = res_DEG, aes(x = baseMean, y = log2FoldChange),
             color = "red", alpha = 1, size = 1.8) +        
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Mean expression",
       y = "log2 fold change") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
