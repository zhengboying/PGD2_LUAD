setwd("D:/workspace")
library(corrplot)
library(reshape2)
library(dplyr)

gene_order = c("SLCO2A1", "PTGDR", "PTGDS", "HPGDS")
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
files = dir("Hotscore_PDL1_spearman_heatmap")

for (file in files) {
  filename = strsplit(file, "_")[[1]][1]
  df = read.table(paste0("./Hotscore_PDL1_spearman_heatmap/", file), header=T, sep=",", check.names=F)
  sub_df = filter(df, Var2 %in% gene_order)
  sub_df = filter(sub_df, Var1 %in% c("PDL1", "HOT_score"))
  cor_matrix <- dcast(sub_df[, c(1,2,3)], Var2 ~ Var1)
  rownames(cor_matrix) = cor_matrix$Var2
  cor_matrix = cor_matrix[gene_order, 2:3]
  
  p_matrix <- dcast(sub_df[, c(1,2,4)], Var2 ~ Var1)
  rownames(p_matrix) = p_matrix$Var2
  p_matrix = p_matrix[gene_order, 2:3]
  pdf(file=paste0("./Hotscore_PDL1_spearman_heatmap/", filename, "_Hotscore_PDL1_spearman_heatmap.pdf"), width=2.5, height=5)
  corrplot(
    as.matrix(cor_matrix),
    method = c('circle'), col = mycol(100),
    outline = 'grey', 
    p.mat = as.matrix(p_matrix),
    sig.level = c(.001, .01, .05),
    insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
    pch.cex = 1.2, #显著性标记大小
    pch.col = 'black' #显著性标记颜色
  )
  dev.off()
}







