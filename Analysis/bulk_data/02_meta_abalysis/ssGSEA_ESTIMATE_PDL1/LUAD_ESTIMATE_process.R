setwd("D:/workspace/")
library(dplyr)
library(ggplot2)
library(stringr)
library(corrplot)
library(psych)
library(estimate)
## 读取TCGA和GEO的数据
data = readRDS("./dat_4cohort.rds")
## 去掉TCGA有正常样本的数据集
data = data[2:5]
## 计算ESTIMATE得分
dir.create("ESTIMATE_score")
write.table(data$GSE30219, "./ESTIMATE_score/GSE30219_raw.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(data$GSE31210, "./ESTIMATE_score/GSE31210_raw.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(data$GSE37745, "./ESTIMATE_score/GSE37745_raw.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(data$TCGA_T, "./ESTIMATE_score/TCGA_T_raw.txt", col.names=T, row.names=T, quote=F, sep="\t")
#############################################################################
input_file_dir = "./ESTIMATE_score/GSE30219_raw.txt"
output_file_dir = "./ESTIMATE_score/GSE30219_raw.gct"
output_estimate = "./ESTIMATE_score/GSE30219_ESTIMATE.gct"
filterCommonGenes(input.f=input_file_dir, output.f=output_file_dir, id="GeneSymbol")
estimateScore(input.ds=output_file_dir, output.ds=output_estimate, platform="affymetrix")
#############################################################################
input_file_dir = "./ESTIMATE_score/GSE31210_raw.txt"
output_file_dir = "./ESTIMATE_score/GSE31210_raw.gct"
output_estimate = "./ESTIMATE_score/GSE31210_ESTIMATE.gct"
filterCommonGenes(input.f=input_file_dir, output.f=output_file_dir, id="GeneSymbol")
estimateScore(input.ds=output_file_dir, output.ds=output_estimate, platform="affymetrix")
#############################################################################
input_file_dir = "./ESTIMATE_score/GSE37745_raw.txt"
output_file_dir = "./ESTIMATE_score/GSE37745_raw.gct"
output_estimate = "./ESTIMATE_score/GSE37745_ESTIMATE.gct"
filterCommonGenes(input.f=input_file_dir, output.f=output_file_dir, id="GeneSymbol")
estimateScore(input.ds=output_file_dir, output.ds=output_estimate, platform="affymetrix")
#############################################################################
input_file_dir = "./ESTIMATE_score/TCGA_T_raw.txt"
output_file_dir = "./ESTIMATE_score/TCGA_T_raw.gct"
output_estimate = "./ESTIMATE_score/TCGA_T_ESTIMATE.gct"
filterCommonGenes(input.f=input_file_dir, output.f=output_file_dir, id="GeneSymbol")
estimateScore(input.ds=output_file_dir, output.ds=output_estimate, platform="affymetrix")
#############################################################################
## 计算spearman相关性，画热图
ESTIMATE_dir = dir("ESTIMATE_score")
ESTIMATE_files = ESTIMATE_dir[c(grep("_ESTIMATE\\.gct", ESTIMATE_dir))]
Exp_files = dir("select_gene_exp")

mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
dir.create("ESTIMATE_spearman_heatmap")
for (i in 1:length(Exp_files)) {
  Exp_df = read.table(paste0("./select_gene_exp/", names(data)[i], "_exp.txt"), header=T, row.names=1, check.names=F, sep="\t")
  ESTIMATE_df = read.table(paste0("./ESTIMATE_score/", names(data)[i], "_ESTIMATE.gct"), header=T, row.names=1, check.names=F, sep="\t")
  ESTIMATE_df = as.data.frame(t(ESTIMATE_df))
  
  spearman_correlation = corr.test(ESTIMATE_df, Exp_df, method="spearman")
  pdf(file=paste0("./ESTIMATE_spearman_heatmap/", names(data)[i], "_spearman_heatmap.pdf"), width=5, height=5)
  corrplot(
    spearman_correlation$r,
    method = c('circle'), col = mycol(100),
    outline = 'grey', 
    p.mat = spearman_correlation$p,
    sig.level = c(.001, .01, .05),
    insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
    pch.cex = 1.2, #显著性标记大小
    pch.col = 'black' #显著性标记颜色
  )
  dev.off()
}







