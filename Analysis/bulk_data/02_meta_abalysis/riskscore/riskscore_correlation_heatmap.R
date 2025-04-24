·setwd("D:/workspace")
library(dplyr)
library(psych)
library(ggplot2)
library(corrplot)
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE31210[, "risk_score"], row.names=rownames(riskscore_data$GSE31210))
ESTIMATE_data = read.table("./ESTIMATE_score/GSE31210_ESTIMATE.gct", header=T, row.names=1, sep="\t")
ESTIMATE_data = as.data.frame(t(ESTIMATE_data))
spearman_correlation = corr.test(ESTIMATE_data, riskscore_data, method="spearman")
pdf("./ESTIMATE_score/GSE31210_riskscore_heatmap.pdf", width=4, height=4)
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
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE37745[, "risk_score"], row.names=rownames(riskscore_data$GSE37745))
ESTIMATE_data = read.table("./ESTIMATE_score/GSE37745_ESTIMATE.gct", header=T, row.names=1, sep="\t")
ESTIMATE_data = as.data.frame(t(ESTIMATE_data))
spearman_correlation = corr.test(ESTIMATE_data, riskscore_data, method="spearman")
pdf("./ESTIMATE_score/GSE37745_riskscore_heatmap.pdf", width=4, height=4)
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
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$TCGA[, "risk_score"], row.names=rownames(riskscore_data$TCGA))
ESTIMATE_data = read.table("./ESTIMATE_score/TCGA_T_ESTIMATE.gct", header=T, row.names=1, sep="\t")
ESTIMATE_data = as.data.frame(t(ESTIMATE_data))
spearman_correlation = corr.test(ESTIMATE_data, riskscore_data, method="spearman")
pdf("./ESTIMATE_score/TCGA_riskscore_heatmap.pdf", width=4, height=4)
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
#####################################################################################################
cell_order = c("PTGDR+ NK1", "PTGDR+ NK2", "PTGDR- NK", "DNAJB1+PTGDR+ NK", "KLRC2+PTGDR+ NK", "PTGDS+PTGDR+ NK", "CD4+ T", "CD8+ T1", "CD8+ T2", "CD8+ T3", "Endo", "Fibro", "Mast", "B", "Myeloid", "pDC", "Epi", "Oligo")
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE31210[, "risk_score"], row.names=rownames(riskscore_data$GSE31210))
ssGSEA_data = read.table("./ssGSEA_score/GSE31210_ssGSEA_score.txt", header=T, row.names=1, sep="\t", check.names=F)
spearman_correlation = corr.test(ssGSEA_data, riskscore_data, method="spearman")
r_matrix = as.data.frame(spearman_correlation$r)
p_matrix = as.data.frame(spearman_correlation$p)
r_matrix = data.frame(risk_score=r_matrix[cell_order, ], row.names=cell_order)
p_matrix = data.frame(risk_score=p_matrix[cell_order, ], row.names=cell_order)
pdf("./ssGSEA_score/GSE31210_riskscore_heatmap.pdf", width=12, height=12)
corrplot(
  as.matrix(r_matrix),
  method = c('circle'), col = mycol(100),
  outline = 'grey', 
  p.mat = as.matrix(p_matrix),
  sig.level = c(.001, .01, .05),
  insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
  pch.cex = 1.2, #显著性标记大小
  pch.col = 'black' #显著性标记颜色
)
dev.off()
#####################################################################################################
cell_order = c("PTGDR+ NK1", "PTGDR+ NK2", "PTGDR- NK", "DNAJB1+PTGDR+ NK", "KLRC2+PTGDR+ NK", "PTGDS+PTGDR+ NK", "CD4+ T", "CD8+ T1", "CD8+ T2", "CD8+ T3", "Endo", "Fibro", "Mast", "B", "Myeloid", "pDC", "Epi", "Oligo")
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE37745[, "risk_score"], row.names=rownames(riskscore_data$GSE37745))
ssGSEA_data = read.table("./ssGSEA_score/GSE37745_ssGSEA_score.txt", header=T, row.names=1, sep="\t", check.names=F)
spearman_correlation = corr.test(ssGSEA_data, riskscore_data, method="spearman")
r_matrix = as.data.frame(spearman_correlation$r)
p_matrix = as.data.frame(spearman_correlation$p)
r_matrix = data.frame(risk_score=r_matrix[cell_order, ], row.names=cell_order)
p_matrix = data.frame(risk_score=p_matrix[cell_order, ], row.names=cell_order)
pdf("./ssGSEA_score/GSE37745_riskscore_heatmap.pdf", width=12, height=12)
corrplot(
  as.matrix(r_matrix),
  method = c('circle'), col = mycol(100),
  outline = 'grey', 
  p.mat = as.matrix(p_matrix),
  sig.level = c(.001, .01, .05),
  insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
  pch.cex = 1.2, #显著性标记大小
  pch.col = 'black' #显著性标记颜色
)
dev.off()
#####################################################################################################
cell_order = c("PTGDR+ NK1", "PTGDR+ NK2", "PTGDR- NK", "DNAJB1+PTGDR+ NK", "KLRC2+PTGDR+ NK", "PTGDS+PTGDR+ NK", "CD4+ T", "CD8+ T1", "CD8+ T2", "CD8+ T3", "Endo", "Fibro", "Mast", "B", "Myeloid", "pDC", "Epi", "Oligo")
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$TCGA[, "risk_score"], row.names=rownames(riskscore_data$TCGA))
ssGSEA_data = read.table("./ssGSEA_score/TCGA_T_ssGSEA_score.txt", header=T, row.names=1, sep="\t", check.names=F)
spearman_correlation = corr.test(ssGSEA_data, riskscore_data, method="spearman")
r_matrix = as.data.frame(spearman_correlation$r)
p_matrix = as.data.frame(spearman_correlation$p)
r_matrix = data.frame(risk_score=r_matrix[cell_order, ], row.names=cell_order)
p_matrix = data.frame(risk_score=p_matrix[cell_order, ], row.names=cell_order)
pdf("./ssGSEA_score/TCGA_riskscore_heatmap.pdf", width=12, height=12)
corrplot(
  as.matrix(r_matrix),
  method = c('circle'), col = mycol(100),
  outline = 'grey', 
  p.mat = as.matrix(p_matrix),
  sig.level = c(.001, .01, .05),
  insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
  pch.cex = 1.2, #显著性标记大小
  pch.col = 'black' #显著性标记颜色
)
dev.off()
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE31210[, "risk_score"], row.names=rownames(riskscore_data$GSE31210))
HOTscore_data = read.table("./HOTscore_PDL1/hotcold_GSE31210.tsv", header=T, row.names=1, sep="\t", check.names=F)
HOTscore_data = HOTscore_data[, c("PDL1", "HOT_score")]
spearman_correlation = corr.test(HOTscore_data, riskscore_data, method="spearman")
pdf("./HOTscore_PDL1/GSE31210_riskscore_heatmap.pdf", width=3, height=3)
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
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$GSE37745[, "risk_score"], row.names=rownames(riskscore_data$GSE37745))
HOTscore_data = read.table("./HOTscore_PDL1/hotcold_GSE37745.tsv", header=T, row.names=1, sep="\t", check.names=F)
HOTscore_data = HOTscore_data[, c("PDL1", "HOT_score")]
spearman_correlation = corr.test(HOTscore_data, riskscore_data, method="spearman")
pdf("./HOTscore_PDL1/GSE37745_riskscore_heatmap.pdf", width=3, height=3)
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
#####################################################################################################
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
riskscore_data = readRDS("./model_score.rds")
riskscore_data = data.frame(risk_score=riskscore_data$TCGA[, "risk_score"], row.names=rownames(riskscore_data$TCGA))
HOTscore_data = read.table("./HOTscore_PDL1/hotcold_TCGA.tsv", header=T, row.names=1, sep="\t", check.names=F)
HOTscore_data = HOTscore_data[, c("PDL1", "HOT_score")]
spearman_correlation = corr.test(HOTscore_data, riskscore_data, method="spearman")
pdf("./HOTscore_PDL1/TCGA_riskscore_heatmap.pdf", width=3, height=3)
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

