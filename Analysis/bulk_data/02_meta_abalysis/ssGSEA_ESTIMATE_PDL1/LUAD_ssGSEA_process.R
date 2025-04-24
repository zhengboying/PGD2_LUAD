setwd("D:/workspace/")
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(corrplot)
library(GSVA)
library(psych)
Marker_dir = dir("Marker")
Marker_dir_sig = Marker_dir[c(grep("_pvalue_0.05_log2FC_", Marker_dir))]
## 去掉两个unknown
Marker_dir_sig = Marker_dir_sig[1:18]
## for循环获得每个文件内top20基因，再做成list
Marker_list = list()
for (file in Marker_dir_sig) {
  cell_name = str_split(file, "_")[[1]][1]
  Marker_df = read.table(paste0("./Marker/", file), header=T, sep=",", check.names=F)
  Marker_df = na.omit(Marker_df)
  Marker_df = Marker_df[1:20, ]
  Marker_list[cell_name] = Marker_df["GeneName"]
}
###########################################################################################################
## 读取TCGA和GEO的数据
data = readRDS("./dat_4cohort.rds")
## 去掉TCGA有正常样本的数据集
data = data[2:5]
## for循环计算每个数据集对细胞marker基因的ssGSEA富集得分
dir.create("ssGSEA_score")
ssGSEA_list = list()
for (i in 1:length(data)) {
  ssGSEA_score = as.data.frame(t(gsva(expr=as.matrix(data[i][[1]]), gset.idx.list=Marker_list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)))
  write.table(ssGSEA_score, paste0("./ssGSEA_score/", names(data)[i], "_ssGSEA_score.txt"), sep="\t", quote=F, col.names=T, row.names=T)
  ssGSEA_list[[names(data)[i]]] = ssGSEA_score
}
###########################################################################################################
## 提取表达矩阵内感兴趣基因的表达值
gene = c("SLCO2A1", "PTGDR", "PTGDS", "HPGDS")
dir.create("select_gene_exp")
data_gene_list = list()
for (i in 1:length(data)) {
  sub_df = as.data.frame(t(data[i][[1]][gene, ]))
  write.table(sub_df, paste0("./select_gene_exp/", names(data)[i], "_exp.txt"), sep="\t", quote=F, col.names=T, row.names=T)
  data_gene_list[[names(data)[i]]] = sub_df
}
###########################################################################################################
## 计算感兴趣基因表达值和ssGSEA富集得分的相关性
cell_order = c("PTGDR+ NK1", "PTGDR+ NK2", "PTGDR- NK", "DNAJB1+PTGDR+ NK", "KLRC2+PTGDR+ NK", "PTGDS+PTGDR+ NK", "CD4+ T", "CD8+ T1", "CD8+ T2", "CD8+ T3", "Endo", "Fibro", "Mast", "B", "Myeloid", "pDC", "Epi", "Oligo")
mycol <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)
dir.create("ssGSEA_spearman_heatmap")
for (i in 1:length(ssGSEA_list)) {
  spearman_correlation = corr.test(ssGSEA_list[i][[1]], data_gene_list[i][[1]], method="spearman")
  pdf(file=paste0("./ssGSEA_spearman_heatmap/", names(ssGSEA_list)[i], "_spearman_heatmap.pdf"), width=8, height=8)
  corrplot(
    spearman_correlation$r[cell_order, ],
    method = c('circle'), col = mycol(100),
    outline = 'grey', 
    p.mat = spearman_correlation$p[cell_order, ],
    sig.level = c(.001, .01, .05),
    insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
    pch.cex = 1.2, #显著性标记大小
    pch.col = 'black' #显著性标记颜色
  )
  dev.off()
}

