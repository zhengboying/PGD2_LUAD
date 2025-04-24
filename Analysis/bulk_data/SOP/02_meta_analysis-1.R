# 00.Init -----

{
  rm(list = ls())
  options(stringsAsFactors = FALSE)
  gc()
  options(future.globals.maxSize = 2000 * 1024^2)
  #options(bitmapType = "cairo")
  set.seed(16)
  num.cores = 12
  OUTPUTPDF  = FALSE
}


# 01.Setup -----
  
project_path =  "/Users/zhengboying/Documents/project_demo/Analysis/bulk_data/"
setwd(project_path)


suppressMessages({
  #library(METAFlux)
  library(pROC)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(ggprism)
  library(ggsignif)
  library(GSVA)
  library(BiocParallel)
  library(corrplot)
  library(plotROC)
  
  library(openxlsx)
  library(survival)
  library(survminer)
})

{
  source("../../Function/Prognosis_function.R")
  
  InputPath <- c("./Rds/")
  
}

# Input Load -----
dat<-readRDS(paste0(InputPath,"dat_4cohort.rds"))
combined_data<-readRDS(paste0(InputPath ,"combined_data_4cohort.rds"))
cln<-readRDS(paste0(InputPath,"cln.rds"))





## Step1：expr -----

ResultPath = "./01_plot/"
##画图
## TN 区别
{
  markers<-c("HPGDS","PTGDS","PTGDR","SLCO2A1")
  df_long<-as.data.frame(t(dat$TCGA[markers,]))
  
  df_long <- df_long %>%
    mutate(group = gsub(".*_(.*)_.*", "\\1", rownames(.)))  # 提取下划线之间的字段
  
  df_long<-melt(df_long)
  ##增加一列group，取行名中下划线之间的字段
  
  p1<-ggplot(df_long,#
             aes(x =group, y = value, fill = group)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c("N"= "#187C3A",
                                 "T" = "#FFA500"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("T", "N")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  
  ggsave(p2, file = paste0(ResultPath, "/gene_expr.pdf"),
         bg = "white", width = 8, height = 8)
}


##进展相关图
{
  markers<-c("HPGDS","PTGDS","PTGDR","SLCO2A1")
  df_long<-as.data.frame(t(dat$TCGA_T[markers,]))
  
  df_long <- df_long %>%
    mutate(ID = gsub("^(.*_).*$", "\\1", rownames(.)))  # 提取下划线之间的字段
  ##删除所有“_T_”
  df_long$ID <- gsub("_T_", "", df_long$ID)
  
  ##
  clin = "/Users/zhengboying/Documents/project_demo/Inputdata/bulk_data/TCGA/clin.xlsx"
  clin <- openxlsx::read.xlsx(clin)
  row.names(clin) <- clin$ID
  
  ## bind
  df_long <- merge(df_long, clin, by = "ID")
  df_long<-melt(df_long)
  
  ##stage=NA 过滤
  df_long<-df_long[!is.na(df_long$Stage),]
  ##if stage ==I II III Early Stage ,IV Advanced Stage
  df_long$Stage.1<-ifelse(df_long$Stage %in% c("Stage I","Stage II","Stage III"),"Early Stage","Advanced Stage")
  
  
  p1<-ggplot(df_long,#
             aes(x =Stage.1, y = value, fill = Stage.1)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c(#"Stage I"= "#187C3A",
      "Early Stage" = "orange",
      "Advanced Stage"= "red"
      
      #scale_fill_manual(values = c("Stage I"= "#187C3A",
      #                             "Stage II" = "yellow",
      #                             "Stage III"= "orange",
      #                             "Stage IV" = "red"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = t.test, 
                size = 0.2, 
                comparisons = list(c("Early Stage", "Advanced Stage")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  ##不显著变化
  ggsave(p2, file = paste0(ResultPath, "/gene_expr_stage.1.pdf"),
         bg = "white", width = 8, height = 8)
}

## Step2:HOT COLD 分析 -----
{
  refdat = "../../DATABASE/tumor_hot_cold/hotColdGenes.xlsx"
  geneSets <- openxlsx::read.xlsx(refdat)
  geneSets <- list(HOT_score=unique(geneSets$SYMBOL))
  ResultPath = "./02_meta_abalysis/"
  
}

##根据基因集的 GSVA 分析结果将样本分类为 "Hot" 或 "Cold"
{
  gsvaParam<-gsvaParam(
    as.matrix(dat$TCGA_T),
    geneSets,
    kcdf = "Gaussian"
  )
  res0 <- gsva(gsvaParam,
               verbose = T,
               BPPARAM = MulticoreParam(workers = 4))
  
  res <- as.data.frame(t(res0))
  res$Phenotype <- sapply(res$HOT_score,function(i){
    if (i>0 & i<1){
      z <- "Hot"
    } else if(i>-1 & i<0){
      z <- "Cold"
    } else {
      z <- ""
    }
    return(z)
  })
  ###PD-L1（程序性死亡配体1）是编码由 CD274 基因产生的蛋白质的通用名称
  res$PDL1 <- as.numeric(dat$TCGA_T[row.names(dat$TCGA_T)=="CD274",])
  res$HPGDS <- as.numeric(dat$TCGA_T[row.names(dat$TCGA_T)=="HPGDS",])
  res$PTGDS <- as.numeric(dat$TCGA_T[row.names(dat$TCGA_T)=="PTGDS",])
  res$PTGDR <- as.numeric(dat$TCGA_T[row.names(dat$TCGA_T)=="PTGDR",])
  res$SLCO2A1 <- as.numeric(dat$TCGA_T[row.names(dat$TCGA_T)=="SLCO2A1",])
  
  ##写出
  write.table(res, file = paste0(ResultPath,"hotcold_","TCGA",".tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}
resTCGA<-res


##dat$GSE37745
{
  gsvaParam<-gsvaParam(
    as.matrix(dat$GSE37745),
    geneSets,
    kcdf = "Gaussian"
  )
  res0 <- gsva(gsvaParam,
               verbose = T,
               BPPARAM = MulticoreParam(workers = 4))
  
  res <- as.data.frame(t(res0))
  res$Phenotype <- sapply(res$HOT_score,function(i){
    if (i>0 & i<1){
      z <- "Hot"
    } else if(i>-1 & i<0){
      z <- "Cold"
    } else {
      z <- ""
    }
    return(z)
  })
  ###PD-L1（程序性死亡配体1）是编码由 CD274 基因产生的蛋白质的通用名称
  res$PDL1 <- as.numeric(dat$GSE37745[row.names(dat$GSE37745)=="CD274",])
  res$HPGDS <- as.numeric(dat$GSE37745[row.names(dat$GSE37745)=="HPGDS",])
  res$PTGDS <- as.numeric(dat$GSE37745[row.names(dat$GSE37745)=="PTGDS",])
  res$PTGDR <- as.numeric(dat$GSE37745[row.names(dat$GSE37745)=="PTGDR",])
  res$SLCO2A1 <- as.numeric(dat$GSE37745[row.names(dat$GSE37745)=="SLCO2A1",])
  ##写出
  write.table(res, file = paste0(ResultPath,"hotcold_","GSE37745",".tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
  
  
}


resGSE37745<-res


##dat$GSE31210
{
  gsvaParam<-gsvaParam(
    as.matrix(dat$GSE31210),
    geneSets,
    kcdf = "Gaussian"
  )
  res0 <- gsva(gsvaParam,
               verbose = T,
               BPPARAM = MulticoreParam(workers = 4))
  
  res <- as.data.frame(t(res0))
  res$Phenotype <- sapply(res$HOT_score,function(i){
    if (i>0 & i<1){
      z <- "Hot"
    } else if(i>-1 & i<0){
      z <- "Cold"
    } else {
      z <- ""
    }
    return(z)
  })
  
  ###PD-L1（程序性死亡配体1）是编码由 CD274 基因产生的蛋白质的通用名称
  res$PDL1 <- as.numeric(dat$GSE31210[row.names(dat$GSE31210)=="CD274",])
  res$HPGDS <- as.numeric(dat$GSE31210[row.names(dat$GSE31210)=="HPGDS",])
  res$PTGDS <- as.numeric(dat$GSE31210[row.names(dat$GSE31210)=="PTGDS",])
  res$PTGDR <- as.numeric(dat$GSE31210[row.names(dat$GSE31210)=="PTGDR",])
  res$SLCO2A1 <- as.numeric(dat$GSE31210[row.names(dat$GSE31210)=="SLCO2A1",])
  ##写出
  write.table(res, file = paste0(ResultPath,"hotcold_","GSE31210",".tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}

resGSE31210<-res


###dat$GSE30219
{
  gsvaParam<-gsvaParam(
    as.matrix(dat$GSE30219),
    geneSets,
    kcdf = "Gaussian"
  )
  res0 <- gsva(gsvaParam,
               verbose = T,
               BPPARAM = MulticoreParam(workers = 4))
  
  res <- as.data.frame(t(res0))
  res$Phenotype <- sapply(res$HOT_score,function(i){
    if (i>0 & i<1){
      z <- "Hot"
    } else if(i>-1 & i<0){
      z <- "Cold"
    } else {
      z <- ""
    }
    return(z)
  })
  ###PD-L1（程序性死亡配体1）是编码由 CD274 基因产生的蛋白质的通用名称
  res$PDL1 <- as.numeric(dat$GSE30219[row.names(dat$GSE30219)=="CD274",])
  res$HPGDS <- as.numeric(dat$GSE30219[row.names(dat$GSE30219)=="HPGDS",])
  res$PTGDS <- as.numeric(dat$GSE30219[row.names(dat$GSE30219)=="PTGDS",])
  res$PTGDR <- as.numeric(dat$GSE30219[row.names(dat$GSE30219)=="PTGDR",])
  res$SLCO2A1 <- as.numeric(dat$GSE30219[row.names(dat$GSE30219)=="SLCO2A1",])
  
  
  ##写出
  write.table(res, file = paste0(ResultPath,"hotcold_","GSE30219",".tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}

resGSE30219<-res



## PLOT

{ 
  df_long<-melt(resTCGA)
  
  p1<-ggplot(df_long,#
             aes(x =Phenotype, y = value, fill = Phenotype)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c(#"Stage I"= "#187C3A",
      "Hot" = "orange",
      "Cold"= "blue"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = t.test, 
                size = 0.2, 
                comparisons = list(c("Hot", "Cold")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  ##显著变化
  ggsave(p2, file = paste0(ResultPath, "/gene_expr_cold_hot.pdf"),
         bg = "white", width = 8, height = 8)
  
}

{ 
  df_long<-melt(resGSE30219)
  
  p1<-ggplot(df_long,#
             aes(x =Phenotype, y = value, fill = Phenotype)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c(#"Stage I"= "#187C3A",
      "Hot" = "orange",
      "Cold"= "blue"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = t.test, 
                size = 0.2, 
                comparisons = list(c("Hot", "Cold")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  ##显著变化
  ggsave(p2, file = paste0(ResultPath, "/gene_expr_cold_hot_resGSE30219.pdf"),
         bg = "white", width = 8, height = 8)
  
}


{ 
  df_long<-melt(resGSE31210)
  
  p1<-ggplot(df_long,#
             aes(x =Phenotype, y = value, fill = Phenotype)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c(#"Stage I"= "#187C3A",
      "Hot" = "orange",
      "Cold"= "blue"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = t.test, 
                size = 0.2, 
                comparisons = list(c("Hot", "Cold")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  ##显著变化
  ggsave(p2, file = paste0(ResultPath, "/gene_expr_cold_hot_resGSE31210.pdf"),
         bg = "white", width = 8, height = 8)
  
}

{ 
  df_long<-melt(resGSE37745)
  
  p1<-ggplot(df_long,#
             aes(x =Phenotype, y = value, fill = Phenotype)) +
    geom_violin(trim = F, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    #scale_x_discrete(limits = order) +
    theme_prism() +
    theme(axis.text.x = element_text( hjust = 1)) +
    #xlab("Stage") +
    ylab("expr") +
    scale_fill_manual(values = c(#"Stage I"= "#187C3A",
      "Hot" = "orange",
      "Cold"= "blue"
    ))+
    #spit wrap celltype,一行一个
    facet_wrap(~variable,scales = "free_y")+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p1
  
  ##添加显著性
  p2<-p1 +
    geom_signif(test = t.test, 
                size = 0.2, 
                comparisons = list(c("Hot", "Cold")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p2
  ##显著变化
  ggsave(p2, file = paste0(ResultPath, "/gene_expr_cold_hot_resGSE37745.pdf"),
         bg = "white", width = 8, height = 8)
  
}


##cor analysis

source("../../Function/diff_test.R")
####Calculate the Spearman correlation between the HOT_score and the expression of each gene
{
  data <- resTCGA[, c("HOT_score", "PDL1", "HPGDS", "PTGDS", "PTGDR", "SLCO2A1")]
  
  cor_test(data =data, 
           name="TCGA",
           method = "spearman")
  
  #resGSE37745
  data <- resGSE37745[, c("HOT_score", "PDL1", "HPGDS", "PTGDS", "PTGDR", "SLCO2A1")]
  
  cor_test(data =data, 
           name="GSE37745",
           method = "spearman")
  
  #resGSE30219
  data <- resGSE30219[, c("HOT_score", "PDL1", "HPGDS", "PTGDS", "PTGDR", "SLCO2A1")]
  
  cor_test(data =data, 
           name="GSE30219",
           method = "spearman")
  
  
  #resGSE31210
  data <- resGSE31210[, c("HOT_score", "PDL1", "HPGDS", "PTGDS", "PTGDR", "SLCO2A1")]
  
  cor_test(data =data, 
           name="GSE31210",
           method = "spearman")
}

####cal auc
{
  roc <- pROC::roc(resTCGA$Phenotype, resTCGA$PTGDR)
  
  ##画曲线
  
  p<-ggplot(resTCGA, aes(d =Phenotype , m = PTGDR)) + 
    geom_roc() + 
    theme_bw() + 
    labs(title = "ROC", x = "1-Specificity", y = "Sensitivity")
  
  resTCGA$Dataset <- "TCGA"
  resGSE37745$Dataset <- "GSE37745"
  resGSE31210$Dataset <- "GSE31210"
  resGSE30219$Dataset <- "GSE30219"
  
  
  combined_data <- rbind(resTCGA, resGSE37745, resGSE31210, resGSE30219)
  
  
  ##PTGDR
  p <- ggplot(combined_data, aes(d = Phenotype, m = PTGDR, color = Dataset)) +
    geom_roc() +
    theme_bw() +
    labs(title = "PTGDR ROC Curves for Multiple Datasets", x = "1-Specificity", y = "Sensitivity") +
    scale_color_brewer(palette = "Set1")# 设置颜色
  
  auc_table<-data.frame()
  for(dataset in unique(combined_data$Dataset)){
    auc_data <- combined_data %>%
      filter(Dataset == dataset)
    # 在每个数据集上计算 AUC
    roc_curve <- roc(auc_data$Phenotype, auc_data$PTGDR)
    auc_value <- auc(roc_curve)  # 计算 AUC
    auc_table<-rbind(auc_table,data.frame(dataset,auc_value))
  }
  
  auc_table$gene<-"PTGDR"
  
  write.csv(auc_table, file = paste0(ResultPath,"PTGDR_auc_table.csv"),  quote = F, row.names = F)
  
  ggsave(p, file = paste0(ResultPath, "/PTGDR_ROC.pdf"),
         bg = "white", width = 6, height = 6)
  
  
  
  
  ##PTGDS
  p <- ggplot(combined_data, aes(d = Phenotype, m = PTGDS, color = Dataset)) +
    geom_roc() +
    theme_bw() +
    labs(title = "PTGDS ROC Curves for Multiple Datasets", x = "1-Specificity", y = "Sensitivity") +
    scale_color_brewer(palette = "Set1")# 设置颜色
  
  
  auc_table<-data.frame()
  for(dataset in unique(combined_data$Dataset)){
    auc_data <- combined_data %>%
      filter(Dataset == dataset)
    # 在每个数据集上计算 AUC
    roc_curve <- roc(auc_data$Phenotype, auc_data$PTGDS)
    auc_value <- auc(roc_curve)  # 计算 AUC
    auc_table<-rbind(auc_table,data.frame(dataset,auc_value))
  }
  
  auc_table$gene<-"PTGDS"
  
  write.csv(auc_table, file = paste0(ResultPath,"PTGDS_auc_table.csv"),  quote = F, row.names = F)
  
  ggsave(p, file = paste0(ResultPath, "/PTGDS_ROC.pdf"),
         bg = "white", width = 6, height = 6)
  
  
  
  p <- ggplot(combined_data, aes(d = Phenotype, m = SLCO2A1, color = Dataset)) +
    geom_roc() +
    theme_bw() +
    labs(title = "SLCO2A1 ROC Curves for Multiple Datasets", x = "1-Specificity", y = "Sensitivity") +
    scale_color_brewer(palette = "Set1")# 设置颜色
  
  auc_table<-data.frame()
  for(dataset in unique(combined_data$Dataset)){
    auc_data <- combined_data %>%
      filter(Dataset == dataset)
    # 在每个数据集上计算 AUC
    roc_curve <- roc(auc_data$Phenotype, auc_data$SLCO2A1)
    auc_value <- auc(roc_curve)  # 计算 AUC
    auc_table<-rbind(auc_table,data.frame(dataset,auc_value))
  }
  
  auc_table$gene<-"SLCO2A1"
  
  write.csv(auc_table, file = paste0(ResultPath,"SLCO2A1_auc_table.csv"),  quote = F, row.names = F)
  
  ggsave(p, file = paste0(ResultPath, "/SLCO2A1_ROC.pdf"),
         bg = "white", width = 6, height = 6)
  
  
  
  
  p <- ggplot(combined_data, aes(d = Phenotype, m = HPGDS, color = Dataset)) +
    geom_roc() +
    theme_bw() +
    labs(title = "HPGDS ROC Curves for Multiple Datasets", x = "1-Specificity", y = "Sensitivity") +
    scale_color_brewer(palette = "Set1")# 设置颜色
  
  
  auc_table<-data.frame()
  for(dataset in unique(combined_data$Dataset)){
    auc_data <- combined_data %>%
      filter(Dataset == dataset)
    # 在每个数据集上计算 AUC
    roc_curve <- roc(auc_data$Phenotype, auc_data$HPGDS)
    auc_value <- auc(roc_curve)  # 计算 AUC
    auc_table<-rbind(auc_table,data.frame(dataset,auc_value))
  }
  
  auc_table$gene<-"HPGDS"
  
  write.csv(auc_table, file = paste0(ResultPath,"HPGDS_auc_table.csv"),  quote = F, row.names = F)
  
  ggsave(p, file = paste0(ResultPath, "/HPGDS_ROC.pdf"),
         bg = "white", width = 6, height = 6)
}


## 保存数据
saveRDS(combined_data, file = "./Rds/combined_data_4cohort.rds")

# estmate ssgsea 见“肺腺癌代谢物_ssGSEA_ESTIMATE_PDL1_相关性热图” by chenliuxin




