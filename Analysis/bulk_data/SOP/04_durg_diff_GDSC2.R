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

library(ggplot2)
# 01.Setup -----

project_path =  "/Users/zhengboying/Documents/project_demo/Analysis/bulk_data/"
setwd(project_path)
InputPath <- c("./Rds/")
model_score<-readRDS(paste0(InputPath,"model_score.rds"))

##GDSC2 198

## diff test for 药物敏感性
{
  source("../../Function/diff_test.R")
  
  DrugPredictions<-read.csv("./03_calcPhenotype/TGCA_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  
  #DrugPredictions<-read.csv("./03_calcPhenotype/TCGA_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  rownames(DrugPredictions) <- gsub("_T_LUAD","",rownames(DrugPredictions))
  
  #densityplot FOR DrugPredictions$Ponatinib_155(),极端值太多 不符合正态，只能用非参数检验
  ggplot(DrugPredictions, aes(x = JQ1_2172)) +
    geom_density(fill = "blue", alpha = 0.5) + # 你可以调整填充颜色和透明度
    labs(title = "Density Plot of JQ1_2172", x = "JQ1_2172", y = "Density") # 
  
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  
  
  
  #group_table = combined_data[combined_data$Dataset == "TCGA",]
  group_table = model_score$TCGA
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  #group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > median(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > median(group_table$PTGDR),"High","Low")
  
  diffDrugs_risksocre_level <- perform_diff_test(df=DrugPredictions, 
                                                 group_table = group_table, 
                                                 group_name = "risk_score_level",
                                                 case="High",
                                                 control="Low",
                                                 diff_test_method = "w.test",##目前仅支持t_test,Wilcoxon test(w.test),auc test
                                                 p.adjust.method = "BH",
                                                 remove_sum0_features = FALSE
  )
  diffDrugs_risksocre_level_sig <- diffDrugs_risksocre_level[diffDrugs_risksocre_level$p.adjust < 0.05,]
  dim(diffDrugs_risksocre_level_sig)
  tcga_diffDrugs<-diffDrugs_risksocre_level_sig
  write.csv(tcga_diffDrugs, file = "./03_calcPhenotype/GDSC2_tcga_diffDrugs.csv", row.names = T,quote = F)
  
  
}


{
  
  
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE37745_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  group_table = model_score$GSE37745
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  #group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > median(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > median(group_table$PTGDR),"High","Low")
  
  
  diffDrugs_risksocre_level <- perform_diff_test(df=DrugPredictions, 
                                                 group_table = group_table, 
                                                 group_name = "risk_score_level",
                                                 case="High",
                                                 control="Low",
                                                 diff_test_method = "w.test",##目前仅支持t_test,Wilcoxon test(w.test),auc test
                                                 p.adjust.method = "BH",
                                                 remove_sum0_features = FALSE
  )
  diffDrugs_risksocre_level_sig <- diffDrugs_risksocre_level[diffDrugs_risksocre_level$p.adjust < 0.05,]
  dim(diffDrugs_risksocre_level_sig)
  GSE37745_diffDrugs<-diffDrugs_risksocre_level_sig
  write.csv(GSE37745_diffDrugs, file = "./03_calcPhenotype/GDSC2_GSE37745_diffDrugs.csv", row.names = T,quote = F)
  
}


{
  source("../../Function/diff_test.R")
  
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE31210_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  group_table = model_score$GSE31210
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  #group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > median(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > median(group_table$PTGDR),"High","Low")
  
  
  #
  diffDrugs_risksocre_level <- perform_diff_test(df=DrugPredictions, 
                                                 group_table = group_table, 
                                                 group_name = "risk_score_level",
                                                 case="High",
                                                 control="Low",
                                                 diff_test_method = "w.test",##目前仅支持t_test,Wilcoxon test(w.test),auc test
                                                 p.adjust.method = "BH",
                                                 remove_sum0_features = FALSE
  )
  
  diffDrugs_risksocre_level_sig <- diffDrugs_risksocre_level[diffDrugs_risksocre_level$p.adjust < 0.05,]
  dim(diffDrugs_risksocre_level_sig)
  GSE31210_diffDrugs<-diffDrugs_risksocre_level_sig
  write.csv(GSE31210_diffDrugs, file = "./03_calcPhenotype/GDSC2_GSE31210_diffDrugs.csv", row.names = T,quote = F)
}

as.data.frame(rownames(tcga_diffDrugs))
as.data.frame(rownames(GSE31210_diffDrugs))
as.data.frame(rownames(GSE37745_diffDrugs))

##tcga_diffDrugs ,GSE31210_diffDrugs,GSE37745_diffDrugs 交集
overlap_drugs<-Reduce(intersect, list(rownames(tcga_diffDrugs),rownames(GSE31210_diffDrugs),rownames(GSE37745_diffDrugs)))

#提取交集的表格并合并
tcga_diffDrugs<-tcga_diffDrugs[rownames(tcga_diffDrugs) %in% overlap_drugs,]
GSE37745_diffDrugs<-GSE37745_diffDrugs[rownames(GSE37745_diffDrugs) %in% overlap_drugs,]
GSE31210_diffDrugs<-GSE31210_diffDrugs[rownames(GSE31210_diffDrugs) %in% overlap_drugs,]

#合并rbind
diffDrugs<-cbind(tcga_diffDrugs,GSE37745_diffDrugs,GSE31210_diffDrugs)
names(diffDrugs)[1:6] <- paste0("TCGA_", names(diffDrugs)[1:6])
names(diffDrugs)[7:12] <- paste0("GSE37745_", names(diffDrugs)[7:12])
names(diffDrugs)[13:18] <- paste0("GSE31210_", names(diffDrugs)[13:18])
#tcga_diffDrugs<-read.csv("./03_calcPhenotype/GDSC2_tcga_diffDrugs.csv",row.names = 1)
#GSE37745_diffDrugs<-read.csv("./03_calcPhenotype/GDSC2_GSE37745_diffDrugs.csv",row.names = 1)
#GSE31210_diffDrugs<-read.csv("./03_calcPhenotype/GDSC2_GSE31210_diffDrugs.csv",row.names = 1)




df<-diffDrugs
all_gt_0 <- df$TCGA_log2FC > 0 & df$GSE37745_log2FC > 0 & df$GSE31210_log2FC > 0

# 创建一个逻辑向量，表示所有三列都小于0的行
all_lt_0 <- df$TCGA_log2FC < 0 & df$GSE37745_log2FC < 0 & df$GSE31210_log2FC < 0

# 合并这两个逻辑向量，选择同时满足所有三列都大于0或都小于0的行
selected_rows <- df[all_gt_0 | all_lt_0, ]

# 查看选中的行
print(selected_rows)

write.csv(selected_rows, file = "./03_calcPhenotype/GDSC2_diffDrugs_overlap.csv", row.names = T,quote = F)

selected_drugs<-rownames(selected_rows)

library("pheatmap")


annot_cols1<-tcga_diffDrugs[rownames(tcga_diffDrugs) %in% selected_drugs,]
annot_cols2<-GSE37745_diffDrugs[rownames(GSE37745_diffDrugs) %in% selected_drugs,]
annot_cols3<-GSE31210_diffDrugs[rownames(GSE31210_diffDrugs) %in% selected_drugs,]
annot_cols1<-annot_cols1[order(-annot_cols1$p.adjust),] 
#heatmap P<0.05
{
  DrugPredictions<-read.csv("./03_calcPhenotype/TCGA_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,selected_drugs]
  rownames(DrugPredictions) <- gsub("_T_LUAD","",rownames(DrugPredictions))
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$TCGA
  #group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  
  group_table_sorted <- group_table[order(-group_table$risk_score),] 
  DrugPredictions_sorted <- DrugPredictions[rownames(group_table_sorted), ]
  
  ##去掉含有NA的行
  DrugPredictions_sorted <- DrugPredictions_sorted[complete.cases(DrugPredictions_sorted),]
  
  #group_table_sorted去掉不在DrugPredictions_sorted中的行
  group_table_sorted <- group_table_sorted[rownames(group_table_sorted) %in% rownames(DrugPredictions_sorted),]
  
  
  DrugPredictions_sorted <- DrugPredictions_sorted[, rownames(annot_cols1)]
  
  colnames(DrugPredictions_sorted)
  rownames(annot_cols1)
  
  
  p.adjust_colors <- colorRampPalette(c("#2E7025", "#4EBF3F","#D5E1D3"))(20)
  p<-pheatmap(log(DrugPredictions_sorted),cluster_rows = F,cluster_cols = F,main = "DrugPredictions",
              fontsize = 8,show_rownames = F,
              scale = "column",  show_colnames = T,
              #color
              color = colorRampPalette(c("#00008B","blue", "white", "red","#993366"))(200),
              annotation_row = group_table_sorted[,c("risk_score","risk_score_level")],
              annotation_col = annot_cols1[,c("log2FC","p.adjust")],
              annotation_colors = list(risk_score_level = c("High" = "#C23930", "Low" = "#39537A"),
                                       p.adjust = p.adjust_colors[order(-annot_cols1$p.adjust)] )
  )
  
  ggsave("./03_calcPhenotype/GDSC2_heatmap-TCGA.pdf", p, width = 5, height = 3.5, units = "in")
}





##GSE 37745

{
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE37745_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,selected_drugs]
  
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$GSE37745
  
  
  group_table_sorted <- group_table[order(-group_table$risk_score),] 
  DrugPredictions_sorted <- DrugPredictions[rownames(group_table_sorted), ]
  
  ##去掉含有NA的行
  DrugPredictions_sorted <- DrugPredictions_sorted[complete.cases(DrugPredictions_sorted),]
  
  #group_table_sorted去掉不在DrugPredictions_sorted中的行
  group_table_sorted <- group_table_sorted[rownames(group_table_sorted) %in% rownames(DrugPredictions_sorted),]
  
  
  
  # 使用match函数获取排序后的索引(tcga保持一致)
  
  DrugPredictions_sorted <- DrugPredictions_sorted[, rownames(annot_cols1)]
  
  colnames(DrugPredictions_sorted)
  rownames(annot_cols1)
  
  
  p.adjust_colors <- colorRampPalette(c("#2E7025", "#4EBF3F","#D5E1D3"))(20)
  p1<-pheatmap(log(DrugPredictions_sorted),cluster_rows = F,cluster_cols = F,main = "DrugPredictions",
               fontsize = 8,show_rownames = F,
               scale = "column",  show_colnames = T,
               #color
               color = colorRampPalette(c("#00008B","blue", "white", "red","#993366"))(200),
               annotation_row = group_table_sorted[,c("risk_score","risk_score_level")],
               annotation_col = annot_cols2[,c("log2FC","p.adjust")],
               annotation_colors = list(risk_score_level = c("High" = "#C23930", "Low" = "#39537A"),
                                        p.adjust = p.adjust_colors[order(-annot_cols2$p.adjust)] )
  )
  
  ggsave("./03_calcPhenotype/GDSC2_heatmap-GSE37745.pdf", p1, width = 5, height = 2.5, units = "in")
}

##GSE31210
{
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE31210_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,selected_drugs]
  
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$GSE31210
  
  
  
  group_table_sorted <- group_table[order(-group_table$risk_score),] 
  DrugPredictions_sorted <- DrugPredictions[rownames(group_table_sorted), ]
  
  ##去掉含有NA的行
  DrugPredictions_sorted <- DrugPredictions_sorted[complete.cases(DrugPredictions_sorted),]
  
  #group_table_sorted去掉不在DrugPredictions_sorted中的行
  group_table_sorted <- group_table_sorted[rownames(group_table_sorted) %in% rownames(DrugPredictions_sorted),]
  
  
  
  # 使用match函数获取排序后的索引(tcga保持一致)
  
  DrugPredictions_sorted <- DrugPredictions_sorted[, rownames(annot_cols1)]
  
  colnames(DrugPredictions_sorted)
  rownames(annot_cols1)
  
  
  p.adjust_colors <- colorRampPalette(c("#2E7025", "#4EBF3F","#D5E1D3"))(20)
  p2<-pheatmap(log(DrugPredictions_sorted),
               cluster_rows = F,cluster_cols = F,main = "DrugPredictions",
               fontsize = 8,show_rownames = F,
               scale = "column",  show_colnames = T,
               #color
               color = colorRampPalette(c("#00008B","blue", "white", "red","#993366"))(200),
               annotation_row = group_table_sorted[,c("risk_score","risk_score_level")],
               annotation_col = annot_cols3[,c("log2FC","p.adjust")],
               annotation_colors = list(risk_score_level = c("High" = "#C23930", "Low" = "#39537A"),
                                        p.adjust = p.adjust_colors[order(-annot_cols3$p.adjust)] )
  )
  
  ggsave("./03_calcPhenotype/GDSC2_heatmap-GSE31210.pdf", p2, width = 5, height = 3.5, units = "in")
}



