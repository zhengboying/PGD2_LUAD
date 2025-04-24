

## diff test for 药物敏感性
{
  source("../../Function/diff_test.R")
  
  #DrugPredictions<-read.csv("./03_calcPhenotype/TGCA_GDSC2Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  
  DrugPredictions<-read.csv("./03_calcPhenotype/TCGA_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  rownames(DrugPredictions) <- gsub("_T_LUAD","",rownames(DrugPredictions))
  
  #densityplot FOR DrugPredictions$Ponatinib_155(),极端值太多 不符合正态，只能用非参数检验
  ggplot(DrugPredictions, aes(x = Ponatinib_155)) +
    geom_density(fill = "blue", alpha = 0.5) + # 你可以调整填充颜色和透明度
    labs(title = "Density Plot of Ponatinib_155", x = "Ponatinib_155", y = "Density") # 
  
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  
  
  
  #group_table = combined_data[combined_data$Dataset == "TCGA",]
  group_table = model_score$TCGA
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
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
  
}


{
  
  
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE37745_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  group_table = model_score$GSE37745
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
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
  
}


{
  source("../../Function/diff_test.R")
  
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE31210_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  
  DrugPredictions<-as.data.frame(t(DrugPredictions))
  
  group_table = model_score$GSE31210
  
  ##split group_table,if PTGDS > mean of PTGDS,then PTGDS_level = "High",else PTGDS_level = "Low"
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > mean(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > mean(group_table$PTGDR),"High","Low")
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  #group_table$PTGDS_level <- ifelse(group_table$PTGDS > median(group_table$PTGDS),"High","Low")
  #group_table$PTGDR_level <- ifelse(group_table$PTGDR > median(group_table$PTGDR),"High","Low")
  
  
  #
  diffDrugs_risksocre_level <- perform_diff_test(df=DrugPredictions, 
                                                 group_table = group_table, 
                                                 group_name = "risk_score_level",
                                                 case="High",
                                                 control="Low",
                                                 diff_test_method = "t.test",##目前仅支持t_test,Wilcoxon test(w.test),auc test
                                                 p.adjust.method = "BH",
                                                 remove_sum0_features = FALSE
  )
  
  diffDrugs_risksocre_level_sig <- diffDrugs_risksocre_level[diffDrugs_risksocre_level$p.adjust < 0.05,]
  dim(diffDrugs_risksocre_level_sig)
  GSE31210_diffDrugs<-diffDrugs_risksocre_level_sig
}

as.data.frame(rownames(tcga_diffDrugs))
as.data.frame(rownames(GSE31210_diffDrugs))
as.data.frame(rownames(GSE37745_diffDrugs))


overlap_drugs <- c(
  "WH.4.023_56",
  "WZ.1.84_59",
  "Ponatinib_155",
  "Enzastaurin_229",
  "Cabozantinib_249",
  "Sphingosine.Kinase.1.Inhibitor.II_408",
  "Alisertib_431",
  "Veliparib_1018",
  "Lenalidomide_1020",
  "SL0101_1039",
  "PD0325901_1060",
  "Selumetinib_1062",
  "HG.5.88.01_1143",
  "KU.60019_1185",
  "UNC1215_1262",
  "RAF_9304_1453",
  "Pictilisib_1527"
)


overlap_drugs_1 <- c(
  "WH.4.023_56",
  "WZ.1.84_59",
  "Ponatinib_155",
  #"Enzastaurin_229",
  "Cabozantinib_249",
  #"Sphingosine.Kinase.1.Inhibitor.II_408",
  "Alisertib_431",
  #"Veliparib_1018",
  #"Lenalidomide_1020",
  #"SL0101_1039",
  "PD0325901_1060",
  "Selumetinib_1062",
  "HG.5.88.01_1143",
  #"KU.60019_1185",
  #"UNC1215_1262",
  "RAF_9304_1453",
  "Pictilisib_1527"
)


overlap_drugs_2<-c("Z.LLNle.CHO_45", "Dasatinib_51", "WH.4.023_56", "WZ.1.84_59", "Pyrimethamine_71", 
                   "TGX221_94", "CP466722_152", "Ponatinib_155", "JNK.9L_157", "HG6.64.1_159", 
                   "Enzastaurin_229", "GSK429286A_230", "Idelalisib_238", "Cabozantinib_249", 
                   "CP724714_255", "NSC.207895_269", "Selisistat_341", "Brivanib..BMS.540215_376", 
                   "Alisertib_431", "CD532_449", "LDN.193189_478", "FTY.720_546", "Refametinib_1014", 
                   "Olaparib_1017", "Veliparib_1018", "Lenalidomide_1020", "Axitinib_1021", "VX.702_1028", 
                   "Vismodegib_1033", "SL0101_1039", "Mirin_1048", "Selumetinib_1062", "AZD6482_1066", 
                   "EHT.1864_1069", "HG.5.88.01_1143", "TW.37_1149", "ZG.10_1161", "KU.60019_1185", 
                   "SB505124_1194", "Trametinib_1372", "Temozolomide_1375", "Doxorubicin_1386", 
                   "AZD8186_1444", "AZD1208_1449", "Selumetinib_1498", "Refametinib_1526", "Pictilisib_1527")



overlap_drugs_3<-c("Z.LLNle.CHO_45", "Dasatinib_51", "WH.4.023_56", "WZ.1.84_59", "TGX221_94", 
                   "CP466722_152", "Ponatinib_155", "HG6.64.1_159", "Idelalisib_238", 
                   "Cabozantinib_249", "Alisertib_431", "LDN.193189_478", "FTY.720_546", 
                   "Refametinib_1014", "Selumetinib_1062", "HG.5.88.01_1143", "ZG.10_1161", 
                   "Trametinib_1372", "Selumetinib_1498", "Refametinib_1526", "Pictilisib_1527")



overlap_drugs_4<- c("Pictilisib_1527", "FTY.720_546", "WH.4.023_56", 
                    "LDN.193189_478", "Dasatinib_51", "HG6.64.1_159", "Z.LLNle.CHO_45")
#取diffDrugs_risksocre_level_sig  在overlap_drugs中的行


write.csv(tcga_diffDrugs, file = "./03_calcPhenotype/tcga_diffDrugs.csv", row.names = T,quote = F)
write.csv(GSE37745_diffDrugs, file = "./03_calcPhenotype/GSE37745_diffDrugs.csv", row.names = T,quote = F)
write.csv(GSE31210_diffDrugs, file = "./03_calcPhenotype/GSE31210_diffDrugs.csv", row.names = T,quote = F)


tcga_diffDrugs<-read.csv("./03_calcPhenotype/tcga_diffDrugs.csv",row.names = 1)
GSE37745_diffDrugs<-read.csv("./03_calcPhenotype/GSE37745_diffDrugs.csv",row.names = 1)
GSE31210_diffDrugs<-read.csv("./03_calcPhenotype/GSE31210_diffDrugs.csv",row.names = 1)

##检查foldchange一致性
annot_cols1<-tcga_diffDrugs[rownames(tcga_diffDrugs) %in% overlap_drugs_3,]
annot_cols2<-GSE37745_diffDrugs[rownames(GSE37745_diffDrugs) %in% overlap_drugs_3,]
annot_cols3<-GSE31210_diffDrugs[rownames(GSE31210_diffDrugs) %in% overlap_drugs_3,]
annot_cols1<-annot_cols1[order(-annot_cols1$p.adjust),] 



library("pheatmap")
#heatmap P<0.05
{
  DrugPredictions<-read.csv("./03_calcPhenotype/TCGA_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,overlap_drugs_3]
  rownames(DrugPredictions) <- gsub("_T_LUAD","",rownames(DrugPredictions))
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$TCGA
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  
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
  
  ggsave("./03_calcPhenotype/heatmap-TCGA.pdf", p, width = 10, height = 3.5, units = "in")
}





##GSE 37745

{
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE37745_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,overlap_drugs_3]
  
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$GSE37745
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  
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
  
  ggsave("./03_calcPhenotype/heatmap-GSE37745.pdf", p1, width = 10, height = 3.5, units = "in")
}

##GSE31210
{
  DrugPredictions<-read.csv("./03_calcPhenotype/GSE31210_GDSC1Data_calcPhenotype_Output/calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
  DrugPredictions<-as.data.frame(DrugPredictions)
  DrugPredictions<-DrugPredictions[,overlap_drugs_3]
  
  
  ##group_table "risk_score"从大到小排序
  group_table<-model_score$GSE31210
  group_table$risk_score_level <- ifelse(group_table$risk_score > mean(group_table$risk_score),"High","Low")
  
  
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
  
  ggsave("./03_calcPhenotype/heatmap-GSE31210.pdf", p2, width = 10, height = 3.5, units = "in")
}

grob_p <- grid.grabExpr(p)
grob_p1 <- grid.grabExpr(p1)
grob_p2 <- grid.grabExpr(p2)

# 使用 grid.arrange 来排列这些 grob 对象
p_all<-grid.arrange(grob_p, grob_p1, grob_p2, ncol = 1, nrow = 3)

ggsave("./03_calcPhenotype/heatmap-all.pdf", p_all, width = 10, height = 20, units = "in")

#ggplot(DrugPredictions, aes(x = ZG.10_1161)) +
#  geom_density(fill = "blue", alpha = 0.5) + # 你可以调整填充颜色和透明度
#  labs(title = "Density Plot of Ponatinib_155", x = "ZG.10_1161", y = "Density") # 