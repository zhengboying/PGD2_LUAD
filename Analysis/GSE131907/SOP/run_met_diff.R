#diff-test for all cell type and NT on mets
{
  library(ComplexHeatmap)
  library(dplyr)
  library(tidyr)
  project_path =  "/Users/zhengboying/Documents/DMPJ0001-PGD2-byzheng/Analysis/GSE131907/"
  setwd(project_path)
  source("../../Function/diff_test.R")
  mets_db1<-'../../Database/metabolite_annotation_HMDB_summary.tsv' #from metbost_db
  mets_db2<-'../../Database/met_sen_Dec-4-2024.tsv' #from metbost_db
  
}



## load data
{
  ##database  
  mets_db1<-read.csv(mets_db1,check.names = F,sep='\t')
  mets_db2<-read.csv(mets_db2,check.names = F,sep='\t')
  
  annot_heatmap<-mets_db1[,c("HMDB_ID","super_class","sub_class","class","metabolite")]
    ##提取HMDB_ID 和 metabolite
    mets_db1<-mets_db1[,c("HMDB_ID","metabolite")]
    length(unique(mets_db1$HMDB_ID))##1241
    mets_db2<-mets_db2[,c("HMDB_ID","standard_metName")]
    length(unique(mets_db2$HMDB_ID))##147
    mets_with_senser<-unique(mets_db2$HMDB_ID)
    mets_db<-mets_db1
    mets_db$senser<-NA
    ##如果mets_db1$HMDB_ID在mets_with_senser，
    ##则mets_db$senser=TRUE
    mets_db$senser[mets_db$HMDB_ID %in% mets_with_senser]<-TRUE
    mets_db$senser[is.na(mets_db$senser)]<-FALSE
    ##
    
    ##met
    #met_mat<-read.csv('./03_mebocost/met_mat.csv',row.names=1,check.names = F)
    #met_mat[1:4,1:6]
    ## 保存
    #saveRDS(met_mat,"./Rds/Mets/met_mat.rds")
    met_mat<-readRDS("./Rds/Mets/met_mat.rds")
    met_mat[1:4,1:6]
    dim(met_mat)
    ##annot infor
    group_table<-read.csv('./03_mebocost/cell_ann.csv',row.names=1,check.names = F)
    
}




##filtered
{
  ## 过滤整行为0的值
  #df <- met_mat[rowSums(met_mat) > 0, ]
}

##执行差异检验
{
  res_all<-perform_pairwise_tests(df = met_mat, group_table =group_table, group_name = "celltype", 
                                  diff_test_method = "w.test", p.adjust.method = "BH") 
}

##生成dataframe
{
  cells<-unique(group_table$celltype)
  
  res_all_df<-res_all
  res_all_df <- lapply(cells, function(i) {
    # 为每个数据框添加一个新列'group'，包含数据框的名称
    df<-data.frame()
    df <- res_all[[i]]
    df$group <- names(res_all)[i]
    return(df)
  })
  
  # 使用do.call和rbind将所有数据框合并成一个大的数据框
  res_all_df <- do.call(rbind,res_all_df)
  
}  

##注释名称和是否存在sensers
{
  res_all_df_mets<- merge(res_all_df,mets_db, by = "HMDB_ID", all.x = TRUE)
  write.csv(res_all_df_mets, "./03.1_met_diff/res_all_celltype.csv") 
  ##read
}  



##生成mean_case宽表
{
  mean_table <- res_all_df_mets %>%
    pivot_wider(names_from = group, values_from = mean_case,
                id_cols = "HMDB_ID") %>%
    ungroup() 
  
  write.csv(mean_table,"./03.1_met_diff/res_mets_mean.csv",row.names = F,quote = F)
}



mean_table<-read.csv("./03.1_met_diff/res_mets_mean.csv",check.names = F,sep=',',row.names=1)
VariableFeatures200<-readRDS("./Rds/Mets/2.metobj_VariableFeatures_200.Rds")
res_all_df_mets<-read.csv("./03.1_met_diff/res_all_celltype.csv",check.names = F,sep=',',row.names=1)





## 筛选mean>0.01."log2FC"> 0,

{
  res_all_df_high<-res_all_df_mets[res_all_df_mets[,"mean_case"] > 0.01,]
  res_all_df_high<-res_all_df_high[res_all_df_high[,"log2FC"] > 0,]
  ##去掉NA
  res_all_df_high <-res_all_df_high[!is.na(res_all_df_high$log2FC),]
  ##去掉Inf
  res_all_df_high<-res_all_df_high[!is.infinite(res_all_df_high$log2FC),]
  
  dim(res_all_df_high)
  table(res_all_df_high$group)
}


##细胞上调UP差异火山图

{
  ##load
  sapply(c('reshape2','pheatmap','ggrepel','dplyr','RColorBrewer','ggprism','ggforce','ggnewscale'),require,character.only=TRUE)
  
  ##annte
  {
    plot_dat<-res_all_df_high
    ##增加一列label，up,down,unchange
    plot_dat$label<-ifelse(plot_dat$p.adjust>0.001,"Unchange","NA")
    ##log2FC>0 && p.adjust<0.01
    plot_dat$label[plot_dat$log2FC>0 & plot_dat$p.adjust<0.001]<-"Up"
    #plot_dat$label[plot_dat$log2FC<0 & plot_dat$p.adjust<0.001]<-"Down"

    #label_data<-plot_dat%>%group_by(group)%>%arrange(group,desc(log2FC))%>%do(head(.,n=2))
  }
  
  ##plot
  {
    p1 <- ggplot(plot_dat)+
      # 抖动散点：
      geom_jitter(aes(group, log2FC, color = label),
                  size=1.5, width = 0.4)+
      # 分组方块：
      geom_tile(aes(group, 0, fill = group),
                height=0.5,
                color = "white",
                alpha = 0.5,
                show.legend = F,
                width=0.85) +
      # 文字：
      geom_text(data = plot_dat[!duplicated(plot_dat$group), ], 
                aes(group, 0, label = group),
                size =5,
                color ="black") +
      # 基因标签：
      #geom_text_repel(
      #  data = label_data,
      #  aes(group, log2FC, label = metabolite),
      #  size=4, max.overlaps = 100
      #) +
      xlab("Group")+
      ylab("log2FoldChange")+
      # 颜色模式
      scale_fill_manual(values = c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2")))+
      scale_color_manual(values = c('gray','#C00C06'))+
      theme_prism()+
      theme(axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "top")
    p1
  }
  
  ##Add label
  {
    ## commscore>1 prob>0.1,p<0.01
    over_mets<-c("Cholesterol", "Choline", "gamma-Aminobutyric acid", "D-Mannose", "Riboflavin", "L-Cysteine",
                 "L-Glutamine", "Leukotriene B4", "Prostaglandin E2", "Prostaglandin D2",
                 "all-trans-Retinoic acid","9-cis-Retinoic acid","13-cis-Retinoic acid",
                 "Sphingosine 1-phosphate","Vitamin A","Adenine","Histamine","Glycerol","Iron",
                 "Uridine diphosphate glucose")
    
    ## top10_mets<- plot_dat$log2FC top 5 
    top20_mets<-plot_dat%>%group_by(group)%>%arrange(group,desc(log2FC))%>%do(head(.,n=20))
    top20_mets<-top20_mets$metabolite
    
    ##
    marker_mets<-intersect(over_mets,top20_mets)
    ##[1] "Cholesterol"      "Choline"          "D-Mannose"        "Riboflavin"       "L-Cysteine"       "Prostaglandin D2"
    label_data <- plot_dat %>% filter(metabolite %in% marker_mets)
    
    
    
    
    p2<-p1 + new_scale_fill()+
      geom_mark_circle(data = label_data,
                       aes(group, log2FC, fill = metabolite, label = metabolite),
                       expand= unit(1.5,"mm"),
                       show.legend = T,
                       label.fontsize = 12,
                       con.cap = 1)
    p2
    ggsave("./03.1_met_diff/valcano_padjust.pdf", p2,height = 5.5, width = 12)
    
    
    ggsave("./03.1_met_diff/plot_dat.pdf",p1,width=10,height=5,units="in")
  }
}



##提取HD_ID res_all_df_high$metabolite %in% marker_mets
{
  marker_mets_ID<-res_all_df_high[res_all_df_high$metabolite %in% marker_mets,]
  marker_mets_ID<-unique(marker_mets_ID$HMDB_ID)
  
  ##
  mean_table_plot <- mean_table[rownames(mean_table) %in% marker_mets_ID, ]
  dim(mean_table_plot)
  
  annot_heatmap1<-annot_heatmap[annot_heatmap$HMDB_ID %in% rownames(mean_table_plot),]
  
  rownames(annot_heatmap1) <- annot_heatmap1$HMDB_ID
  annot_heatmap1$HMDB_ID <- NULL
  annot_heatmap1$sub_class <- NULL
  annot_heatmap1$class <- NULL
  annot_heatmap1$metabolite <- NULL
  ha = rowAnnotation(df = annot_heatmap1)
  
  ####注释行
  ht<-Heatmap(mean_table_plot,
              row_names_side = "right",  
              show_row_names = T, 
              show_column_names = T,
              right_annotation  = ha,
              split = annot_heatmap1$super_class)
  
  #
  pdf(file = "./03.1_met_diff/maker_mets_heatmap.pdf", width =8, height =6)
  draw(ht)
  dev.off()  
  
  
}






##提取行名Prostaglandin D2的结果
#HMDB0001403
{
  PGD2_table <- t(met_mat[rownames(met_mat) == "HMDB0001403", ])
  
  ##与annot信息合并
  merged_df <- merge(PGD2_table, group_table, by = "row.names")
  head(merged_df)
  #write
  write.csv(merged_df,"./03.1_met_diff/PGD2_table.csv",row.names = F,quote = F)
}

##画vlnplot图并计算显著性
library(ggsignif)
#Olig Mast pDC Fibro Myeloid Endo NK

## 画location 组的vlnplot

## 比较N/T/M
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro","Oligo"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
  
  
  ## 横坐标顺序
  location_order<- c("Normal Lung","Tumor Lung","Metastasis")#,"Normal LN","Metastatic LN","MPE") 
  ## 画图
  p3<-ggplot(df,
             aes(x = location, y = HMDB0001403, fill = location)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = location_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Location") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c(
                                 "Normal Lung" = "#008000",
                                 "Tumor Lung"= "#B41D23",
                                 "Metastasis"="#FFA500"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p3
  
  ##添加显著性
  p4<-p3 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("Tumor Lung", "Metastasis"),
                                   c("Normal Lung", "Tumor Lung")),
                vjust = 0.3) +
  theme(legend.position = "none")
  p4
}  
ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_NTM.pdf", p4,height = 8, width = 4)

## 比所有
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro","Myeloid", "Endo", "NK","Oligo"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  #df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
  
  
  ## 横坐标顺序
  location_order<- c("Normal Lung","Tumor Lung","Metastasis","Normal LN","Metastatic LN","MPE") 
  ## 画图
  p5<-ggplot(df,#
             aes(x = location, y = HMDB0001403, fill = location)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = location_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Location") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c("Metastatic LN" = "#AF689B",
                                 "MPE"= "#B7B577",
                                 "Normal LN" = "#11E693",
                                 "Normal Lung" = "#008000",
                                 "Tumor Lung"= "#B41D23",
                                 "Metastasis"="#FFA500"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p5
  
  
  ##添加显著性
  p6<-p5 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("Tumor Lung", "Metastasis"),
                                   c("Normal Lung", "Tumor Lung"),
                                   c("Normal LN", "Metastatic LN")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p6
}  
ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_all.pdf", p6,height = 12, width = 4)


## 比进展
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro","Oligo"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  df<-subset(df, location %in% c("Tumor Lung"))
  
  
  ## 横坐标顺序
  stage.1_order<- c("early stage","advanced stage") 
  stage.2_order<- c("I","II","III","IV") 
  ## 画图
  p7<-ggplot(df,#
             aes(x =stage.1, y = HMDB0001403, fill = stage.1)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = stage.1_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #xlab("Stage") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c("early stage"= "#B7B577",
                                 "advanced stage"= "#B41D23"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p7
  
  
  ##添加显著性
  p8<-p7 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("early stage", "advanced stage")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p8
}  
ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_stage1.pdf", p8,height = 7, width = 4)



## 比进展(所有)
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro","Myeloid", "Endo", "NK","Oligo"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  df<-subset(df, location %in% c("Tumor Lung"))
  
  ## 横坐标顺序
  stage.1_order<- c("early stage","advanced stage") 
  stage.2_order<- c("I","II","III","IV") 
  ## 画图
  p9<-ggplot(df,#
             aes(x =stage.2, y = HMDB0001403, fill = stage.2)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = stage.2_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #xlab("Stage") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c("I"= "#B7B577",
                                 "II"= "#9FA500",
                                 "III" = "#FFA500",
                                 "IV"= "#B41D23"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p9
  
  ##添加显著性
  p10<-p9 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("I", "II"),
                                   c("II", "III"),
                                   c("III", "IV")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p10
}

ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_stage2.pdf", p10,height = 12, width = 4)


## 比原位癌细胞分化程度
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro","Myeloid", "Endo", "NK","Oligo"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  df<-subset(df, location %in% c("Tumor Lung"))
  # 去掉df$Differentiated_state ==na
  df<-subset(df, Differentiated_state != "na")
  
  ## 横坐标顺序
  Differentiated_state_order<- c("WD","MD","PD") 
  ## 画图
  p11<-ggplot(df,#
             aes(x =Differentiated_state, y = HMDB0001403, fill = Differentiated_state)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = Differentiated_state_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #xlab("Stage") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c("WD"= "#B7B577",
                                 "MD" = "#FFA500",
                                 "PD"= "#B41D23"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p11
  
  ##添加显著性
  p12<-p11 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("WD", "MD"),
                                   c("MD", "PD"),
                                   c("PD", "WD")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p12
}

ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_DIFF_ALL.pdf", p12,height = 11, width = 4)


## 比原位癌细胞分化程度
{
  #提取数据
  df<-subset(merged_df, celltype %in% c( "Mast", "pDC", "Fibro"))
  #提取"Normal Lung","Tumor Lung","Metastasis"
  df<-subset(df, location %in% c("Tumor Lung"))
  # 去掉df$Differentiated_state ==na
  df<-subset(df, Differentiated_state != "na")
  
  ## 横坐标顺序
  Differentiated_state_order<- c("WD","MD","PD") 
  ## 画图
  p13<-ggplot(df,#
              aes(x =Differentiated_state, y = HMDB0001403, fill = Differentiated_state)) +
    geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_discrete(limits = Differentiated_state_order) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #xlab("Stage") +
    ylab("Prostaglandin D2") +
    scale_fill_manual(values = c("WD"= "#B7B577",
                                 "MD" = "#FFA500",
                                 "PD"= "#B41D23"))+
    #spit wrap celltype,一行一个
    facet_wrap(~celltype, ncol = 1)+
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
  p13
  
  ##添加显著性
  p14<-p13 +
    geom_signif(test = wilcox.test, 
                size = 0.2, 
                comparisons = list(c("WD", "MD"),
                                   c("MD", "PD"),
                                   c("PD", "WD")),
                vjust = 0.3) +
    theme(legend.position = "none")
  p14
}

ggsave("./03.1_met_diff/Prostaglandin_D2_vlnplot_DIFF.pdf", p14,height = 6, width = 4)


























———————————以下暂时没做————————————

##提取细胞间均值方差最大化合物
res_var<-as.data.frame(apply(res_mean_mat,1,var))
colnames(res_var)<-"var"
res_var$metabolite<-rownames(res_var)
write.csv(res_var,"res_var.csv",row.names = F,quote = F)
#提取前20
top20_var <- head(res_var[order(-res_var$var), ], 20)
top20_var <- top20_var$metabolite
Heatmap(res_mean_mat[top20_var,], 
        row_names_side = "right",  
        show_row_names = TRUE, 
        show_column_names = TRUE)

#提取前100
top100_var <- head(res_var[order(-res_var$var), ], 100)
top100_var <- top100_var$metabolite

Heatmap(res_mean_mat[top100_var,], 
        row_names_side = "right",  
        show_row_names = TRUE, 
        show_column_names = TRUE)




##提取热图注释file
annot_heatmap<-mets_db1[,c("super_class","sub_class","class","metabolite")]

##group_table$cell_type 没有Unkown 


#focus_mets<-c("Choline", "L-Glutamine", "Folic acid", "Cholesterol", "Iron", "9-cis-Retinoic acid", 
#              "all-trans-Retinoic acid", "13-cis-Retinoic acid")

##mebocost 代谢物均值计算结果，上个项目已经核查过，本数据集略过
#avg_met<-read.csv('../avg_met.csv',row.names=1,check.names = F)
# 将 avg_met 数据框与 mets_db 数据框合并，基于 HMDB_ID 列
#avg_met$rowname <- rownames(avg_met)
#avg_met <- merge(avg_met, mets_db, by.x = "rowname", by.y = "HMDB_ID", all.x = TRUE)
##metabolite设置为行名
#rownames(avg_met) <- avg_met$metabolite
# 删除多余的rowname列
#avg_met$rowname <- NULL
# 删除多余的metabolite列
#avg_met$metabolite <- NULL
#avg_met<-as.matrix(avg_met)
#Heatmap(avg_met[focus_mets,], 
#        row_names_side = "right",  
#        show_row_names = TRUE, 
#        show_column_names = TRUE)









##提取均值矩阵(res_mean_mat与avg_met完全一致)

##仅保留列名中有mean_case的列
res_mean_mat<-res_df[,grep("mean_case",colnames(res_df))]

#删除列名中的mean_case
colnames(res_mean_mat)<-gsub("mean_case","",colnames(res_mean_mat))

#3提取padj
res_padj_mat<-res_df[,grep("p.adjust",colnames(res_df))]
colnames(res_padj_mat)<-gsub("p.adjust","",colnames(res_padj_mat))

write.csv(res_mean_mat,"res_mean_mat.csv",row.names = T,quote = T)
write.csv(res_padj_mat,"res_padj_mat.csv",row.names = T,quote = T)

#画图

all_mets<-rownames(res_mean_mat)

annot_heatmap1<-annot_heatmap[annot_heatmap$metabolite %in% all_mets,]

rownames(annot_heatmap1) <- annot_heatmap1$metabolite
annot_heatmap1$metabolite <- NULL
annot_heatmap1$sub_class <- NULL

ha = rowAnnotation(df = annot_heatmap1)


Heatmap(res_mean_mat[all_mets,], 
        row_names_side = "right",  
        show_row_names = TRUE, 
        show_column_names = TRUE)
##注释行
Heatmap(res_mean_mat[all_mets,], 
        row_names_side = "right",  
        show_row_names = TRUE, 
        show_column_names = TRUE,
        right_annotation  = ha,
        split = annot_heatmap1$super_class)
#手动导出图片





