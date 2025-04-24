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

###############SOP_pre_mRNA_Multiple.R#####################

# 01.Setup -----

project_path =  "/Users/zhengboying/Documents/DMPJ0001-PGD2-byzheng/Analysis/GSE131907/"
setwd(project_path)


suppressMessages({
  library(Seurat)
  library(ggplot2)
})

# 02.Input Load -----

suppressMessages({
  library(reshape2)
})

{ 
  # param1：ProjectName
  ProjectName <- "GSE131907"
  # param2：sample_map
  sample_map <- c(
    "LUNG_N01"="GSM3827114",
    "LUNG_N06"="GSM3827115",
    "LUNG_N08"="GSM3827116",
    "LUNG_N09"="GSM3827117",
    "LUNG_N18"="GSM3827118",
    "LUNG_N19"="GSM3827119",
    "LUNG_N20"="GSM3827120",
    "LUNG_N28"="GSM3827121",
    "LUNG_N30"="GSM3827122",
    "LUNG_N31"="GSM3827123",
    "LUNG_N34"="GSM3827124",
    "LUNG_T06"="GSM3827125",
    "LUNG_T08"="GSM3827126",
    "LUNG_T09"="GSM3827127",
    "LUNG_T18"="GSM3827128",
    "LUNG_T19"="GSM3827129",
    "LUNG_T20"="GSM3827130",
    "LUNG_T25"="GSM3827131",
    "LUNG_T28"="GSM3827132",
    "LUNG_T30"="GSM3827133",
    "LUNG_T31"="GSM3827134",
    "LUNG_T34"="GSM3827135",
    "EBUS_06"="GSM3827136",##有转移的原位样本tL/B
    "EBUS_28"="GSM3827137",##有转移的原位样本tL/B
    #"EBUS_49"="GSM3827138",##有转移的原位样本tL/B,有红细胞污染
    "BRONCHO_58"="GSM3827139",##有转移的原位样本tL/B
    "EBUS_10"="GSM3827140",## Metastatic LN
    "BRONCHO_11"="GSM3827141", ## Metastatic LN
    "EBUS_12"="GSM3827142",## Metastatic LN
    "EBUS_13"="GSM3827143",## Metastatic LN
    "EBUS_15"="GSM3827144",## Metastatic LN
    "EBUS_19"="GSM3827145",## Metastatic LN
    "EBUS_51"="GSM3827146",## Metastatic LN
    "LN_01"="GSM3827147", ##Normal LN
    "LN_02"="GSM3827148", ##Normal LN
    "LN_03"="GSM3827149",
    "LN_04"="GSM3827150",
    "LN_05"="GSM3827151",
    "LN_06"="GSM3827152",
    "LN_07"="GSM3827153",
    "LN_08"="GSM3827154",
    "LN_11"="GSM3827155",
    "LN_12"="GSM3827156",## Normal LN
    "EFFUSION_06"="GSM3827157", ##MPE
    "EFFUSION_11"="GSM3827158", ##MPE
    "EFFUSION_12"="GSM3827159", ## MPE
    "EFFUSION_13"="GSM3827160", ## MPE
    "EFFUSION_64"="GSM3827161", ## MPE
    "NS_02"="GSM3827162",##脑转移
    "NS_03"="GSM3827163",
    "NS_04"="GSM3827164",
    "NS_06"="GSM3827165",
    "NS_07"="GSM3827166",
    "NS_12"="GSM3827167",
    "NS_13"="GSM3827168",
    "NS_16"="GSM3827169",
    "NS_17"="GSM3827170",
    "NS_19"="GSM3827171"##脑转移
  )
  SampleName <- basename(sample_map) 
  print(SampleName)
  
  # param3：SampleSpecies
  SampleSpecies <- "Human"
}  
source("../../Function/PlotTheme.R")

# 03.Load Data -----
if (TRUE){
  #scObject<-readRDS("./Rds/Seurat/4.1_filtercells_nFeature3000.Rds")
  ## readRDS
  scObject<-readRDS("./Rds/Results/1.celltype.Rds")
}

## Create result path
ResultPath <- c("./ploting_results/")
#dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)

#table(scObject$celltype)

#B    Endo     Epi   Fibro    Mast Myeloid      NK   Oligo     pDC       T 
#21903    2107   26144    4430    3485   36294   24489     655     861   56674 



# 04. Plot -----
## Fig1: FeaturePlot
{
  FigPath <- c("./ploting_results/Fig1")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  source("../../Function/ClusterStat.R")
  source("../../Function/ClusterStatPlot.R")
  
  ##Fig1A
  {
    celltype<-c("0" ="T",
                "1" = "Myeloid",##
                "2" = "Epi",
                "3" = "NK",##
                "4" = "B",
                "5" = "Myeloid",
                "6" = "Fibro",
                "7" = "Mast",
                "8" = "Endo",
                "9" = "pDC",
                "10" = "Epi",
                "11" = "Oligo",
                "12" = "T",
                "13" = "T",
                "14" = "T")
    
    col_plates_cluster<-c("#AF689B","#B7B577","#C87452","#32CD32","#00BFFF","#FFD700","#FFA500","#DFABCF","grey","#11E693")
    celltype<-as.data.frame(celltype)
    celltype$RNA_snn_res.0.1Harmony <-rownames(celltype)
    ## Create celltype dataframe (all information for each cell)
    celltype_df <- data.frame(CELL = rownames(scObject@meta.data), RNA_snn_res.0.1Harmony = scObject@meta.data$RNA_snn_res.0.1Harmony)
    celltype_df <- merge(celltype_df, celltype, by = "RNA_snn_res.0.1Harmony", all.x = TRUE)
    row.names(celltype_df)<-celltype_df$CELL
    celltype_df$CELL<-NULL
    head(celltype_df)
    
    #Add metadata to Seurat object
    scObject <- AddMetaData(
      object = scObject,
      metadata = celltype_df
    )
    
    Idents(scObject) <- "celltype"
    
    p1<-DimPlot(scObject, reduction = "UMAP", label = FALSE, pt.size = 0.05,
                cols = col_plates_cluster, group.by = "celltype",raster=FALSE)
    ggsave(p1, file = paste0(FigPath, "/Fig1A_celltype_umap.pdf"),
           bg = "white", width = 6, height = 5)

  }

  ## Save
  # saveRDS(scObject, file = "./Rds/Results/1.celltype.Rds")
  
  ##Fig1B
  {
    Idents(scObject) <- "celltype"
    ##filter badmarkers
    suppressMessages({
      library(RColorBrewer)
      library(tidyverse)
      library(ComplexHeatmap)
    })
    source("../../Function/Plot.R")
    SampleSpecies = "Human"
    CanonicalMarker<-read.csv("../../Database/CellType.annotdata/ref_cell_type.csv",header = T)
    
    #去掉CellTye =neutrophils,Pericytes,mregDCs,Lam,"RTMs",CMono"
    #CanonicalMarker<-CanonicalMarker[!CanonicalMarker$CellType %in% c("Pericytes","DC"),]#
    ##删掉一些不怎么表达的基因
    CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CMA1","KRT1","IL1RAPL1","PF4","FOLR2","FABP3",
                                                                    "PLA2G2D","SH2D1B","RGS5","IGLL1","AMPD1","PPBP","TBX21",
                                                                    "TNFRSF4","BATF","FOXP3","CTLA4","LAIR2","EPCAP"),]
    
    CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CD34","CDH5","CLEC9A","CLEC10A","CD1C","IL1RL1",
                                                                    "CTSG","HDC","CPA3","S100A8","S100A9","COL1A1","COL3A1",
                                                                    "SDC1","PAX5","IGLL5","CD38","IGHD","LAMP3", "ACTA2","TAGLN",
                                                                    "CSF3R","MKI67","JCHAIN","CD19","PECAM1","PTPRC","GATA2","FCGR3A",
                                                                    "KLRF1","CLDN11","SSR4","XBP1","CD4","CD8A","CD8B"
    ),]
    
    if (str_to_title(SampleSpecies) == "Human") {
      col_plates_init<-c("NK" = "grey","T"="#90E663","Leukocyte" = "#39A093",
                         "B"="#90E663","Plasma"="#11E693","Oligo" = "brown",
                         "Mast" = "#63B3E6","Myeloid" = "#13F4FF","Monocyte" = "#78C4FF", 
                         "pDC" = "#13F4FF","DC" = "#78C4FF","Neutrophil" = "pink","Pericyte"= "white",
                         "Epithelium" = "#AB49F0", "Platelet" = "black","Megakaryocyte" = "orange",
                         "Endothelial"="#F5802C","MSCs" = "#E5202C","Stromal" = "#F5E02C","Fibroblast" = "#DEF02C"
      )
      ht<-anno_hematmap(obj = scObject,
                        markers = CanonicalMarker,
                        col_plates = col_plates_init,
                        cell_cluster = 9,#not used
                        group = "celltype",
                        row_names_fontsize = 8,
                        seed = 12)
      ##save
      pdf(file = paste0(FigPath, "/Fig1B_CanonicalMarker_DotPlot.pdf"), width =6, height =6)
      draw(ht)
      dev.off()  
    }
    
    ##FigS2B_FeaturePlot
    CanonicalMarker_1<-c("CD79A",#B
                         "CLDN5","VWF", #Endothelial
                         "EPCAM","CDH1",#Epithelium
                         "DCN",#fibroblast
                         "MOG",#Oligo
                         "TPSAB1","KIT",#Mast
                         "LYZ","CD68",#Mye
                         "KLRB1","NKG7",#NK
                         "LILRA4", #pDC
                         "CD3D","CD3E")#T
    
    p1 <- FeaturePlot(scObject, features = CanonicalMarker_1, reduction = "UMAP",
                      raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
                      pt.size = 0.5, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    
    ggsave(p1, file = paste0(FigPath, "/FigS2B_CanonicalMarker_FeaturePlot.png"),
           bg = "white", width = 8, height = 23)
  }
  
  ##Fig1C
  {
    Idents(scObject) <- "celltype"
    OUTPUTPDF = TRUE
    ##location
    {
      ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                       RowType = "celltype", ColType = "location")
      write.csv(ClusterStatResult, file = paste0(FigPath, "/Fig1B-1_celltype_Sample_Stat.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Plot base on number
      ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                      RowType = "celltype", ColType = "Sample",
                      PlotPrefixName = "Fig1B-1_celltype",
                      RowTypeAngle = FALSE, ColTypeAngle = TRUE,
                      ResultPath = paste0(FigPath,"/"))
      
      
      ## Calculate percent
      ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                       RowType = "celltype", ColType = "location")
      write.csv(ClusterStatResult, file = paste0(FigPath, "/Fig1B-1_celltype_Sample_Stat_Percent.csv"),
                quote = TRUE, row.names = FALSE)
    }
  }
  
  ##Fig1D
  {  
    ## stage.1 Tumor Lung
    { 
      temp_obj <- subset(scObject, subset = location == "Tumor Lung")
      table(temp_obj@meta.data$location)
      #47157 cells
      ClusterStatResult <- ClusterStat(object = temp_obj, percent = FALSE,
                                       RowType = "celltype", ColType = "stage.1")
      write.csv(ClusterStatResult, file = paste0(FigPath, "/Fig1B-2_celltype_Sample_Stat.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Plot base on number
      ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                      RowType = "celltype", ColType = "Sample",
                      PlotPrefixName = "Fig1B-2_celltype",
                      RowTypeAngle = FALSE, ColTypeAngle = TRUE,
                      ResultPath = paste0(FigPath,"/"))
      
      
      ## Calculate percent
      ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                       RowType = "celltype", ColType = "stage.1")
      write.csv(ClusterStatResult, file = paste0(FigPath, "/Fig1B-2_celltype_Sample_Stat_Percent.csv"),
                quote = TRUE, row.names = FALSE)
    }

  }
  
  ## save
  #saveRDS(temp_obj, file = "./Rds/Results/1.tumorlung.Rds")
  
  ##Fig1E/F
  {
    Idents(scObject) <- "celltype"
    temp_obj <- subset(scObject, subset = location %in% c("Tumor Lung","Normal Lung"))
    table(temp_obj@meta.data$location)
    
    
    p1<-DimPlot(temp_obj, reduction = "UMAP", label = FALSE, pt.size = 0.005, alpha = 0.1,
                cols = c("#B7B577","#AF689B"), group.by = "location",raster=FALSE)+NoLegend()
    p1
    ggsave(p1, file = paste0(FigPath, "/Fig1E_celltype_umap_location.pdf"),
           bg = "white", width = 5, height = 5)
    #F
    temp_obj <- subset(scObject, subset = location == "Tumor Lung")
    table(temp_obj@meta.data$location)
    
    #col_plates_cluster<-c("#AF689B","#B7B577","#C87452","#32CD32","#00BFFF","#FFD700","#FFA500","#DFABCF","grey","#11E693")
    
    p1<-DimPlot(temp_obj, reduction = "UMAP", label = FALSE, pt.size = 0.005,alpha = 0.1,
                cols = c("#AF689B","#B7B577"), group.by = "stage.1",raster=FALSE)+NoLegend()
    p1
    ggsave(p1, file = paste0(FigPath, "/Fig1F_celltype_umap_stage.1.pdf"),
           bg = "white", width = 5, height = 5)
    
  }
}

## Fig2: MetsPlot
{
  FigPath <- c("./ploting_results/Fig2")
  #dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  
  ##Fig2A(TOP200)
  {
    MetsObject<-readRDS("./Rds/Mets/1.metobj.Rds")
    Idents(MetsObject) <- "celltype"
    col_plates_cluster<-c("grey","#00BFFF","#32CD32","#B7B577","#DFABCF","#FFA500","#FFD700","#11E693","#C87452","#AF689B")
  
    p1<-DimPlot(MetsObject, reduction = "umap", label = FALSE, pt.size = 0.05,
              cols = col_plates_cluster, group.by = "celltype",raster=FALSE)
    p1
    ggsave(p1, file = paste0(FigPath, "/Fig2A_umap_mets.pdf"),
          bg = "white", width = 6, height = 5)))
  }
  
  ## Fig2B/C/D/E
  {
    #saveRDS(VariableFeatures(MetsObject),"./Rds/Mets/2.metobj_VariableFeatures_200.Rds")
    file.copy("./03_mebocost/details/count_dot_plot.pdf", paste0(FigPath, "/Fig2B_count_dot_plot.pdf"))
    file.copy("./03.1_met_diff/valcano_padjust.pdf", paste0(FigPath, "/Fig2C_valcano_padjust.pdf"))
    file.copy("./03.1_met_diff/maker_mets_heatmap.pdf", paste0(FigPath, "/Fig2D_maker_mets_heatmap.pdf"))
    file.copy("./03_mebocost/details/commu_dotmap_6mets.pdf", paste0(FigPath, "/Fig2E_commu_dotmap_6mets.pdf"))
  }

}

## Fig3: PGD2
{
  FigPath <- c("./ploting_results/Fig3")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  
  ##Fig3A/B/C/D
  {
    file.copy("./03_mebocost/details/FlowPlot_PGD2-1.pdf", paste0(FigPath, "/Fig3A_FlowPlot_PGD2-1.pdf"))
    file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_NTM.pdf", paste0(FigPath, "/Fig3B_Prostaglandin_D2_vlnplot_NTM.pdf"))
    file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_stage1.pdf", paste0(FigPath, "/Fig3C_Prostaglandin_D2_vlnplot_stage1.pdf"))
    file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_DIFF.pdf", paste0(FigPath, "/Fig3D_Prostaglandin_D2_vlnplot_DIFF.pdf"))
  }
  
  ## Fig3E
  {
    scObject<-readRDS("./Rds/Results/1.celltype.Rds")
    
    ## 
    p1<-FeaturePlot(scObject, features = c("PTGDS","HPGDS"), reduction = "UMAP",
                raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
                pt.size = 0.5, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ##
    ggsave(p1, file = paste0(FigPath, "/Fig3E_FeaturePlot.pdf"),
           bg = "white", width = 12, height = 5)
    ggsave(p1, file = paste0(FigPath, "/Fig3E_FeaturePlot.png"),
           bg = "white", width = 12, height = 5)
  }
  
   ## fig3G
  {
    ##scObject<-readRDS("./Rds/Results/1.celltype.Rds")
    
    ## 
    p1<-FeaturePlot(scObject, features = c("PTGDR","SLCO2A1"), reduction = "UMAP",
                    raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
                    pt.size = 0.5, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ##
    ggsave(p1, file = paste0(FigPath, "/Fig3G_FeaturePlot.pdf"),
           bg = "white", width = 12, height = 5)
    ggsave(p1, file = paste0(FigPath, "/Fig3G_FeaturePlot.png"),
           bg = "white", width = 12, height = 5)
  }
  
  
  {
    #Idents(scObject)<-"celltype"
 ## 不画图，提取scObject$RNA$data以下基因("PTGDR","PTGDR2","SLCO2A1","HPGDS","PTGDS")的表达矩阵到新的matrx
    #matrix<-scObject@assays$RNA$data[c("PTGDR","PTGDR2","SLCO2A1","HPGDS","PTGDS"),]
    #df<-as.data.frame(t(matrix))
    #merged_df <- merge(scObject@meta.data, df, by = "row.names")
    
    #write
    #write.csv(merged_df,paste0(FigPath, "/genes_table.csv"),row.names = F,quote = F)
    merged_df<-read.csv(paste0(FigPath, "/genes_table.csv"))
    
  }
  
  ## 计算各组均值
  {
    #计算PTGDS的NK中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanNK<-df %>% group_by(location) %>% summarise(NK.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Fibro中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Fibro"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanFibro<-df %>% group_by(location) %>% summarise(Fibro.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Oligo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Oligo"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanOligo<-df %>% group_by(location) %>% summarise(Oligo.PTGDS=mean(PTGDS))
    
    
    ## 计算PTGDS pDC中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("pDC"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanpDC<-df %>% group_by(location) %>% summarise(pDC.PTGDS=mean(PTGDS))
    
    
    
    
    
    
    #计算PTGDR的T中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("T"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanT<-df %>% group_by(location) %>% summarise(T.PTGDR=mean(PTGDR))
    
    
    #计算PTGDR的nk中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanNK_<-df %>% group_by(location) %>% summarise(NK.PTGDR=mean(PTGDR))
    
    
    ## 计算SLCO2A1的Endo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Endo"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanEndo<-df %>% group_by(location) %>% summarise(Endo.SLCO2A1=mean(SLCO2A1))
    
    
    
    #计算PTGDS的SLCO2A1的Epi中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Epi"))
    df<-subset(df, location %in% c("Normal Lung","Tumor Lung","Metastasis"))
    meanEpi<-df %>% group_by(location) %>% summarise(Epi.SLCO2A1=mean(SLCO2A1))
    
    
    
    
    #将多个均值合并
    mean<-merge(meanNK,meanFibro,by="location")
    mean<-merge(mean,meanOligo,by="location",all=TRUE)
    mean<-merge(mean,meanpDC,by="location")
    mean<-merge(mean,meanT,by="location")
    mean<-merge(mean,meanNK_,by="location")
    mean<-merge(mean,meanEndo,by="location")
    mean<-merge(mean,meanEpi,by="location")
    
    row.names(mean)<-mean$location
    
    mean1<-mean
    
    
    ht<-pheatmap::pheatmap(mean[,2:9],
                           cluster_rows =F,
                           scale = "column",
                           cluster_cols = F,show_rownames = T,
                           show_colnames = T ) #
    
    ggsave(paste0(FigPath, "/Fig3F_heatmap.pdf"),ht,width=10,height=3.5,units="in")
    
    
    
    FigPath <- c("./ploting_results/Fig3")
    
    
    
    
    
    
    
    #计算PTGDS的NK中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanNK<-df %>% group_by(stage.1) %>% summarise(NK.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Fibro中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Fibro"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanFibro<-df %>% group_by(stage.1) %>% summarise(Fibro.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Oligo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Oligo"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanOligo<-df %>% group_by(stage.1) %>% summarise(Oligo.PTGDS=mean(PTGDS))
    
    
    ## 计算PTGDS pDC中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("pDC"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanpDC<-df %>% group_by(stage.1) %>% summarise(pDC.PTGDS=mean(PTGDS))
    
    
    
    
    
    
    #计算PTGDR的T中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("T"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanT<-df %>% group_by(stage.1) %>% summarise(T.PTGDR=mean(PTGDR))
    
    
    #计算PTGDR的nk中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanNK_<-df %>% group_by(stage.1) %>% summarise(NK.PTGDR=mean(PTGDR))
    
    
    ## 计算SLCO2A1的Endo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Endo"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanEndo<-df %>% group_by(stage.1) %>% summarise(Endo.SLCO2A1=mean(SLCO2A1))
    
    
    
    #计算PTGDS的SLCO2A1的Epi中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Epi"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanEpi<-df %>% group_by(stage.1) %>% summarise(Epi.SLCO2A1=mean(SLCO2A1))
    
    
    
    
    #将多个均值合并
    mean<-merge(meanNK,meanFibro,by="stage.1")
    mean<-merge(mean,meanOligo,by="stage.1",all=TRUE)
    mean<-merge(mean,meanpDC,by="stage.1")
    mean<-merge(mean,meanT,by="stage.1")
    mean<-merge(mean,meanNK_,by="stage.1")
    mean<-merge(mean,meanEndo,by="stage.1")
    mean<-merge(mean,meanEpi,by="stage.1")
    
    
    row.names(mean)<-mean$stage.1
    
    mean2<-mean
    
    ht1<-pheatmap::pheatmap(mean[,2:9],
                            cluster_rows =F,
                            scale = "column",
                            cluster_cols = F,show_rownames = T,
                            show_colnames = T ) #
    
    ggsave(paste0(FigPath, "/Fig3F_heatmap-1.pdf"),ht,width=10,height=3.5,units="in")
    
    #画折线图以stage
    
    
    
    #计算PTGDS的NK中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanNK<-df %>% group_by(Differentiated_state) %>% summarise(NK.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Fibro中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Fibro"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanFibro<-df %>% group_by(Differentiated_state) %>% summarise(Fibro.PTGDS=mean(PTGDS))
    
    
    #计算PTGDS的Oligo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Oligo"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanOligo<-df %>% group_by(Differentiated_state) %>% summarise(Oligo.PTGDS=mean(PTGDS))
    
    
    ## 计算PTGDS pDC中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("pDC"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanpDC<-df %>% group_by(Differentiated_state) %>% summarise(pDC.PTGDS=mean(PTGDS))
    
    
    
    
    
    
    #计算PTGDR的T中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("T"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanT<-df %>% group_by(Differentiated_state) %>% summarise(T.PTGDR=mean(PTGDR))
    
    
    #计算PTGDR的nk中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("NK"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanNK_<-df %>% group_by(Differentiated_state) %>% summarise(NK.PTGDR=mean(PTGDR))
    
    
    ## 计算SLCO2A1的Endo中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Endo"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanEndo<-df %>% group_by(Differentiated_state) %>% summarise(Endo.SLCO2A1=mean(SLCO2A1))
    
    
    
    #计算PTGDS的SLCO2A1的Epi中的 按照location的均值
    df<-subset(merged_df, celltype %in% c("Epi"))
    df<-subset(df, location %in% c("Tumor Lung"))
    meanEpi<-df %>% group_by(Differentiated_state) %>% summarise(Epi.SLCO2A1=mean(SLCO2A1))
    
    
    
    
    #将多个均值合并
    mean<-merge(meanNK,meanFibro,by="Differentiated_state")
    mean<-merge(mean,meanOligo,by="Differentiated_state",all=TRUE)
    mean<-merge(mean,meanpDC,by="Differentiated_state")
    mean<-merge(mean,meanT,by="Differentiated_state")
    mean<-merge(mean,meanNK_,by="Differentiated_state")
    mean<-merge(mean,meanEndo,by="Differentiated_state")
    mean<-merge(mean,meanEpi,by="Differentiated_state")
    
    
    row.names(mean)<-mean$Differentiated_state
    
    #去掉Differentiated_state=na
    mean<-subset(mean, Differentiated_state != "na")
  }
  ##Fig3F
  {
    #画折线图以stage
    
    meantable<-rbind(mean1[,2:9],mean2[,2:9],mean[,2:9])
    
    write.csv(meantable, file = paste0(FigPath, "/Fig3F_heatmap-3.tsv"),  quote = F, row.names = F)
    
    
    meantable<-t(meantable)
    
    desired_row_order <- c("Fibro.PTGDS", "pDC.PTGDS", "NK.PTGDS", "Oligo.PTGDS", 
                           "NK.PTGDR", "T.PTGDR", "Endo.SLCO2A1", "Epi.SLCO2A1")
    
    desired_col_order <- c("Normal Lung", "Tumor Lung", "Metastasis", "early stage", 
                           "advanced stage", "WD", "MD", "PD")
    
    # 根据自定义顺序重新排序数据框或矩阵的列
    meantable_reordered <- meantable[desired_row_order,]
    meantable_reordered <- meantable_reordered[,desired_col_order]
    
    ## 
    ht3<-pheatmap::pheatmap(meantable_reordered,
                            cluster_rows =F,
                            scale = "row",
                            cluster_cols = F,show_rownames = T,
                            color = colorRampPalette(c("#00008B", "#0000CD", "#FFFFFF", "#FF0000", "#8B0000"))(30),
                            show_colnames = T ) #
    
    ggsave(paste0(FigPath, "/Fig3F_heatmap.pdf"),ht3,width=5,height=10,units="in")
    
  }
  
  
}


## Fig4: Function in cell
{
  FigPath <- c("./ploting_results/Fig4")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  library(AUCell)
  library(ggplot2)
  library(ggpubr)
  library(ggthemes)
  library(tidyr)
  scObject<-readRDS("./Rds/Results/1.celltype.Rds")
  Idents(scObject) <- "celltype"
  NKcell<-subset(scObject, idents = c("NK"))
  scObject<-NKcell
  
  {
    cytotoxicity_sets<-c("GZMA", "GZMB", "GZMH", "GZMM", "GZMK", "GNLY", "PRF1", "CTSW") 
    inflammatory_sets<-c("CCL2", "CCL3", "CCL4", "CCL5", "CXCL10", "CXCL9", "IL1B", "IL6","IL7", "IL15","IL18")
    general_stress_sets<-c("BAG3", "CALU", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HIF1A", "HSP90AA1", "HSP90AB1",
                           "HSP90B1", "HSPA1A", "HSPA1B", "HSPA6", "HSPB1", "HSPH1", "IER2", "JUN", "JUNB", "NFKBIA",
                           "NFKBIZ","RGS2", "SLC2A3", "SOCS3", "UBC", "ZFAND2A", "ZFP36", "ZFP36L1")
    
    
    ##
    cells_rankings <- AUCell_buildRankings(scObject@assays$RNA$data,  plotStats=TRUE,  nCores=4)
    
    ##aucells
    AUCcells_cytotoxicity <- AUCell_calcAUC(cytotoxicity_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs1 <- t(as.data.frame(getAUC(AUCcells_cytotoxicity)))
    ##改写列名
    colnames(aucs1) <- "cytotoxicity_score"
    AUCcells_inflammatory <- AUCell_calcAUC(inflammatory_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs2 <- t(as.data.frame(getAUC(AUCcells_inflammatory)))
    colnames(aucs2) <- "inflammatory_score"
    AUCcells_general_stress <- AUCell_calcAUC(general_stress_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs3 <- t(as.data.frame(getAUC(AUCcells_general_stress)))
    colnames(aucs3) <- "general_stress_score"
    
    aucs <- cbind(aucs1,aucs2,aucs3)
    
    ##
    extracted_rows <- scObject@assays$RNA$data[rownames(scObject@assays$RNA$data) %in% c("PTGDS","PTGDR"), ]
    df<-as.data.frame(t(extracted_rows))
    
    ##PTGDR>0,+
    df$PTGDR_level <- ifelse(df$PTGDR>0, "+", "-")
    df$PTGDS_level <- ifelse(df$PTGDS>0, "+", "-")
    
    ##
    merge_df <- cbind(aucs,df)
    
    write.csv(merge_df, file = paste0(FigPath, "/Fig4_NK_AUCell_score.csv"), quote = TRUE, row.names = FALSE)
    
    df_long <- gather(merge_df, key = "variable", 
                      value = "value", "cytotoxicity_score",
                      "inflammatory_score", "general_stress_score")
    
    
    
    
    ##
    scObject<-Tcell
    cells_rankings <- AUCell_buildRankings(scObject@assays$RNA$data,  plotStats=TRUE,  nCores=4)
    
    ##aucells
    AUCcells_cytotoxicity <- AUCell_calcAUC(cytotoxicity_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs1 <- t(as.data.frame(getAUC(AUCcells_cytotoxicity)))
    ##改写列名
    colnames(aucs1) <- "cytotoxicity_score"
    AUCcells_inflammatory <- AUCell_calcAUC(inflammatory_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs2 <- t(as.data.frame(getAUC(AUCcells_inflammatory)))
    colnames(aucs2) <- "inflammatory_score"
    AUCcells_general_stress <- AUCell_calcAUC(general_stress_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs3 <- t(as.data.frame(getAUC(AUCcells_general_stress)))
    colnames(aucs3) <- "general_stress_score"
    
    aucs <- cbind(aucs1,aucs2,aucs3)
    
    ##
    extracted_rows <- scObject@assays$RNA$data[rownames(scObject@assays$RNA$data) %in% c("PTGDS","PTGDR"), ]
    df<-as.data.frame(t(extracted_rows))
    df$PTGDS <-NULL
    ##PTGDR>0,+
    df$PTGDR_level <- ifelse(df$PTGDR>0, "+", "-")
    #df$PTGDS_level <- ifelse(df$PTGDS>0, "+", "-")
    
    ##
    merge_df <- cbind(aucs,df)
    
    write.csv(merge_df, file = paste0(FigPath, "/Fig4_T_AUCell_score.csv"), quote = TRUE, row.names = FALSE)
    merge_df<-read.csv(paste0(FigPath, "/Fig4_T_AUCell_score.csv"))
    df_long <- gather(merge_df, key = "variable", 
                      value = "value", "cytotoxicity_score",
                      "inflammatory_score", "general_stress_score")
    
  }
  ##Fig4C
  {
    
    p1<-ggplot(merge_df, aes(x = PTGDR, y = PTGDS)) +
      geom_point(aes(), size = 3, alpha = 0.4) +
      geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
      #theme_prism() +
      xlab("PTGDR") +
      ylab("PTGDS") +
      theme(legend.position = "none")
    ##增加文本
    #geom_text(aes(label = "cor=0.0645, p-value < 2.2e-16"), 
    #           x = Inf, y = Inf, vjust = 1.5, hjust = 1.5)
    
    
    cor.test(df$PTGDR, df$PTGDS, method = "pearson")
    
    ggsave(p1, file = paste0(FigPath, "/Fig4A_PTGDR_PTGDScor.pdf"),
           bg = "white", width = 8, height = 8)
  }
  
  ##Fig4D-NK
  {
    p1<-ggplot(df_long,#
               aes(x =PTGDS_level, y = value, fill = PTGDS_level)) +
      geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white") +
      #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
      #scale_x_discrete(limits = order) +
      #theme_prism() +
      theme(axis.text.x = element_text( hjust = 1)) +
      #xlab("Stage") +
      ylab("AUCell_score") +
      scale_fill_manual(values = c("+"= "#187C3A",
                                   "-" = "#FFA500"
      ))+
      #spit wrap celltype,一行一个
      facet_wrap(~variable, ncol = 1)+
      stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
    p1
    
    ##添加显著性
    p2<-p1 +
      geom_signif(test = wilcox.test, 
                  size = 0.2, 
                  comparisons = list(c("+", "-")),
                  vjust = 0.3) +
      theme(legend.position = "none")
    p2
    
    p3<-ggplot(df_long,#
               aes(x =PTGDR_level, y = value, fill = PTGDR_level)) +
      geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white") +
      #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
      #scale_x_discrete(limits = order) +
      theme_prism() +
      theme(axis.text.x = element_text( hjust = 1)) +
      #xlab("Stage") +
      ylab("AUCell_score") +
      scale_fill_manual(values = c("+"= "#187C3A",
                                   "-" = "#FFA500"
      ))+
      #spit wrap celltype,一行一个
      facet_wrap(~variable, ncol = 1)+
      stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
    p3
    
    ##添加显著性
    p4<-p3 +
      geom_signif(test = wilcox.test, 
                  size = 0.2, 
                  comparisons = list(c("+", "-")),
                  vjust = 0.3) +
      theme(legend.position = "none")
    p4
    
    p5<-p2+p4
    
    ##
    ggsave(p5, file = paste0(FigPath, "/Fig4B_PTGDS_PTGDR_AUCell_score.pdf"),
           bg = "white", width = 10, height = 10)
  }
  
  ##Fig4d-T
  {
    {
      p1<-ggplot(df_long,#
                 aes(x =PTGDR_level, y = value, fill = PTGDR_level)) +
        geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
        geom_boxplot(width = 0.1, fill = "white") +
        #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
        #scale_x_discrete(limits = order) +
        #theme_prism() +
        theme(axis.text.x = element_text( hjust = 1)) +
        #xlab("Stage") +
        ylab("AUCell_score") +
        scale_fill_manual(values = c("+"= "#187C3A",
                                     "-" = "#FFA500"
        ))+
        #spit wrap celltype,一行一个
        facet_wrap(~variable, ncol = 1)+
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
      p1
      
      ##添加显著性
      p2<-p1 +
        geom_signif(test = t.test, 
                    size = 0.2, 
                    comparisons = list(c("+", "-")),
                    vjust = 0.3) +
        theme(legend.position = "none")
      p2
      
      
      ##
      ggsave(p2, file = paste0(FigPath, "/Fig4G_T_PTGDR_AUCell_score.pdf"),
             bg = "white", width = 4, height = 10)
    }
  }
  
  ## Fig4C
  {
  ##scObject<-readRDS("./Rds/Results/1.celltype.Rds")
  
  ##   
  la<-c( "LDHB", "LDHA") 
  p1<-FeaturePlot(scObject, features = la, reduction = "UMAP",
                  raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
                  pt.size = 0.5, max.cutoff = "q99") &
    PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
  ##
  ggsave(p1, file = paste0(FigPath, "/Fig4C_FeaturePlot.pdf"),
         bg = "white", width = 12, height = 5)
  ggsave(p1, file = paste0(FigPath, "/Fig4C_FeaturePlot.png"),
         bg = "white", width = 12, height = 5)
  
  
  ##FigS6J,,la的senser有些有表达，不是很好
  
  {
    FigPath <- c("./ploting_results/FigS6")
    
    la_senser<-c( "HCAR1", "SLC16A1","SLC16A7","SLC16A3","SLC2A1") 
    p1<-FeaturePlot(scObject, features = la_senser, reduction = "UMAP",split.by = "location",
                    raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =5,
                    pt.size = 0.5, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    p1
    #
    ggsave(p1, file = paste0(FigPath, "/FigS6J_FeaturePlot.png"),
           bg = "white", width = 20, height = 5)
  }
  }
  
  
  ## Fig4D
  {
    library(AUCell)
    scObject<-readRDS("./Rds/Results/1.celltype.Rds")
    Idents(scObject) <- "celltype"
    Endocell<-subset(scObject, idents = c("Endo"))
    saveRds(Endocell,"./Rds/Results/2.Endo.Rds")
    scObject<-Endocell
    #la<-c("LDHAL6A", "LDHB","LDHC", "LDHA", "LDHAL6B")
    
    #HMDB0001311 D-Lactic acid
    #HMDB0000190		L-Lactic acid
    la<-c( "LDHB", "LDHA")
    cells_rankings <- AUCell_buildRankings(scObject@assays$RNA$data,  plotStats=TRUE,  nCores=4)
    AUCcells<- AUCell_calcAUC(la, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
    aucs <- t(as.data.frame(getAUC(AUCcells)))
    colnames(aucs) <- "la_score"
    
    ##
    extracted_rows <- scObject@assays$RNA$data[rownames(scObject@assays$RNA$data) %in% c("SLCO2A1","CD8A"), ]
    df<-as.data.frame(t(extracted_rows))
    df$SLCO2A1_level <- ifelse(df$SLCO2A1>0, "+", "-")
    df$CD8A<-NULL
    
    ##
    merge_df <- cbind(df,aucs)
    
    ##save
    write.csv(merge_df, file = paste0(FigPath, "/Fig4C_la.scoretable.csv"), quote = TRUE, row.names = FALSE)
    
    df<-merge_df
    
    ##
    ## 画图
    p1<-ggplot(df,#
               aes(x =SLCO2A1_level, y = la_score, fill = SLCO2A1_level)) +
      geom_violin(trim = TRUE, scale = "width",color= NA,alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white") +
      #geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
      #scale_x_discrete(limits = order) +
      #theme_prism() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      #xlab("Stage") +
      ylab("la_score") +
      scale_fill_manual(values = c("+"= "#B7B577",
                                   "-" = "#FFA500"
      ))+
      #spit wrap celltype,一行一个
      #facet_wrap(~celltype, ncol = 1)+
      stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") 
    p1
    
    ##添加显著性
    p2<-p1 +
      geom_signif(test = wilcox.test, 
                  size = 0.2, 
                  comparisons = list(c("+", "-")),
                  vjust = 0.3) +
      theme(legend.position = "none")
    p2
    
    ggsave(p2, file = paste0(FigPath, "/Fig4C_La.score.pdf"),
           bg = "white", width = 8, height = 8)
    
  }
  
}

##Fig5
{
  scObject<-readRDS("./Rds/Results/2.NK.Rds")
  
  
}
  

  

  
  
  
  
  







## FigS1: ClusterPlot
{
  Idents(scObject) <- "RNA_snn_res.0.1Harmony"
  FigPath <- c("./ploting_results/FigS1")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  file.copy("./5.reCluster/HarmonyIntegration/Phase_umap.png", paste0(FigPath,"/FigS1A_Phase_umap.png"))
  file.copy("./5.reCluster/HarmonyIntegration/Sample_umap.png",paste0(FigPath,"/FigS1B_Sample_umap.png"))
  file.copy("./5.reCluster/HarmonyIntegration/res.0.1_umap.png", paste0(FigPath,"/FigS1C_res.0.1_umap.png"))
  file.copy("./5.reCluster/HarmonyIntegration/res.0.1_Sample_Stat_Percent.png",paste0(FigPath,"/FigS1D_res.0.1_Sample_Stat_Percent.png"))
}


## FigS2: FeaturePlot
{
  Idents(scObject) <- "RNA_snn_res.0.1Harmony"
  FigPath <- c("./ploting_results/FigS2")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  FigS2A_raw<-file.copy("./5.reCluster/RNA_snn_res.0.1Harmony_CanonicalMarkerDotPlot.png", paste0(FigPath, "/RNA_snn_res.0.1Harmony_CanonicalMarkerDotPlot.png"))
  
  ##filter badmarkers
  suppressMessages({
    library(RColorBrewer)
    library(tidyverse)
    library(ComplexHeatmap)
  })
  source("../../Function/Plot.R")
  SampleSpecies = "Human"
  CanonicalMarker<-read.csv("../../Database/CellType.annotdata/ref_cell_type.csv",header = T)
  
  #去掉CellTye =neutrophils,Pericytes,mregDCs,Lam,"RTMs",CMono"
  #CanonicalMarker<-CanonicalMarker[!CanonicalMarker$CellType %in% c("Pericytes","DC"),]#
  ##删掉一些不怎么表达的基因
  CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CMA1","KRT1","IL1RAPL1","PF4","FOLR2","FABP3",
                                                      "PLA2G2D","SH2D1B","RGS5","IGLL1","AMPD1","PPBP","TBX21",
                                                      "TNFRSF4","BATF","FOXP3","CTLA4","LAIR2","EPCAP"),]
  
  CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CD34","CDH5","CLEC9A","CLEC10A","CD1C","IL1RL1",
                                                                  "CTSG","HDC","CPA3","S100A8","S100A9","COL1A1","COL3A1",
                                                                  "SDC1","PAX5","IGLL5","CD38","IGHD","LAMP3", "ACTA2","TAGLN",
                                                                  "CSF3R","MKI67"
                                                                  ),]
    
  if (str_to_title(SampleSpecies) == "Human") {
    col_plates_init<-c("NK" = "grey","T"="#90E663","Leukocyte" = "#39A093",
                       "B"="#90E663","Plasma"="#11E693","Oligo" = "brown",
                       "Mast" = "#63B3E6","Myeloid" = "#13F4FF","Monocyte" = "#78C4FF", 
                       "pDC" = "#13F4FF","DC" = "#78C4FF","Neutrophil" = "pink","Pericyte"= "white",
                       "Epithelium" = "#AB49F0", "Platelet" = "black","Megakaryocyte" = "orange",
                       "Endothelial"="#F5802C","MSCs" = "#E5202C","Stromal" = "#F5E02C","Fibroblast" = "#DEF02C"
    )
    ResNubmer=0.1
    Idents(scObject)<-paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
    ht<-anno_hematmap(obj = scObject,
                      markers = CanonicalMarker,
                      col_plates = col_plates_init,
                      cell_cluster = 9,#not used
                      group = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                      row_names_fontsize = 8,
                      seed = 12)
    ##save
    pdf(file = paste0(FigPath, "/FigS2A_CanonicalMarker_DotPlot.pdf"), width =7, height = 9)
    draw(ht)
    dev.off()  
  }
  
  ##FigS2B_FeaturePlot
  CanonicalMarker_1<-c("CD79A",#B
                       "CLDN5","VWF", #Endothelial
                       "EPCAM","CDH1",#Epithelium
                       "DCN",#fibroblast
                       "MOG",#Oligo
                       "TPSAB1","KIT",#Mast
                       "LYZ","CD68",#Mye
                       "KLRB1","NKG7",#NK
                       "LILRA4", #pDC
                       "CD3D","CD3E")#T
  
      p1 <- FeaturePlot(scObject, features = CanonicalMarker_1, reduction = "UMAP",
                        raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
                        pt.size = 0.5, max.cutoff = "q99") &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      
      ggsave(p1, file = paste0(FigPath, "/FigS2B_CanonicalMarker_FeaturePlot.png"),
             bg = "white", width = 8, height = 23)
    }
    


## FigS3: ClusterStatPlot
{ 
  Idents(scObject) <- "celltype"
  FigPath <- c("./ploting_results/FigS3")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  source("../../Function/ClusterStat.R")
  source("../../Function/ClusterStatPlot.R")
  
  
  ## FigS3A/FigS3B
  {
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = "celltype", ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(FigPath, "/Fig1B_celltype_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(PlotData = ClusterStatResult,PDF = TRUE,
                    RowType = "celltype", ColType = "Sample",
                    PlotPrefixName = "FigS3_celltype",
                    RowTypeAngle = TRUE, ColTypeAngle = TRUE,
                    ResultPath = paste0(FigPath,"/"))
    
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = "celltype", ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(FigPath, "/FigS3_celltype_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }
}



## FigS4
{ 
  FigPath <- c("./ploting_results/FigS4")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  
  file.copy("./03_mebocost/details/mebocost_eventnum.pdf", paste0(FigPath, "/FigS4A_mebocost_eventnum.pdf"))
  file.copy("./03_mebocost/details/FlowPlot.pdf", paste0(FigPath, "/FigS4B_FlowPlot.pdf"))
  file.copy("./03_mebocost/details/commu_dotmap.pdf", paste0(FigPath, "/Fig4C_commu_dotmap.pdf"))
  
}



## FigS5
{ 
  FigPath <- c("./ploting_results/FigS5")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  
  file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_all.pdf", paste0(FigPath, "/FigS5A_Prostaglandin_D2_vlnplot_all.pdf"))
  file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_stage2.pdf", paste0(FigPath, "/FigS5B_Prostaglandin_D2_vlnplot_stage2.pdf"))
  file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_DIFF_ALL.pdf", paste0(FigPath, "/Fig5C_Prostaglandin_D2_vlnplot_DIFF_ALL.pdf"))
  
  ##FigS5D
  p1<-FeaturePlot(scObject, features = c("PTGDR2"), reduction = "UMAP",
                  raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                  pt.size = 0.5, max.cutoff = "q99") &
    PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
  ##
  ggsave(p1, file = paste0(FigPath, "/FigS5D_FeaturePlot.pdf"),
         bg = "white", width = 5.5, height = 5)
  ggsave(p1, file = paste0(FigPath, "/FigS5D_FeaturePlot.png"),
         bg = "white", width = 5.5, height = 5)
  
}



## FigS6
{ 
  FigPath <- c("./ploting_results/FigS6")
  dir.create(FigPath, showWarnings = FALSE, recursive = TRUE)
  
  #file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_all.pdf", paste0(FigPath, "/FigS5A_Prostaglandin_D2_vlnplot_all.pdf"))
  #file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_stage2.pdf", paste0(FigPath, "/FigS5B_Prostaglandin_D2_vlnplot_stage2.pdf"))
  #file.copy("./03.1_met_diff/Prostaglandin_D2_vlnplot_DIFF_ALL.pdf", paste0(FigPath, "/Fig5C_Prostaglandin_D2_vlnplot_DIFF_ALL.pdf"))
  
}



}




{
  Idents(scObject) <- "celltype"
  NKcell<-subset(scObject, idents = c("NK"))
  
  Tcell<-subset(scObject, idents = c("T"))
  Endocell<-subset(scObject, idents = c("Endo"))
  
  
  saveRDS(NKcell,"./Rds/Results/2.NK.Rds")
  
  saveRDS(Tcell,"./Rds/Results/2.T.Rds")
  
  saveRDS(Endocell,"./Rds/Results/2.Endo.Rds")
}







VlnPlot(scObject, features = c("PTGDR","PTGDR2","SLCO2A1","HPGDS","PTGDS"), pt.size = 0) &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())

VlnPlot(subset(scObject,idents = c("pDC","Fibro")), features = c("PTGDS"), pt.size = 0,split.by = "location") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


FeaturePlot(scObject, features = c("NCAM1","FCGR3A"), reduction = "UMAP",
            raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
            pt.size = 0.5, max.cutoff = "q99") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())



FeaturePlot(temp_obj, features = c("PTGDS"), reduction = "UMAP",split.by = "stage.2",
            raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,
            pt.size = 0.5, max.cutoff = "q99") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())

VlnPlot(temp_obj, features = c("PTGDS"), pt.size = 0,split.by = "stage.1") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


VlnPlot(subset(scObject,idents = c("pDC","Fibro","NK")), features = c("PTGDS"), pt.size = 0,split.by = "location") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


VlnPlot(subset(temp_obj,idents = c("Mast")), features = c("HPGDS"), pt.size = 0,split.by = "stage.2") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


  

