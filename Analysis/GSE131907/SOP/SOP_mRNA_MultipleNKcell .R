# 00.Init -----

{
  rm(list = ls())
  options(stringsAsFactors = FALSE)
  gc()
  options(future.globals.maxSize = 2000 * 1024^2)
  #options(bitmapType = "cairo")
  set.seed(16)
  num.cores = 12
  OUTPUTPDF = FALSE
  pre_process = FALSE
  SampleInfoTable = TRUE
  NoIntegration = FALSE #多样本一般需要整合，默认为FALSE
  HarmonyIntegration = TRUE
  RunSingleR = TRUE
  Test_nFeature = FALSE
}

###############SOP_pre_mRNA_Multiple.R#####################

# 01.Setup -----
suppressMessages({
  library(Seurat)
  library(ggplot2)
})

## Setup
{
  project_path =  "/Users/zhengboying/Documents/DMPJ0001-PGD2-byzheng/Analysis/GSE131907/"
  setwd(project_path)
  source("../../Function/PlotTheme.R")
}


# 02.Input Load -----
suppressMessages({
  library(reshape2)
})

## Setup
{
  # param1：ProjectName
  ProjectName <- "GSE131907"
  
  # param2：SampleSpecies
  SampleSpecies <- "Human"
  
  # param3：Input DateType
  DateType = "MEX"
  
  # param4：Inputdata path(h5 or CellRanger formats)
  CellRangerPath = paste0(project_path, "../../Inputdata/")
  h5Path = "/"
  
  # param5：sample_map
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
}  



# 03.Basic QC -----
suppressMessages({
  library(stringr)
})

## Setup
{
  source("../../Function/CellStat.R")
  ResultPath <- c("./NK_Cluster/") 
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  scObject<-readRDS("./Rds/Results/2.NK.Rds")
  scObject<-readRDS("./Rds/Results/5.NK_subtype.Rds")
}

############################################SOP_mRNA_Multiple.R#########################################



# 06.Basic analysis -----

## Load Object


suppressMessages({
  library(stringr)
})

## Setup
{
  ## FindVariableFeatures Parameters

  FindVariableFeatures_nfeatures <- c(2000)

  
  FindVariableFeatures_selection.method = "vst"
  
  ## ScaleData Parameters
  ScaleData_vars.to.regress = c("nCount_RNA", "nFeature_RNA", "CC.Difference")
  
  ## PCA Parameters
  PCA_npcs = 30
}

## 06.1 CellCycle -----
{
  
  if (str_to_title(SampleSpecies) == "Human") {
    load("../../Database/Human_CellCycleGene.rda")
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    load("../../Database/Mouse_CellCycleGene.rda")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  scObject <- NormalizeData(scObject, assay = "RNA", verbose = FALSE)
  scObject <- CellCycleScoring(scObject, set.ident = TRUE,
                               s.features = rownames(scObject)[rownames(scObject) %in% CellCycleGene$s.genes],
                               g2m.features = rownames(scObject)[rownames(scObject) %in% CellCycleGene$g2m.genes])
  
  ## CC.Difference is the difference between the S and G2M scores
  scObject$CC.Difference <- scObject$S.Score - scObject$G2M.Score
  
  Idents(scObject) <- "Phase"
  Idents(scObject) <- factor(Idents(scObject), levels = sort(levels(scObject)))
  scObject$Phase <- scObject@active.ident
}

## 06.2 Basic Seurat -----





for (n in FindVariableFeatures_nfeatures){
  ## FindVariableFeatures
  scObject <- FindVariableFeatures(scObject, assay = "RNA", selection.method = FindVariableFeatures_selection.method,
                                   nfeatures = n, verbose = FALSE)
  
  ## Remove bad genes
  if (str_to_title(SampleSpecies) == "Human") {
    VariableFeatures(scObject) <- grep(pattern = "^MT-", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^RP[SL]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^HB[ABD]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^IG[HKL]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "AC233755.1", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^TR[ABDG]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    VariableFeatures(scObject) <- grep(pattern = "^mt-", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Rp[sl]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Hb[abd]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Ig[hkl]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Tr[abdg]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  ## ScaleData
  scObject <- ScaleData(scObject, assay = "RNA", features =VariableFeatures(scObject),
                        vars.to.regress = ScaleData_vars.to.regress,
                        verbose = FALSE)
  
  ## PCA
  scObject <- RunPCA(scObject, assay = "RNA", npcs = PCA_npcs,
                     seed.use = 16, verbose = FALSE)
  
  ## Save PCA result
  saveRDS(scObject, file = paste0("./Rds/Results/3.NK_PCA_nFeature",n,".Rds"))
  }
}


# 07. Cluster -----
if (TRUE){
  scObject<-readRDS("./Rds/Results/3.NK_PCA_nFeature2000.Rds")
}
suppressMessages({
  library(harmony)
  library(dplyr)
  library(reshape2)
})

## Setup
{
  source("../../Function/ClusterStat.R")
  source("../../Function/ClusterStatPlot.R")
  source("../../Function/PlotTheme.R")
  ResultPath <- c("./NK_Cluster/")
  
  ## NoIntegration and HarmonyIntegration parameters
  dims = 1:30
  RNA_snn_res = c( 0.1,0.2,0.4,0.6)
  TSNE = FALSE
  UMAP = TRUE
  min.dist = 0.3
  
  ## Load data
  nFeature = 2000
  #scObject<-readRDS(paste0("./Rds/Seurat/3.PCA_nFeature",nFeature,".Rds"))
}


## 07.2 HarmonyIntegration -----
if (HarmonyIntegration == TRUE){
  
  ## Create result path
  dir.create(paste0(ResultPath,"HarmonyIntegration/"), showWarnings = FALSE, recursive = TRUE)
  
  ##
  scObject <- RunHarmony(scObject, assay.use = "RNA", group.by.vars = "Sample",
                         theta = 1, lambda = 0.5, max.iter = 20,
                         early_stop = -Inf, verbose = FALSE) # epsilon.harmony>early_stop
  
  scObject <- FindNeighbors(scObject, reduction = "harmony", dims = dims,
                            verbose = FALSE)
  if(UMAP == TRUE){
  scObject <- RunUMAP(scObject, assay = "RNA", reduction = "harmony", dims = dims,
                      n.components = 2, n.neighbors = 40, min.dist = min.dist,
                      reduction.name = "UMAP", reduction.key = "UMAP_",
                      seed.use = 16, verbose = FALSE)
  }
  
  if(TSNE == TRUE){
    scObject <- RunTSNE(scObject, assay = "RNA", reduction = "harmony", dims = dims,
                        dim.embed = 2,
                        reduction.name = "TSNE", reduction.key = "TSNE_",
                        seed.use = 16, verbose = FALSE)
  }

  


  
  
  ## Stat and plot
  for(i in RNA_snn_res) {
    scObject <- FindClusters(scObject, resolution = i, random.seed = 16, 
                             algorithm = 1, n.start = 10, n.iter = 10,verbose = FALSE)
    
    scObject@meta.data[[paste0("RNA_snn_res.", i)]] <- NULL
    scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony")]] <- scObject$seurat_clusters
    scObject$seurat_clusters <- NULL
    
    ## Plot
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "UMAP", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", "res.", sprintf("%0.1f", i), "_umap.png"),
             bg = "white", width = 8, height = 7)
    }
    
    if (TSNE == TRUE) {
      p1 <- DimPlot(scObject, reduction = "TSNE", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }

  }
  
}

## Save step 07 result
{
  saveRDS(scObject, file = paste0("./Rds/Results/4.NK_umap_tsne_nFeature",nFeature,".Rds"))
}


# 08. ResNubmer selection -----
if(FALSE){
  scObject<-readRDS(paste0("./Rds/Results/4.NK_umap_tsne_nFeature",nFeature,".Rds"))
}


## Setup
{
  source("../../Function/ResChoose.R")
  
  ## Create result path
  ResultPath <- c("./NK_Cluster/")
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  
  ## ResNubmer selection parameters
  RNA_snn_res = RNA_snn_res # c(0.1, 0.2, 0.4, 0.6) defult is 07.Cluster
  VarNumber = 15 ## 总群数量一般15个clusters足够
  KeyWord = "Harmony" #RAW
}

{
    ##
    ResNubmer <- ResChoose(object = scObject, KeyWord = KeyWord, VarNumber = VarNumber,
                           ReductionNumber = RNA_snn_res)
  message(paste0("When selecting close to ", VarNumber, " clusters, the corresponding resolution is ", ResNubmer))
}


# 09.CanonicalMarker -----

suppressMessages({
  library(Nebulosa)##plot_density
  library(stringr)
  library(RColorBrewer)
  library(tidyverse)
  library(ComplexHeatmap)
})

## Setup
{
  source("../../Function/Plot.R")
  ResultPath <- c("./NK_Cluster/")
  FeaturePlot = TRUE
  DotPlot = TRUE ## surpport Human only.
  plot_density = FALSE ## plot_density耗时较久，非交付需求可选择FALSE
  ResNubmer = 0.4 ## from step 08
  reduction = "UMAP" ##tSNE,UMAP,umap
}  

## Load CanonicalMarker  
{
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)  
  ##
  if (str_to_title(SampleSpecies) == "Human") {
    ## all markers,which need to be filtered after check results,
    CanonicalMarker<-read.csv("../../Database/CellType.annotdata/ref_nk_cell_type.csv",header = T)
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    source("../../Database/Mouse_CanonicalMarker.R")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
}  

## 09.1 CanonicalMarker DotPlot -----
{
  if (str_to_title(SampleSpecies) == "Human") {
    col_plates_init<-c("NK" = "grey","NK1"="#90E663","NK2" = "#39A093",
                       "B"="#90E663","Plasma"="#11E693","PGD2" = "brown",
                       "Mast" = "#63B3E6","Myeloid" = "#13F4FF","Monocyte" = "#78C4FF", 
                       "pDC" = "#13F4FF","DC" = "#78C4FF","Neutrophil" = "pink","Pericyte"= "white",
                       "Epithelium" = "#AB49F0", "Platelet" = "black","Megakaryocyte" = "orange",
                       "Endothelial"="#F5802C","MSCs" = "#E5202C","Stromal" = "#F5E02C","Fibroblast" = "#DEF02C")
    
    
    CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("RGS1","GZMH","CCL3","CREM","IL7R"),]
    
    
    Idents(scObject)<-paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
    ht<-anno_hematmap(obj = scObject,
                      markers = CanonicalMarker,
                      col_plates = col_plates_init,
                      cell_cluster = 9,#not used
                      group = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                      row_names_fontsize = 8,
                      seed = 12)
    ##save
    png(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot.png"), width =12, height = 16, units = "in", res = 600)
    draw(ht)
    dev.off()  
    pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot.pdf"), width =8, height = 6)
    draw(ht)
    dev.off()
    
    ## Filtered markers
    if(TRUE){
      #去掉CellTye =neutrophils,Pericytes,mregDCs,Lam,"RTMs",CMono"
      #CanonicalMarker<-CanonicalMarker[!CanonicalMarker$CellType %in% c("Pericytes","DC"),]#
      ##删掉一些不怎么表达的基因
      
      celltype<-c("0" ="KLRC2+/PTGDR+ NK",
                  "1" = "PTGDS+/PTGDR+ NK",##
                  "2" = "PTGDR- NK",
                  "3" = "PTGDR- NK",##
                  "4" = "PTGDR+ NK1",
                  "5" = "KLRC2+/PTGDR+ NK",
                  "6" = "PTGDR- NK",
                  "7" = "KLRC2+/PTGDR+ NK",
                  "8" = "DNAJB1+/PTGDR+ NK",
                  "9" = "PTGDR+ NK2",
                  "10" = "KLRC2+/PTGDR+ NK",
                  "11" = "KLRC2+/PTGDR+ NK",
                  "12" = "Unknown",
                  "13" = "KLRC2+/PTGDR+ NK")
                  
      
      col_plates_cluster<-c("#AF689B","#B7B577","#C87452","#32CD32","#00BFFF","#FFD700","#FFA500","#DFABCF","#11E693")
      celltype<-as.data.frame(celltype)
      celltype$RNA_snn_res.0.4Harmony <-rownames(celltype)
      
      ## Create celltype dataframe (all information for each cell)
      celltype_df <- data.frame(CELL = rownames(scObject@meta.data), RNA_snn_res.0.4Harmony = scObject@meta.data$RNA_snn_res.0.4Harmony)
      celltype_df <- merge(celltype_df, celltype, by = "RNA_snn_res.0.4Harmony", all.x = TRUE)
      row.names(celltype_df)<-celltype_df$CELL
      celltype_df$CELL<-NULL
      head(celltype_df)
      
      #Add metadata to Seurat object
      scObject <- AddMetaData(
        object = scObject,
        metadata = celltype_df
      )
      
      Idents(scObject) <- "celltype"
      
      CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("KLRF1","NKG7"),]
      
      ht<-anno_hematmap(obj = scObject,
                        markers = CanonicalMarker,
                        col_plates = col_plates_init,
                        cell_cluster = 9,#not used
                        group = "celltype",
                        row_names_fontsize = 8,
                        seed = 12)
      ##save
      png(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.png"), width =12, height = 16, units = "in", res = 600)
      draw(ht)
      dev.off()  
      pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.pdf"), width =6, height = 3.5)
      draw(ht)
      dev.off()
      
    }
  }
}



p1<-DimPlot(scObject, reduction = "UMAP", label = FALSE, pt.size = 0.05,
            cols = col_plates_cluster, group.by = "celltype",raster=FALSE)
ggsave(p1, file = paste0(ResultPath, "/celltype_umap.pdf"),
       bg = "white", width = 6, height = 5)


## 09.2 CanonicalMarker FeaturePlot ----
{
  ##
  if (str_to_title(SampleSpecies) == "Human") {
  CanonicalMarker<-c("EPCAM","krt19","CDH1",
                    "PTPRC","MKI67",
                    "LYZ","CD68",
                    "CD3D")
  
  
  }
  
  
  CanonicalMarker<-c(#"NCAM1","FCGR3A",
                "IL7R","GZMH","CCL3","CREM","RGS1",
                "IL32","KLRC2","ZNF90","MKI67","NFKBIA","DNAJB1","NR4A3")
  
  
  ##FEATURE PLOT
  {
    p1<-FeaturePlot(scObject, features = CanonicalMarker, reduction = "UMAP",
                    raster = FALSE, order = TRUE, cols = c("black", "#00ff00"),alpha = 0.4,
                    pt.size = 1, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    
    p1
    
    ggsave(p1, file = paste0(ResultPath, "_FeaturePlot_Density.png"),
           bg = "white", width = 14, height = 10)
    
    
    
    p1<-FeaturePlot(scObject, features = c("KLRC2","CX3CR1","DNAJB1","NR4A3"), reduction = "UMAP",
                    raster = FALSE, order = TRUE, cols = c("black", "#00ff00"),alpha = 0.4,
                    pt.size = 2, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    
    p1
    
    ggsave(p1, file = paste0(ResultPath, "4GENE_FeaturePlot_Density.png"),
           bg = "white", width = 10, height = 9)
    
    
    p1<-FeaturePlot(scObject, features = c("PTGDR","PTGDS","FCGR3A","NCAM1"), reduction = "UMAP",
                    raster = FALSE, order = TRUE, cols = c("black", "#00ff00"),alpha = 0.4,
                    pt.size = 2, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    
    ggsave(p1, file = paste0(ResultPath, "/FeaturePlot_PTGDR_PTGDS_FCGR3A_NCAM1.png"),
           bg = "white", width = 10, height = 9)
  }
  
  
 
  
  for (i in CanonicalMarker) {
    try({
      #p1 <- FeaturePlot(scObject, features = i, reduction = reduction,
         #               raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
        #                pt.size = 0.5, max.cutoff = "q99") &
       # PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      #ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot.png"),
      #       bg = "white", width = 8, height = 7)
      if (plot_density == TRUE){
        p1 <- plot_density(scObject, features = i, reduction = reduction, size = 1) &
          PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
        ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot_Density.png"),
               bg = "white", width = 8, height = 7)
      }
    })
  }
}



ResultPath <- c("./NK_Cluster/")





saveRDS(scObject, file = paste0("./Rds/Results/5.NK_subtype.Rds"))



##更新density plot


library(Nebulosa)
CanonicalMarker<-c(#"NCAM1","FCGR3A",
  "IL7R","GZMH","CCL3","CREM","RGS1",
  "IL32","KLRC2","ZNF90","MKI67","NFKBIA","DNAJB1","NR4A3")


features <- c("FCGR3A","NCAM1","PTGDR","PTGDS")
features1 <- c("KLRC2","CX3CR1","DNAJB1","NR4A3")
p1 <- plot_density(scObject, features = features, 
                   reduction = "UMAP", method = c( "wkde"), size = 1) &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
ggsave(p1, file = paste0(ResultPath, "FIG5B_FeaturePlot_Density.png"),
       bg = "white", width = 12, height = 9)

p2 <- plot_density(scObject, features = features1, 
                   reduction = "UMAP", method = c( "wkde"), size = 1) &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
ggsave(p2, file = paste0(ResultPath, "FIG5D_FeaturePlot_Density.png"),
       bg = "white", width = 12, height = 9)


p3 <- plot_density(scObject, features = CanonicalMarker, 
                   reduction = "UMAP", method = c( "wkde"), size = 1) &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
ggsave(p3, file = paste0(ResultPath, "FIGs7c_FeaturePlot_Density.png"),
       bg = "white", width = 16, height = 12)




