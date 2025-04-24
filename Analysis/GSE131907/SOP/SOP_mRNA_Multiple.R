# 00.Init -----

{
  rm(list = ls())
  options(stringsAsFactors = FALSE)
  gc()
  options(future.globals.maxSize = 2000 * 1024^2)
  # options(bitmapType = "cairo")
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
  
if(pre_process == TRUE){  
  # Input data：expression data
  # Extract the required results in cellranger
  # LibraryType has to one form the list: "mRNA", "Flex", "mRNA_TCR", "mRNA_TCR_BCR".
  source("../../Function/LoadmRNA.R")
  dir.create("./Rds/Seurat/", showWarnings = FALSE, recursive = TRUE)
  
  ## DateType has to one form the list: "MEX, "h5", "MGI".
  scObject <- LoadmRNA(SampleName = SampleName, DateType = DateType,
                       CellRangerPath = CellRangerPath,
                       Project = ProjectName,
                       h5Path =  h5Path,
                       min.cells = 3, min.features = 1,
                       names.field = 1, names.delim = "-")
  
}


# 03.Basic QC -----
suppressMessages({
  library(stringr)
})

## Setup
{
  source("../../Function/CellStat.R")
  ResultPath <- c("./1.QC/") 
}

if(pre_process == TRUE){
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)

  ## Calculate the percentage of mt/rib/hgb/GenesPerUMI
  if (str_to_title(SampleSpecies) == "Human") {
    scObject <- PercentageFeatureSet(scObject, pattern = "^MT-", col.name = "percent.mt")
    scObject <- PercentageFeatureSet(scObject, pattern = "^RP[SL]", col.name = "percent.rib")
    scObject <- PercentageFeatureSet(scObject, pattern = "^HB[ABD]", col.name = "percent.hgb")
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    scObject <- PercentageFeatureSet(scObject, pattern = "^mt-", col.name = "percent.mt")
    scObject <- PercentageFeatureSet(scObject, pattern = "^Rp[sl]", col.name = "percent.rib")
    scObject <- PercentageFeatureSet(scObject, pattern = "^Hb[abd]", col.name = "percent.hgb")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  scObject$GenesPerUMI <- log10(scObject$nFeature_RNA)/log10(scObject$nCount_RNA)
  
  ## Basic_Plot
  for (i in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rib", "percent.hgb", "GenesPerUMI")) {
    p1 <- VlnPlot(scObject, features = i, pt.size = 0.05, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Basic_", i, "_ViolinPlot_point.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
    
    p1 <- VlnPlot(scObject, features = i, pt.size = 0, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Basic_", i, "_ViolinPlot_nopoint.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
  }
  
  {
    p1 <- VlnPlot(scObject, features = "percent.mt", pt.size = 0.05, y.max = 30, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Basic_percent.mt_ViolinPlot_point_ymax30.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
    
    p1 <- VlnPlot(scObject, features = "percent.mt", pt.size = 0, y.max = 30, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Basic_percent.mt_ViolinPlot_nopoint_ymax30.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
  }
  
  {
    p1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5,
                         shuffle = TRUE, raster = FALSE, seed = 16) &
      scale_color_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath, "Basic_nCount_RNA-nFeature_RNA_Correlation.png"),
           bg = "white", width = 8, height = 7)
    
    p1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5,
                         shuffle = TRUE, raster = FALSE, seed = 16) &
      scale_color_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath, "Basic_nCount_RNA-percent.mt_Correlation.png"),
           bg = "white", width = 8, height = 7)
  }
  
  ## Cell information statistics
  Basic_scObjectStat <- CellStat(object = scObject, Type = "Sample", assay.use = "RNA")
  colnames(Basic_scObjectStat) <- c("Sample", paste0("Basic_", colnames(Basic_scObjectStat)[2:6]))
  print(Basic_scObjectStat)

}

## Save the basic object
if(pre_process == TRUE){
  saveRDS(scObject, file = "./Rds/Seurat/1.Basic.Rds")
}

############################################SOP_mRNA_Multiple.R#########################################

## Load Basic Seurat object
if(pre_process == FALSE){
  scObject<-readRDS("./Rds/Seurat/1.Basic.Rds")
  ## remove sample (optional)
  # scObject<-subset(scObject, Sample != "GSM3827138")
}


# 04.FilterLowQuality -----
suppressMessages({
  library(stringr)
  library(DoubletFinder)
  ResultPath <- c("./1.QC/") 
})

## Setup
{
  ResultPath <- "./1.QC/"
  source("../../Function/CellFilter.R")
  source("../../Function/DoubletFinder_k.R")
  source("../../Function/CellStat.R")
  
  ## parameters
  DoubletFinderPlot = TRUE 
  DoubletFinderPlotRds = FALSE
}

## 04.1 FilterLowQuality -----
{
  ## confidence interval 90% 1.645, 95% 1.96, 99% 2.575
  scObject <- CellFilter(object = scObject, Type = "Sample",
                         MinCount = 500, MinFeature = 200,MinGenesPerUMI = 0.8,
                         Sigma = 2.575, SigmaTop = TRUE, SigmaBottom = TRUE,
                         MinPercent.mt = NULL, MaxPercent.mt = 10,
                         MaxPercent.hgb = 1, MaxPercent.rib = 60
                         )
}
   
## 04.2 FilterDoublet -----
{
  ## Create result path
  dir.create(paste0(ResultPath, "DoubletFinder/"), showWarnings = FALSE, recursive = TRUE)
  dir.create("./Rds/DoubletFinder/", showWarnings = FALSE, recursive = TRUE)
  
  ## Get DoubletFinder result(拆开分析啦)
  
  temp_scObject<-subset(scObject, idents = SampleName[1:10])
  scObject_DoubletFinder_10 <- DoubletFinder_k(object = temp_scObject, Sample = "Sample", PCs = 1:20,
                                            DoubletRatePer = 0.005, num.cores = num.cores,
                                            plot.output = DoubletFinderPlot, rds.output = DoubletFinderPlotRds,
                                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                                            rds.dir = "./Rds/DoubletFinder/")
  saveRDS(scObject_DoubletFinder_10, file = "./Rds/DoubletFinder/DoubletFinder_10.Rds")
  
  temp_scObject<-subset(scObject, idents = SampleName[11:20])
  scObject_DoubletFinder_20 <- DoubletFinder_k(object = temp_scObject, Sample = "Sample", PCs = 1:20,
                                            DoubletRatePer = 0.005, num.cores = num.cores,
                                            plot.output = DoubletFinderPlot, rds.output = DoubletFinderPlotRds,
                                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                                            rds.dir = "./Rds/DoubletFinder/")
  saveRDS(scObject_DoubletFinder_20, file = "./Rds/DoubletFinder/DoubletFinder_20.Rds")
  
  temp_scObject<-subset(scObject, idents = SampleName[21:30])
  scObject_DoubletFinder_30 <- DoubletFinder_k(object = temp_scObject, Sample = "Sample", PCs = 1:20,
                                            DoubletRatePer = 0.005, num.cores = num.cores,
                                            plot.output = DoubletFinderPlot, rds.output = DoubletFinderPlotRds,
                                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                                            rds.dir = "./Rds/DoubletFinder/")
  saveRDS(scObject_DoubletFinder_30, file = "./Rds/DoubletFinder/DoubletFinder_30.Rds")
  
  temp_scObject<-subset(scObject, idents = SampleName[31:40])
  scObject_DoubletFinder_40 <- DoubletFinder_k(object = temp_scObject, Sample = "Sample", PCs = 1:20,
                                            DoubletRatePer = 0.005, num.cores = num.cores,
                                            plot.output = DoubletFinderPlot, rds.output = DoubletFinderPlotRds,
                                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                                            rds.dir = "./Rds/DoubletFinder/")
  saveRDS(scObject_DoubletFinder_40, file = "./Rds/DoubletFinder/DoubletFinder_40.Rds")
 
  temp_scObject<-subset(scObject, idents = SampleName[41:57])
  scObject_DoubletFinder_57 <- DoubletFinder_k(object = temp_scObject, Sample = "Sample", PCs = 1:20,
                                            DoubletRatePer = 0.005, num.cores = num.cores,
                                            plot.output = DoubletFinderPlot, rds.output = DoubletFinderPlotRds,
                                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                                            rds.dir = "./Rds/DoubletFinder/") 
  saveRDS(scObject_DoubletFinder_57, file = "./Rds/DoubletFinder/DoubletFinder_57.Rds")
  
  scObject_DoubletFinder <- rbind(scObject_DoubletFinder_10, scObject_DoubletFinder_20, 
                                  scObject_DoubletFinder_30, scObject_DoubletFinder_40, scObject_DoubletFinder_57)
  saveRDS(scObject_DoubletFinder, file = "./Rds/DoubletFinder/DoubletFinder.Rds")
  
  ## Remove Doublet
  {
    scObject <- AddMetaData(scObject, scObject_DoubletFinder)
    Idents(scObject) <- "DoubletFinder"
    scObject <- subset(scObject, idents = c("Singlet"))
    scObject$DoubletFinder <- scObject@active.ident
    Idents(scObject) <- "Sample"
    scObject
  }
}

## 04.3 Recreate SeuratObject -----
{
  ## Filtered min.features==0 features
  scObject <- CreateSeuratObject(counts = scObject@assays$RNA$counts,
                                 min.cells = 3, min.features = 1,
                                 meta.data = FetchData(scObject, vars = c("Sample")))
  Idents(scObject) <- "Sample"
  scObject$Sample <- scObject@active.ident
}

## 04.4 FilteredPlot -----
{
  Idents(scObject) <- "Sample"
  
  ## Calculate the percentage of percent.mt/percent.rib/percent.hgb/GenesPerUMI
  if (str_to_title(SampleSpecies) == "Human") {
    scObject <- PercentageFeatureSet(scObject, pattern = "^MT-", col.name = "percent.mt")
    scObject <- PercentageFeatureSet(scObject, pattern = "^RP[SL]", col.name = "percent.rib")
    scObject <- PercentageFeatureSet(scObject, pattern = "^HB[ABD]", col.name = "percent.hgb")
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    scObject <- PercentageFeatureSet(scObject, pattern = "^mt-", col.name = "percent.mt")
    scObject <- PercentageFeatureSet(scObject, pattern = "^Rp[sl]", col.name = "percent.rib")
    scObject <- PercentageFeatureSet(scObject, pattern = "^Hb[abd]", col.name = "percent.hgb")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  scObject$GenesPerUMI <- log10(scObject$nFeature_RNA)/log10(scObject$nCount_RNA)
  
  ## FilteredPlot
  for (i in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rib", "percent.hgb", "GenesPerUMI")) {
    p1 <- VlnPlot(scObject, features = i, pt.size = 0.05, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Filter_", i, "_ViolinPlot_point.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
    
    p1 <- VlnPlot(scObject, features = i, pt.size = 0, raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "Filter_", i, "_ViolinPlot_nopoint.png"),
           bg = "white", width = (0.6*nlevels(scObject) + 2), height = 7)
  }
  
  {
    p1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5,
                         shuffle = TRUE, raster = FALSE, seed = 16) &
      scale_color_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath, "Filter_nCount_RNA-nFeature_RNA_Correlation.png"),
           bg = "white", width = 8, height = 7)
    
    p1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5,
                         shuffle = TRUE, raster = FALSE, seed = 16) &
      scale_color_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath, "Filter_nCount_RNA-percent.mt_Correlation.png"),
           bg = "white", width = 8, height = 7)
  }
  
  ## Cell information statistics
  Filter_scObjectStat <- CellStat(object = scObject, Type = "Sample", assay.use = "RNA")
  colnames(Filter_scObjectStat) <- c("Sample", paste0("Filter_", colnames(Filter_scObjectStat)[2:6]))
  
  scObjectStat <- merge(x = Basic_scObjectStat, y = Filter_scObjectStat, by = "Sample")
  scObjectStat$percent <- as.numeric(scObjectStat$Filter_CellNumber)/as.numeric(scObjectStat$Basic_CellNumber)
  scObjectStat$percent <- scales::label_percent(accuracy = 0.01)(scObjectStat$percent)
  
  scObjectStat <- scObjectStat[, c(1, 2, 7, 12, 3, 8, 4, 9, 5, 10, 6, 11)]
  rownames(scObjectStat) <- scObjectStat$Sample
  scObjectStat <- scObjectStat[c(SampleName, "Total"), ]
  write.csv(scObjectStat, file = paste0(ResultPath, "CellStat.csv"),
            quote = TRUE, row.names = FALSE)
}


# 05.AddSampleInfo -----
suppressMessages({
  library(dplyr)
  library(stringr)
  library(reshape2)
})

## 05.1 Check SampleInfo -----
if (SampleInfoTable == TRUE){
  SampleInfo <- read.csv("../../meta.data.csv", header = TRUE)
  head(SampleInfo)
  ## Extracted SampleInfo$ID %in% SampleName Information
  SampleInfo <- SampleInfo[SampleInfo$ID %in% SampleName, ]
  
  ## Check all SampleName with SampleInfo
  if (dim(SampleInfo)[1] == length(SampleName)) {
    print(SampleName)
    message(paste0("The meta.data.csv file could provide all ", dim(SampleInfo)[1] ," Sample Info."))
  }else{
    ## print SampleName which is missing in SampleInfo
    print(paste(SampleName[!SampleName %in% SampleInfo$ID]," is missing.", sep = "\n"))
    stop(paste0("The meta.data.csv file could provide ", length(SampleName)-dim(SampleInfo)[1] ," Sample Info."))
  }
}

## 05.2 Add SampleNInfo to SeuratObject -----
if (SampleInfoTable == TRUE){
  ## Create CellSampleInfo dataframe (all information for each cell)
  CellSampleInfo <- data.frame(CELL = rownames(scObject@meta.data), Sample = scObject@meta.data$Sample)
  
  CellSampleInfo <- merge(CellSampleInfo, SampleInfo, by.x = "Sample",by.y = "ID", all.x = TRUE)
  row.names(CellSampleInfo)<-CellSampleInfo$CELL
  CellSampleInfo$CELL<-NULL
  head(CellSampleInfo)
  
  #Add metadata to Seurat object
  scObject <- AddMetaData(
    object = scObject,
    metadata = CellSampleInfo
  )
  
  Idents(scObject) <- "Sample"
}

## Save step 04 and 05 result
## if step 05 is True, the Filter.Rds should be saved with SampleInfo in object meta.data
{
  saveRDS(scObject, file = "./Rds/Seurat/2.Filter.Rds")
}


# 06.Basic analysis -----

## Load Object
if(FALSE){
  scObject<-readRDS("./Rds/Seurat/2.Filter.Rds")
  #scObject <- subset(scObject, cells = sample(colnames(scObject), 12786))
}

suppressMessages({
  library(stringr)
})

## Setup
{
  ## FindVariableFeatures Parameters
  if (Test_nFeature == TRUE){
    FindVariableFeatures_nfeatures <- c(1500,2000,3000,5000)
  }else{
    FindVariableFeatures_nfeatures <- c(2000)
  }
  
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
                                   nfeatures = FindVariableFeatures_nfeatures, verbose = FALSE)
  
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
  scObject <- ScaleData(scObject, assay = "RNA", features = VariableFeatures(scObject),
                        vars.to.regress = ScaleData_vars.to.regress,
                        verbose = FALSE)
  
  ## PCA
  scObject <- RunPCA(scObject, assay = "RNA", npcs = PCA_npcs,
                     seed.use = 16, verbose = FALSE)
  
  ## Save PCA result
  saveRDS(scObject, file = paste0("./Rds/Seurat/3.PCA_nFeature",n,".Rds"))
  }
}


# 07. Cluster -----
if (TRUE){
  scObject<-readRDS("./Rds/Seurat/3.PCA_nFeature2000.Rds")
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
  ResultPath <- c("./2.Cluster/")
  
  ## NoIntegration and HarmonyIntegration parameters
  dims = 1:30
  RNA_snn_res = c( 0.1)
  TSNE = FALSE
  UMAP = TRUE
  min.dist = 0.5
  
  ## Load data
  nFeature = 3000
  #scObject<-readRDS(paste0("./Rds/Seurat/3.PCA_nFeature",nFeature,".Rds"))
}

## 07.1 NoIntegration -----
if (NoIntegration == TRUE){
  ## Create result path
  dir.create(paste0(ResultPath,"NoIntegration/"), showWarnings = FALSE, recursive = TRUE)
  
  ## RUN UMAP and tSNE
  scObject <- FindNeighbors(scObject, reduction = "pca", dims = dims,
                            verbose = FALSE)
  if(UMAP == TRUE){
    scObject <- RunUMAP(scObject, assay = "RNA", reduction = "pca", dims = dims,
                        n.components = 2, n.neighbors = 40, min.dist = min.dist,
                        reduction.name = "umap", reduction.key = "umap_",
                        seed.use = 16, verbose = FALSE)
  }
  
  if(TSNE == TRUE){
    scObject <- RunTSNE(scObject, assay = "RNA", reduction = "pca", dims = dims,
                        dim.embed = 2,
                        reduction.name = "tSNE", reduction.key = "tSNE_",
                        seed.use = 16, verbose = FALSE)
  }
  
  
  ## Plot mix umap and tsne
  for (i in c("Sample", "Phase")) {
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "umap", group.by = i,
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/",i, "_umap.png"),
             bg = "white", width = 8, height = 7)
    }
    
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "tSNE", group.by = i,
                  label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                  pt.size = 0.5, label.size = 5, seed = 16) &
      labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
      theme(legend.title = element_blank())
    if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
    ggsave(p1, file = paste0(ResultPath, "NoIntegration/",i, "_tsne.png"),
           bg = "white", width = 8, height = 7)
    }
  }
  
  ##
  for(i in RNA_snn_res) {
    scObject <- FindClusters(scObject, resolution = i, random.seed = 16, verbose = FALSE)
    
    scObject@meta.data[[paste0("RNA_snn_res.", i)]] <- NULL
    scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW")]] <- scObject$seurat_clusters
    scObject$seurat_clusters <- NULL
    
    ## Plot
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "umap", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/","res.", sprintf("%0.1f", i), "_umap.png"),
             bg = "white", width = 8, height = 7)
    
    }
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "tSNE", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"),
                  label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                  pt.size = 0.5, label.size = 5, seed = 16) &
      scale_color_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath,"NoIntegration/", "res.", sprintf("%0.1f", i), "_tsne.png"),
           bg = "white", width = 8, height = 7)
    }
    
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "NoIntegration/", "res.", sprintf("%0.1f", i), "_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                    RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample",
                    PlotPrefixName = paste0("res.", sprintf("%0.1f", i)),
                    RowTypeAngle = FALSE, ColTypeAngle = TRUE)
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath,"NoIntegration/", "res.", sprintf("%0.1f", i), "_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }
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

  
  ## Plot mix umap and tsne
  
  for (i in c("Sample", "Phase")) {
    if(UMAP == TRUE){
        p1 <- DimPlot(scObject, reduction = "UMAP", group.by = i,
                      label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                      pt.size = 0.5, label.size = 5, seed = 16) &
          labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.title = element_blank())
        if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
        ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", i, "_umap.png"),
               bg = "white", width = 8, height = 7)
    }
    
    if (TSNE == TRUE) {
       p1 <- DimPlot(scObject, reduction = "TSNE", group.by = i,
                    label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", i, "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }

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
    
    
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                    RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample",
                    PlotPrefixName = paste0("res.", sprintf("%0.1f", i)),
                    RowTypeAngle = FALSE, ColTypeAngle = TRUE,
                    ResultPath = paste0(ResultPath, "HarmonyIntegration/"))
    
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }
  
}

## Save step 07 result
{
  saveRDS(scObject, file = paste0("./Rds/Seurat/4.umap_tsne_nFeature",nFeature,".Rds"))
}


# 08. ResNubmer selection -----
if(FALSE){
  scObject<-readRDS(paste0("./Rds/Seurat/4.umap_tsne_nFeature",nFeature,".Rds"))
}


## Setup
{
  source("../../Function/ResChoose.R")
  
  ## Create result path
  ResultPath <- c("./3.res.Marker/")
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
  ResultPath <- c("./3.res.Marker/")
  FeaturePlot = TRUE
  DotPlot = TRUE ## surpport Human only.
  plot_density = FALSE ## plot_density耗时较久，非交付需求可选择FALSE
  ResNubmer = 0.1 ## from step 08
  reduction = "UMAP" ##tSNE,UMAP,umap
}  

## Load CanonicalMarker  
{
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)  
  ##
  if (str_to_title(SampleSpecies) == "Human") {
    ## all markers,which need to be filtered after check results,
    CanonicalMarker<-read.csv("../../Database/CellType.annotdata/ref_cell_type.csv",header = T)
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    source("../../Database/Mouse_CanonicalMarker.R")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
}  

## 09.1 CanonicalMarker DotPlot -----
{
  if (str_to_title(SampleSpecies) == "Human") {
    col_plates_init<-c("NK" = "grey","T"="#90E663","Leukocyte" = "#39A093",
                       "B"="#90E663","Plasma"="#11E693","Oligo" = "brown",
                       "Mast" = "#63B3E6","Myeloid" = "#13F4FF","Monocyte" = "#78C4FF", 
                       "pDC" = "#13F4FF","DC" = "#78C4FF","Neutrophil" = "pink","Pericyte"= "white",
                       "Epithelium" = "#AB49F0", "Platelet" = "black","Megakaryocyte" = "orange",
                       "Endothelial"="#F5802C","MSCs" = "#E5202C","Stromal" = "#F5E02C","Fibroblast" = "#DEF02C")
    
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
    pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot.pdf"), width =12, height = 16)
    draw(ht)
    dev.off()
    
    ## Filtered markers
    if(TRUE){
      #去掉CellTye =neutrophils,Pericytes,mregDCs,Lam,"RTMs",CMono"
      #CanonicalMarker<-CanonicalMarker[!CanonicalMarker$CellType %in% c("Pericytes","DC"),]#
      ##删掉一些不怎么表达的基因
      CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CMA1","KRT1","IL1RAPL1","PF4","FOLR2","FABP3",
                                                                      "PLA2G2D","SH2D1B","RGS5","IGLL1","AMPD1","PPBP","TBX21",
                                                                      "TNFRSF4","BATF","FOXP3","CTLA4","LAIR2","EPCAP"),]
      
      CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CD34","CDH5","CLEC9A","CLEC10A","CD1C","IL1RL1",
                                                                      "CTSG","HDC","CPA3","S100A8","S100A9","COL1A1","COL3A1"),]
      
      ht<-anno_hematmap(obj = scObject,
                        markers = CanonicalMarker,
                        col_plates = col_plates_init,
                        cell_cluster = 9,#not used
                        group = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                        row_names_fontsize = 8,
                        seed = 12)
      ##save
      png(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.png"), width =12, height = 16, units = "in", res = 600)
      draw(ht)
      dev.off()  
      pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.pdf"), width =12, height = 16)
      draw(ht)
      dev.off()
      
    }
  }
}

## 09.2 CanonicalMarker FeaturePlot ----
{
  ##
  if (str_to_title(SampleSpecies) == "Human") {
  CanonicalMarker<-c("EPCAM","krt19","CDH1",
                    "PTPRC","MKI67",
                    "LYZ","CD68",
                    "CD3D")
  }
  
  for (i in CanonicalMarker) {
    try({
      p1 <- FeaturePlot(scObject, features = i, reduction = reduction,
                        raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                        pt.size = 0.5, max.cutoff = "q99") &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot.png"),
             bg = "white", width = 8, height = 7)
      if (plot_density == TRUE){
        p1 <- plot_density(scObject, features = i, reduction = reduction, size = 1) &
          PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
        ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot_Density.png"),
               bg = "white", width = 8, height = 7)
      }
    })
  }
}







# 10. SingleR -----
suppressMessages({
  library(SingleR)
  library(stringr)
  library(dplyr)
  library(reshape2)
})

## Setup
if(RunSingleR == TRUE){
  source("../../Function/ResChoose.R")
  source("../../Function/ClusterStat.R")
  source("../../Function/ClusterStatPlot.R")
  
  ## Create result path
  ResultPath <- c("./4.SingleR/")
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  SampleSpecies = SampleSpecies
  ResNubmer = ResNubmer## from step 08
  UMAP = FALSE
  UMAP_reduction = "UMAP" ## "UMAP" is from afer HarmonyIntegration,"umap" is from NoIntegration 
  TSNE = TRUE
  TSNE_reduction = "TSNE" ## "TSNE" is from afer HarmonyIntegration,"tSNE" is from NoIntegration 
}

## 10.1 run SingleR
if(RunSingleR == TRUE){
  ## Load ref_SingleR data
  if (str_to_title(SampleSpecies) == "Human") {
    load("../../Database/human_singleR.rda")
    AnnoResult <- SingleR(test = LayerData(scObject, assay = "RNA", layer = "counts"),
                          ref = human_singleR, labels =human_singleR$label.main,
                          clusters = scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")]])
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    load("../../Database/mouse_singleR.rda")
    AnnoResult <- SingleR(test = LayerData(scObject, assay = "RNA", layer = "counts"),
                          ref = mouse_singleR, labels =mouse_singleR$label.main,
                          clusters = scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")]])
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }

  ## 
  NewCluster <- AnnoResult@listData$labels
  names(NewCluster) <- levels(scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")]])
  
  Idents(scObject) <- paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
  scObject <- RenameIdents(scObject, NewCluster)
  scObject$singleR <- scObject@active.ident
  
  ## Dimplot
  {
    ## UMAP
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = UMAP_reduction, group.by = "singleR",
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = TRUE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "singleR_umap.png"),
             bg = "white", width = 9, height = 7)
      
      if(OUTPUTPDF == TRUE){
        ggsave(p1, file = paste0(ResultPath, "singleR_umap.pdf"),
               bg = "white", width = 9, height = 7)
      }
    }
    
    ## TSNE
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = TSNE_reduction, group.by = "singleR",
                  label = TRUE, shuffle = TRUE, raster = FALSE, repel = TRUE,
                  pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "singleR_tsne.png"),
           bg = "white", width = 9, height = 7)
      
      if(OUTPUTPDF == TRUE){
        ggsave(p1, file = paste0(ResultPath, "singleR_tsne.pdf"),
            bg = "white", width = 9, height = 7)
      }
    }
    
  } 
  
  ## SingleRStat
  {
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = "singleR", ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "singleR_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(object = scObject, PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                    RowType = "singleR", ColType = "Sample", PlotPrefixName = "singleR",
                    RowTypeAngle = TRUE, ColTypeAngle = TRUE,
                    ResultPath = ResultPath)
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = "singleR", ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "singleR_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }  
  
}

## Save SingleR results 
if(RunSingleR == TRUE){
  saveRDS(scObject, file = paste0("./Rds/SingleR/5.SingleR_nFeature",nFeature,".Rds"))
  
}