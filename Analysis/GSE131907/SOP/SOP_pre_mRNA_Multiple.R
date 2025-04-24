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
  pre_process = TRUE
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
