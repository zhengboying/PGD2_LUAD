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

############################################SOP_mRNA_Multiple.R#########################################

## Load subset Seurat object
{
  scObject<-readRDS("./Rds/Results/2.celltype_sub6000.Rds")
}


# 01.Extracted Infor -----
suppressMessages({
  library(stringr)
  library(dplyr)
  library(stringr)
  library(reshape2)

})




# 05.AddSampleInfo -----



scObject<-readRDS("./Rds/Results/5.NK_subtype.Rds")
celltype_NK <- scObject@meta.data %>% select(celltype)

scObject<-readRDS("./Rds/Results/4.T_umap_tsne_nFeature2000.Rds")
celltype_T <- scObject@meta.data %>% select(celltype)

subtype<-rbind(celltype_NK,celltype_T) 
saveRDS(file = "./Rds/Results/6.celltype_NK_T.Rds")




scObject<-readRDS("./Rds/Results/2.celltype_sub6000.Rds")

subtype$ID<-rownames(subtype)

subtype_1 <- subtype[subtype$ID %in% row.names(scObject@meta.data), ]

subtype_1$ID<-NULL
改写colname
colnames(subtype_1)<-c("subcelltype")

scObject <- AddMetaData(
  object = scObject,
  metadata = subtype_1
)


#if meta.data$subcelltype is NA, meta.data$subcelltype = meta.data$celltype
scObject@meta.data$subcelltype[is.na(scObject@meta.data$subcelltype)] <- 
  scObject@meta.data$celltype[is.na(scObject@meta.data$subcelltype)]    

##scObject@meta.data$subcelltype ==T, 替换为unknow
scObject@meta.data$subcelltype[scObject@meta.data$subcelltype =="Unknown"] <- "UnknownT" 

scObject@meta.data$subcelltype[scObject@meta.data$subcelltype =="PTGDS+/PTGDR+ NK"] <- "PTGDS+PTGDR+ NK" 
scObject@meta.data$subcelltype[scObject@meta.data$subcelltype =="KLRC2+/PTGDR+ NK"] <- "KLRC2+PTGDR+ NK" 
scObject@meta.data$subcelltype[scObject@meta.data$subcelltype =="DNAJB1+/PTGDR+ NK"] <- "DNAJB1+PTGDR+ NK" 



saveRDS(scObject, file = "./Rds/Results/7.celltype_sub6000_TYPE.Rds")



## find markers
{
  source("../../Function/CalculateMarker.R")
  ResultPath <- "./Rds/Results/Marker/"
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  
  TopNumber <- 20
  AnnoList <- read.table("../../Database/Human_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
}


TopGene <- c()
Idents(scObject) <- "subcelltype"
for (i in levels(scObject)) {
  CalculateResult <- CalculateMarker(object = scObject,
                                     ClusterName = "subcelltype", ClusterName.1 = i,
                                     min.pct = 0.01, logfc.threshold = 0.01,
                                     max.cells.per.ident = 3000, only.pos = TRUE,
                                     annotation = AnnoList)
  write.csv(CalculateResult, file = paste0(ResultPath, i, "_Marker.csv"),
            quote = TRUE, row.names = FALSE)
  
  TopGene <- c(TopGene, top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName)
  
  CalculateResult <- filter(CalculateResult, p_val < 0.05)
  CalculateResult <- filter(CalculateResult, avg_log2FC >= 0.1)
  write.csv(CalculateResult, file = paste0(ResultPath, i, "_Marker_pvalue_0.05_log2FC_0.1.csv"),
            quote = TRUE, row.names = FALSE)
  
}
  
  
  
