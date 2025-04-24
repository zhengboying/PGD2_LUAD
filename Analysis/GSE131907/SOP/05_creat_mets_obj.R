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
  # param3：SampleSpecies
  SampleSpecies <- "Human"
}  
source("../../Function/PlotTheme.R")

# 03.Load Data -----
if (TRUE){
  mets_mat<-readRDS("./Rds/Mets/met_mat.rds")
  meta_data<-read.csv("./03_mebocost/cell_ann.csv",row.names = 1)
}

##obj
{
  scObject <- CreateSeuratObject(counts = mets_mat,
                               min.cells = 3, min.features = 1,
                               meta.data = meta_data)

Idents(scObject) <- "celltype"
scObject$celltype <- scObject@active.ident
}


# 06.Basic analysis -----

suppressMessages({
  library(stringr)
})

## Setup
{
  ## FindVariableFeatures Parameters
  FindVariableFeatures_nfeatures <- c(200)
  
  FindVariableFeatures_selection.method = "vst"
  ## PCA Parameters
  PCA_npcs = 30
}



## 06.2 Basic Seurat -----
{ 
  ## 
  scObject@assays$RNA$data <-  scObject@assays$RNA$counts
  scObject <- FindVariableFeatures(scObject, assay = "RNA", selection.method = FindVariableFeatures_selection.method,
                                 nfeatures = FindVariableFeatures_nfeatures, verbose = FALSE)

  ## ScaleData,默认做过
  
  scObject <- ScaleData(scObject, assay = "RNA", features = VariableFeatures(scObject),
                      verbose = FALSE)
  
  ## PCA
  scObject <- RunPCA(scObject, assay = "RNA", npcs = PCA_npcs,
                     seed.use = 16, verbose = FALSE)
}


# 07. Cluster -----

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
  ResultPath <- c("./03_mebocost/")
  
  ## NoIntegration and HarmonyIntegration parameters
  dims = 1:30
  RNA_snn_res = c( 0.1)
  TSNE = FALSE
  UMAP = TRUE
  min.dist = 0.5
  

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
}

  
  ##
  for(i in RNA_snn_res) {
    scObject <- FindClusters(scObject, resolution = i, random.seed = 16, verbose = FALSE)
    
    scObject@meta.data[[paste0("RNA_snn_res.", i)]] <- NULL
    scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW")]] <- scObject$seurat_clusters
    scObject$seurat_clusters <- NULL
    
    ## Plot
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "umap", group.by = "celltype",
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/","res.", sprintf("%0.1f", i), "_umap.png"),
             bg = "white", width = 8, height = 7)
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/","res.", sprintf("%0.1f", i), "_umap.pdf"),
             bg = "white", width = 8, height = 7)
      
    }
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "tSNE", group.by = "celltype",
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath,"NoIntegration/", "res.", sprintf("%0.1f", i), "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }
    
  }}

## 07.2 HarmonyIntegration# 效果不好
}

## Save step 07 result
{
  Idents(scObject)<-"celltype"
  saveRDS(scObject, file = paste0("./Rds/Mets/1.metobj.Rds"))
}






##HMDB0001403 = Prostaglandin D2
## test
FeaturePlot(scObject, features = c("HMDB0001403"), reduction = "umap",split.by = "location",
            raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),ncol =2,slot = "counts",
            pt.size = 0.5, max.cutoff = "q99") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())

VlnPlot(subset(scObject,idents = c("pDC","Mast","Fibro","Endo","Oligo")), features = c("HMDB0001403"), pt.size = 0,split.by = "location",slot = "counts") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


VlnPlot(subset(scObject,idents = c("pDC","Fibro")), features = c("HMDB0001403"), pt.size = 0.1,split.by = "location",slot = "counts") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())



test<-subset(scObject,idents = c("pDC","Fibro"),subset = location %in% c("Normal Lung","Tumor Lung") )
             

VlnPlot(test, features = c("HMDB0001403"), pt.size = 0,split.by = "location",slot = "counts") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())


VlnPlot(subset(test,subset = location == "Tumor Lung"), features = c("HMDB0001403"), pt.size = 0,split.by = "stage.1",slot = "counts") &
  PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
