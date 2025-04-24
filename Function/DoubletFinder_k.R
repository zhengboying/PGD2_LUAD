
DoubletFinder_k <- function(object = scObject, Sample = "orig.ident", PCs = 1:20,
                            DoubletRate = NULL, DoubletRatePer = 0.008, num.cores = 1,
                            plot.output = TRUE, rds.output = TRUE,
                            plot.dir = paste0(ResultPath, "DoubletFinder/"),
                            rds.dir = paste0("./Rds/DoubletFinder/")) {
  if (!is(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }
  
  if (!(is.null(DoubletRate) | is.numeric(DoubletRate))) {
    stop("DoubletRate must be NULL or numeric.")
  }
  
  # split object by Sample
  object_tmp <- SplitObject(object, split.by = Sample)
  
  ## run Seurat
  for(i in names(object_tmp)) {
    DF_tmp <- object_tmp[[i]]
    
    DF_tmp <- NormalizeData(DF_tmp, verbose = FALSE)
    DF_tmp <- FindVariableFeatures(DF_tmp, selection.method = "vst", verbose = FALSE)
    DF_tmp <- ScaleData(DF_tmp, verbose = FALSE)
    DF_tmp <- RunPCA(DF_tmp, seed.use = 6, verbose = FALSE)
    DF_tmp <- RunUMAP(DF_tmp, reduction = "pca", dims = PCs,
                      reduction.name = "umap", reduction.key = "umap_",
                      seed.use = 6, verbose = FALSE)
    
    ## pK Identification(no ground-truth)
    ##(11/21/2023) Made compatible with Seurat v5 and removed '_v3' flag from relevant function names
    ## pN ~ This defines the number of generated artificial doublets, expressed as a proportion 
    ## of the merged real-artificial data. Default is set to 25%, based on observation that 
    ## DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).
    ## pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion 
    ## of the merged real-artificial data. No default is set, as pK should be adjusted for each 
    ## scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.
    ## nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. 
    ## This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, 
    ## and adjusted according to the estimated proportion of homotypic doublets.
    
    sweep.list <- paramSweep(DF_tmp, PCs = PCs, sct = FALSE, num.cores = num.cores)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    ## 10X support: 0.8% per 1000 cells.
    if (is.null(DoubletRate)) {
      nExp_poi <- ceiling((DoubletRatePer * (ncol(DF_tmp) / 1000)) * ncol(DF_tmp))
    } else {
      nExp_poi <- ceiling(DoubletRate * ncol(DF_tmp))
    }
    
    ## Run DoubletFinder
    DF_tmp <- doubletFinder(DF_tmp, PCs = PCs, pN = 0.25, pK = mpK, nExp = nExp_poi[1],
                               reuse.pANN = FALSE, sct = FALSE)
    
    DF_tmp$Doublet <- DF_tmp@meta.data[paste0("DF.classifications_0.25_", mpK, "_", nExp_poi[1])]
    DF_tmp$Doublet <- factor(DF_tmp$Doublet, levels = c("Doublet", "Singlet"))
    DF_tmp$DoubletPoint <- DF_tmp@meta.data[paste0("pANN_0.25_", mpK, "_", nExp_poi[1])]
    
    if (isTRUE(plot.output)) {
      p1 <- DimPlot(DF_tmp, reduction = "umap", group.by = "Doublet",
                    pt.size = 0.75, raster = FALSE, cols = c("red", "black")) &
        labs(title = NULL) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0("./", plot.dir, "/", i, "_DoubletFinder.png"),
             bg = "white", width = 7.5, height = 6)
      
      p1 <- FeaturePlot(DF_tmp, reduction = "umap", features = "DoubletPoint",
                        pt.size = 0.75, raster = FALSE) &
        labs(title = NULL) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0("./", plot.dir, "/", i, "_DoubletFinderScore.png"),
             bg = "white", width = 7, height = 6)
    }
    
    if (isTRUE(rds.output)) {
      saveRDS(DF_tmp, file = paste0("./", rds.dir, "/", i, ".Rds"))
    }
    
    if(i == names(object_tmp)[1]) {
      scObject_DF <- FetchData(DF_tmp, vars = c("Doublet", "DoubletPoint"))
    } else {
      scObject_DF <- rbind(scObject_DF, FetchData(DF_tmp, vars = c("Doublet", "DoubletPoint")))
    }
  }
  
  colnames(scObject_DF) <- c("DoubletFinder", "DoubletFinderScore")
  return(scObject_DF)
}