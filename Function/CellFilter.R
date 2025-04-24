
CellFilter <- function(object = scObject, Type = "Sample",
                       MinCount = 500, MinFeature = 200,MinGenesPerUMI = 0.8,
                       Sigma = 3, SigmaTop = TRUE, SigmaBottom = TRUE,
                       MinPercent.mt = NULL, MaxPercent.mt = NULL,
                       MaxPercent.hgb = NULL, MaxPercent.rib = NULL) {
  if (!is(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }
  
  #Idents(object) <- Type
  
  ## 
  ## percent.mt
  if(is.numeric(MinPercent.mt)) {
    object <- subset(object, subset = percent.mt >= MinPercent.mt)
  }
  
  if(is.numeric(MaxPercent.mt)) {
    object <- subset(object, subset = percent.mt < MaxPercent.mt)
  }
  
  ## percent.hgb
  if(is.numeric(MaxPercent.hgb)) {
    object <- subset(object, subset = percent.hgb < MaxPercent.hgb)
  }
  
  ## percent.rib
  if(is.numeric(MaxPercent.rib)) {
    object <- subset(object, subset = percent.rib < MaxPercent.rib)
  }
  
  ## filter nCount_RNA nFeature_RNA with sd
  NeedCell <- lapply(levels(object), function(i) {
    ## subset one sample
    object_tmp <- subset(object, idents = i)
    
    ## nCount
    nCountTop <- 10^(mean(log10(object_tmp$nCount_RNA)) + Sigma * sd(log10(object_tmp$nCount_RNA)))
    nCountBottom <- 10^(mean(log10(object_tmp$nCount_RNA)) - Sigma * sd(log10(object_tmp$nCount_RNA)))
    
    ## nFeature
    nFeatureTop <- 10^(mean(log10(object_tmp$nFeature_RNA)) + Sigma * sd(log10(object_tmp$nFeature_RNA)))
    nFeatureBottom <- 10^(mean(log10(object_tmp$nFeature_RNA)) - Sigma * sd(log10(object_tmp$nFeature_RNA)))
    
    ## subset cell not within the confidence interval
    if (SigmaTop) {
      object_tmp <- subset(object_tmp, subset = (nCount_RNA < nCountTop &
                                                   nFeature_RNA < nFeatureTop))
    }
    
    if (SigmaBottom) {
      object_tmp <- subset(object_tmp, subset = (nCount_RNA > nCountBottom &
                                                   nFeature_RNA > nFeatureBottom))
    }
    
    return(Cells(object_tmp))
  })
  
  object <- subset(object, cells = unlist(NeedCell))
  
  ## filter low nCount_RNA
  if(is.numeric(MinCount)) {
    object <- subset(object, subset = nCount_RNA > MinCount)
  }
  
  ## filter low GenesPerUMI
  if(is.numeric(MinGenesPerUMI)) {
    object <- subset(object, subset = GenesPerUMI > MinGenesPerUMI)
  }
  ## filter low nFeature_RNA
  if(is.numeric(MinFeature)) {
    object <- subset(object, subset = nFeature_RNA > MinFeature)
  }
  
  return(object)
}
