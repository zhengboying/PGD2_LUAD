
ResChoose <- function(object = scObject, KeyWord = NULL, VarNumber = NULL,
                      ReductionNumber = NULL) {
  if (!is(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }
  
  if (is.null(KeyWord)) {
    stop('argument "KeyWord" is missing.')
  }
  
  if (is.null(VarNumber)) {
    stop('argument "VarNumber" is missing.')
  }
  
  if (is.null(ReductionNumber)) {
    stop('argument "ReductionNumber" is missing.')
  }
  
  AllMetaName <- colnames(object@meta.data)
  FilterMetaName <- AllMetaName[stringr::str_detect(AllMetaName, KeyWord)]
  
  if (length(FilterMetaName) == 0) {
    stop("Please use true KeyWord.")
  }
  
  ## Specify the FilterMetaName of the cluster whose level is closest to the specified VarNumber
  Meta <- FetchData(object, vars = FilterMetaName)
  LevelNumber <- apply(Meta, 2, FUN = function(x) {x <- nlevels(as.factor(x))})
  LevelNumber <- abs(LevelNumber - VarNumber)
  ResNubmer <- ReductionNumber[LevelNumber == min(LevelNumber)][1]
  
  return(as.numeric(ResNubmer))
}
