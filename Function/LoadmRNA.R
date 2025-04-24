# Define a function to load 10X scRNA-seq data
LoadmRNA <- function(SampleName = SampleName, DateType = NULL,
                     CellRangerPath = paste0(project_path, "/Inputdata/"),
                     Project = "demo",
                     MEXPath = "/",
                     h5Path = "/outs/filtered_feature_bc_matrix.h5",
                     min.cells = 1, min.features = 1,
                     names.field = 1, names.delim = "-") {
  # Check if DateType input is valid and not null
  if (is.null(DateType) || !(DateType %in% c("MEX", "h5", "MGI"))) {
    stop("DateType must be one of the following: MEX, h5, MGI")
  }
  
  # Load data for each sample and create a Seurat object
  message(paste0(Sys.time(), " Load start."))
  object <- lapply(SampleName, function(i) {
    # Load data depending on the DateType input
    if (DateType == "MEX") {
      object <- Read10X(paste0(CellRangerPath, "/", i, MEXPath),gene.column = 2,strip.suffix = TRUE)#Remove trailing "-1" if present in all cell barcodes.
    }
    
    if (DateType == "h5") {
      object <- Read10X_h5(paste0(CellRangerPath, "/", i, h5Path))
    }
    
    if (DateType == "MGI") {
      object <- Read10X(paste0(CellRangerPath, "/", i, MEXPath), gene.column = 1)
    }
    
    # Create a Seurat object using the loaded data
    object <- CreateSeuratObject(counts = object, assay = "RNA",
                                 min.cells = min.cells, min.features = min.features,
                                 names.field = names.field, names.delim = names.delim)
    
    # Modify the cell names in the Seurat object
    #object <- RenameCells(object, new.names = gsub("\\-1", "", Cells(object)))
    object <- RenameCells(object, add.cell.id = i)
    object@meta.data$Sample<-i
    # Print message indicating completion of loading for each SampleName
    message(paste0(Sys.time(), " Load ", i, " finish."))
    
    return(object)
  })
  
  # If there are more than one SampleName, merge the Seurat objects for each SampleName
  if (length(SampleName) > 1) {
    message(paste0(Sys.time(), " Start merge all samples."))
    #object <- merge(object[[1]], y = object[-1])
    object<- merge(object[[1]], y = object[2:length(object)], add.cell.ids = NULL, project = Project)
    object <- JoinLayers(object)
    Idents(object)<-object@meta.data$Sample
    message(paste0(Sys.time(), " Samples merge completed."))
  }
  
  # If there is only one SampleName, assign the Seurat object to the object variable
  if (length(SampleName) == 1) {
    object <- object[[1]]
  }
  print(object)
  return(object)
}
