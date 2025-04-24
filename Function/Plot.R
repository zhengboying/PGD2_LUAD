anno_hematmap<-function(obj = obj2,
                        markers = mye_markers,
                        col_plates = col_plates,
                        group = "sub_celltype",#any column in meta.data
                        cell_cluster = 3,#any
                        row_names_fontsize = 5,
                        seed = 12){
  ##process markers
  print(paste0("Markers number: ",length(markers$Marker)))
  print(paste0("Including cell type: ",unique(markers$CellType)))
  #提取markers表达情况
  Idents(obj)<-group
  p_dotplot <- DotPlot(object = obj, features = unique(markers$Marker,idents = group))
  df<-p_dotplot$data
  ##提取表达scale均值
  exp_mat <- df %>%  
    dplyr::select(features.plot, id, avg.exp.scaled) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%  
    as.data.frame()                      
  row.names(exp_mat) <- exp_mat$features.plot
  exp_mat["features.plot"] <-NULL
  head(exp_mat)
  #删除整行NaN的行
  exp_mat <- exp_mat[complete.cases(exp_mat), ]
  dim(exp_mat)
  #提取细胞比例
  percent_mat <- df %>%  
    dplyr::select(features.plot, id, pct.exp) %>% 
    pivot_wider(names_from = id, values_from = pct.exp) %>%  
    as.data.frame()           
  row.names(percent_mat) <- percent_mat$features.plot
  percent_mat["features.plot"] <-NULL
  #回填所有NA为0
  percent_mat[is.na(percent_mat)] <- 0
  head(percent_mat)
  dim(percent_mat)
  percent_mat <- percent_mat[!apply(percent_mat, 1, function(x) all(x == 0)), ]
  dim(percent_mat)
  #plot
  #col_fun = circlize::colorRamp2(c(-0.6, 0, 2), brewer.pal(11, "PiYG")[c(11,8,1)])
  col_fun = circlize::colorRamp2(c(min(exp_mat), 0, max(exp_mat)), c("#424da7","#ffffff","#dd2b19"))
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h,gp = gpar(col = NA, fill = NA))      
    grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit(2, "mm")),              
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
  annotation_df <- markers[match(rownames(exp_mat), markers$Marker),]
  row_ha <- rowAnnotation(
    celltype = anno_simple(annotation_df$CellType, 
                           col = col_plates))# 不需要 show_legend 参数
  lgd_list = list(
    Legend(labels = c(0,0.25,0.5,0.75,1), title = "Percentage of Cell Type",          
           graphics = list(function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(2, "mm"),
                                                            gp = gpar(fill = "black")),            
                           function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"), 
                                                            gp = gpar(fill = "black")),            
                           function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"), 
                                                            gp = gpar(fill = "black")),            
                           function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"), 
                                                            gp = gpar(fill = "black")),            
                           function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"), 
                                                            gp = gpar(fill = "black")))))
  set.seed(seed)  #
  ha<-Heatmap(exp_mat,        
              heatmap_legend_param = list(title = "Scaled Expression level"),
              col = col_fun, rect_gp = gpar(type = "none"),cell_fun = cell_fun,       
              row_names_gp = gpar(fontsize = row_names_fontsize),
              #row_km = cell_cluster,
              row_split = annotation_df$CellType,
              border = "black",
              left_annotation = row_ha,
              cluster_rows = FALSE, 
              cluster_columns = TRUE,
              show_parent_dend_line = F)
  ht<-draw(ha, annotation_legend_list = lgd_list,heatmap_legend_side = "right")        
  return(ht)          
}





#' @title Umap_plot_with_anno
#' 
#' @description This function plotting the cell by umaps and annotate the cell type.
#' 
#' @param obj obj after cellAnnotPipe, or readRDS(paste0(workdir,"/01_singlecell_seurat/final_obj_cellAnno.rds"))
#' @return p
#' @examples
# p<-Umap_plot_with_anno(obj)
#

Umap_plot_with_anno<-function(obj){
                    umaps <- as.data.frame(obj@reductions$umap@cell.embeddings)
                    celltype <- Idents(obj)
                    umaps <- cbind.data.frame(umaps,celltype)
                    labs <- umaps %>% group_by(celltype) %>% summarise(UMAP_1 = median(umap_1),UMAP_2 = median(umap_2))
                    labs <- as.data.frame(labs)
                    pal <- palettes(category = c("random","default")[1],
                                    methods = c("continuous","discrete")[2],
                                    X = NULL,
                                    seeds = 0,
                                    n = length(unique(umaps$celltype)),
                                    Alpha = 1,
                                    show_col = F)

                    p <- ggplot(data = umaps,aes(x = umap_1,y = umap_2,color = celltype)) + 
                                 geom_point(size = 0.01,alpha = 1) +
                                 #画圈的，一般不用
                                 # stat_ellipse(aes(color = celltype,fill = celltype),level = .95,geom = "polygon",linetype = 1,show.legend = F,alpha = .4) +
                                 theme_void() +
                      theme_dr(xlength = 0.2,ylength = 0.2,arrow = grid::arrow(length = unit(0.1,"inches"),ends = "last",type="closed")) +
                      geom_text(data = labs,aes(x = UMAP_1,y = UMAP_2,label = celltype),color = "black",size = 2) +
                      scale_color_manual(values = pal) +
                      scale_fill_manual(values = pal) +
                      theme(
                           panel.grid = element_blank(),
                           plot.margin = margin(5.5,15,5.5,5.5),
                           legend.key.size = unit(0.2,"inches")
                      ) + labs(color = "Cluster") + 
                      guides(shape = guide_legend(override.aes = list(size = 8)))
                    return(p)
}




#' @title Featureplot_pipe
#' 
#' @description This function plotting the cell markers in every cell type.
#' 
#' @param cell_type seurat object from Obj_create.
#' @param obj obj after cellAnnotPipe, or readRDS(paste0(workdir,"/01_singlecell_seurat/final_obj_cellAnno.rds"))
#' @param method 归一化方法，包含seurat开发的ScaleData和SCTransform。
#' @param FindVariableFeatures.data 默认是"ALL"，所有样本一起找高变基因，这时候的高变基因包含样本间和样本内部细胞间，当样本间处理差异较大时，或者分群结果校正后仍有较为严重的批次，考虑使用代表性的sample,专注于细胞类型间高变基因
#' @return A seurat object containing the loaded data.
#' @examples
# Featureplot_pipe(obj)
#
Featureplot_pipe<-function(obj,obj.anno){
  dir.create("./01_singlecell_seurat/markers_plot")
  cell_type<-unique(obj@meta.data$cell_type)
  if("Unknown" %in% cell_type){
    cell_type<-cell_type[-which(cell_type=="Unknown")]
    print("Unknown cell type removed!")
  } else{
    print("No Unknown cell type!")
  }
  ##cell_type 元素中的空格以_代替，+ () 删除生成names
  names(cell_type)<-gsub(" ","_",cell_type)
  names(cell_type)<-gsub("\\(","",names(cell_type))
  names(cell_type)<-gsub("\\)","",names(cell_type))
  names(cell_type)<-gsub("\\+","",names(cell_type))
  print(cell_type)
  for (i in 1:length(cell_type)){
    print(paste0("Ploting ", cell_type[i]))
    markers<- na.omit(obj.anno$gene.set$Symbol[obj.anno$gene.set$cell_name == cell_type[i]])
    print(paste0(cell_type[i]," markers number:", length(markers)))
    print(markers)
    p <- FeaturePlot(obj, features = markers, ncol = 4,label = F,label.size = 1)
    nrows<-ceiling(length(markers)/4)
    ggsave(paste0(workdir,"/01_singlecell_seurat/markers_plot/",names(cell_type[i]),".jpg"),plot = p,dpi = 600,units = "in",width = 14,height = 1.5*nrows)
    #ggsave(paste0(workdir,"/01_singlecell_seurat/markers_plot/",names(cell_type[i]),".pdf"),plot = p,dpi = 600,units = "in",width = 14,height = 12)
  }
  ##plot Tumor cell markers
  tumor_markers <- c('MUC16','EPCAM','MET','EGFR','SMARCA4','ROS1','ALK','RET','NAPSA','NKX2-1','CEACAM5','MUC1','ALDH1A1')
  p <- FeaturePlot(obj, features = tumor_markers, ncol = 4, label = F,label.size = 1)
  nrows<-ceiling(length(tumor_markers)/4)
  ggsave(paste0(workdir,"/01_singlecell_seurat/markers_plot/cancerMarker.jpg"),plot = p,dpi = 600,units = "in",width = 14,height = 1.5*nrows)
  #ggsave(paste0(workdir,"/01_singlecell_seurat/markers_plot/cancerMarker.pdf"),plot = p,dpi = 600,units = "in",width = 14,height = 12)
  print("all the results saving in ./01_singlecell_seurat/markers_plot!!")
  
}




#' @title totalCNV_plot
#' 
#' @description This function plotting total cnv for eachcell, group by group_column.
#' 
#' @param data input data, a data.frame with columns "cell_type" ,"orig.ident","seurat_clusters" and "cnv_score"
#' @param group_column violin plot group by group_column
#' @examples
# totalCNV_plot(totalCNV_table, "cell_type")
#
totalCNV_plot<- function(totalCNV_table, 
					group_column = "cell_type",
					outputdir = "./02_infercnv"##without "/"
					) {
  # 检查 totalCNV_table 是否为数据框
  if (!is.data.frame(totalCNV_table)) {
    stop("totalCNV_table must be a data frame")
  }
  
  # 检查 group_column 是否为字符串
  if (!is.character(group_column)) {
    stop("The 'group_column' parameter must be a string (e.g., 'cell_type').")
  }
  
  # 检查 group_column 列是否存在于 totalCNV_table 中
  if (!(group_column %in% colnames(totalCNV_table))) {
    stop(paste("Column", group_column, "does not exist in totalCNV_table"))
  }
  n_color <- n_distinct(totalCNV_table[[group_column]])
  color_v <- palettes(category = c("random", "default")[1],
                      methods = c("continuous", "discrete")[2],
                      X = NULL,
                      seeds = 0,
                      n = n_color,
                      Alpha = 1,
                      show_col = F)
  totalCNV_table[[group_column]] <- factor(totalCNV_table[[group_column]])
  p <- ggplot(totalCNV_table, aes(x = !!sym(group_column), y = cnv_score, fill = !!sym(group_column))) +
   geom_violin() +
   scale_fill_manual(values = color_v) +
   theme_bw() +
   theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16, color = "black")
  )
  
  ggsave(paste0(outputdir, "/CNV_by_",group_column,".jpg"), plot = p, width = 10, height = 6, units = "in")
}




#' @title DoHeatmapPlot
DoHeatmapPlot <- function(object, groupBy, features) {
  require(ComplexHeatmap)
  # 获取绘图数据
  plot_data <- SeuratObject::FetchData(object = object,
                                       vars = c(features, groupBy), 
                                       slot = 'counts') %>% 
    dplyr::mutate(across(.cols = where(is.numeric), .fns = ~ log2(.x + 1))) %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    dplyr::arrange(group) %T>% 
    assign(x = 'clusterInfo', value = .$group, envir = .GlobalEnv) %>% 
    dplyr::select(-group) %>% 
    t()
}
theme_dr <- function(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")) {
  ggplot2::theme(
    # 这里可以设置你的自定义主题参数
    # panel.grid.major = element_line(size = xlength),
    
    # 使用grid::arrow来自定义箭头样式
    # legend.key = element_rect(fill = "blue"),
    # 你需要根据实际需要来设置这些参数
  )
}
