#
{
library(dplyr)
library(data.table)
library(tools)
library(Seurat)
}



args <- commandArgs(trailingOnly = TRUE)
workdir = args[1]
setwd(paste0(workdir,"03_mebocost"))



obj<-readRDS("../Rds/Results/1.celltype.Rds")
table(obj@meta.data$celltype)
#提取mebocost 需要的inputfile，超过5000的cell取5000
Idents(obj) <- "celltype"
tmp1<-subset(x = obj,idents = c("Endo","Fibro","Mast","Oligo","pDC"))
tmp2<-subset(x = obj,idents = c("B","Epi","Myeloid","NK","T"), downsample = 6000)

sub_obj <- merge(x = tmp1, y = tmp2)
sub_obj <- JoinLayers(sub_obj)

#cells <-row.names(obj@meta.data[obj@meta.data$orig.ident == "GSM7157674",])
#subset_obj <- subset(obj, cells = cells)


table(sub_obj@meta.data$celltype)

saveRDS(sub_obj, file = "../Rds/Results/2.celltype_sub6000.Rds")

write.csv(sub_obj@meta.data, file = "cell_ann.csv",quote = F)
exp_mat<-as.matrix(sub_obj@assays$RNA$data)#取scale之后的数据
exp_mat[1:4,1:4]
fwrite(exp_mat, "exp_mat.csv.gz",row.names = TRUE)

print("finished extract data with downsample 6000！！！")





