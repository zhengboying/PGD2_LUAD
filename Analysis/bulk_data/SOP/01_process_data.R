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


# 01.Setup -----
  
project_path =  "/Users/zhengboying/Documents/project_demo/Analysis/bulk_data/"
setwd(project_path)


suppressMessages({
  library(data.table)
  library(reshape2)
  library(dplyr)
  library(openxlsx)
})

{
  source(paste0(scripts,"diff_test.R"))
  
  indir = "/Users/zhengboying/Documents/project_demo/Inputdata/bulk_data"
  ResultPath <- c("./01_plot/")
  
}

# 02.Input Load -----


##input data
file.symlink(from = "/Users/zhengboying/Documents/data/LC_bulk_data/TCGA" , to = indir)
file.symlink(from = "/Users/zhengboying/Documents/data/LC_bulk_data/GSE37745" , to = indir)
file.symlink(from = "/Users/zhengboying/Documents/data/LC_bulk_data/GSE31210" , to = indir)
file.symlink(from = "/Users/zhengboying/Documents/data/LC_bulk_data/GSE30219" , to = indir)



# 03. data processing -----

dat<-c()
datL <- c("TCGA","GSE30219","GSE31210","GSE37745")

## expr_TPM
dat$TCGA <- fread(paste0(indir,"/","TCGA","/tpm_final.txt"))
dat$GSE30219 <- fread(paste0(indir,"/","GSE30219","/tpm_final.txt"))
dat$GSE31210 <- fread(paste0(indir,"/","GSE31210","/tpm_final.txt"))
dat$GSE37745 <- fread(paste0(indir,"/","GSE37745","/tpm_final.txt"))

## Gene列作为行名
for (i in datL){
  dat[[i]] <- as.data.frame(dat[[i]])
  rownames(dat[[i]]) <- dat[[i]]$Gene
  dat[[i]] <- subset(dat[[i]], select = -c(Gene,Type))
}
##提取一组TCGA中只含肿瘤的样本
dat$TCGA_T <- dat$TCGA[,grepl("_T_",names(dat$TCGA))]

saveRDS(dat, file = "./Rds/dat_4cohort.rds")

## zscore
dat_zscore <- c()
dat_zscore$TCGA_T <- as.data.frame(sapply(dat$TCGA_T , function(x) (x - mean(x)) / sd(x)))
dat_zscore$TCGA <- as.data.frame(sapply(dat$TCGA , function(x) (x - mean(x)) / sd(x)))
dat_zscore$GSE30219 <- as.data.frame(sapply(dat$GSE30219, function(x) (x - mean(x)) / sd(x)))
dat_zscore$GSE31210 <- as.data.frame(sapply(dat$GSE31210, function(x) (x - mean(x)) / sd(x)))
dat_zscore$GSE37745 <- as.data.frame(sapply(dat$GSE37745, function(x) (x - mean(x)) / sd(x)))

saveRDS(dat_zscore, file = "./Rds/dat_zscore_4cohort.rds")


## Clinical data
cln<-c()
cln$TCGA <- openxlsx::read.xlsx("../../Inputdata/bulk_data/TCGA/clin.xlsx")
row.names(cln$TCGA) <- cln$TCGA$ID
cln$GSE30219 <- openxlsx::read.xlsx("../../Inputdata/bulk_data/GSE30219/clin.xlsx")
row.names(cln$GSE30219) <- cln$GSE30219$ID
cln$GSE31210 <- openxlsx::read.xlsx("../../Inputdata/bulk_data/GSE31210/clin.xlsx")
row.names(cln$GSE31210) <- cln$GSE31210$ID
cln$GSE37745 <- openxlsx::read.xlsx("../../Inputdata/bulk_data/GSE37745/clin.xlsx")
row.names(cln$GSE37745) <- cln$GSE37745$ID

saveRDS(cln, file = "./Rds/cln.rds")

