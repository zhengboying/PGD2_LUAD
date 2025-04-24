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
  #library(METAFlux)
  library(pROC)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(ggprism)
  library(ggsignif)
  library(GSVA)
  library(BiocParallel)
  library(corrplot)
  library(plotROC)
  library(survAUC)
  library(openxlsx)
  library(survival)
  library(survminer)
})

{
  source("../../Function/Prognosis_function.R")
  
  InputPath <- c("./Rds/")
  
}

# Input Load -----
dat<-readRDS(paste0(InputPath,"dat_4cohort.rds"))
combined_data<-readRDS(paste0(InputPath ,"combined_data_4cohort.rds"))
cln<-readRDS(paste0(InputPath,"cln.rds"))
dat_zscore<-readRDS(paste0(InputPath,"dat_zscore_4cohort.rds"))
model_score<-readRDS(paste0(InputPath,"model_score.rds"))




##读取../../Input/



## Step1：oncoPredict::calcPhenotype  -----
{
ResultPath <- c("./03_calcPhenotype/")
dir.create(ResultPath, showWarnings = FALSE)
file.symlink(from = "/Users/zhengboying/Documents/data/DATABASE/oncoPredict" , to = "../../Database")
source("../../Function/oncoPredict.R")
#TrainingData Information
#GDSC1Data: 367 drugs, 958 training samples
#GDSC2Data:198 drugs, 805 training samples
#CTRP2Data: 545 drugs, 829 training samples (提供TPM,RPKM 默认用TPM)

#GDSC1Data
{
  df<-as.matrix(dat$TCGA_T)
  
  calcPhenotype_pipeline(TrainingData = "GDSC1Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/TCGA_GDSC1Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE30219)
  
  calcPhenotype_pipeline(TrainingData = "GDSC1Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE30219_GDSC1Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE31210)
  
  calcPhenotype_pipeline(TrainingData = "GDSC1Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE31210_GDSC1Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE37745)
  
  calcPhenotype_pipeline(TrainingData = "GDSC1Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE37745_GDSC1Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  
}

#GDSC2Data （new）
{
  df<-as.matrix(dat$TCGA_T)
  
  calcPhenotype_pipeline(TrainingData = "GDSC2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/TGCA_GDSC2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE30219)
  
  calcPhenotype_pipeline(TrainingData = "GDSC2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE30219_GDSC2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE31210)
  
  calcPhenotype_pipeline(TrainingData = "GDSC2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE31210_GDSC2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE37745)
  
  calcPhenotype_pipeline(TrainingData = "GDSC2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE37745_GDSC2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  
  
  
  
}


#CTRP2Data(too big no run)
{
  df<-as.matrix(dat$TCGA_T)
  
  calcPhenotype_pipeline(TrainingData = "CTRP2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/TGCA_CTRP2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE30219)
  
  calcPhenotype_pipeline(TrainingData = "CTRP2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE30219_CTRP2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE31210)
  
  calcPhenotype_pipeline(TrainingData = "CTRP2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE31210_CTRP2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  df<-as.matrix(dat$GSE37745)
  
  calcPhenotype_pipeline(TrainingData = "CTRP2Data",##GDSC1Data, GDSC2Data, CTRP2Data
                         Outputpath = paste0(project_path,"/03_calcPhenotype/GSE37745_CTRP2Data_calcPhenotype_Output/"),
                         Inputdata = df,
                         TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                         batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                         
                         cc=FALSE,##correlation coefficients for biomarker discovery
                         rsq=FALSE,##output the R^2 values for the data you train on from true and predicted values
                         
                         #pcr must be FALSE if cc is TRUE
                         pcr = FALSE,
                         report_pc = FALSE,
                         percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
  )
  
  
}





## Step2：sur -----
### single gene survival analysis
model_score<-c()
##TCGA
{
  
  
  
  ## make data
  {
    expr_table <- combined_data[combined_data$Dataset == "TCGA",]
    #expr_table <- expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1")]
    #expr_table <- dat_zscore$TCGA_T[rownames(dat_zscore$TCGA_T) %in% c("PTGDS", "PTGDR", "HPGDS", "SLCO2A1"), ]
    #expr_table<-t(expr_table)
    ##因为density发现数据分布差异是明显的批次效应
    #expr_table <- as.data.frame(scale(expr_table))
   
    #REMOVER  string _T_LUAD in rowname
    rownames(expr_table) <- gsub("_T_LUAD","",rownames(expr_table))
    #expr_table = log2(expr_table + 0.001)
    
    clin <- cln$TCGA
    ##colname time  to OS
    colnames(clin)[colnames(clin) == "time"] <- "OS"
    ##  Status == Alive ,0,  Status == Dead,1
    clin$Status <- ifelse(clin$Status == "Alive",0,1)
    head(clin)
    rownames(clin)<-clin$ID
    
    data <- cbind(clin[,c('OS', 'Status')],expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1"),drop = FALSE])
  }
  
  ## cox
  
  ##单基因
  {
    FUN_surv <- function(i, surv_expr) {  
    Sur <- Surv(surv_expr$OS, surv_expr$Status)  
    cox <- coxph(Sur ~ surv_expr[, i], data = surv_expr)  
    coxSummary <- summary(cox)  
    gene <- names(surv_expr)[i]  
    HR <- coxSummary$coefficients[, "exp(coef)"]  
    lower <- coxSummary$conf.int[, 3]  
    upper <- coxSummary$conf.int[, 4]  
    PValue <- round(coxSummary$coefficients[, 5], 6)  
    res <- data.frame(gene = gene, HR = HR, lower.95 = lower, upper.95 = upper, pvalue = PValue)  
    return(res)
  }
    l <- pbapply::pblapply(3:ncol(data), FUN = FUN_surv, surv_expr = data)
    Univar <- do.call(rbind, l)
    
    write.csv(Univar, file = "./04_sur/Univar_result_OS.csv", row.names = F)
  }
  #plot
  {library(forestplot)
    data <- read.csv("./04_sur/Univar_result_OS.csv")
    head(data)
    #2,3,4,5列取4位有效数字
    data$HR <- round(data$HR, 4)
    data$lower.95 <- round(data$lower.95, 4)
    data$upper.95 <- round(data$upper.95, 4)
    data$pvalue <- round(data$pvalue, 4)
    
    data1 <- data %>% dplyr::filter(HR > 1)
    data2 <- data %>% dplyr::filter(HR < 1)
    hrtable <- rbind(  c("Increase in Hazard", NA, NA, NA, NA, NA),
                       data1,  c("Reduce in Hazard", NA, NA, NA, NA, NA), 
                       data2)
    # 将要在图中展示的文本
    tabletext <- cbind(  
      c("Gene", "Increase in Hazard", data1$gene, "Reduce in Hazard", data2$gene),  
      c("P-value", NA, data1$pvalue, NA, data2$pvalue), 
      c("Hazard Ratio", NA, data1$`HR`, NA, data2$`HR`))
    which(tabletext[, 1] == "Reduction in Hazard")
    nrow(tabletext) + 1
    quantile(data$lower)
    quantile(data$upper)
    pdf("./04_sur/Univar_forestplot.pdf", width =6, height = 6) 
    
    
    ##所有tabletext数值取4位有效数字
    
    p<-forestplot(  
      labeltext = tabletext, 
      # 图中的文本  
      mean = c(NA, as.numeric(hrtable$HR)), 
      # HR  
      lower = c(NA, as.numeric(hrtable$lower)), 
      upper = c(NA, as.numeric(hrtable$upper)),
      graph.pos = 4, 
      # 图在表中的列位置  graphwidth = unit(.4, "npc"), 
      # 图在表中的宽度比例  fn.ci_norm = "fpDrawDiamondCI", 
      # box类型选择钻石
      col = fpColors(box = "#CD534CFF",lines="#0073C2FF", zero = "#0073C2FF"),
      #box颜色  
      boxsize = 0.7, 
      # box大小根据样本量设置  
      lwd.ci = 3, 
      ci.vertices.height = 0.125,  
      ci.vertices = T, 
      # 置信区间用线宽、高、型  
      zero = 1, 
      # zero线横坐标  
      lwd.zero = 2,  
      xticks = seq(0.4, 1.2, 0.2), 
      # 横坐标刻度根据需要可随意设置  
      lwd.xaxis = 2, # X轴线宽  
      xlab = "Hazard Ratio",  
      hrzl_lines = list("1" = gpar(lwd = 2, col = "black"), 
                        # 第二行顶部加黑实线    
                        "2" = gpar(lwd = 1, col = "grey50", lty = 2),
                        # 第二行顶部加灰色虚线    
                        "3" = gpar(lwd = 1, col = "grey50", lty = 2), 
                        "4" = gpar(lwd = 1, col = "grey50", lty = 2),
                        # 第二行顶部加灰色虚线    
                        "5" = gpar(lwd = 1, col = "grey50", lty = 2), 
                        # 第九行顶部加灰色虚线    
                        "8" = gpar(lwd = 2, col = "black")  ),
      # 最后一行底部加黑线，""中数字为nrow(tabletext) + 1  
      lineheight = unit(0.8, "cm"),   
      txt_gp = fpTxtGp(label = gpar(cex = 1), 
                       # 各种字体大小设置    
                       ticks = gpar(cex = 1.25),    
                       xlab = gpar(cex = 1.85),    
                       title = gpar(cex = 1.5)),  
      clip = c(0, 2),  
      colgap = unit(0.15, "cm"),   
      mar = unit(rep(1.25, times = 4), "cm"), 
      #  new_page = F # 是否新页
    )
    p
    
    dev.off()
  }
  
 
  {
    ##
    fit_cox <- coxph(Surv(OS, Status) ~ PTGDR+HPGDS+PTGDS, data = data)
    summary(fit_cox)
    #coxph(formula = Surv(OS, Status) ~ PTGDR + HPGDS + PTGDS, data = data)
    
    #n= 504, number of events= 181 
    
    #coef exp(coef) se(coef)      z Pr(>|z|)  
    #PTGDR -0.09648   0.90803  0.14135 -0.683   0.4949  
    #HPGDS -0.03077   0.96970  0.01663 -1.850   0.0643 .
    #PTGDS -0.07928   0.92378  0.04139 -1.916   0.0554 .
    # ---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    
    #exp(coef) exp(-coef) lower .95 upper .95
    #PTGDR    0.9080      1.101    0.6883     1.198
    #HPGDS    0.9697      1.031    0.9386     1.002
    #PTGDS    0.9238      1.083    0.8518     1.002
    
    #Concordance= 0.611  (se = 0.023 )
    #Likelihood ratio test= 13.07  on 3 df,   p=0.004
    #Wald test            = 11.5  on 3 df,   p=0.009
    #Score (logrank) test = 11.47  on 3 df,   p=0.009
    
    
    ##resuts
    ##library(car)
    vif(fit_cox)
    #PTGDR    HPGDS    PTGDS 可接受的共线性
    #1.154256 1.045590 1.200671

    
    
    #
    #fit <- glm(Status ~ PTGDS + PTGDR, family=binomial(link = "logit"), data = data_TCGA)
    data$risk_score <- with(data, coef(fit_cox)[1] * PTGDR + coef(fit_cox)[2] * HPGDS, coef(fit_cox)[3] * PTGDS)
  }
  model_score$TCGA<-data
  
  #多因素cox
  {
    data <- cbind(clin[,c('OS','Status','Stage', 'T','N','M','Gender','Age')],model_score$TCGA[,c("risk_score"),drop = FALSE])
    
    
    # Stage
    data[which(data$Stage %in% c('Stage I', 'Stage IA', 'Stage IB')),'Stage'] <- '1'
    data[which(data$Stage %in% c('Stage II', 'Stage IIA', 'Stage IIB')),'Stage'] <- '2'
    data[which(data$Stage %in% c('Stage IIIA', 'Stage IIIB')),'Stage'] <- '3'
    data[which(data$Stage %in% c('Stage IV')),'Stage'] <- '4'
    data$Stage <- as.numeric(data$Stage)
    # T
    data[which(data$T %in% c('T1', 'T1a', 'T1b')),'T'] <- '1'
    data[which(data$T %in% c('T2', 'T2a', 'T2b')),'T'] <- '2'
    data[which(data$T %in% c('T3')),'T'] <- '3'
    data[which(data$T %in% c('T4')),'T'] <- '4'
    data[which(data$T %in% c('T5')),'T'] <- '5'
    data$T <- as.numeric(data$T)
    
    # N
    data[which(data$N %in% c('N0')),'N'] <- '0'
    data[which(data$N %in% c('N1')),'N'] <- '1'
    data[which(data$N %in% c('N2')),'N'] <- '2'
    data[which(data$N %in% c('N3')),'N'] <- '3'
    data[which(data$N %in% c('NX')),'N'] <- '4'
    data$N <- as.numeric(data$N)
    
    # M
    data[which(data$M %in% c('M0')),'M'] <- '0'
    data[which(data$M %in% c('M1', 'M1a', 'M1b')),'M'] <- '1'
    data[which(data$M %in% c('MX')),'M'] <- '2'
    data$M <- as.numeric(data$M)
    
    # Female=2;Male=1
    data$Gender <- ifelse(data$Gender == 'male',1,2)
    data$Gender <- as.numeric(data$Gender)
    data$Age <- as.numeric(data$Age)
    
    
    #多因素cox
    multiCox=coxph(Surv(OS, Status)~., data = data)
    
    coxSummary <- summary(multiCox)  
    gene <- names(data)[-c(1,2)]
    HR <- coxSummary$coefficients[, "exp(coef)"]  
    lower <- coxSummary$conf.int[, 3]  
    upper <- coxSummary$conf.int[, 4]  
    PValue <- round(coxSummary$coefficients[, 5], 6)  
    multiTab <- data.frame(gene = gene, HR = HR, lower.95 = lower, upper.95 = upper, pvalue = PValue) 
    write.csv(multiTab, file = "./04_sur/multiCox.csv", row.names = F)
    
    {library(forestplot)
      data <- read.csv("./04_sur/multiCox.csv")
      head(data)
      #2,3,4,5列取4位有效数字
      data$HR <- round(data$HR, 4)
      data$lower.95 <- round(data$lower.95, 4)
      data$upper.95 <- round(data$upper.95, 4)
      data$pvalue <- round(data$pvalue, 4)
      
      data1 <- data %>% dplyr::filter(HR > 1)
      data2 <- data %>% dplyr::filter(HR < 1)
      hrtable <- rbind(  c("Increase in Hazard", NA, NA, NA, NA, NA),
                         data1,  c("Reduce in Hazard", NA, NA, NA, NA, NA), 
                         data2)
      # 将要在图中展示的文本
      tabletext <- cbind(  
        c("Gene", "Increase in Hazard", data1$gene, "Reduce in Hazard", data2$gene),  
        c("P-value", NA, data1$pvalue, NA, data2$pvalue), 
        c("Hazard Ratio", NA, data1$`HR`, NA, data2$`HR`))
      which(tabletext[, 1] == "Reduction in Hazard")
      nrow(tabletext) + 1
      quantile(data$lower)
      quantile(data$upper)
      pdf("./04_sur/mcox_forestplot.pdf", width =6, height = 6) 
      
      
      ##所有tabletext数值取4位有效数字
      
      p<-forestplot(  
        labeltext = tabletext, 
        # 图中的文本  
        mean = c(NA, as.numeric(hrtable$HR)), 
        # HR  
        lower = c(NA, as.numeric(hrtable$lower)), 
        upper = c(NA, as.numeric(hrtable$upper)),
        graph.pos = 4, 
        # 图在表中的列位置  graphwidth = unit(.4, "npc"), 
        # 图在表中的宽度比例  fn.ci_norm = "fpDrawDiamondCI", 
        # box类型选择钻石
        col = fpColors(box = "#CD534CFF",lines="#0073C2FF", zero = "#0073C2FF"),
        #box颜色  
        boxsize = 0.7, 
        # box大小根据样本量设置  
        lwd.ci = 3, 
        ci.vertices.height = 0.125,  
        ci.vertices = T, 
        # 置信区间用线宽、高、型  
        zero = 1, 
        # zero线横坐标  
        lwd.zero = 2,  
        xticks = seq(0, 8, 1), 
        # 横坐标刻度根据需要可随意设置  
        lwd.xaxis = 2, # X轴线宽  
        xlab = "Hazard Ratio",  
        hrzl_lines = list("1" = gpar(lwd = 2, col = "black"), 
                          # 第二行顶部加黑实线    
                          "2" = gpar(lwd = 1, col = "grey50", lty = 2),
                          # 第二行顶部加灰色虚线    
                          #"3" = gpar(lwd = 1, col = "grey50", lty = 2), 
                          #"4" = gpar(lwd = 1, col = "grey50", lty = 2),
                          # 第二行顶部加灰色虚线    
                          "9" = gpar(lwd = 1, col = "grey50", lty = 2), 
                          # 第九行顶部加灰色虚线    
                          "10" = gpar(lwd = 2, col = "black")  ),
        # 最后一行底部加黑线，""中数字为nrow(tabletext) + 1  
        lineheight = unit(0.8, "cm"),   
        txt_gp = fpTxtGp(label = gpar(cex = 1), 
                         # 各种字体大小设置    
                         ticks = gpar(cex = 1.25),    
                         xlab = gpar(cex = 1.85),    
                         title = gpar(cex = 1.5)),  
        clip = c(0, 2),  
        colgap = unit(0.15, "cm"),   
        mar = unit(rep(1.25, times = 4), "cm"), 
        #  new_page = F # 是否新页
      )
      p
      
      dev.off()
    }
    
    
  }
  
  #sur
  {
    get_km(model_score$TCGA,methods="sur",output = './04_sur', data.name = 'TCGA_LUAD')
    ##TCGA_LUAD cutoffs: -0.178541104534964
    #n events median 0.95LCL 0.95UCL
    #group=High 244    109   1194     995    1528
    #group=Low  260     72   1778    1499    2617
    median_time_high_risk <- surv_median(surv_fit, group = "high")
    
    ##高低风险组按照-0.178541104534964分组
    model_score$TCGA$risk_score_level <- ifelse(model_score$TCGA$risk_score > -0.178541104534964,"High","Low")
  } 
  

  
}


##GSE31210

{
  {
    expr_table <- combined_data[combined_data$Dataset == "GSE31210",]
    
    #expr_table <- expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1")]
    #expr_table = log2(expr_table + 0.001)
    
    
    #expr_table <- dat_zscore$GSE31210[rownames(dat_zscore$GSE31210) %in% c("PTGDS", "PTGDR", "HPGDS", "SLCO2A1"), ]
    #expr_table<-t(expr_table)
    #expr_table <- as.data.frame(scale(expr_table))
    
    clin <- cln$GSE31210
    colnames(clin)[colnames(clin) == "OS"] <- "Status"
    ##colname time  to OS
    colnames(clin)[colnames(clin) == "days.before.death"] <- "OS"
    # numeric
    clin$OS <- as.numeric(clin$OS)
    rownames(clin)<-clin$geo_accession
    head(clin)
    data <- cbind(clin[,c('OS', 'Status')], expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1"),drop = FALSE])
  }
  
  ## cox
  {
    
   
    #fit <- glm(Status ~ PTGDS + PTGDR, family=binomial(link = "logit"), data = data_TCGA)
    data$risk_score <- with(data, coef(fit_cox)[1] * PTGDR + coef(fit_cox)[2] * HPGDS, coef(fit_cox)[3] * PTGDS)
    #data$risk_score_1 <- with(data, coef(fit_cox_mean)[1] * PTGmean + coef(fit_cox_mean)[2] * HPGDS)
    #data$risk_score_2 <- with(data, coef(fit_cox_product)[1] * PTGproduct + coef(fit_cox_product)[2] * HPGDS)

  }
  model_score$GSE31210<-data
  
  #sur
  {
    get_km(model_score$GSE31210,methods="sur",output = './04_sur', data.name = 'GSE31210')
    #GSE31210 7cutoffs: -0.84423637337437
    #Call: survfit(formula = Surv(OS, Status) ~ group, data = osdf)
    
    #n events median 0.95LCL 0.95UCL
    #group=High 193     32     NA      NA      NA
    #group=Low   33      3     NA    2661      NA
    model_score$GSE31210$risk_score_level <- ifelse(model_score$GSE31210$risk_score > -0.84423637337437,"High","Low")
    
    
     }


##GSE37745
{
  {
    expr_table <- combined_data[combined_data$Dataset == "GSE37745",]
    
    #expr_table <- expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1")]
    #expr_table = log2(expr_table + 0.001)
    
    #expr_table <- dat_zscore$GSE37745[rownames(dat_zscore$GSE37745) %in% c("PTGDS", "PTGDR", "HPGDS", "SLCO2A1"), ]
    #expr_table <- GSE37745[rownames(GSE37745) %in% c("PTGDS", "PTGDR", "HPGDS", "SLCO2A1"), ]
    #expr_table<-t(expr_table)
    #expr_table <- as.data.frame(scale(expr_table))
    
    clin <- cln$GSE37745
    colnames(clin)[colnames(clin) == "OS"] <- "Status"
    ##colname time  to OS
    colnames(clin)[colnames(clin) == "days.to.determined.death"] <- "OS"
    clin$OS <- as.numeric(clin$OS)
    rownames(clin)<-clin$geo_accession
    head(clin)
    
    
    data <- cbind(clin[,c('OS', 'Status')], expr_table[,c("PTGDS","PTGDR","HPGDS","SLCO2A1"),drop = FALSE])
  }
  
  ## cox
  {
    
    data$risk_score <- with(data, coef(fit_cox)[1] * PTGDR + coef(fit_cox)[2] * HPGDS, coef(fit_cox)[3] * PTGDS)
    
  }
  model_score$GSE37745<-data
  
  #sur
  {
    get_km(model_score$GSE37745,methods="sur",output = './04_sur',cutoffs= NULL,data.name = 'GSE37745')
    
    #GSE37745 7cutoffs: -0.643924755201412
    #Call: survfit(formula = Surv(OS, Status) ~ group, data = osdf)
    
    #n events median 0.95LCL 0.95UCL
    #group=High 83     65   1205    1000    2060
    #group=Low  23     12   3545    1463      NA
    
    model_score$GSE37745$risk_score_level <- ifelse(model_score$GSE37745$risk_score > -0.643924755201412,"High","Low")
  }
  
  
}

  model_score$TCGA$Dataset <- "TCGA"
  model_score$GSE37745$Dataset <- "GSE37745"
  model_score$GSE31210$Dataset <- "GSE31210"
  model_score$GSE30219$Dataset <- "GSE30219"
saveRDS(model_score, file = "./Rds/model_score.rds")

combined_data_1 <- rbind(model_score$TCGA, model_score$GSE37745, model_score$GSE31210)

rownames(combined_data) <- gsub("_T_LUAD","",rownames(combined_data))
combined_data <- combined_data[combined_data$Dataset != "GSE30219",]
{
  #makedata combined_data_1, combined_data 
  combined_all<-cbind(combined_data_1, combined_data[,c("HOT_score","Phenotype","PDL1")])
  
  saveRDS(combined_all, file = "./Rds/combined_all.rds")
  write.csv(combined_all, file = "./04_sur/combined_all.csv", row.names = F)
  
  
  ##画曲线
  

  library(plotROC)
  

  p <- ggplot(combined_all, aes(d = Phenotype, m = -risk_score, color = Dataset)) +
    geom_roc() +
    theme_bw() +
    labs(title = "Risk score ROC Curves for Multiple Datasets", x = "1-Specificity", y = "Sensitivity") +
    scale_color_brewer(palette = "Set1")# 设置颜色
  
  auc_table<-data.frame()
  for(dataset in unique(combined_all$Dataset)){
    auc_data <- combined_all %>%
      filter(Dataset == dataset)
    # 在每个数据集上计算 AUC
    roc_curve <- roc(auc_data$Phenotype, auc_data$risk_score)
    auc_value <- auc(roc_curve)  # 计算 AUC
    auc_table<-rbind(auc_table,data.frame(dataset,auc_value))
  }
  
  #dataset auc_value
  #1     TCGA 0.7000789
  #2 GSE37745 0.6385027
  #3 GSE31210 0.6568581
  
  
  ggsave(p, file = paste0( "./03_calcPhenotype/hot_riskscore_ROC.pdf"),
         bg = "white", width = 7, height = 6)
}


##KM 时间依赖性ROC-
cutoff=365
cutoff1=1095
cutoff2=1828


library(survivalROC)
pdf("./04_sur/tcga_KM_ROC_plot.pdf", width = 6, height = 6) # 设置PDF文件的宽度和高度（单位
ROC <- survivalROC(Stime = model_score$TCGA$OS,                   
                   status = model_score$TCGA$Status,                   
                   marker = model_score$TCGA$risk_score,                   
                   predict.time = cutoff,                   
                   method = 'KM')

plot(ROC$FP,ROC$TP,type = 'l',col='#FFFFCC',
     xlim = c(0,1),ylim = c(0,1),lwd=3,     
     xlab = 'Fales Postitve Rate',ylab = 'True Positive Rate')
abline(0,1,col='gray',lty=2)
aucText1=paste0("1 years"," (AUC=",sprintf("%.3f",ROC$AUC),")")

ROC1 <- survivalROC(Stime = model_score$TCGA$OS,                   
                    status = model_score$TCGA$Status,                   
                    marker = model_score$TCGA$risk_score,                   
                    predict.time = cutoff1,                   
                    method = 'KM')
lines(ROC1$FP, ROC1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
      col="#EFB24C",lwd = 2)
aucText2=paste0("3 years"," (AUC=",sprintf("%.3f",ROC1$AUC),")")


ROC2 <- survivalROC(Stime = model_score$TCGA$OS,                   
                    status = model_score$TCGA$Status,                   
                    marker = model_score$TCGA$risk_score,                   
                    predict.time = cutoff2,                   
                    method = 'KM')

lines(ROC2$FP, ROC2$TP, type="l", xlim=c(0,1), ylim=c(0,1),
      col="red",lwd = 2)
aucText3=paste0("5 years"," (AUC=",sprintf("%.3f",ROC2$AUC),")")
legend("bottomright", c(aucText1,aucText2,aucText3),       lwd=2,bty="n",col=c("#FFFFCC","#EFB24C","#E31A1C"),cex=1.2)


dev.off()