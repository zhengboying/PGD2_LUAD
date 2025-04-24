#file.symlink(from = "/Users/zhengboying/Documents/data/DATABASE/oncoPredict" , to = "../../Database")
# input : A matrix with rownames() as genes and colnames() as samples.
#df<-"../../Database/oncoPredict/TrainingData/"
#df=as.matrix(read.table(Inputdata, header=TRUE, row.names=1))
calcPhenotype_pipeline <- function(TrainingData = "GDSC1Data",##GDSC1Data, GDSC2Data, CTRP2Data
                                   Outputpath = paste0(project_path,"/03_calcPhenotype/GDSC1Data_calcPhenotype_Output/"),
                                   Inputdata = df,
                                   TrainingDatapath = '/Users/zhengboying/Documents/project_demo/Database/oncoPredict/TrainingData/',
                                   batchCorrect = "eb",##eb(ComBat), qn, standardize, none
                                   cc=TRUE,##correlation coefficients for biomarker discovery
                                   rsq=TRUE,##output the R^2 values for the data you train on from true and predicted values
                                   
                                   #pcr must be FALSE if cc is TRUE
                                   pcr = FALSE,
                                   report_pc = FALSE,
                                   percent = 80 ##percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
                                   
){
  library(oncoPredict)
  set.seed(12345)
  dir.create(Outputpath, showWarnings = FALSE)
  ##parameters
  {
    
    ##load traning data
    {
      if(TrainingData == "GDSC1Data"){
        #Read GDSC training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
        trainingExprData=readRDS(paste0(TrainingDatapath,"GDSC1_Expr.rds"))
        dim(trainingExprData) #17419 958 
        #Read GDSC1 response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
        trainingPtype = readRDS(paste0(TrainingDatapath,"GDSC1_Res.rds"))
        dim(trainingPtype) #958 367 For GDSC1
        #trainingPtype<-trainingPtype[,1:2] #Just 2 drugs for the vignette. 
        
        #e^IC50 since the IC50s are actual ln values/log transformed already
        trainingPtype<-exp(trainingPtype) 
      }
      else if(TrainingData == "GDSC2Data"){
        #Read GDSC training expression data. rownames() are genes and colnames() are samples.
        trainingExprData=readRDS(file=paste0(TrainingDatapath,'GDSC2_Expr.rds'))
        dim(trainingExprData) #17419 805
        #Read GDSC2 response data. rownames() are samples, colnames() are drugs. 
        trainingPtype = readRDS(file=paste0(TrainingDatapath,'GDSC2_Res.rds'))
        dim(trainingPtype) #805 198
        
        #e^IC50 since the IC50s are actual ln values/log transformed already
        trainingPtype<-exp(trainingPtype) 
      }
      else if(TrainingData == "CTRP2Data"){
        #Read CTRP training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
        trainingExprData = readRDS(file = "CTRP2_Expr_TPM.rds")
        dim(trainingExprData) #51847 829
        #Read CTRP training response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
        trainingPtype = readRDS(file = "CTRP2_Res.rds")
        dim(trainingPtype) #829 545 
        
      }else{
        stop("Please select a training dataset from GDSC1Data, GDSC2Data, or CTRP2Data.")
      }
      
      
    }
    
    ## calcPhenotype parameters
    {
      #batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
      #"eb" is good to use when you use microarray training data to build models on microarray testing data.
      #"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
      batchCorrect<-"eb"
      
      #Determine whether or not to power transform the phenotype data.
      #Default is TRUE.
      powerTransformPhenotype<-TRUE
      
      #Determine percentage of low varying genes to remove.
      #Default is 0.2 (seemingly arbitrary).
      removeLowVaryingGenes<-0.2
      
      #Determine method to remove low varying genes.
      #Options are 'homogenizeData' and 'rawData'
      #homogenizeData is likely better if there is ComBat batch correction, raw data was used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data.
      removeLowVaringGenesFrom<-"homogenizeData"
      
      #Determine the minimum number of training samples required to train on.
      #Note: this shouldn't be an issue if you train using GDSC or CTRP because there are many samples in both training datasets.
      #10, I believe, is arbitrary and testing could be done to get a better number.
      minNumSamples=10
      
      #Determine how you would like to deal with duplicate gene IDs.
      #Sometimes based on how you clean the data, there shouldn't be any duplicates to deal with.
      #Options are -1 for ask user, 1 for summarize by mean, and 2 for disregard duplicates
      selection<- 1
      
      #Determine if you'd like to print outputs.
      #Default is TRUE.
      printOutput=TRUE
      
      #Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
      #Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
      pcr=pcr
      
      #Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
      report_pc = report_pc
      
      #Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
      #These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
      cc=cc
      
      #Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
      #These values represent the percentage in which the optimal model accounts for the variance in the training data.
      #Options are 'TRUE' and 'FALSE'.
      rsq=rsq
      
      #Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
      percent=percent
    }
    ## run calcPhenotype
    {
      savedir<-setwd(Outputpath)
      dim(Inputdata)
      calcPhenotype(trainingExprData=trainingExprData,
                    trainingPtype=trainingPtype,
                    testExprData=Inputdata,
                    batchCorrect=batchCorrect,
                    powerTransformPhenotype=powerTransformPhenotype,
                    removeLowVaryingGenes=removeLowVaryingGenes,
                    minNumSamples=minNumSamples,
                    selection=selection,
                    printOutput=printOutput,
                    pcr=pcr,
                    removeLowVaringGenesFrom=removeLowVaringGenesFrom,
                    report_pc=report_pc,
                    cc=cc,
                    percent=percent,
                    rsq=rsq)
      
      #If pcr is performed, you can view a drug's first two principal components (and so on) using the code below. 
      #View(load('./calcPhenotype_Output/Vinblastine_1004.RData'))
      #View(pcs[,1,1]) #The first pc. 
      #View(pcs[,1,2]) #The second pc.
    }
  }
}