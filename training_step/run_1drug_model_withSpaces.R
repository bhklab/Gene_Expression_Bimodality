library(caret)
library(Biobase)
#.libPaths("/mnt/work1/users/bhklab/Rlib")
library(PharmacoGx)
library(mRMRe)

#.libPaths("/mnt/work1/users/home2/wbaalawi/R/x86_64-pc-linux-gnu-library/3.5")
library(rlobico)
library(mCI)


cat("\014")
rm(list = ls())
options(warn = FALSE)
options(LD_LIBRARY_PATH="/Applications/CPLEX_Studio_Community201/cplex/bin/x86-64_osx")
seed = 12345




args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Needed arguments: drug dataset K_CV numOfSolutions outputPath dictionary data", call.=FALSE)
}

drug <- args[1]
dataSet <- args[2]
k_CV = as.numeric(args[3])
numOfSolutions = as.numeric(args[4])
outputPath <- args[5]



sharedPath <- try(system("echo $AZ_BATCH_NODE_SHARED_DIR/data", intern = TRUE))
setwd(paste(sharedPath,"/AzureRun",sep = ""))

DrugDictionary <- readRDS(args[6])

dataObj <- args[7]

drug <- DrugDictionary[drug]
print(getwd())

source("getPredictions_LOBICO.R")
source("crossValidate_Lobico_updated_Rversion.R")
source("resample_Cindex_PValue.R")
source("lobico.R")
source("showformula.R")
source("getsolution.R")
source("runLobicoRWrapper.R")


load(dataObj)

dataSet_PSet <- readRDS(paste(dataSet,".rds",sep = ""))

newSet <- finalSetOfGenes

AAC <- summarizeSensitivityProfiles(dataSet_PSet,sensitivity.measure = "aac_recomputed",fill.missing = F)
expressionData = ccle2_binary[,finalSetOfGenes]

A <- tryCatch(crossValidate_Lobico_Outer_final_ensemble(AAC,expressionData = expressionData,drug = drug,newSet = newSet,numOfSolutions = numOfSolutions,k_CV = k_CV), error=function(e) return(e))

saveRDS(A,file = paste(outputPath,"/",drug,"_",dataSet,"_CV",k_CV,"NS",numOfSolutions,".rda",sep = ""))






