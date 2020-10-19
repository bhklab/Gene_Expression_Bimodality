gCSI <- readRDS("~/Projects/final_psets/gCSI_2018.rds")
GDSC1000 <- readRDS("~/Projects/final_psets/GDSCv2.rds")

AUC_gCSI <- summarizeSensitivityProfiles(gCSI,sensitivity.measure = "aac_recomputed",fill.missing = F)
AUC_GDSC1000 <- summarizeSensitivityProfiles(GDSC1000,sensitivity.measure = "aac_recomputed",fill.missing = F)

commonDrugsGDSC_gCSI <- Reduce(intersect,x = list(rownames(AUC_gCSI)
                                                  ,rownames(AUC_GDSC1000)
                                                  ,rownames(output_ensemble_CTRPv2_5CV3NS_evaluation)[output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"]>=mciTh]))



#######################################
# validate on GCSI


gCSIRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(gCSI,mDataType = "rnaseq",fill.missing = F)))
gCSITissues <- gCSI@cell[rownames(gCSIRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(gCSITissues$tissueid=="lung")
gCSITissues_WithoutHAEMATO <- gCSITissues[ibx,,drop=F]
gCSIRnaSeq <- gCSIRnaSeq_ALL[rownames(gCSITissues_WithoutHAEMATO),]

gcsi <- gCSIRnaSeq[,finalSetOfGenes]
gcsi_binary <- getBinaryValues(gcsi,cutoffs_ccle[finalSetOfGenes])


drugs_common <- intersect(rownames(AUC_gCSI),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_gCSI <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=gcsi_binary[,finalSetOfGenes],AUCs_Drugs=AUC_gCSI)




names(results_final_common_gCSI) <- final_drugs_common

results_final_common_evaluation_gCSI <- lapply(results_final_common_gCSI, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_gCSI <- do.call(rbind,results_final_common_evaluation_gCSI)


data <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_evaluation[rownames(results_final_common_evaluation_gCSI),"mCI"]
                   ,"gCSI_lung"=results_final_common_evaluation_gCSI[,"mCI"],"drug"=rownames(results_final_common_evaluation_gCSI),stringsAsFactors=F)


saveRDS(data,"globalCTRPv2_lungGCSI_predictions.rds")


###############################################################
###############################################################
###############################################################
###############################################################

load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.8_LUNG.RData",verbose = T)

gCSIRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(gCSI,mDataType = "rnaseq",fill.missing = F)))
gCSITissues <- gCSI@cell[rownames(gCSIRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(gCSITissues$tissueid=="lung")
gCSITissues_WithoutHAEMATO <- gCSITissues[ibx,,drop=F]
gCSIRnaSeq <- gCSIRnaSeq_ALL[rownames(gCSITissues_WithoutHAEMATO),]

gcsi <- gCSIRnaSeq[,finalSetOfGenes]
gcsi_binary <- getBinaryValues(gcsi,cutoffs_ccle2[finalSetOfGenes])



#######################################
# validate on GCSI

drugs_common_lung <- intersect(rownames(AUC_gCSI),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common_lung,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common_lung,])[ibx]


results_final_common_gCSI <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG_final,experssion=gcsi_binary[,finalSetOfGenes],AUCs_Drugs=AUC_gCSI)




names(results_final_common_gCSI) <- final_drugs_common

results_final_common_evaluation_gCSI <- lapply(results_final_common_gCSI, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_gCSI <- do.call(rbind,results_final_common_evaluation_gCSI)


data <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation_gCSI),"mCI"]
                   ,"gCSI_lung_specific"=results_final_common_evaluation_gCSI[,"mCI"],"drug"=rownames(results_final_common_evaluation_gCSI),stringsAsFactors=F)



data_panCancer <- readRDS("globalCTRPv2_lungGCSI_predictions.rds")

common_2 <- intersect( rownames(data),rownames(data_panCancer))

data_final = data.frame("drugs" = common_2, "LungSpecificRules"=data[common_2,"gCSI_lung_specific"]
                        ,"PanCancerRules"=data_panCancer[common_2,"gCSI_lung"],stringsAsFactors = F)

library(reshape2)


data_final1 <- melt(data_final,id.vars ="drugs",variable.name = "RulesType")
data_final1 <- data_final1[!is.na(data_final1$value),]
data_final1$RulesType <- gsub("Rules","",data_final1$RulesType)

pdf("lungSpecificVSpanCancerRules_gCSI.pdf",height = 5)
ggplot(data_final1,aes(fill=RulesType, y=value, x=drugs)) +
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Validating on gCSI") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


############################################
# validate on GDSC

GDSCRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(GDSC1000,mDataType = "rnaseq",fill.missing = F)))
GDSCTissues <- GDSC1000@cell[rownames(GDSCRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(GDSCTissues$tissueid=="lung")
GDSCTissues_WithoutHAEMATO <- GDSCTissues[ibx,,drop=F]
GDSCRnaSeq <- GDSCRnaSeq_ALL[rownames(GDSCTissues_WithoutHAEMATO),]

GDSC <- GDSCRnaSeq[,finalSetOfGenes]
GDSC_binary <- getBinaryValues(GDSC,cutoffs_ccle2[finalSetOfGenes])



#######################################
# validate on GDSC

drugs_common_lung <- intersect(rownames(AUC_GDSC1000),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common_lung,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common_lung,])[ibx]


results_final_common_GDSC <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG_final,experssion=GDSC_binary[,finalSetOfGenes],AUCs_Drugs=AUC_GDSC1000)




names(results_final_common_GDSC) <- final_drugs_common

results_final_common_evaluation_GDSC <- lapply(results_final_common_GDSC, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_GDSC <- do.call(rbind,results_final_common_evaluation_GDSC)


data <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation_GDSC),"mCI"]
                   ,"GDSC_lung_specific"=results_final_common_evaluation_GDSC[,"mCI"],"drug"=rownames(results_final_common_evaluation_GDSC),stringsAsFactors=F)


saveRDS(data,"lungCTRPv2_lungGDSC_predictions.rds")



############################################
# validate on GDSC (Pan cancer)
load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.75.RData")


GDSCRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(GDSC1000,mDataType = "rnaseq",fill.missing = F)))
GDSCTissues <- GDSC1000@cell[rownames(GDSCRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(GDSCTissues$tissueid=="lung")
GDSCTissues_WithoutHAEMATO <- GDSCTissues[ibx,,drop=F]
GDSCRnaSeq <- GDSCRnaSeq_ALL[rownames(GDSCTissues_WithoutHAEMATO),]

GDSC <- GDSCRnaSeq[,finalSetOfGenes]
GDSC_binary <- getBinaryValues(GDSC,cutoffs_ccle[finalSetOfGenes])



#######################################
# validate on GDSC

drugs_common <- intersect(rownames(AUC_GDSC1000),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_GDSC <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=GDSC_binary[,finalSetOfGenes],AUCs_Drugs=AUC_GDSC1000)




names(results_final_common_GDSC) <- final_drugs_common

results_final_common_evaluation_GDSC <- lapply(results_final_common_GDSC, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_GDSC <- do.call(rbind,results_final_common_evaluation_GDSC)


data_panCancer <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_evaluation[rownames(results_final_common_evaluation_GDSC),"mCI"]
                   ,"GDSC_lung_specific"=results_final_common_evaluation_GDSC[,"mCI"],"drug"=rownames(results_final_common_evaluation_GDSC),stringsAsFactors=F)


data_lung <- readRDS("lungCTRPv2_lungGDSC_predictions.rds")


common_2 <- intersect( rownames(data_lung),rownames(data_panCancer))

data_final = data.frame("drugs" = common_2, "LungSpecificRules"=data_lung[common_2,"GDSC_lung_specific"]
                        ,"PanCancerRules"=data_panCancer[common_2,"GDSC_lung_specific"],stringsAsFactors = F)

library(reshape2)

data_final1 <- melt(data_final,id.vars ="drugs",variable.name = "RulesType")
data_final1 <- data_final1[!is.na(data_final1$value),]
data_final1$RulesType <- gsub("Rules","",data_final1$RulesType)

pdf("lungSpecificVSpanCancerRules_GDSCv2.pdf",height = 5)
ggplot(data_final1,aes(fill=RulesType, y=value, x=drugs)) +
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Validating on GDSCv2") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


########################################################
########################################################
########################################################


load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.75.RData",verbose = T)

ccleLung <- intersect(rownames(CCLE@cell)[CCLE@cell$tissueid == "lung"],rownames(ccle2_binary))

drugs_common <- intersect(rownames(output_ensemble_CTRPv2_5CV3NS_evaluation),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_CTRPv2_PanCancer_lung <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=ccle2_binary[ccleLung,],AUCs_Drugs=AAC_CTRPv2)




names(results_final_common_CTRPv2_PanCancer_lung) <- final_drugs_common

results_final_common_evaluation_CTRPv2_PanCancer_lung <- lapply(results_final_common_CTRPv2_PanCancer_lung, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_CTRPv2_PanCancer_lung <- do.call(rbind,results_final_common_evaluation_CTRPv2_PanCancer_lung)


data_panCancer <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_evaluation[rownames(results_final_common_evaluation_CTRPv2_PanCancer_lung),"mCI"]
                             ,"CTRPv2_lung_PanCancer"=results_final_common_evaluation_CTRPv2_PanCancer_lung[,"mCI"],"drug"=rownames(results_final_common_evaluation_CTRPv2_PanCancer_lung),stringsAsFactors=F)

saveRDS(data_panCancer,"panCancerCTRPv2_lungCTRPv2_predictions.rds")


###################################################
###################################################
###################################################



load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.8_LUNG.RData",verbose = T)

ccleLung <- intersect(rownames(CCLE@cell)[CCLE@cell$tissueid == "lung"],rownames(ccle2_binary))

#drugs_common <- intersect(rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))
#ibx <- which(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,"mCI"]>mciTh)
#final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,])[ibx]
final_drugs_common <- unionDrugs

results_final_common_CTRPv2_lung_specfic <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  print(x)
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  if(!is.list(output_ensemble_CTRPv2_5CV3NS_LUNG[[x]][[2]]))
    return(cbind("Vote"=NA,"Obs"=NA))
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    
    SolMat <- output[[x]][[3]][[y]][[2]]
    if(is.null(SolMat))
      return(NA)
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  
  if(is.na(final_results[[1]]))
    return(cbind("Vote"=NA,"Obs"=NA))
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG,experssion=ccle2_binary[ccleLung,],AUCs_Drugs=AAC_CTRPv2)




names(results_final_common_CTRPv2_lung_specfic) <- final_drugs_common

results_final_common_evaluation_CTRPv2_lung_specific <- lapply(results_final_common_CTRPv2_lung_specfic, function(x){
  
  
  if(is.na(x[1,1]))
    results <- c("CI"=NA,"CI.pval"=NA,"mCI"=NA,"mCI.pval"=NA)
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation_CTRPv2_lung_specific <- do.call(rbind,results_final_common_evaluation_CTRPv2_lung_specific)

common_1 <- intersect(rownames(results_final_common_evaluation_CTRPv2_lung_specific),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))

data <- data.frame("CTRPv2_trainging"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[common_1,"mCI"]
                             ,"CTRPv2_lung_Specific"=results_final_common_evaluation_CTRPv2_lung_specific[common_1,"mCI"],"drug"=common_1,stringsAsFactors=F)

saveRDS(data,"lungCTRPv2_lungCTRPv2_predictions_all.rds")

data_panCancer <- readRDS("panCancerCTRPv2_lungCTRPv2_predictions.rds")


common_2 <- intersect(rownames(data),rownames(data_panCancer))
unionDrugs <- union(rownames(data),rownames(data_panCancer))
saveRDS(unionDrugs,"unionDrugs_PancancerVsLungSpecific.rds")
LungSpecificOnly <- setdiff(rownames(data),rownames(data_panCancer))
PanCancerOnly <- setdiff(rownames(data_panCancer),rownames(data))

library(VennDiagram)
pdf("CTRPv2_LungVSpanCancerRules_vennDiagram.pdf")
draw.pairwise.venn(area1 = dim(data_panCancer)[1], area2 = dim(data)[1]
                   , cross.area = length(common_2)
                   , category = c("Pan-cancer","lung-specific")
                   ,ext.text = T,margin = c(0.1,0.1,0.1,0.1)
                   )
dev.off()

data_final = data.frame("drugs" = common_2, "LungSpecific"=data[common_2,"CTRPv2_lung_Specific"]
                        ,"PanCancer"=data_panCancer[common_2,"CTRPv2_lung_PanCancer"],stringsAsFactors = F)

library(reshape2)

data_final1 <- melt(data_final,id.vars =c("drugs"),variable.name = "RulesType")
data_final1 <- data_final1[!is.na(data_final1$value),]
data_final1$RulesType <- gsub("Rules","",data_final1$RulesType)


colnames(data_final)

pdf("lungSpecificVSpanCancerRules_CTRPv2.pdf",width = 5,height = 5)
ggplot(data_final,aes(LungSpecific,PanCancer)) + 
  geom_point() +
  geom_abline(intercept =0 , slope = 1,col="lightgray",lty=2) +
  ylim(c(0.6,1)) + xlim(c(0.6,1)) + 
  xlab("Lung-specific rules [rCI]") + ylab("Pan-cancer rules [rCI]") +
  ggtitle("Lung-specific vs Pan-cancer rules [Lung - CTRPv2]") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


colnames(data_final1)[3] <- "rCI"
data_final1$RulesType <- factor(data_final1$RulesType)
ggplot(data_final1,aes(RulesType,rCI)) +
  geom_boxplot() #+
  
  






