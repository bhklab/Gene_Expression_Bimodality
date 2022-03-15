crossValidate_Lobico_Outer_final_ensemble <- function(AUCs_Drugs, expressionData, drug, k_CV=10, newSet,numOfSolutions=3) {
  
  
  AUCs_Drug <- AUCs_Drugs[drug,,drop=F]
  AUCs_Drug <- AUCs_Drug[,!is.na(AUCs_Drug)]
  cellLines_drugTested <- names(AUCs_Drug)
  cellLines_RNA <- rownames(expressionData)
  
  commonCellLines <- intersect(cellLines_drugTested,cellLines_RNA)
  
  expressionData_sub <- expressionData[commonCellLines,newSet,drop=F]
  AUCs_Drug <- AUCs_Drug[commonCellLines,drop=F]
  
  
  Data_LOBICO <- data.frame('Cell Lines'=commonCellLines ,'drug_AUC'=AUCs_Drug[commonCellLines,drop=F],expressionData_sub[commonCellLines,,drop=F],stringsAsFactors = F)
  # Data_LOBICO <- as.data.frame(Data_LOBICO,stringsAsFactors=F)
  #  lapply(2:dim(Data_LOBICO)[2],function(r,Data_LOBICO){
  #    class(Data_LOBICO[,r]) <<- "numeric"
  #  },Data_LOBICO=Data_LOBICO)
  formulaOut <- list()
  output_Final <- NULL
  
  output_Final <- NULL
  set.seed(1245)
  flds <- createFolds(Data_LOBICO[,"drug_AUC"], k = k_CV, list = F, returnTrain = F)
  
  formulaOut <- list()
  
  
  Sol_out <- list()
  
  Pred_Output <- matrix(data = NA,nrow = 1,ncol = numOfSolutions+2)
  start.time = proc.time()
  for (i in 1:k_CV){
    print(paste(">>>>>>>>>>>>>>>>>>%%%",i,"%%%<<<<<<<<<<<<<<<<"))
    
    formulaOut[[i]] <- ""
    testInd <- which(flds==i)
    trainInd <- which(flds!=i)
    trainData <- Data_LOBICO[trainInd,2:dim(Data_LOBICO)[2],drop=F]
    testData <- Data_LOBICO[testInd,2:dim(Data_LOBICO)[2],drop=F]
    
    
    dd <- mRMR.data(data = trainData)
    test <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = numOfSolutions, feature_count = 30, method="bootstrap")
    
    mrmrOut <- lapply(1:numOfSolutions, function(p){
      formulaOut1 <- ""
      newSetTest <- c(apply(solutions(test)[[1]][,p,drop=F], 2, function(t, y) { return(y[t]) }, y=featureNames(dd)))
      namesIbx <- match(newSetTest,newSet)
      names(newSetTest) <- names(newSet)[namesIbx]
      
      ibx <- which(newSetTest == colnames(trainData)[1])
      if(length(ibx)!=0){
        newSetTest <- newSetTest[-ibx]
      }
      print(">>>>>>>>>>>>>>>>>>>>>> Inner CV start <<<<<<<<<<<<<<<<<<<<<")
      print(paste(">>>>>> mRMRe",p,"<<<<<<<"))
      A <- crossValidate_Lobico(AUCs_Drugs = trainData[,1,drop=F],expressionData = trainData[,2:dim(trainData)[2],drop=F],k_CV = k_CV,drug = drug,newSet = newSetTest)
      print(">>>>>>>>>>>>>>>>>>>>>> Inner CV End <<<<<<<<<<<<<<<<<<<<<")
      
      train_CV_output_drug <- A[[1]]
      K <- NULL
      M <- NULL
      Pvalue_annot <- NULL
      
      ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
      if(length(ibx)==1){
        K = train_CV_output_drug[ibx,"K"]
        M = train_CV_output_drug[ibx,"M"]
        Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
      }else{
        K = train_CV_output_drug[ibx[1],"K"]
        M = train_CV_output_drug[ibx[1],"M"]
        Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
      }
      
      Data_LOBICO_training <- data.frame('Cell Lines'=rownames(trainData) ,'drug_AUC'=trainData[,1],trainData[,newSetTest,drop=F],stringsAsFactors = F)
      
      #   source("~/Downloads/R-LOBICO/lobico.R")
      
      SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2],drop=F]),K = K,M = M)
      SolMat <- round(SolMat)
      
      colnames(SolMat) <- newSetTest
      
      
      return(list(A,SolMat,K,M,newSetTest))
      #      return(newSetTest)
    })
    
    print(">>>>>>>>>>>>>>>>>>>>>> End of modeling <<<<<<<<<<<<<<<<<<<<<")
    outputs_solutions <- lapply(mrmrOut,function(x){
      
      SolMat <- x[[2]]
      K <- x[[3]]
      M <- x[[4]]
      finalGeneSet <- x[[5]]
      
      if(all(SolMat==0)){
        output <- cbind(rep("0",length(as.numeric(testData[,"drug_AUC"]))),testData[,"drug_AUC"])
        return(output[,1])
      }
      
      #      write.table(showFormula(SolMat,K,M,newSet),file = "fluvastatin.CV.txt",quote = F,append = T)
      output <- getPredictions_LOBICO_outer(SolMat,testData[,c("drug_AUC",finalGeneSet),drop=F],K,M)
      
      return(output[,1])
      
    })
    print(">>>>>>>>>>>>>>>>>>>>>> End of Testing <<<<<<<<<<<<<<<<<<<<<")
    outputs_solutions <- do.call(cbind,outputs_solutions)
    
    outputs_solutions <- cbind(outputs_solutions,"Vote"=apply(outputs_solutions,1,function(x){
      round(mean(x))
    }))
    
    
    formulas_solutions <- lapply(mrmrOut,function(x){
      
      SolMat <- x[[2]]
      K <- x[[3]]
      M <- x[[4]]
      finalGeneSet <- x[[5]]
      
      return(showformula(SolMat,K,M,finalGeneSet))
      
    })
    
    
    formulaOut[[i]] <- formulas_solutions
    
    output <- cbind(outputs_solutions,"Obs"=testData[rownames(outputs_solutions),c("drug_AUC")])
    
    expressedCellLines <- rownames(output)[output[,"Vote"]==1]
    notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
    
    #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
    CI <- paired.concordance.index(observations = output[c(notExpressedCellLines,expressedCellLines),"Obs"],predictions = output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                   ,delta.pred = 0,delta.obs = 0,outx = F)
    
    mCI <- paired.concordance.index(observations = output[c(notExpressedCellLines,expressedCellLines),"Obs"],predictions = output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                    ,delta.pred = 0,delta.obs = 0.2,outx = F)
    
    
    Sol_out[[i]] <- list(CI["cindex"],CI["p.value"],mCI["cindex"],mCI["p.value"],mrmrOut)
    Pred_Output <- rbind(Pred_Output,output)
    
    
  }
  
  endtime <- proc.time() - start.time
  
  print(endtime)
  
  print(">>>>>>>>>>>>>>>>>>>>>> End of CV <<<<<<<<<<<<<<<<<<<<<")
  
  
  class(Pred_Output) <- "numeric"
  
  Pred_Output <- Pred_Output[-1,]
  
  expressedCellLines <- rownames(Pred_Output)[Pred_Output[,"Vote"]==1]
  notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,"Vote"]==0]
  
  mrmrOut_All <- NULL
  formulas_solutions <- NULL
  
  if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
    
    output_Final <- rbind(output_Final,c('CI'=0.5,'mCI'=0.5,'Pvalue.CI'=1,'Pvalue.mCI'=1,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = "",
                                         'CI.balanced'=0.5,'Pvalue.CI.balanced'=1,'Pvalue.balanced'=1,'Equal.N'=0,
                                         'mCI.balanced'=0.5,'Pvalue.mCI.balanced'=1))
    
    
  }else{
    
    
    CI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                   ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),"Obs"]
                                   ,outx = F,delta.pred = 0,delta.obs = 0)[c("cindex","p.value")]
    
    mCI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                    ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),"Obs"],delta.pred = 0,delta.obs = 0.2
                                    ,outx = F)[c("cindex","p.value")]
    
    Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,"Obs"], Pred_Output[expressedCellLines,"Obs"],alternative = "greater")$p.value
    Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,"Obs"], Pred_Output[expressedCellLines,"Obs"],alternative = "less")$p.value
    
    
    Pvalue_annot <- NULL
    
    if( Pvalue_greater < Pvalue_less){
      Pvalue_annot <- "greater"
    }else{
      Pvalue_annot <- "less"
    }
    
    data = Pred_Output
    balanced_CI <- resample_Cindex_PValue(Pred_Output[,"Vote"],Pred_Output[,"Obs"],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
    balanced_mCI <- resample_Cindex_PValue(Pred_Output[,"Vote"],Pred_Output[,"Obs"],permutation = 100,Pvalue_annot = Pvalue_annot
                                           ,delta.pred = 0,delta.obs = 0.2,outx=F)
    
    
    
    Data_LOBICO_All <- Data_LOBICO[,2:dim(Data_LOBICO)[2]]
    
    dd <- mRMR.data(data = Data_LOBICO_All)
    test <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = numOfSolutions, feature_count = 30,method="bootstrap")
     print("######## getting the final model ########")
    mrmrOut_All <- lapply(1:numOfSolutions, function(p){
      formulaOut1 <- ""
      newSetTest <- c(apply(solutions(test)[[1]][,p,drop=F], 2, function(t, y) { return(y[t]) }, y=featureNames(dd)))
      namesIbx <- match(newSetTest,newSet)
      names(newSetTest) <- names(newSet)[namesIbx]
      
      ibx <- which(newSetTest == colnames(Data_LOBICO_All)[1])
      if(length(ibx)!=0){
        newSetTest <- newSetTest[-ibx]
      }
      
      A <- crossValidate_Lobico(AUCs_Drugs = Data_LOBICO_All[,1,drop=F],expressionData = Data_LOBICO_All[,2:dim(Data_LOBICO_All)[2]],k_CV = k_CV, drug = drug,newSet = newSetTest)
      train_CV_output_drug <- A[[1]]
      K <- NULL
      M <- NULL
      Pvalue_annot <- NULL
      
      ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
      if(length(ibx)==1){
        K = train_CV_output_drug[ibx,"K"]
        M = train_CV_output_drug[ibx,"M"]
        Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
      }else{
        K = train_CV_output_drug[ibx[1],"K"]
        M = train_CV_output_drug[ibx[1],"M"]
        Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
      }
      
      Data_LOBICO_training <- data.frame('Cell Lines'=rownames(Data_LOBICO_All) ,'drug_AUC'=Data_LOBICO_All[,1],Data_LOBICO_All[,newSetTest],stringsAsFactors = F)
      
      #   source("~/Downloads/R-LOBICO/lobico.R")
      
      SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
      SolMat <- round(SolMat)
      
      colnames(SolMat) <- newSetTest
      
      
      return(list(A,SolMat,K,M,newSetTest))
      #      return(newSetTest)
    })
    
    
    formulas_solutions <- lapply(mrmrOut_All,function(x){
      
      SolMat <- x[[2]]
      K <- x[[3]]
      M <- x[[4]]
      finalGeneSet <- x[[5]]
      
      return(showformula(SolMat,K,M,finalGeneSet))
      
    })
    
    
    output_Final <- rbind(output_Final,c('CI'=CI[["cindex"]],'mCI'=mCI[["cindex"]],'Pvalue.CI'=CI[["p.value"]],'Pvalue.mCI'=mCI[["p.value"]],'Pvalue'=ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less),'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                                         'CI.balanced'=balanced_CI$Cindex,'Pvalue.CI.balanced'=balanced_CI$Pvalue,'Pvalue.balanced'=balanced_CI$Pvalue,'Equal.N'=balanced_CI$N,
                                         'mCI.balanced'=balanced_mCI$Cindex,'Pvalue.mCI.balanced'=balanced_mCI$PvalueCI))
    # write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
    #     par(font=2)
    #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
    #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
    print("######## Obtained he final model ########")
  }
  
  print("######## Done ########")
  return(list(output_Final,formulas_solutions,mrmrOut_All,Pred_Output))
}



crossValidate_Lobico_Outer_final_single <- function(AUCs_Drugs, expressionData, drug, k_CV=10, newSet,numOfSolutions=3) {
  
  
  
  
  AUCs_Drug <- AUCs_Drugs[drug,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  cellLines_drugTested <- names(AUCs_Drug)
  cellLines_RNA <- rownames(expressionData)
  
  commonCellLines <- intersect(cellLines_drugTested,cellLines_RNA)
  
  expressionData_sub <- expressionData[commonCellLines,newSet]
  AUCs_Drug <- AUCs_Drug[commonCellLines]
  
  
  Data_LOBICO <- data.frame('Cell Lines'=commonCellLines ,'drug_AUC'=AUCs_Drug[commonCellLines],expressionData_sub[commonCellLines,],stringsAsFactors = F)
  # Data_LOBICO <- as.data.frame(Data_LOBICO,stringsAsFactors=F)
  #  lapply(2:dim(Data_LOBICO)[2],function(r,Data_LOBICO){
  #    class(Data_LOBICO[,r]) <<- "numeric"
  #  },Data_LOBICO=Data_LOBICO)
  formulaOut <- list()
  output_Final <- NULL
  
  output_Final <- NULL
  set.seed(1245)
  flds <- createFolds(Data_LOBICO[,"drug_AUC"], k = k_CV, list = F, returnTrain = F)
  
  formulaOut <- list()
  
  
  Sol_out <- list()
  
  Pred_Output <- matrix(data = NA,nrow = 1,ncol = 2)
  start.time = proc.time()
  for (i in 1:k_CV){
    print(paste(">>>>>>>>>>>>>>>>>>",i,"<<<<<<<<<<<<<<<<"))
    
    formulaOut[[i]] <- ""
    testInd <- which(flds==i)
    trainInd <- which(flds!=i)
    trainData <- Data_LOBICO[trainInd,2:dim(Data_LOBICO)[2]]
    testData <- Data_LOBICO[testInd,2:dim(Data_LOBICO)[2]]
    
    
    dd <- mRMR.data(data = trainData)
    test <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = numOfSolutions, feature_count = 10,method="bootstrap")
    
    mrmrOut <- lapply(1:numOfSolutions, function(p){
      formulaOut1 <- ""
      newSetTest <- c(apply(solutions(test)[[1]][,p,drop=F], 2, function(t, y) { return(y[t]) }, y=featureNames(dd)))
      namesIbx <- match(newSetTest,newSet)
      names(newSetTest) <- names(newSet)[namesIbx]
      
      ibx <- which(newSetTest == colnames(trainData)[1])
      if(length(ibx)!=0){
        newSetTest <- newSetTest[-ibx]
      }
      # print(">>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<")
      # print(i)
      # print(">>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<")
      # A <- crossValidate_Lobico(AUCs_Drugs = trainData[,1,drop=F],expressionData = trainData[,2:dim(trainData)[2]],drug = drug,newSet = newSetTest)
      # print(paste(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",i,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"))
      # train_CV_output_drug <- A[[1]]
      # K <- NULL
      # M <- NULL
      # Pvalue_annot <- NULL
      #
      # ibx <- which(train_CV_output_drug[,"CI"]==max(train_CV_output_drug[,"CI"]))
      # if(length(ibx)==1){
      #   K = train_CV_output_drug[ibx,"K"]
      #   M = train_CV_output_drug[ibx,"M"]
      #   Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
      # }else{
      #   K = train_CV_output_drug[ibx[1],"K"]
      #   M = train_CV_output_drug[ibx[1],"M"]
      #   Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
      # }
      #
      # Data_LOBICO_training <- data.frame('Cell Lines'=rownames(trainData) ,'drug_AUC'=trainData[,1],trainData[,newSetTest],stringsAsFactors = F)
      #
      # #   source("~/Downloads/R-LOBICO/lobico.R")
      #
      # SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
      # SolMat <- round(SolMat)
      #
      # colnames(SolMat) <- newSetTest
      #
      #
      #return(list(A,SolMat,K,M,newSetTest))
      return(newSetTest)
    })
    
    newSetTest <- unique(unlist(mrmrOut))
    
    A <- crossValidate_Lobico(AUCs_Drugs = trainData[,1,drop=F],expressionData = trainData[,newSetTest],drug = drug,newSet = newSetTest)
    print(paste(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",i,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"))
    train_CV_output_drug <- A[[1]]
    K <- NULL
    M <- NULL
    Pvalue_annot <- NULL
    
    ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
    if(length(ibx)==1){
      K = train_CV_output_drug[ibx,"K"]
      M = train_CV_output_drug[ibx,"M"]
      Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
    }else{
      K = train_CV_output_drug[ibx[1],"K"]
      M = train_CV_output_drug[ibx[1],"M"]
      Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
    }
    
    Data_LOBICO_training <- data.frame('Cell Lines'=rownames(trainData) ,'drug_AUC'=trainData[,1],trainData[,newSetTest],stringsAsFactors = F)
    
    #   source("~/Downloads/R-LOBICO/lobico.R")
    
    SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
    SolMat <- round(SolMat)
    
    colnames(SolMat) <- newSetTest
    
    if(all(SolMat==0)){
      output <- cbind(rep("0",length(as.numeric(testData[,"drug_AUC"]))),testData[,"drug_AUC"])
      return(output[,1])
    }
    finalGeneSet <- newSetTest
    output <- getPredictions_LOBICO_outer(SolMat,testData[,c("drug_AUC",finalGeneSet)],K,M)
    colnames(output) <- c("Vote","Obs")
    # SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
    #
    # outputs_solutions <- lapply(mrmrOut,function(x){
    #
    #   SolMat <- x[[2]]
    #   K <- x[[3]]
    #   M <- x[[4]]
    #   finalGeneSet <- x[[5]]
    #
    #   if(all(SolMat==0)){
    #     output <- cbind(rep("0",length(as.numeric(testData[,"drug_AUC"]))),testData[,"drug_AUC"])
    #     return(output[,1])
    #   }
    #
    #   #      write.table(showFormula(SolMat,K,M,newSet),file = "fluvastatin.CV.txt",quote = F,append = T)
    #   output <- getPredictions_LOBICO_outer(SolMat,testData[,c("drug_AUC",finalGeneSet)],K,M)
    #
    #   return(output[,1])
    #
    # })
    #
    # outputs_solutions <- do.call(cbind,outputs_solutions)
    #
    # outputs_solutions <- cbind(outputs_solutions,"Vote"=apply(outputs_solutions,1,function(x){
    #   round(mean(x))
    # }))
    #
    #
    # formulas_solutions <- lapply(mrmrOut,function(x){
    #
    #   SolMat <- x[[2]]
    #   K <- x[[3]]
    #   M <- x[[4]]
    #   finalGeneSet <- x[[5]]
    #
    #   return(showformula(SolMat,K,M,finalGeneSet))
    #
    # })
    #
    #
    # formulaOut[[i]] <- formulas_solutions
    #
    # output <- cbind(outputs_solutions,"Obs"=testData[rownames(outputs_solutions),c("drug_AUC")])
    
    expressedCellLines <- rownames(output)[output[,"Vote"]==1]
    notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
    
    #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
    CI <- paired.concordance.index(observations = output[c(notExpressedCellLines,expressedCellLines),"Obs"],predictions = output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                   ,delta.pred = 0,delta.obs = 0,outx = F)
    
    mCI <- paired.concordance.index(observations = output[c(notExpressedCellLines,expressedCellLines),"Obs"],predictions = output[c(notExpressedCellLines,expressedCellLines),"Vote"]
                                    ,delta.pred = 0,delta.obs = 0.2,outx = F)
    
    
    Sol_out[[i]] <- list(CI["cindex"],CI["p.value"],mCI["cindex"],mCI["p.value"],K,M,mrmrOut)
    Pred_Output <- rbind(Pred_Output,output)
    
    
  }
  
  endtime <- proc.time() - start.time
  
  print(endtime)
  
  
  
  
  class(Pred_Output) <- "numeric"
  
  Pred_Output <- Pred_Output[-1,]
  
  expressedCellLines <- rownames(Pred_Output)[Pred_Output[,"Vote"]==1]
  notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,"Vote"]==0]
  
  if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
    
    output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=0.5,'mCI'=0.5,'Pvalue.CI'=1,'Pvalue.mCI'=1,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = "",
                                         'CI.balanced'=0.5,'Pvalue.CI.balanced'=1,'Pvalue.balanced'=1,'Equal.N'=0,
                                         'mCI.balanced'=0.5,'Pvalue.mCI.balanced'=1))
    next
  }
  
  
  CI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),1]
                                 ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,outx = F,delta.pred = 0,delta.obs = 0)[c("cindex","p.value")]
  
  mCI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),1]
                                  ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),2],delta.pred = 0,delta.obs = 0.2
                                  ,outx = F)[c("cindex","p.value")]
  
  Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "greater")$p.value
  Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "less")$p.value
  
  
  Pvalue_annot <- NULL
  
  if( Pvalue_greater < Pvalue_less){
    Pvalue_annot <- "greater"
  }else{
    Pvalue_annot <- "less"
  }
  
  data = Pred_Output
  balanced_CI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
  balanced_mCI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot
                                         ,delta.pred = 0,delta.obs = 0.2,outx=F)
  
  output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=CI[["cindex"]],'mCI'=mCI[["cindex"]],'Pvalue.CI'=CI[["p.value"]],'Pvalue.mCI'=mCI[["p.value"]],'Pvalue'=ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less),'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                                       'CI.balanced'=balanced_CI$Cindex,'Pvalue.CI.balanced'=balanced_CI$Pvalue,'Pvalue.balanced'=balanced_CI$Pvalue,'Equal.N'=balanced_CI$N,
                                       'mCI.balanced'=balanced_mCI$Cindex,'Pvalue.mCI.balanced'=balanced_mCI$PvalueCI))
  # write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
  #     par(font=2)
  #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
  #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  
  
  return(list(output_Final,formulaOut,formulaOut_Final_all))
}


crossValidate_Lobico <- function(AUCs_Drugs, expressionData, drug, k_CV=10, newSet) {
  
  AUCs_Drug <- AUCs_Drugs[,1]
  names(AUCs_Drug) <- rownames(AUCs_Drugs)
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  cellLines_drugTested <- names(AUCs_Drug)
  cellLines_RNA <- rownames(expressionData)
  
  commonCellLines <- intersect(cellLines_drugTested,cellLines_RNA)
  print(length(expressionData))
  expressionData_sub <- expressionData[commonCellLines,newSet,drop=F]
  AUCs_Drug <- AUCs_Drug[commonCellLines]
  
  
  Data_LOBICO <- data.frame('Cell Lines'=commonCellLines ,'drug_AUC'=AUCs_Drug[commonCellLines],expressionData_sub[commonCellLines,],stringsAsFactors = F)
  formulaOut <- list()
  output_Final <- NULL
  
  # if(sum(as.numeric(Data_LOBICO[,"drug_AUC"])>0.2)==0){
  #   output_Final <- rbind(output_Final,c('K'=0,'M'=0,'Cindex'=0,'Pvalue'="greater",'N1'=length(Data_LOBICO[,"drug_AUC"]), 'N2'=0,'Pvalue_annot' = 1,
  #                                        'Cindex.balanced'=0,'Pvalue.balanced'=1,'Equal.N'=0, 'time.taken'=0))
  #
  #   formulaOut[1] <- "no senstive cell lines"
  #   return(list(output_Final,formulaOut))
  # }
  
  output_Final <- NULL
  set.seed(12435)
  flds <- createFolds(Data_LOBICO[,"drug_AUC"], k = k_CV, list = F, returnTrain = F)
  
  #  sequenceOfKM <-matrix(data = c(1,1,1,2,1,3,1,4,2,2,2,1,3,1,4,1,2,3,3,2,3,3,2,4,4,2,3,4,4,4,5,1,5,2,1,5,2,5),nrow = 19, ncol = 2,byrow = T)
  sequenceOfKM <-matrix(data = c(1,1,1,2,1,3,1,4,2,2,2,1,3,1,4,1),nrow = 8, ncol = 2,byrow = T)#,3,3,2,4,4,2,3,4,4,4,5,1,5,2,1,5,2,5),nrow = 19, ncol = 2,byrow = T)
#    sequenceOfKM <-matrix(data = c(1,1),nrow = 1, ncol = 2,byrow = T)
  
  formulaOut <- list()
  
  # sequenceOfKM <-matrix(data = c(1,1,1,2),nrow = 2, ncol = 2,byrow = T)
  #  formulaOut <- list()
  for (j in 1:dim(sequenceOfKM)[1]){
    
    
    start.time5 <- proc.time()
    
    K <- sequenceOfKM[j,1]
    M <- sequenceOfKM[j,2]
    
    
    Sol_out <- list()
    
    #trainData <- Data_LOBICO
    #testData <- Data_LOBICO
    
    Pred_Output <- matrix(data = NA,nrow = 1,ncol = 2)
    
    
    formulaOut[[j]] <- ""
    for (r in 1:k_CV){
      
      testInd <- which(flds==r)
      trainInd <- which(flds!=r)
      trainData <- Data_LOBICO[trainInd,]
      testData <- Data_LOBICO[testInd,]
      
      
      #source("~/Downloads/R-LOBICO/runLobicoRWrapper.R")
      
      SolMat <- runLobicoRWrapper(drugPheno = as.numeric(trainData[,"drug_AUC"]),mat = (trainData[,3:dim(trainData)[2]]),K = K,M = M)
      SolMat <- round(SolMat)
      colnames(SolMat) <- newSet
      if(all(SolMat==0)){
        formulaOut[j] <- "no solution for training"
        output <- cbind(rep("0",length(as.numeric(testData[,"drug_AUC"]))),testData[,"drug_AUC"])
        rownames(output) <- rownames(testData)
        Pred_Output <- rbind(Pred_Output,output)
        next
      }
      
      formulaOut[[j]] <- c(formulaOut[[j]],showformula(SolMat,K,M,newSet)[[1]])
      
      print("&&&&&&&&&&&&&")
      print(showformula(SolMat,K,M,newSet))
      print(j)
      print(r)
      print("&&&&&&&&&&&&&")
      output <- getPredictions_LOBICO(SolMat,testData,K,M)
      Pred_Output <- rbind(Pred_Output,output)
      
      
    }
    
    
    class(Pred_Output) <- "numeric"
    Pred_Output <- Pred_Output[c(-1),]
    
    
    
    
    expressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==1]
    notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==0]
    
    if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
      end.time5 <- proc.time()
      time.taken5 <- end.time5 - start.time5
      output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=0.5,'mCI'=0.5,'Pvalue.CI'=1,'Pvalue.mCI'=1,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = "",
                                           'CI.balanced'=0.5,'Pvalue.CI.balanced'=1,'Pvalue.balanced'=1,'Equal.N'=0,
                                           'mCI.balanced'=0.5,'Pvalue.mCI.balanced'=1))
      next
    }
    
    #integrCindex <- Hmisc::rcorr.cens(S=Pred_Output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
    
    CI <- paired.concordance.index(predictions = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),1])
                                   ,observations = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),2])
                                   ,outx = F,delta.pred = 0,delta.obs = 0)[c("cindex","p.value")]
    
    mCI <- paired.concordance.index(predictions = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),1])
                                    ,observations = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),2]),delta.pred = 0,delta.obs = 0.2
                                    ,outx = F)[c("cindex","p.value")]
    
    Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "greater")$p.value
    Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "less")$p.value
    
    
    Pvalue_annot <- NULL
    
    if( Pvalue_greater < Pvalue_less){
      Pvalue_annot <- "greater"
    }else{
      Pvalue_annot <- "less"
    }
    
    data = Pred_Output
    balanced_CI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
    balanced_mCI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot
                                           ,delta.pred = 0,delta.obs = 0.2,outx=F)
    end.time5 <- proc.time()
    time.taken5 <- (end.time5 - start.time5)["elapsed"]
    
    output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=CI[["cindex"]],'mCI'=mCI[["cindex"]],'Pvalue.CI'=CI[["p.value"]],'Pvalue.mCI'=mCI[["p.value"]],'Pvalue'=ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less),'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                                         'CI.balanced'=balanced_CI$Cindex,'Pvalue.CI.balanced'=balanced_CI$Pvalue,'Pvalue.balanced'=balanced_CI$Pvalue,'Equal.N'=balanced_CI$N,
                                         'mCI.balanced'=balanced_mCI$Cindex,'Pvalue.mCI.balanced'=balanced_mCI$PvalueCI))
    #    write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
    #     par(font=2)
    #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
    #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
    
    
    
  }
  
  
  return(list(output_Final,formulaOut))
}


crossValidate_Lobico_Outer <- function(AUCs_Drugs, expressionData, drug, k_CV=10, newSet,numOfSolutions=1) {
  
  
  
  
  AUCs_Drug <- AUCs_Drugs[drug,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  cellLines_drugTested <- names(AUCs_Drug)
  cellLines_RNA <- rownames(expressionData)
  
  commonCellLines <- intersect(cellLines_drugTested,cellLines_RNA)
  
  expressionData_sub <- expressionData[commonCellLines,newSet]
  AUCs_Drug <- AUCs_Drug[commonCellLines]
  
  
  Data_LOBICO <- data.frame('Cell Lines'=commonCellLines ,'drug_AUC'=AUCs_Drug[commonCellLines],expressionData_sub[commonCellLines,],stringsAsFactors = F)
  # Data_LOBICO <- as.data.frame(Data_LOBICO,stringsAsFactors=F)
  #  lapply(2:dim(Data_LOBICO)[2],function(r,Data_LOBICO){
  #    class(Data_LOBICO[,r]) <<- "numeric"
  #  },Data_LOBICO=Data_LOBICO)
  formulaOut <- list()
  output_Final <- NULL
  
  output_Final <- NULL
  set.seed(1245)
  flds <- createFolds(Data_LOBICO[,"drug_AUC"], k = k_CV, list = F, returnTrain = F)
  
  formulaOut <- list()
  
  
  Sol_out <- list()
  
  Pred_Output <- matrix(data = NA,nrow = 1,ncol = 2)
  start.time = proc.time()
  for (i in 1:k_CV){
    print(paste(">>>>>>>>>>>>>>>>>>",i,"<<<<<<<<<<<<<<<<"))
    
    formulaOut[[i]] <- ""
    testInd <- which(flds==i)
    trainInd <- which(flds!=i)
    trainData <- Data_LOBICO[trainInd,2:dim(Data_LOBICO)[2]]
    testData <- Data_LOBICO[testInd,2:dim(Data_LOBICO)[2]]
    
    
    dd <- mRMR.data(data = trainData)
    test <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 3, feature_count = 10)
    
    mrmrOut <- lapply(1:numOfSolutions, function(p){
      formulaOut1 <- ""
      newSetTest <- c(apply(solutions(test)[[1]][,p,drop=F], 2, function(t, y) { return(y[t]) }, y=featureNames(dd)))
      namesIbx <- match(newSetTest,newSet)
      names(newSetTest) <- names(newSet)[namesIbx]
      
      ibx <- which(newSetTest == colnames(trainData)[1])
      if(length(ibx)!=0){
        newSetTest <- newSetTest[-ibx]
      }
      print(">>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<")
      print(i)
      print(">>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<")
      A <- crossValidate_Lobico(AUCs_Drugs = trainData[,1,drop=F],expressionData = trainData[,2:dim(trainData)[2]],drug = drug,newSet = newSetTest)
      print(paste(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",i,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"))
      train_CV_output_drug <- A[[1]]
      K <- NULL
      M <- NULL
      Pvalue_annot <- NULL
      
      ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
      if(length(ibx)==1){
        K = train_CV_output_drug[ibx,"K"]
        M = train_CV_output_drug[ibx,"M"]
        Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
      }else{
        K = train_CV_output_drug[ibx[1],"K"]
        M = train_CV_output_drug[ibx[1],"M"]
        Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
      }
      
      Data_LOBICO_training <- data.frame('Cell Lines'=rownames(trainData) ,'drug_AUC'=trainData[,1],trainData[,newSetTest],stringsAsFactors = F)
      
      #   source("~/Downloads/R-LOBICO/lobico.R")
      
      SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
      
      
      
      
      #    write.table(Data_LOBICO_training, file='test5.tsv', quote=F, sep='\t',row.names = FALSE)
      #    evaluate(matlab, paste("[SolMat]=runLobicoWrapper('test5.tsv',",K,",",M,");",sep = ""))
      #     SolMat <- getVariable(matlab, c("SolMat"))$SolMat
      
      colnames(SolMat) <- newSetTest
      
      if(all(SolMat==0)){
        formulaOut1 <- ""
      }else{
        formulaOut1 <- showformula(SolMat,K,M,newSetTest)
      }
      
      
      
      
      return(list(A,formulaOut1,newSetTest))
    })
    
    
    
    #finalGeneSet <- unique(unlist(lapply(mrmrOut, "[[",3)))
    
    allGenes <- unlist(lapply(mrmrOut, "[[",2))
    allGenes <- gsub("[~&|]","",x = allGenes)
    allGenes <- unlist(strsplit(allGenes," "))
    finalGeneSet <- unique(allGenes)
    finalGeneSet <- finalGeneSet[finalGeneSet!=""]
    ibx <- match(finalGeneSet,newSet)
    names(finalGeneSet) <- names(newSet)[ibx]
    
    
    B <- crossValidate_Lobico(AUCs_Drugs = trainData[,1,drop=F],expressionData = trainData[,finalGeneSet],drug = drug,newSet = finalGeneSet)
    
    
    train_CV_output_drug <- B[[1]]
    # train_CV_output_drug <- mrmrOut[[1]][[1]][[1]]
    K <- NULL
    M <- NULL
    Pvalue_annot <- NULL
    
    ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
    if(length(ibx)==1){
      K = train_CV_output_drug[ibx,"K"]
      M = train_CV_output_drug[ibx,"M"]
      Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
    }else{
      K = train_CV_output_drug[ibx[1],"K"]
      M = train_CV_output_drug[ibx[1],"M"]
      Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
    }
    
    Data_LOBICO_training <- cbind('Cell Lines'=rownames(trainData) ,'drug_AUC'=trainData[,1],trainData[,finalGeneSet])
    
    SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,3:dim(Data_LOBICO_training)[2]]),K = K,M = M)
    
    
    
    #    write.table(Data_LOBICO_training, file='test6.tsv', quote=F, sep='\t',row.names = FALSE)
    #   print(drug)
    
    #    evaluate(matlab, paste("[SolMat,labels]=runLobicoWrapper('test6.tsv',",K,",",M,");",sep = ""))
    print(drug)
    #    SolMat <- getVariable(matlab, c("SolMat"))$SolMat
    colnames(SolMat) <- finalGeneSet
    if(all(SolMat==0)){
      formulaOut[[i]] <- "no solution for training"
      output <- cbind(rep("0",length(as.numeric(testData[,"drug_AUC"]))),testData[,"drug_AUC"])
      Sol_out[[i]] <- c(0.5,K,M)
      rownames(output) <- rownames(testData)
      Pred_Output <- rbind(Pred_Output,output)
      next
    }
    
    formulaOut[[i]] <- c(formulaOut[[i]],showformula(SolMat,K,M,finalGeneSet)[[1]])
    #      write.table(showFormula(SolMat,K,M,newSet),file = "fluvastatin.CV.txt",quote = F,append = T)
    output <- getPredictions_LOBICO_outer(SolMat,testData[,c("drug_AUC",finalGeneSet)],K,M)
    
    
    expressedCellLines <- rownames(output)[output[,1]==1]
    notExpressedCellLines <- rownames(output)[output[,1]==0]
    
    #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
    integrCindex <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),2],output[c(notExpressedCellLines,expressedCellLines),1]
                                             ,delta.pred = 0,delta.obs = 0,outx = F)
    
    
    Sol_out[[i]] <- c(integrCindex["cindex"],K,M)
    Pred_Output <- rbind(Pred_Output,output)
    
    
  }
  
  endtime <- proc.time() - start.time
  
  print(endtime)
  
  
  formulaOut_Final <- gsub("\\||&|!|~","",unlist(formulaOut))
  
  finalGeneSet_all <- unlist(lapply(formulaOut_Final, function(q){
    unlist(strsplit(q,split = " "))
  }))
  finalGeneSet_all <- finalGeneSet_all[finalGeneSet_all!=""]
  finalGeneSet_all <- unique(finalGeneSet_all)
  
  ibx <- match(finalGeneSet_all,newSet)
  names(finalGeneSet_all) <- names(newSet)[ibx]
  Data_LOBICO_training <- cbind('drug_AUC'=Data_LOBICO[,2],Data_LOBICO[,(finalGeneSet_all)])
  
  dd <- mRMR.data(data = Data_LOBICO_training)
  test <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 3, feature_count = 10)
  
  mrmrOut1 <- lapply(1:3, function(p){
    formulaOut1 <- ""
    newSetTest <- c(apply(solutions(test)[[1]][,p,drop=F], 2, function(t, y) { return(y[t]) }, y=featureNames(dd)))
    ibx <- match(newSetTest,newSet)
    names(newSetTest) <- names(newSet)[ibx]
    
    ibx <- which(newSetTest == colnames(Data_LOBICO_training)[1])
    if(length(ibx)!=0){
      newSetTest <- newSetTest[-ibx]
    }
    A <- crossValidate_Lobico(AUCs_Drugs = Data_LOBICO_training[,1,drop=F],expressionData = Data_LOBICO_training[,2:dim(Data_LOBICO_training)[2]],drug = drug,newSet = newSetTest)
    print(paste(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",p,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"))
    train_CV_output_drug <- A[[1]]
    K <- NULL
    M <- NULL
    Pvalue_annot <- NULL
    
    ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
    if(length(ibx)==1){
      K = train_CV_output_drug[ibx,"K"]
      M = train_CV_output_drug[ibx,"M"]
      Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
    }else{
      K = train_CV_output_drug[ibx[1],"K"]
      M = train_CV_output_drug[ibx[1],"M"]
      Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
    }
    
    
    SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,newSetTest]),K = K,M = M)
    
    
    
    
    
    colnames(SolMat) <- newSetTest
    
    if(all(SolMat==0)){
      formulaOut1 <- ""
    }else{
      formulaOut1 <- showformula(SolMat,K,M,newSetTest)
    }
    
    
    
    
    return(list(A,formulaOut1,newSetTest))
  })
  
  
  
  #finalGeneSet <- unique(unlist(lapply(mrmrOut, "[[",3)))
  
  
  allGenes <- unlist(lapply(mrmrOut1, "[[",2))
  allGenes <- gsub("[~&|]","",x = allGenes)
  allGenes <- unlist(strsplit(allGenes," "))
  finalGeneSet1 <- unique(allGenes)
  finalGeneSet1 <- finalGeneSet1[finalGeneSet1!=""]
  ibx <- match(finalGeneSet1,newSet)
  names(finalGeneSet1) <- names(newSet)[ibx]
  
  
  # finalGeneSet1 <- unique(unlist(lapply(mrmrOut1, "[[",2)))
  # finalGeneSet1 <- finalGeneSet1[finalGeneSet1!=""]
  #  ibx <- match(finalGeneSet1,newSet)
  # names(finalGeneSet1) <- names(newSet)[ibx]
  
  
  B <- crossValidate_Lobico(AUCs_Drugs = Data_LOBICO_training[,1,drop=F],expressionData = Data_LOBICO_training[,2:dim(Data_LOBICO_training)[2]],drug = drug,newSet = finalGeneSet1)
  
  
  train_CV_output_drug <- B[[1]]
  # train_CV_output_drug <- mrmrOut[[1]][[1]][[1]]
  K <- NULL
  M <- NULL
  Pvalue_annot <- NULL
  
  ibx <- which(train_CV_output_drug[,"mCI.balanced"]==max(train_CV_output_drug[,"mCI.balanced"],na.rm = T))
  if(length(ibx)==1){
    K = train_CV_output_drug[ibx,"K"]
    M = train_CV_output_drug[ibx,"M"]
    Pvalue_annot <- train_CV_output_drug[ibx,"Pvalue_annot"]
  }else{
    K = train_CV_output_drug[ibx[1],"K"]
    M = train_CV_output_drug[ibx[1],"M"]
    Pvalue_annot <- train_CV_output_drug[ibx[1],"Pvalue_annot"]
  }
  
  
  SolMat <- runLobicoRWrapper(drugPheno = as.numeric(Data_LOBICO_training[,"drug_AUC"]),mat = (Data_LOBICO_training[,finalGeneSet1]),K = K,M = M)
  
  colnames(SolMat) <- finalGeneSet1
  
  
  formulaOut_Final_all <- showformula(SolMat,K,M,finalGeneSet1)
  #      write.table(showFormula(SolMat,K,M,newSet),file = "fluvastatin.CV.txt",quote = F,append = T)
  #  output <- getPredictions_LOBICO_outer(SolMat,Data_LOBICO_training[,c("drug_AUC",finalGeneSet1)],K,M)
  
  
  
  #  Pred_Output <- output
  
  
  
  class(Pred_Output) <- "numeric"
  
  Pred_Output <- Pred_Output[-1,]
  
  expressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==1]
  notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==0]
  
  if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
    
    output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=0.5,'mCI'=0.5,'Pvalue.CI'=1,'Pvalue.mCI'=1,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = "",
                                         'CI.balanced'=0.5,'Pvalue.CI.balanced'=1,'Pvalue.balanced'=1,'Equal.N'=0,
                                         'mCI.balanced'=0.5,'Pvalue.mCI.balanced'=1))
    next
  }
  
  
  CI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),1]
                                 ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,outx = F,delta.pred = 0,delta.obs = 0)[c("cindex","p.value")]
  
  mCI <- paired.concordance.index(predictions = Pred_Output[c(notExpressedCellLines,expressedCellLines),1]
                                  ,observations = Pred_Output[c(notExpressedCellLines,expressedCellLines),2],delta.pred = 0,delta.obs = 0.2
                                  ,outx = F)[c("cindex","p.value")]
  
  Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "greater")$p.value
  Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "less")$p.value
  
  
  Pvalue_annot <- NULL
  
  if( Pvalue_greater < Pvalue_less){
    Pvalue_annot <- "greater"
  }else{
    Pvalue_annot <- "less"
  }
  
  data = Pred_Output
  balanced_CI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
  balanced_mCI <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot
                                         ,delta.pred = 0,delta.obs = 0.2,outx=F)
  
  output_Final <- rbind(output_Final,c('K'=K,'M'=M,'CI'=CI[["cindex"]],'mCI'=mCI[["cindex"]],'Pvalue.CI'=CI[["p.value"]],'Pvalue.mCI'=mCI[["p.value"]],'Pvalue'=ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less),'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                                       'CI.balanced'=balanced_CI$Cindex,'Pvalue.CI.balanced'=balanced_CI$Pvalue,'Pvalue.balanced'=balanced_CI$Pvalue,'Equal.N'=balanced_CI$N,
                                       'mCI.balanced'=balanced_mCI$Cindex,'Pvalue.mCI.balanced'=balanced_mCI$PvalueCI))
  # write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
  #     par(font=2)
  #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
  #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  
  
  return(list(output_Final,formulaOut,formulaOut_Final_all))
}


applyLBIOCOrules <- function(AUCs_Drugs_training, expressionData_training, drug, LOBCIO_CV_Out,newSet_training,newSet_final=NULL,AUCs_Drugs_testing, expressionData_testing
                             ,newSet_testing,validate=T) {
  
  AUCs_Drug_training <- AUCs_Drugs_training[drug,]
  AUCs_Drug_training <- AUCs_Drug_training[!is.na(AUCs_Drug_training)]
  cellLines_drugTested_training <- names(AUCs_Drug_training)
  cellLines_RNA_training <- rownames(expressionData_training)
  
  LOBCIO_CV_Out_sub <- LOBCIO_CV_Out[[drug]]
  
  K <- LOBCIO_CV_Out_sub[[1]][1]
  M <- LOBCIO_CV_Out_sub[[1]][2]
  
  newSet_final1 <- NULL
  
  if(!is.null(newSet_final)){
    K="1"
    M="1"
    newSet_final1 <- newSet_final
  }else{
    newSet_final1 <- unlist(strsplit(gsub(pattern = "[!|&]",replacement = "", x = LOBCIO_CV_Out_sub[[3]]),split = " "))
    newSet_final1 <- newSet_final1[which(newSet_final1!="")]
  }
  finalGeneSet <- newSet_training[newSet_final1]
  
  commonCellLines_training <- intersect(cellLines_drugTested_training,cellLines_RNA_training)
  
  expressionData_sub_training <- expressionData_training[commonCellLines_training,finalGeneSet,drop=F]
  AUCs_Drug_training <- AUCs_Drug_training[commonCellLines_training]
  
  
  Data_LOBICO <- cbind('Cell Lines'=commonCellLines_training ,'drug_AUC'=AUCs_Drug_training[commonCellLines_training],expressionData_sub_training[commonCellLines_training,,drop=F])
  Data_LOBICO <- as.data.frame(Data_LOBICO,stringsAsFactors=F)
  lapply(2:dim(Data_LOBICO)[2],function(r,Data_LOBICO){
    class(Data_LOBICO[,r]) <<- "numeric"
  },Data_LOBICO=Data_LOBICO)
  
  start.time5 <- proc.time()
  
  
  write.table(Data_LOBICO, file='test6.tsv', quote=F, sep='\t',row.names = FALSE)
  print(drug)
  evaluate(matlab, paste("[SolMat,labels]=runLobicoWrapper('test6.tsv',",K,",",M,");",sep = ""))
  print(drug)
  SolMat <- getVariable(matlab, c("SolMat"))$SolMat
  colnames(SolMat) <- finalGeneSet
  if(all(SolMat==0)){
    output <- cbind(rep("0",length(as.numeric(Data_LOBICO[,"drug_AUC"]))),Data_LOBICO[,"drug_AUC"])
    rownames(output) <- rownames(Data_LOBICO)
    Pred_Output <- rbind(Pred_Output,output)
    
  }
  
  
  
  output <- getPredictions_LOBICO_outer(SolMat,Data_LOBICO[,c("drug_AUC",finalGeneSet)],K,M)
  Pred_Output <- output
  
  
  
  
  
  class(Pred_Output) <- "numeric"
  
  
  
  expressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==1]
  notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==0]
  
  if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
    end.time5 <- proc.time()
    time.taken5 <- end.time5 - start.time5
    output_Final <- rbind(output_Final,c('K'=K,'M'=M,'Cindex'=0,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),
                                         'Cindex.balanced'=0,'Pvalue.balanced'=1,'Equal.N'=0, 'time.taken'=time.taken5))
    
    return(c("CI.training"=0.5,"Pvalue.training"=1
             ,"CI.testing"=0.5,"Pvalue.testing"=1))
  }
  
  integrCindex <- Hmisc::rcorr.cens(S=Pred_Output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  
  C <- integrCindex["C Index"]
  names(C) <- NULL
  Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "greater")$p.value
  Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "less")$p.value
  
  
  Pvalue_annot <- NULL
  
  if( Pvalue_greater < Pvalue_less){
    Pvalue_annot <- "greater"
    Pvalue <- Pvalue_greater
  }else{
    Pvalue_annot <- "less"
    Pvalue <- Pvalue_less
  }
  
  data = Pred_Output
  balanced <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
  end.time5 <- proc.time()
  time.taken5 <- 0
  
  output_Final <- c('K'=K,'M'=M,'Cindex'=C,'Pvalue'=ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less),'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                    'Cindex.balanced'=balanced$Cindex,'Pvalue.balanced'=balanced$Pvalue,'Equal.N'=balanced$N, 'time.taken'=time.taken5)
  # write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
  #     par(font=2)
  #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
  #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  
  boxplot(list("FALSE"=Pred_Output[notExpressedCellLines,2], "TRUE"=Pred_Output[expressedCellLines,2]),main=paste(drug,"\n",ifelse(is.null(newSet_final),LOBCIO_CV_Out_sub[[3]],names(finalGeneSet))),col=c("indianred2","lightskyblue"))
  #  legend("topright",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  C <- as.numeric(LOBCIO_CV_Out_sub[[1]][3])
  Pvalue <- as.numeric(LOBCIO_CV_Out_sub[[1]][4])
  N1 <- as.numeric(LOBCIO_CV_Out_sub[[1]][5])
  N2 <- as.numeric(LOBCIO_CV_Out_sub[[1]][6])
  Cbalanced <- as.numeric(LOBCIO_CV_Out_sub[[1]][8])
  PvalBalanced <- as.numeric(LOBCIO_CV_Out_sub[[1]][9])
  Ne <- as.numeric(LOBCIO_CV_Out_sub[[1]][10])
  legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",N1),paste("N1:",N2),paste("CI balanced:", sprintf("%.2g",Cbalanced)),paste("P-value balanced:", sprintf("%.1E",PvalBalanced)),paste("Ne:",Ne) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  
  
  
  
  AUCs_Drug_testing <- AUCs_Drugs_testing[drug,]
  AUCs_Drug_testing <- AUCs_Drug_testing[!is.na(AUCs_Drug_testing)]
  cellLines_drugTested_testing <- names(AUCs_Drug_testing)
  cellLines_RNA_testing <- rownames(expressionData_testing)
  commonCellLines_testing <- intersect(cellLines_RNA_testing,cellLines_drugTested_testing)
  names(newSet_testing)[grep("ENSG00000090530",newSet_testing)] <- "LEPREL1"
  names(newSet_testing)[grep("ENSG00000187164",newSet_testing)] <- "KIAA1598"
  finalGeneSet_testing <- newSet_testing[names(finalGeneSet)]
  print(finalGeneSet_testing)
  expressionData_sub_testing <- expressionData_testing[commonCellLines_testing,finalGeneSet_testing,drop=F]
  AUCs_Drug_testing <- AUCs_Drug_testing[commonCellLines_testing]
  
  #commonCellLines_testing <- intersect(commonCellLines_testing,commonCellLines_training)
  
  Data_LOBICO <- cbind('Cell Lines'=commonCellLines_testing ,'drug_AUC'=AUCs_Drug_testing[commonCellLines_testing],expressionData_sub_testing[commonCellLines_testing,finalGeneSet_testing,drop=F])
  Data_LOBICO <- as.data.frame(Data_LOBICO,stringsAsFactors=F)
  lapply(2:dim(Data_LOBICO)[2],function(r,Data_LOBICO){
    class(Data_LOBICO[,r]) <<- "numeric"
  },Data_LOBICO=Data_LOBICO)
  
  
  colnames(SolMat) <- finalGeneSet_testing
  
  start.time5 <- proc.time()
  output <- getPredictions_LOBICO_outer(SolMat,Data_LOBICO[,c("drug_AUC",finalGeneSet_testing)],K,M)
  Pred_Output <- output
  
  
  
  
  
  class(Pred_Output) <- "numeric"
  
  
  
  expressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==1]
  notExpressedCellLines <- rownames(Pred_Output)[Pred_Output[,1]==0]
  
  if( length(expressedCellLines) == 0 || length(notExpressedCellLines) == 0 ){
    end.time5 <- proc.time()
    time.taken5 <- end.time5 - start.time5
    output_Final <- rbind(output_Final,c('K'=K,'M'=M,'Cindex'=0,'Pvalue'=1,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),
                                         'Cindex.balanced'=0,'Pvalue.balanced'=1,'Equal.N'=0, 'time.taken'=time.taken5))
    
    return(c("CI.training"=Cbalanced,"Pvalue.training"=PvalBalanced
             ,"CI.testing"=0.5,"Pvalue.testing"=1))
  }
  
  integrCindex <- Hmisc::rcorr.cens(S=Pred_Output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(Pred_Output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  
  C <- integrCindex["C Index"]
  names(C) <- NULL
  Pvalue_greater <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "greater")$p.value
  Pvalue_less <- wilcox.test(Pred_Output[notExpressedCellLines,2], Pred_Output[expressedCellLines,2],alternative = "less")$p.value
  
  
  Pvalue_annot <- NULL
  
  if( Pvalue_greater < Pvalue_less){
    Pvalue_annot <- "greater"
    Pvalue <- Pvalue_greater
  }else{
    Pvalue_annot <- "less"
    Pvalue <- Pvalue_less
  }
  
  data = Pred_Output
  balanced <- resample_Cindex_PValue(Pred_Output[,1],Pred_Output[,2],permutation = 100,Pvalue_annot = Pvalue_annot,outx=F)
  end.time5 <- proc.time()
  time.taken5 <- 0
  
  Pvalue <- ifelse(Pvalue_annot=="greater",Pvalue_greater,Pvalue_less)
  output_Final <- c('K'=K,'M'=M,'Cindex'=C,'Pvalue'=Pvalue,'N1'=length(notExpressedCellLines), 'N2'=length(expressedCellLines),'Pvalue_annot' = Pvalue_annot,
                    'Cindex.balanced'=balanced$Cindex,'Pvalue.balanced'=balanced$Pvalue,'Equal.N'=balanced$N, 'time.taken'=time.taken5)
  # write.table(output_Final,file = "fluvastatin.CV.txt",quote = F,append = T)
  #     par(font=2)
  #    boxplot(cbind("Not Enriched"=Pred_Output[notExpressedCellLines,2],"Enriched"=Pred_Output[expressedCellLines,2]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
  #    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
  
  
  if(validate==T){
    boxplot(list("FALSE"=Pred_Output[notExpressedCellLines,2], "TRUE"=Pred_Output[expressedCellLines,2]),main=paste(drug,"\n",ifelse(is.null(newSet_final),LOBCIO_CV_Out_sub[[3]],names(finalGeneSet))),col=c("indianred2","lightskyblue"))
    #  legend("topright",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines)),paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex)),paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",length(notExpressedCellLines)) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
    
    
    
    legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("N0:",length(notExpressedCellLines))
                                ,paste("N1:",length(expressedCellLines)),paste("CI balanced:", sprintf("%.2g",balanced$Cindex))
                                ,paste("P-value balanced:", sprintf("%.1E",balanced$Pvalue)),paste("Ne:",balanced$N) ),bty="n",pch = c(NA,NA,15,15,NA,NA,NA),col=c("indianred2","lightskyblue") )
    
  }
  
  
  
  
  
  
  
  
  
  Cbalanced <- as.numeric(LOBCIO_CV_Out_sub[[1]][8])
  PvalBalanced <- as.numeric(LOBCIO_CV_Out_sub[[1]][9])
  
  return(c("CI.training"=Cbalanced,"Pvalue.training"=PvalBalanced
           ,"CI.testing"=balanced$Cindex,"Pvalue.testing"=balanced$Pvalue))
}





