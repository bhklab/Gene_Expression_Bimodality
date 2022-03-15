resample_Cindex_PValue <- function(pred,classes,outx=F,permutation=100,Pvalue_annot,delta.pred = 0,delta.obs = 0){

  

  
  data <- cbind('pred'=as.numeric(pred),'class'=classes)
  
  
  pos <- sum(data[,"pred"]==1)
  neg <- sum(data[,"pred"]==0)
  
  draw = NULL
  if (pos >= neg){
    draw  = 1
  }else{
    draw = 0
  }
  
  set.seed(12345)
  CIs <- list()
  Pvalues <- list()
  pvalues_CIs <- list()

  
  
  
  for( i in 1:(permutation)){
    
    ibx_Samples <- sample(which(data[,"pred"]==draw),size = sum(data[,"pred"]!=draw), replace = F)  
    ibx_Samples_notDraw <- which(data[,"pred"]!=draw)
    
    ibx_Samples_ALL <- c(ibx_Samples,ibx_Samples_notDraw)
    
#    integrCindex <- Hmisc::rcorr.cens(S=data[ibx_Samples_ALL,"class"], x = data[ibx_Samples_ALL,"pred"],outx = outx)
    integrCindex <- paired.concordance.index(observations = data[ibx_Samples_ALL,"class"],predictions = data[ibx_Samples_ALL,"pred"]
                                             ,outx = outx,delta.pred = delta.pred,delta.obs = delta.obs)
    CIs[i] <- integrCindex[["cindex"]]
    pvalues_CIs[i] <- integrCindex[["p.value"]]
    Pvalue <- NULL
    if(draw == 1){
      Pvalue <- wilcox.test(data[ibx_Samples_notDraw,"class"], data[ibx_Samples,"class"],alternative = Pvalue_annot)$p.value
    }else{
      Pvalue <- wilcox.test(data[ibx_Samples,"class"], data[ibx_Samples_notDraw,"class"],alternative = Pvalue_annot)$p.value
    }

    Pvalues[i] <- Pvalue
  }
  return(list('Cindex'=mean(unlist(CIs),na.rm = T),'PvalueCI'=mean(unlist(pvalues_CIs),na.rm = T),'Pvalue'=mean(unlist(Pvalues),na.rm = T),'N'=sum(data[,"pred"]!=draw)))
  
}



resample_Cindex_PValue2 <- function(pred,classes,outx=F,permutation=100,Pvalue_annot,threshold){
  
  
  
  
  data <- cbind('pred'=as.numeric(pred),'class'=classes)
  
  
  pos <- sum(data[,"pred"]>=threshold)
  neg <- sum(data[,"pred"]<threshold)
  
  draw = NULL
  if (pos >= neg){
    draw  = 1
  }else{
    draw = 0
  }
  
  set.seed(12345)
  Cindecies <- list()
  Pvalues <- list()
  for( i in 1:(permutation)){
    
    if(draw == 1){
      ibx_Samples <- sample(which(data[,"pred"]>=threshold),size = sum(data[,"pred"]<threshold), replace = F)  
      ibx_Samples_notDraw <- which(data[,"pred"]<threshold)
      
      ibx_Samples_ALL <- c(ibx_Samples,ibx_Samples_notDraw)
      
      integrCindex <- Hmisc::rcorr.cens(S=data[ibx_Samples_ALL,"class"], x = data[ibx_Samples_ALL,"pred"],outx = outx)
      Cindecies[i] <- integrCindex["C Index"]
      Pvalue <- NULL
      if(draw == 1){
        Pvalue <- wilcox.test(data[ibx_Samples_notDraw,"class"], data[ibx_Samples,"class"],alternative = Pvalue_annot)$p.value
      }else{
        Pvalue <- wilcox.test(data[ibx_Samples,"class"], data[ibx_Samples_notDraw,"class"],alternative = Pvalue_annot)$p.value
      }
      
      Pvalues[i] <- Pvalue
    }else{
      ibx_Samples <- sample(which(data[,"pred"]<threshold),size = sum(data[,"pred"]>=threshold), replace = F)  
      ibx_Samples_notDraw <- which(data[,"pred"]>=threshold)
      
      ibx_Samples_ALL <- c(ibx_Samples,ibx_Samples_notDraw)
      
      integrCindex <- Hmisc::rcorr.cens(S=data[ibx_Samples_ALL,"class"], x = data[ibx_Samples_ALL,"pred"],outx = outx)
      Cindecies[i] <- integrCindex["C Index"]
      Pvalue <- NULL
      if(draw == 1){
        Pvalue <- wilcox.test(data[ibx_Samples_notDraw,"class"], data[ibx_Samples,"class"],alternative = Pvalue_annot)$p.value
      }else{
        Pvalue <- wilcox.test(data[ibx_Samples,"class"], data[ibx_Samples_notDraw,"class"],alternative = Pvalue_annot)$p.value
      }
      
      Pvalues[i] <- Pvalue
    }
  }
  return(list('Cindex'=mean(unlist(Cindecies),na.rm = T),'Pvalue'=mean(unlist(Pvalues),na.rm = T),'N'=sum(data[,"pred"]!=draw)))
  
}



