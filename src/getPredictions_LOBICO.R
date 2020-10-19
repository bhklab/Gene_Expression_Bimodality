getPredictions_LOBICO <- function(SolMat,testData,K,M){
  
  

  
  if (K<=M){
    # rows of SolMat have an "|" between them and "&" between columns
    label_or <- list()
    for(i in 1:dim(SolMat)[1]){
      
      ibx <- which(SolMat[i,]!=0)
      
      L <- apply(testData[,c(3:(dim(SolMat)[2]+2))],1,function(x,ibx){
        labels_and <- list()
        for(j in 1:length(ibx)){
          if(round(SolMat[i,ibx[j]])==1){
            labels_and[j] <- x[ibx[j]]
          }else if(round(SolMat[i,ibx[j]])==-1){
            labels_and[j] <- as.numeric(!as.numeric(x[ibx[j]]))
          }
        }
        
        return(as.numeric(all(as.numeric(labels_and)==1)))
      },ibx=ibx)
      
      label_or[[i]] <- L
      
    }
    label_or <- matrix(unlist(label_or), ncol = dim(testData)[1], byrow = TRUE)
    colnames(label_or) <- rownames(testData)
    
    R <- apply(label_or,2,function(x){
      return(as.numeric(any(x==1)))
    })
    
    return(cbind(R,testData[,"drug_AUC"]))
    
  }else{
    # rows of SolMat have an "&" between them and "|" between columns
    
    label_or <- list()
    for(i in 1:dim(SolMat)[1]){
      
      ibx <- which(SolMat[i,]!=0)
      
      L <- apply(testData[,c(3:(dim(SolMat)[2]+2))],1,function(x,ibx){
        labels_and <- list()
        for(j in 1:length(ibx)){
          if(round(SolMat[i,ibx[j]])==1){
            labels_and[j] <- x[ibx[j]]
          }else if(round(SolMat[i,ibx[j]])==-1){
            labels_and[j] <- as.numeric(!as.numeric(x[ibx[j]]))
          }
        }
        
        return(as.numeric(any(as.numeric(labels_and)==1)))
      },ibx=ibx)
      
      label_or[[i]] <- L
      
    }
    label_or <- matrix(unlist(label_or), ncol = dim(testData)[1], byrow = TRUE)
    colnames(label_or) <- rownames(testData)
    
    R <- apply(label_or,2,function(x){
      return(as.numeric(all(x==1)))
    })
    
    return(cbind(R,testData[,"drug_AUC"]))
    
    
    
  }
  
  
}

getPredictions_LOBICO_outer <- function(SolMat,testData,K,M){
  
  
  
  
  if (K<=M){
    # rows of SolMat have an "|" between them and "&" between columns
    label_or <- list()
    for(i in 1:dim(SolMat)[1]){
      
      ibx <- which(SolMat[i,]!=0)
      
      L <- apply(testData[,c(2:(dim(SolMat)[2]+1)),drop=F],1,function(x,ibx){
        labels_and <- list()
        for(j in 1:length(ibx)){
          if(round(SolMat[i,ibx[j]])==1){
            labels_and[j] <- x[ibx[j]]
          }else if(round(SolMat[i,ibx[j]])==-1){
            labels_and[j] <- as.numeric(!as.numeric(x[ibx[j]]))
          }
        }
        
        return(as.numeric(all(as.numeric(labels_and)==1)))
      },ibx=ibx)
      
      label_or[[i]] <- L
      
    }
    label_or <- matrix(unlist(label_or), ncol = dim(testData)[1], byrow = TRUE)
    colnames(label_or) <- rownames(testData)
    
    R <- apply(label_or,2,function(x){
      return(as.numeric(any(x==1)))
    })
    
    return(cbind(R,testData[,"drug_AUC"]))
    
  }else{
    # rows of SolMat have an "&" between them and "|" between columns
    
    label_or <- list()
    for(i in 1:dim(SolMat)[1]){
      
      ibx <- which(SolMat[i,]!=0)
      
      L <- apply(testData[,c(2:(dim(SolMat)[2]+1)),drop=F],1,function(x,ibx){
        labels_and <- list()
        for(j in 1:length(ibx)){
          if(round(SolMat[i,ibx[j]])==1){
            labels_and[j] <- x[ibx[j]]
          }else if(round(SolMat[i,ibx[j]])==-1){
            labels_and[j] <- as.numeric(!as.numeric(x[ibx[j]]))
          }
        }
        
        return(as.numeric(any(as.numeric(labels_and)==1)))
      },ibx=ibx)
      
      label_or[[i]] <- L
      
    }
    label_or <- matrix(unlist(label_or), ncol = dim(testData)[1], byrow = TRUE)
    colnames(label_or) <- rownames(testData)
    
    R <- apply(label_or,2,function(x){
      return(as.numeric(all(x==1)))
    })
    
    return(cbind(R,testData[,"drug_AUC"]))
    
    
    
  }
  
  
}


getPredictions_LOBICO2 <- function(SolMat,testData,K,M){
  
  L <- apply(X = testData[,c(3:(dim(SolMat)[2]+2))],MARGIN = 1,FUN = function(x, SolMat,K,M){
    
      if(K==1 & M==1){
        str = ""
        pos = which(SolMat!=0)
        if (SolMat[1,pos]==-1){
          str = paste(str, '!',sep = "");
        }
        str = paste(str, x[pos], ' ',sep="")
        str <- eval(parse(text=str))
        return(str)
      }
      
      if (K<=M){
        
        str = ""
        for (k in 1:K){
          pos = which(SolMat[k,]!=0)
          for (p in 1:length(pos)){
            if (SolMat[k,pos[p]]==-1){
              str = paste(str, '!',sep = "");
            }
            str = paste(str, x[pos[p]], ' ',sep="")
            if (p!=length(pos)){
              str = paste(str, ' & ',sep = "")
            }else{
              if (k!=K){
                str = paste(str, '  |  ',sep = "")
              }
            }
          }
        }
        str <- paste("(",str,")")
        str <- gsub(" | ", ") | (", str, fixed = T)
        str <- eval(parse(text=str))
        return(str)
      }else{
        # rows of SolMat have an "&" between them and "|" between columns
        Ktemp = K;
        K = M;
        M = Ktemp;
        
        str = ""
        for (k in 1:K){
          pos = which(SolMat[k,]!=0)
          for (p in 1:length(pos)){
            if (SolMat[k,pos[p]]==-1){
              str = paste(str, '!',sep = "")
            }
            str = paste(str, x[pos[p]], " ",sep = "")
            if (p!=length(pos)){
              str = paste(str, ' | ',sep = "")
            }else{
              if (k!=K){
                str = paste(str, '  &  ',sep = "")
              }
            }
          }
        }
        str <- paste("(",str,")")
        str <- gsub(" & ", ") & (", str, fixed = T)
         str <- eval(parse(text=str))
        return(str)
      }
  }, SolMat=SolMat,K=K,M=M)

  return(cbind(as.numeric(L),testData[,"drug_AUC"]))
  
  
  }


