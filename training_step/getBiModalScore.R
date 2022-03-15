## Function to compute BiModality Index (adopted from genefu)

getBiModalScore <- function(data){
  require(mclust)
  
  data <- data[!is.na(data)]
  if(var(data)==0){
    return(0)
  }
  rr2 <- mclust::Mclust(data = data, modelNames="V", G=2,verbose = F)
  res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
  
  if(is.null(rr2[[1]])) { ## EM algorithm did not converge
    # message("EM didn't converage")
    return(0)
  }
  res[1, ] <- rr2$parameters$mean
  res[2, ] <- rr2$parameters$variance$sigmasq
  res[3, ] <- rr2$parameters$pro
  
  ## bimodality index (BI)
  smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2] + res[2, 1]) / 2)
  bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd
  
  return(bi)
  
}





getBiModalScore_Updated <- function(data,type=c("gaussianMix","gammaGaussianMix","compare"),plot=F,trunc=F,iter=100){
  
  
  require(dplyr)
  
  data <- data[!is.na(data)]
  if(var(data)==0){
    return(list("BIndex"=0,"mix"=NULL))
  }
  
  
  init1 <- c(quantile(data,probs = 0.1),quantile(data,probs = 0.9))
#  

  ratios <- c(0.5,0.5)
  
  if(type=="gaussianMix"){
    init2 <- c(abs(range(data)[1]-range(data)[2])/10,abs(range(data)[1]-range(data)[2])/10)
    mix <- tryCatch(getBestDist(data,init1,init2,ratios,type = "gaussianMix",iter = iter),error = function(e){
      print(e)
      return(NULL)
    })
  } else if( type=="gammaGaussianMix"){
    init2 <- c(abs(range(data)[1]-range(data)[2])/10,abs(range(data)[1]-range(data)[2])/10)
    mix <- tryCatch(getBestDist(data,init1,init2,ratios,type = "gammaGaussianMix",iter = iter),error = function(e){
      print(e)
      return(NULL)
    })
    print(mix)
  }else if( type=="negativeBinomialMix"){
    init2 <- c(quantile(data,probs = 0.1),quantile(data,probs = 0.9))+0.1
    mix <- tryCatch(getBestDist(data,init1,init2,ratios,type = "negativeBinomialMix",iter = iter),error = function(e){
      print(e)
      return(NULL)
    })
  } else if(type == "compare"){
    init2 <- c(abs(range(data)[1]-range(data)[2])/10,abs(range(data)[1]-range(data)[2])/10)
    mix1 <- tryCatch(getBestDist(data,init1,init2,ratios,type = "gaussianMix",iter = iter),error = function(e){
      print(e)
      return(NULL)
    })
    mix2 <- tryCatch(getBestDist(data,init1,init2,ratios,type = "gammaGaussianMix",iter = iter),error = function(e){
      print(e)
      return(error)
    })
  
    if(is.null(mix1) & !is.null(mix2)){
      mix = mix2
    }else if(!is.null(mix1) & is.null(mix2)){
      mix = mix1
    }else if(is.null(mix1) & is.null(mix2)){
      return(list("BIndex"=0,"mix"=NULL))
    }else{
      if(mix1$loglik >= mix2$loglik){
        mix = mix1
      }else{
        mix = mix2
      }
    }
  }
  
  if(is.null(mix)){
    return(list("BIndex"=0,"mix"=NULL))
  }
  res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
  
  
  mix$m.step$var[which(mix$m.step$var<1)] <- 1
  res[1, ] <- mix$m.step$mu
  res[2, ] <- mix$m.step$var
  res[3, ] <- mix$m.step$alpha

  smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2] + res[2, 1]) / 2)
  bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd
  
  
  if(plot==T){
    mu <- mix$m.step$mu[1]
    size <- mu+mu^2/mix$m.step$var[1]
    maxProb <- mix$m.step$alpha[1] * dnorm(x = as.integer(mix$m.step$mu[1]),mean = mix$m.step$mu[1],sd = sqrt(mix$m.step$var[1]))
    print(data.frame(x = data) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", 
                     fill = "white") +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mix$m.step$mu[1], sqrt(mix$m.step$var[1]), 
                                lam = mix$m.step$alpha[1], type=mix$type, maxProb=maxProb,trunc=trunc),
                    colour = "red", lwd = 1.5) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mix$m.step$mu[2], sqrt(mix$m.step$var[2]), 
                                lam = mix$m.step$alpha[2], type=ifelse(mix$type=="negativeBinomialMix","negativeBinomialMix","gaussianMix")),
                    colour = "blue", lwd = 1.5) +
      ylab("Density") +
      xlab("Values") +
      ggtitle("Final MM Fit"))
  }
  
  return(list("BIndex"=bi,"mix"=mix))
}


#' Expectation Step of the EM Algorithm
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point.
#'
#' @param first.vector Vector containing the mean of each component ["gaussianMix"], or the shape of the first component and the mean of the second ["gammaGaussianMix"]
#' @param second.vector Vector containing the standard deviation of each component ["gaussianMix"], or the rate of the first component and the standard deviation of the second ["gammaGaussianMix"]
#' @param alpha.vector Vector containing the mixing weights  of each component
#' @return Named list containing the loglik and posterior.df
e_step <- function(x, first.vector, second.vector, alpha.vector,type=c("gaussianMix","gammaGaussianMix","negativeBinomialMix")) {
  
  if(type=="gaussianMix"){
    comp1.prod <- dnorm(x, first.vector[1], second.vector[1]) * alpha.vector[1]
    
  }else if(type == "gammaGaussianMix"){
    if(sum(x<0)>0 ){
      stop("data contain negative or non-integer values and cannot be used to fit a gamma dist",call. = T)
      print("data contain negative values and cannot be used to fit a gamma dist")
      
    }
    shape = first.vector[1]/(second.vector[1]^2)
    scale = sqrt((second.vector[1]^2)/shape)
    comp1.prod <- dgamma(x, shape = shape, scale = scale) * alpha.vector[1]
    
  }else if(type == "negativeBinomialMix"){
    if(sum(x<0)>0| any(!is.integer(x))){
      stop("data contain negative values and cannot be used to fit a negative binomial dist")
      print("data contain negative values and cannot be used to fit a negative binomial dist")
    }
    mu = first.vector[1]
    size <- mu^2/(second.vector[1]^2-mu)
    if(mu==0){
      size <- 0.1
    }
    comp1.prod <- dnbinom(x, size = size,mu = mu) * alpha.vector[1]
 #   print(paste("mu1:",mu,", size1:",size))
    mu = first.vector[2]
    size <- mu^2/(second.vector[2]^2-mu)
    comp2.prod <- dnbinom(x, size = size,mu = mu) * alpha.vector[2]
#    print(paste("mu2:",mu,", size2:",size))
    
  }
#  print("wohoho3")
  if(type != "negativeBinomialMix"){
    comp2.prod <- dnorm(x, first.vector[2], second.vector[2]) * alpha.vector[2]
    
  }
  
#  print(comp2.prod)
  
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod / sum.of.comps
  comp2.post <- comp2.prod / sum.of.comps
  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post))
}

#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
m_step <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])
  
  comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
  comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)
  
  comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
  comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
  
  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)
  
  list("mu" = c(comp1.mu, comp2.mu),
       "var" = c(comp1.var, comp2.var),
       "alpha" = c(comp1.alpha, comp2.alpha))
}






getBestDist <- function(data,init1,init2,ratios,type,iter=100){
  
  loglik.vector <- NULL
  loglik.diff <- NULL
  for (i in 1:iter) {
    if (i == 1) {
      # Initialization
      e.step <- e_step(x = data,first.vector = init1 , second.vector = init2,alpha.vector = ratios,type = type)
      m.step <- m_step(data, e.step[["posterior.df"]])
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
      
      if(is.infinite(e.step$loglik)){
       # return(list("loglik"=-Inf,"m.step"=m.step,"type"=type,"iter"=i))
      }
      
    } else {
      # Repeat E and M steps till convergence
      e.step <- e_step(data, m.step[["mu"]], sqrt(m.step[["var"]]), 
                       m.step[["alpha"]],type=type)
      if(is.infinite(e.step$loglik) ){
        break
      }
      m.step <- m_step(data, e.step[["posterior.df"]])
    
      loglik.diff <- c(loglik.diff,abs((cur.loglik - e.step[["loglik"]])))
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])
      
      if(is.nan(loglik.diff[length(loglik.diff)]) & loglik.diff[length(loglik.diff)] < 1e-6) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }
  
  return(list("loglik"=loglik.vector[length(loglik.vector)],"m.step"=m.step,"type"=type,"iter"=i))
  
  
}

#' Plot a Mixture Component
#' 
#' @param x Input ata.
#' @param mu Mean of component.
#' @param sigma Standard of component.
#' @param lam Mixture weight of component.
plot_mix_comps <- function(x, mu, sigma, lam,type,maxProb=NULL,trunc=F) {
  if(type=="gaussianMix"){
    postR <- dnorm(x, mu, sigma)
    if(!is.null(maxProb) && trunc){
      ibx <- which(postR >= maxProb)
      if(length(ibx)>0){
        postR[ibx] <- maxProb
      }
      return(postR)
    }
    
    return(lam * postR)

    
  }else if(type=="gammaGaussianMix"){
    shape = mu/(sigma^2)
    scale = sqrt((sigma^2)/shape)
    return(lam * dgamma(x, shape = shape, scale = scale))
  }
  else if(type=="negativeBinomialMix"){
    size <- mu^2/(sigma^2-mu)
    return(lam * dnbinom(x,mu=mu, size = size))
  }
}



getBinaryValues <- function(expressionData, cutoffs){
  binaryExpression <- matrix(data = NA,nrow = dim(expressionData)[1], ncol = length(cutoffs))
  colnames(binaryExpression) <- names(cutoffs)
  rownames(binaryExpression) <- rownames(expressionData)
  
  for( gene in names(cutoffs)){ 
    binaryExpression[,gene] <- ifelse(expressionData[,gene]>cutoffs[gene],1,0)
  }
  
  return(binaryExpression)
}


getKappa <- function(expressionData1, expressionData2,chosenGenes){
  expressionData1 <- expressionData1[,chosenGenes]
  expressionData2 <- expressionData2[,chosenGenes]
  commonSamples <- intersect(rownames(expressionData1),rownames(expressionData2))
  
  kappas <- do.call(rbind,lapply(1:dim(expressionData1)[2], function(x){
    return(c(cohen.kappa(cbind(expressionData1[commonSamples,x],expressionData2[commonSamples,x]))$kappa,length(commonSamples)))
  }))
  
  rownames(kappas) <- chosenGenes
  return(kappas)
  
}



library(mccr)

getMCC <- function(expressionData1, expressionData2,chosenGenes){
  expressionData1 <- expressionData1[,chosenGenes]
  expressionData2 <- expressionData2[,chosenGenes]
  commonSamples <- intersect(rownames(expressionData1),rownames(expressionData2))
  
  corrs <- do.call(rbind,lapply(1:dim(expressionData1)[2], function(x){
    return(c(mccr(expressionData1[commonSamples,x],expressionData2[commonSamples,x]),length(commonSamples)))
  }))
  
  rownames(corrs) <- chosenGenes
  return(corrs)
  
}

getCor <- function(expressionData1, expressionData2,chosenGenes){
  expressionData1 <- expressionData1[,chosenGenes]
  expressionData2 <- expressionData2[,chosenGenes]
  commonSamples <- intersect(rownames(expressionData1),rownames(expressionData2))
  
  corrs <- do.call(rbind,lapply(1:dim(expressionData1)[2], function(x){
    return(c(cor(expressionData1[commonSamples,x],expressionData2[commonSamples,x],method = "s"),length(commonSamples)))
  }))
  
  rownames(corrs) <- chosenGenes
  return(corrs)
  
}


