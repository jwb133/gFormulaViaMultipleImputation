#simulation code for paper 'G-formula for causal inference via multiple imputation'

expit <- function(x) exp(x)/(1+exp(x))

library(mice)

#define a function which performs univariate approximate Bayesian bootstrap imputation
#code follows mice.impute.sample with required modification
mice.impute.uniABB <- function(y,ry,x = NULL, wy = NULL, ...) {
  if (is.null(wy)) {
    wy <- !ry
  }
  yry <- y[ry]
  
  #n=number of observed values
  n <- sum(ry)
  #draw of probability values
  unifDraws <- sort(runif(n-1))
  probDraw <- c(unifDraws,1) - c(0,unifDraws)
  #draw sample and return
  sample(x=yry,size=sum(wy),replace=TRUE,prob=probDraw)
}

#function to perform gformula via MI
gformulaViaMI <- function(obsData, M=100,l0ABB=FALSE,missingData,maxit=5) {
  
  n <- nrow(obsData)
  
  #create blank dataset with treatment indicators set as per desired regime
  syntheticDataBlank <- rbind(obsData, obsData)
  #syntheticDataBlank[,c("l1", "l2", "y")] <- NA
  syntheticDataBlank[,] <- NA
  syntheticDataBlank[1:n,"a0"] <- 0
  syntheticDataBlank[1:n,"a1"] <- 0
  syntheticDataBlank[1:n,"a2"] <- 0
  syntheticDataBlank[(n+1):(2*n),"a0"] <- 1
  syntheticDataBlank[(n+1):(2*n),"a1"] <- 1
  syntheticDataBlank[(n+1):(2*n),"a2"] <- 1
  
  predMat <- make.predictorMatrix(obsData)
  predMat[,] <- lower.tri(predMat)
  
  if (l0ABB==TRUE) {
    #use approximate Bayesian bootstrap for L0 imputation
    methodVal <- c("uniABB",rep("norm",ncol(obsData)-1))
  } else {
    methodVal <- c(rep("norm",ncol(obsData)))
  }
  
  #we repeat until all variance estimates are positive
  miVarEst <- c(0,0,0)
  impEsts <- array(0, dim=c(M,3))
  impVars <- array(0, dim=c(M,3))
  #integer to record attempt number
  attempt <- 1
  
  while(sum(miVarEst<=0)>0) {
  
    if (missingData==TRUE) {
      #impute any missing data first
      
      
      #use mice to impute missing data
      intermediateImps <- mice(obsData, m=M, defaultMethod = c("norm", "logreg", "polyreg", "polr"),
                               printFlag = FALSE, maxit=maxit)

      #now impute potential outcomes
      imputedDatasets <- vector(mode = "list", length = M)
      for (i in 1:M) {
        #impute data in synthetic part using mice, with m=1
        inputData <- rbind(complete(intermediateImps,action=i),
                             syntheticDataBlank)
        
        imps <- mice(data=inputData,method=methodVal, predictorMatrix = predMat, m=1,
                     maxit=1, printFlag = FALSE)

        imputedDatasets[[i]] <- complete(imps, action=1)
        #extract just the synthetic part
        imputedDatasets[[i]] <- imputedDatasets[[i]][(n+1):(3*n),]
        
      }
    } else {
      #setup without any regular missing data
      inputData <- rbind(obsData,syntheticDataBlank)
      imps <- mice(inputData, m=M, method=methodVal, predictorMatrix = predMat, maxit=1, printFlag = FALSE)
    }
    
    for (i in 1:M) {
      if (missingData==TRUE) {
        impData <-  imputedDatasets[[i]]
      } else {
        impData <- complete(imps,action=i)
        #extract just the synthetic part
        impData <- impData[(n+1):(3*n),]
      }
      
      #analyse it
      impEsts[i,1] <- mean(impData$y[impData$a0==0])
      impVars[i,1] <- var(impData$y[impData$a0==0])/n
      impEsts[i,2] <- mean(impData$y[impData$a0==1])
      impVars[i,2] <- var(impData$y[impData$a0==1])/n
      impEsts[i,3] <- impEsts[i,2] - impEsts[i,1]
      impVars[i,3] <- impVars[i,1] + impVars[i,2]
      
    }

    if (attempt==1) {
      finalImpEsts <- impEsts
      finalVars <- impVars
      attempt <- attempt + 1
    } else {
      finalImpEsts <- rbind(finalImpEsts, impEsts)
      finalVars <- rbind(finalVars, impVars)
      attempt <- attempt + 1
    }
    
    Vhat <- colMeans(finalVars)
    Bhat <- diag(var(finalImpEsts))
    miVarEst <- (1+1/nrow(finalImpEsts))*Bhat - Vhat
    
    
  }
    
  list(miEst=colMeans(finalImpEsts),miVarEst=miVarEst,M=nrow(finalImpEsts),Bhat=Bhat,Vhat=Vhat)
}

#function to perform gformula MI simulations
gformulaMISim <- function(nSim=1000,M=100,n=500, l0ABB=FALSE, 
                          missingData=FALSE, missingProp=0.5, progress=TRUE,maxit=5) {

  miEst <- array(0, dim=c(nSim,3))
  miVar <- array(0, dim=c(nSim,3))
  Bhat <- array(0, dim=c(nSim,3))
  Vhat <- array(0, dim=c(nSim,3))
  MRequired <- array(0, dim=c(nSim))
  
  for (sim in 1:nSim) {
    
    if (progress==TRUE) {
      print(sim)
    }
      
    l0 <- rnorm(n)
    a0 <- 1*(runif(n)<expit(l0))
    l1 <- l0+a0+rnorm(n)
    a1 <- 1*(runif(n)<expit(l1+a0))
    l2 <- l1+a1+rnorm(n)
    a2 <- 1*(runif(n)<expit(l2+a1))
    y <- l2+a2+rnorm(n)
    
    obsData <- data.frame(l0=l0,a0=a0,l1=l1,a1=a1,l2=l2,a2=a2,y=y)
    
    if (missingData==TRUE) {
      #make some data missing completely at random
      obsData$l1[runif(n)<missingProp] <- NA
      obsData$a1[runif(n)<missingProp] <- NA
      obsData$l2[runif(n)<missingProp] <- NA
      obsData$a2[runif(n)<missingProp] <- NA
      obsData$y[runif(n)<missingProp] <- NA
    }
    
    gformMIResult <- gformulaViaMI(obsData=obsData,M=M,l0ABB=l0ABB,
                                   missingData=missingData,maxit=maxit)
    
    miEst[sim,] <- gformMIResult$miEst
    miVar[sim,] <- gformMIResult$miVarEst
    MRequired[sim] <- gformMIResult$M
    Bhat[sim,] <- gformMIResult$Bhat
    Vhat[sim,] <- gformMIResult$Vhat
  
  }
  
  list(miEst=miEst, miVar=miVar, MRequired=MRequired, Bhat=Bhat, Vhat=Vhat)

}

#specify number of simulations to use throughout
numSims <- 10000

set.seed(738355)
mVals <- c(5,10,25,50,100)
mSims <- vector("list", length(mVals))
for (i in 1:length(mVals)) {
  print(paste("M=",mVals[i],sep=""))
  mSims[[i]] <- gformulaMISim(nSim=numSims,M=mVals[i],progress=FALSE)
}

resTable <- array(0, dim=c(length(mVals),8))
for (i in 1:length(mVals)) {
  print(i)
  #initial M
  resTable[[i,1]] <- mVals[i]
  #mean treatment effect
  resTable[[i,2]] <- round(mean(mSims[[i]]$miEst[,3]),3) - 3
  #empirical SD
  resTable[[i,3]] <- round(sd(mSims[[i]]$miEst[,3]),3)
  #mean SD estimate
  resTable[[i,4]] <- round(mean(sqrt(mSims[[i]]$miVar[,3])),3)
  #CI coverage
  tdf <- (mSims[[i]]$MRequired-1)*(1-(mSims[[i]]$MRequired*mSims[[i]]$Vhat[,3])/((mSims[[i]]$MRequired+1)*mSims[[i]]$Bhat[,3]))^2
  resTable[[i,5]] <- round(100*mean(1*((mSims[[i]]$miEst[,3]-qt(0.975,df=tdf)*sqrt(mSims[[i]]$miVar[,3])<3) &
                               (mSims[[i]]$miEst[,3]+qt(0.975,df=tdf)*sqrt(mSims[[i]]$miVar[,3])>3))),1)
  resTable[[i,6]] <- round(100*mean(1*((mSims[[i]]$miEst[,3]-1.96*sqrt(mSims[[i]]$miVar[,3])<3) &
                                         (mSims[[i]]$miEst[,3]+1.96*sqrt(mSims[[i]]$miVar[,3])>3))),1)
  #mean number of actual imputations performed
  resTable[[i,7]] <- round(mean(mSims[[i]]$MRequired),1)
  #max number of actual imputations performed
  resTable[[i,8]] <- max(mSims[[i]]$MRequired)
  #degrees of freedom
  #resTable[[i,8]] <- median(tdf)
}

resTable
colnames(resTable) <- c("M", "Bias", "Emp. SE", "Est. SE", "Raghu df 95% CI","Z 95% CI", "Mean M", "Max M")
library(xtable)
xtable(resTable, digits=c(0,0,3,3,3,1,1,1,0))

#approximate Bayesian bootstrap, at M=50
set.seed(738355)
abbSim <- gformulaMISim(nSim=numSims,M=50,progress=FALSE,l0ABB=TRUE)
mean(abbSim$miEst[,3])-3
sd(abbSim$miEst[,3])
mean(abbSim$miVar[,3]^0.5)
tdf <- (abbSim$MRequired-1)*(1-(abbSim$MRequired*abbSim$Vhat[,3])/((abbSim$MRequired+1)*abbSim$Bhat[,3]))^2
#Raghu df
100*mean(1*((abbSim$miEst[,3]-qt(0.975,df=tdf)*sqrt(abbSim$miVar[,3])<3) &
                    (abbSim$miEst[,3]+qt(0.975,df=tdf)*sqrt(abbSim$miVar[,3])>3)))
#N(0,1)
100*mean(1*((abbSim$miEst[,3]-1.96*sqrt(abbSim$miVar[,3])<3) &
              (abbSim$miEst[,3]+1.96*sqrt(abbSim$miVar[,3])>3)))


#now with missing data with different proportions of missing data MCAR
set.seed(738355)
missingProps <- c(0.05,0.10,0.25,0.50)
missingMiceSims <- vector("list", length(missingProps))
for (i in 1:length(missingProps)) {
  print(paste("Missingness proportion=",missingProps[i],sep=""))
  if (missingProps[i]==0.5) {
    #more iterations needed in mice for missing data imputation with 50% missingness
    missingMiceSims[[i]] <- gformulaMISim(nSim=numSims,M=50,missingData=TRUE,missingProp=missingProps[i],
                                          progress=FALSE,maxit=50)  
  } else {
    missingMiceSims[[i]] <- gformulaMISim(nSim=numSims,M=50,missingData=TRUE,missingProp=missingProps[i],
                              progress=FALSE)
  }
}

resTable <- array(0, dim=c(length(missingProps),6))
for (i in 1:length(missingProps)) {
  print(i)
  #missingness proportion
  resTable[[i,1]] <- missingProps[i]
  #mean treatment effect
  resTable[[i,2]] <- round(mean(missingMiceSims[[i]]$miEst[,3]),3) - 3
  #empirical SD
  resTable[[i,3]] <- round(sd(missingMiceSims[[i]]$miEst[,3]),3)
  #mean SD estimate
  resTable[[i,4]] <- round(mean(sqrt(missingMiceSims[[i]]$miVar[,3])),3)
  #CI coverage
  tdf <- (missingMiceSims[[i]]$MRequired-1)*(1-(missingMiceSims[[i]]$MRequired*missingMiceSims[[i]]$Vhat[,3])/((missingMiceSims[[i]]$MRequired+1)*missingMiceSims[[i]]$Bhat[,3]))^2
  resTable[[i,5]] <- round(100*mean(1*((missingMiceSims[[i]]$miEst[,3]-qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$miVar[,3])<3) &
                                         (missingMiceSims[[i]]$miEst[,3]+qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$miVar[,3])>3))),1)
  resTable[[i,6]] <- round(100*mean(1*((missingMiceSims[[i]]$miEst[,3]-1.96*sqrt(missingMiceSims[[i]]$miVar[,3])<3) &
                                         (missingMiceSims[[i]]$miEst[,3]+1.96*sqrt(missingMiceSims[[i]]$miVar[,3])>3))),1)
  #median number of actual imputations performed
  #resTable[[i,6]] <- round(mean(missingMiceSims[[i]]$MRequired),1)
  #max number of actual imputations performed
  #resTable[[i,7]] <- max(missingMiceSims[[i]]$MRequired)
}

resTable
colnames(resTable) <- c("\\pi", "Bias", "Emp. SE", "Est. SE", "Raghu df 95% CI","Z 95% CI")
xtable(resTable, digits=c(0,2,3,3,3,1,1))

