#simulation code for paper 'G-formula for causal inference via multiple imputation'

expit <- function(x) exp(x)/(1+exp(x))

library(mice)
library(gFormulaMI)

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
  
  if (l0ABB==TRUE) {
    #use approximate Bayesian bootstrap for L0 imputation
    methodVal <- c("uniABB",rep("norm",ncol(obsData)-1))
  } else {
    methodVal <- c(rep("norm",ncol(obsData)))
  }
  
  #we repeat until variance estimate is positive
  #save seed so we can restore if we have to increase number of imps
  startSeed <- .Random.seed
  currentM <- 0
  miVarEst <- 0
  
  while(miVarEst<=0) {
    
    #add an additional M imputations to what was used previously
    currentM <- currentM + M
    set.seed(startSeed)
  
    if (missingData==TRUE) {
      #use mice to impute missing data
      intermediateImps <- mice(obsData, m=currentM, defaultMethod = c("norm", "logreg", "polyreg", "polr"),
                               printFlag = FALSE, maxit=maxit)

      imps <- gFormulaImpute(intermediateImps,M=currentM,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)),
                             method=methodVal)
      
    } else {
      #setup without any regular missing data
      imps <- gFormulaImpute(obsData,M=currentM,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)),
                              method=methodVal)
    }
    
    #analyse imputed datasets
    fits <- with(imps, lm(y~factor(regime)))
    pooled <- syntheticPool(fits)
    miVarEst <- pooled[2,4]
  }
    
  list(miEst=pooled[2,1],miVarEst=pooled[2,4],M=currentM,Bhat=pooled[2,3],Vhat=pooled[2,2])
}

#function to perform gformula MI simulations
gformulaMISim <- function(nSim=1000,M=100,n=500, l0ABB=FALSE, 
                          missingData=FALSE, missingProp=0.5, progress=TRUE,maxit=5) {

  miEst <- array(0, dim=nSim)
  miVar <- array(0, dim=nSim)
  Bhat <- array(0, dim=nSim)
  Vhat <- array(0, dim=nSim)
  MRequired <- array(0, dim=nSim)
  
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
    
    miEst[sim] <- gformMIResult$miEst
    miVar[sim] <- gformMIResult$miVarEst
    MRequired[sim] <- gformMIResult$M
    Bhat[sim] <- gformMIResult$Bhat
    Vhat[sim] <- gformMIResult$Vhat
  
  }
  
  list(miEst=miEst, miVar=miVar, MRequired=MRequired, Bhat=Bhat, Vhat=Vhat)

}

#specify number of simulations to use throughout
numSims <- 10

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
  resTable[[i,2]] <- round(mean(mSims[[i]]$miEst),3) - 3
  #empirical SD
  resTable[[i,3]] <- round(sd(mSims[[i]]$miEst),3)
  #mean SD estimate
  resTable[[i,4]] <- round(mean(sqrt(mSims[[i]]$miVar)),3)
  #CI coverage
  tdf <- (mSims[[i]]$MRequired-1)*(1-(mSims[[i]]$MRequired*mSims[[i]]$Vhat)/((mSims[[i]]$MRequired+1)*mSims[[i]]$Bhat))^2
  resTable[[i,5]] <- round(100*mean(1*((mSims[[i]]$miEst-qt(0.975,df=tdf)*sqrt(mSims[[i]]$miVar)<3) &
                               (mSims[[i]]$miEst+qt(0.975,df=tdf)*sqrt(mSims[[i]]$miVar)>3))),1)
  resTable[[i,6]] <- round(100*mean(1*((mSims[[i]]$miEst-1.96*sqrt(mSims[[i]]$miVar)<3) &
                                         (mSims[[i]]$miEst+1.96*sqrt(mSims[[i]]$miVar)>3))),1)
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
mean(abbSim$miEst)-3
sd(abbSim$miEst)
mean(abbSim$miVar^0.5)
tdf <- (abbSim$MRequired-1)*(1-(abbSim$MRequired*abbSim$Vhat)/((abbSim$MRequired+1)*abbSim$Bhat))^2
#Raghu df
100*mean(1*((abbSim$miEst-qt(0.975,df=tdf)*sqrt(abbSim$miVar)<3) &
                    (abbSim$miEst+qt(0.975,df=tdf)*sqrt(abbSim$miVar)>3)))
#N(0,1)
100*mean(1*((abbSim$miEst-1.96*sqrt(abbSim$miVar)<3) &
              (abbSim$miEst+1.96*sqrt(abbSim$miVar)>3)))


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
  resTable[[i,2]] <- round(mean(missingMiceSims[[i]]$miEst),3) - 3
  #empirical SD
  resTable[[i,3]] <- round(sd(missingMiceSims[[i]]$miEst),3)
  #mean SD estimate
  resTable[[i,4]] <- round(mean(sqrt(missingMiceSims[[i]]$miVar)),3)
  #CI coverage
  tdf <- (missingMiceSims[[i]]$MRequired-1)*(1-(missingMiceSims[[i]]$MRequired*missingMiceSims[[i]]$Vhat)/((missingMiceSims[[i]]$MRequired+1)*missingMiceSims[[i]]$Bhat))^2
  resTable[[i,5]] <- round(100*mean(1*((missingMiceSims[[i]]$miEst-qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$miVar)<3) &
                                         (missingMiceSims[[i]]$miEst+qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$miVar)>3))),1)
  resTable[[i,6]] <- round(100*mean(1*((missingMiceSims[[i]]$miEst-1.96*sqrt(missingMiceSims[[i]]$miVar)<3) &
                                         (missingMiceSims[[i]]$miEst+1.96*sqrt(missingMiceSims[[i]]$miVar)>3))),1)
  #median number of actual imputations performed
  #resTable[[i,6]] <- round(mean(missingMiceSims[[i]]$MRequired),1)
  #max number of actual imputations performed
  #resTable[[i,7]] <- max(missingMiceSims[[i]]$MRequired)
}

resTable
colnames(resTable) <- c("\\pi", "Bias", "Emp. SE", "Est. SE", "Raghu df 95% CI","Z 95% CI")
xtable(resTable, digits=c(0,2,3,3,3,1,1))

