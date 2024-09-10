#simulation code for paper 'G-formula for causal inference via multiple imputation'

expit <- function(x) exp(x)/(1+exp(x))

library(mice)
library(gFormulaMI)
library(gfoRmula)

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

#function to perform gformula via MI, increasing M or nSyn if necessary to get a positive variance estimate
gformulaViaMI <- function(obsData, M=100,nSynMultiplier=1,increaseM=TRUE,
                          l0ABB=FALSE,missingData,maxit=5) {
  
  n <- nrow(obsData)
  
  if (l0ABB==TRUE) {
    #use approximate Bayesian bootstrap for L0 imputation
    methodVal <- c("uniABB",rep("norm",ncol(obsData)-1))
  } else {
    methodVal <- c(rep("norm",ncol(obsData)))
  }
  
  #we repeat until variance estimate is positive
  #save seed so we can restore if we have to increase number of imps or nSyn
  startSeed <- .Random.seed
  
  if (increaseM==TRUE) {
    currentM <- 0
    currentnSyn <- n*nSynMultiplier
  } else {
    #increase nSyn
    currentM <- M
    currentnSyn <- 0
  }
  
  miVarEst <- 0
  
  while(miVarEst<=0) {
    
    .Random.seed <- startSeed
    
    if (increaseM==TRUE) {
      #add an additional M imputations to what was used previously
      currentM <- currentM + M
    } else {
      #increase nSyn / nSim
      currentnSyn <- currentnSyn + n*nSynMultiplier
    }

    if (missingData==TRUE) {
      #use mice to impute missing data
      intermediateImps <- mice(obsData, m=currentM, defaultMethod = c("norm", "logreg", "polyreg", "polr"),
                               printFlag = FALSE, maxit=maxit)
    } else {
      #if no missing data, set this object to the observed data
      intermediateImps <- obsData
    }

    #perform G-formula imputation
    imps <- gFormulaImpute(intermediateImps,M=currentM,nSim=currentnSyn,
                             trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)),
                             method=methodVal, silent=TRUE)
    
    #analyse imputed datasets
    fits <- with(imps, lm(y~factor(regime)))
    try({pooled <- syntheticPool(fits)})
    if (exists("pooled")==TRUE) {
      miVarEst <- pooled[2,4]
    } else {
      miVarEst <- 0
    }
  }
    
  list(miEst=pooled[2,1],miVarEst=pooled[2,4],M=currentM,nSyn=currentnSyn,
       Bhat=pooled[2,3],Vhat=pooled[2,2])
}

#function to perform gformula using gfoRmula
gfoRmulaRun <- function(obsData, nsimul=NULL) {
  
  #first need to reshape to long format
  longData <- reshape(obsData, direction="long", varying=list(c("a0","a1","a2"),
                                                              c("l0","l1","l2")),
                      v.names=c("a","l"),
                      timevar="time",
                      times=c(0,1,2))
  longData <- longData[order(longData$id,longData$t),]
  #set y to missing at times 0 and 1
  longData$y[longData$t<2] <- NA
  
  id <- 'id'
  time_name <- 'time'
  covnames <- c('l', 'a')
  outcome_name <- 'y'
  covtypes <- c('normal', 'binary')
  histories <- c(lagged)
  histvars <- list(c('a', 'l'))
  covparams <- list(covmodels = c(l ~ lag1_a + lag1_l + lag2_a + lag2_l + factor(time),
                                  a ~ lag1_a + lag1_l + lag2_a + lag2_l + factor(time)))
  ymodel <- y ~ a + lag1_a + lag2_a + l + lag1_l + lag2_l
  intvars <- list('a', 'a')
  interventions <- list(list(c(static, rep(0, 3))),
                        list(c(static, rep(1, 3))))
  int_descript <- c('Never treat', 'Always treat')
  
  #save seed so we can restore random number state afterwards
  startSeed <- .Random.seed
  gform_cont_eof <- gformula_continuous_eof(obs_data = data.table::as.data.table(longData),
                                            id = id,
                                            time_name = time_name,
                                            covnames = covnames,
                                            outcome_name = outcome_name,
                                            covtypes = covtypes,
                                            covparams = covparams, ymodel = ymodel,
                                            intvars = intvars,
                                            interventions = interventions,
                                            int_descript = int_descript,
                                            histories = histories, histvars = histvars,
                                            nsimul = nsimul,
                                            seed=sample.int(2^30, size = 1),
                                            model_fits=TRUE)
  .Random.seed <- startSeed
  
  list(gfoRmulaEst = gform_cont_eof$result[3,4]-gform_cont_eof$result[2,4])
}

simData <- function(n, missingData=FALSE, missingProp=0.5) {
  
  #simulate data
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
  
  obsData
}

#function to perform gformula MI simulations
gformulaMISim <- function(nSim=1000,M=100,n=500,l0ABB=FALSE,nSynMultiplier,increaseM=TRUE,
                          missingData=FALSE, missingProp=0.5, progress=TRUE,maxit=5) {

  #set up lists/arrays to store results
  resultList <- list(initialM=M,
                      initialk=nSynMultiplier,
                      est=array(0, dim=nSim),
                      var=array(0, dim=nSim),
                      Bhat=array(0, dim=nSim),
                      Vhat=array(0, dim=nSim),
                      finalM=array(0, dim=nSim),
                      finalnSyn=array(0, dim=nSim))
  #run simulations
  for (sim in 1:nSim) {
    
    if (progress==TRUE) {
      print(sim)
    }
    
    obsData <- simData(n=n,missingData=missingData,missingProp=missingProp)
    
    #gFormula via MI, increasing M if necessary
    gformMIResult <- gformulaViaMI(obsData=obsData,M=M,increaseM=increaseM,l0ABB=l0ABB,
                                   nSynMultiplier=nSynMultiplier,
                                   missingData=missingData,maxit=maxit)
    resultList$est[sim] <- gformMIResult$miEst
    resultList$var[sim] <- gformMIResult$miVarEst
    resultList$Bhat[sim] <- gformMIResult$Bhat
    resultList$Vhat[sim] <- gformMIResult$Vhat
    resultList$finalM[sim] <- gformMIResult$M
    resultList$finalnSyn[sim] <- gformMIResult$nSyn
    
  }
  
  resultList

}

#function to perform gfoRmula simulations
gfoRmulaSim <- function(nSim=1000,n=500,nsimul=500,missingData=FALSE,
                        missingProp=0.5, progress=TRUE) {
  
  #set up lists/arrays to store results
  resultList <- list(est=array(0, dim=nSim))
  
  #run simulations
  for (sim in 1:nSim) {
    
    if (progress==TRUE) {
      print(sim)
    }
    
    obsData <- simData(n=n,missingData=missingData,missingProp=missingProp)
    
    #gfoRmula package  
    gfoRmulaResult <- gfoRmulaRun(obsData=obsData, nsimul=nsimul)
    resultList$est[sim] <- as.numeric(gfoRmulaResult$gfoRmulaEst)
    
  }
  
  resultList
  
}

#specify number of simulations to use throughout
numSims <- 10000
setwd("D:/Git/gFormulaViaMultipleImputationPaper")

#table 2 simulations - no missing data, varying initial M and n_syn
set.seed(738355)
mVals <- c(5,10,25,50,100)
kVals <- c(1,2,5,10)
noMissingDataSims <- vector("list", length(mVals)*length(kVals))
for (i in 1:length(kVals)) {
  for (j in 1:length(mVals)) {
    noMissingDataSims[[(i-1)*length(mVals)+j]] <- gformulaMISim(nSim=numSims,M=mVals[j],
                                                                nSynMultiplier=kVals[i],progress=FALSE)
  }
}
save(noMissingDataSims, file="table2Results.RData")

# table 3 gfoRmula, no missing data
set.seed(8946565)
nsimulVals <- c(500,1000,2500,5000)
gfoRmulaSims <- vector("list", length(nsimulVals))
for (i in 1:length(nsimulVals)) {
  print(i)
  gfoRmulaSims[[i]] <- gfoRmulaSim(nSim=numSims,nsimul = nsimulVals[i],progress=FALSE)
}
save(nsimulVals, gfoRmulaSims, file="table3Results.RData")

# table 4 simulations - no missing data, n_obs=100, and Approximate Bayesian bootstrap
set.seed(791235091)
mVals <- c(5,10,25,50,100)
kVals <- c(1,2,5,10)
noMissingDataSimsnObsSmall <- vector("list", length(mVals)*length(kVals))
for (i in 1:length(kVals)) {
  for (j in 1:length(mVals)) {
    noMissingDataSimsnObsSmall[[(i-1)*length(mVals)+j]] <- gformulaMISim(n=100,nSim=numSims,M=mVals[j],
                                                                         nSynMultiplier=kVals[i],progress=FALSE)
  }
}

#approximate Bayesian bootstrap, at M=50
set.seed(738355)
abbSim <- gformulaMISim(nSim=numSims,M=50,progress=FALSE,l0ABB=TRUE,nSynMultiplier=1)

save(noMissingDataSimsnObsSmall, abbSim, file="table4Results.RData")

#table 5 simulations - with missing data, MCAR
set.seed(738355)
missingProps <- c(0.05,0.10,0.25,0.50)
missingMiceSims <- vector("list", length(missingProps))
for (i in 1:length(missingProps)) {
  print(paste("Missingness proportion=",missingProps[i],sep=""))
  if (missingProps[i]==0.5) {
    #more iterations needed in mice for missing data imputation with 50% missingness
    missingMiceSims[[i]] <- gformulaMISim(nSim=numSims,M=50,missingData=TRUE,missingProp=missingProps[i],
                                          progress=FALSE,maxit=50,nSynMultiplier=1)  
  } else {
    missingMiceSims[[i]] <- gformulaMISim(nSim=numSims,M=50,missingData=TRUE,missingProp=missingProps[i],
                                          progress=FALSE,nSynMultiplier=1)
  }
}
save(missingProps, missingMiceSims, file="table5Results.RData")