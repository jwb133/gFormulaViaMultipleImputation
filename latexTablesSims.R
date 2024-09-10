# construct latex for tables


#function to print latex results table for tables 2 and 4
constructResTable <- function(resultsSet) {
  
  resTable <- array(0, dim=c(length(resultsSet),8))
  for (i in 1:length(resultsSet)) {
    #initial M
    resTable[[i,1]] <- resultsSet[[i]]$initialM
    #mean treatment effect (true value is 3)
    resTable[[i,2]] <- round(mean(resultsSet[[i]]$est),3) - 3
    #empirical SD
    resTable[[i,3]] <- round(sd(resultsSet[[i]]$est),3)
    #mean SD estimate
    resTable[[i,4]] <- round(mean(sqrt(resultsSet[[i]]$var)),3)
    #CI coverage
    tdf <- (resultsSet[[i]]$finalM-1)*(1-(resultsSet[[i]]$finalM*resultsSet[[i]]$Vhat)/((resultsSet[[i]]$finalM+1)*resultsSet[[i]]$Bhat))^2
    resTable[[i,5]] <- round(100*mean(1*((resultsSet[[i]]$est-qt(0.975,df=tdf)*sqrt(resultsSet[[i]]$var)<3) &
                                           (resultsSet[[i]]$est+qt(0.975,df=tdf)*sqrt(resultsSet[[i]]$var)>3))),1)
    resTable[[i,6]] <- round(100*mean(1*((resultsSet[[i]]$est-1.96*sqrt(resultsSet[[i]]$var)<3) &
                                           (resultsSet[[i]]$est+1.96*sqrt(resultsSet[[i]]$var)>3))),1)
    #mean number of actual imputations performed
    resTable[[i,7]] <- round(mean(resultsSet[[i]]$finalM),1)
    #max number of actual imputations performed
    resTable[[i,8]] <- max(resultsSet[[i]]$finalM)
  }
  
  colnames(resTable) <- c("M", "Bias", "Emp. SE", "Est. SE", 
                          "Raghu df 95% CI","Z 95% CI",
                          "Mean M", "Max M")
  resTable
}

# table 2
load("table2Results.RData")
constructResTable(noMissingDataSims)

library(xtable)
print(xtable(constructResTable(noMissingDataSims), 
             digits=c(0,0,3,3,3,1,1,1,0)),
      include.rownames = FALSE)

# table 3 gfoRmula table
load("table3Results.RData")
gfoRmulaResTable <- array(0, dim=c(length(gfoRmulaSims),3))
for (i in 1:length(gfoRmulaSims)) {
  #nsimul
  gfoRmulaResTable[[i,1]] <- nsimulVals[i]
  #mean treatment effect (true value is 3)
  gfoRmulaResTable[[i,2]] <- round(mean(gfoRmulaSims[[i]]$est),3) - 3
  #empirical SD
  gfoRmulaResTable[[i,3]] <- round(sd(gfoRmulaSims[[i]]$est),3)
}
gfoRmulaResTable


# table 4 simulations - no missing data, n_obs=100, and Approximate Bayesian bootstrap
load("table4Results.RData")
print(xtable::xtable(constructResTable(noMissingDataSimsnObsSmall), 
                     digits=c(0,0,3,3,3,1,1,1,0)),
      include.rownames = FALSE)

# approximate Bayesian bootstrap
mean(abbSim$est)-3
sd(abbSim$est)
mean(abbSim$var^0.5)
tdf <- (abbSim$finalM-1)*(1-(abbSim$finalM*abbSim$Vhat)/((abbSim$finalM+1)*abbSim$Bhat))^2
#Raghu df
100*mean(1*((abbSim$est-qt(0.975,df=tdf)*sqrt(abbSim$var)<3) &
              (abbSim$est+qt(0.975,df=tdf)*sqrt(abbSim$var)>3)))
#N(0,1)
100*mean(1*((abbSim$est-1.96*sqrt(abbSim$var)<3) &
              (abbSim$est+1.96*sqrt(abbSim$var)>3)))


#table 5 simulations - with missing data, MCAR
load("table5Results.RData")
constructResTable(missingMiceSims)

resTable <- array(0, dim=c(length(missingProps),6))
for (i in 1:length(missingProps)) {
  print(i)
  #missingness proportion
  resTable[[i,1]] <- missingProps[i]
  #mean treatment effect
  resTable[[i,2]] <- round(mean(missingMiceSims[[i]]$est),3) - 3
  #empirical SD
  resTable[[i,3]] <- round(sd(missingMiceSims[[i]]$est),3)
  #mean SD estimate
  resTable[[i,4]] <- round(mean(sqrt(missingMiceSims[[i]]$var)),3)
  #CI coverage
  tdf <- (missingMiceSims[[i]]$finalM-1)*(1-(missingMiceSims[[i]]$finalM*missingMiceSims[[i]]$Vhat)/((missingMiceSims[[i]]$finalM+1)*missingMiceSims[[i]]$Bhat))^2
  resTable[[i,5]] <- round(100*mean(1*((missingMiceSims[[i]]$est-qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$var)<3) &
                                         (missingMiceSims[[i]]$est+qt(0.975,df=tdf)*sqrt(missingMiceSims[[i]]$var)>3))),1)
  resTable[[i,6]] <- round(100*mean(1*((missingMiceSims[[i]]$est-1.96*sqrt(missingMiceSims[[i]]$var)<3) &
                                         (missingMiceSims[[i]]$est+1.96*sqrt(missingMiceSims[[i]]$var)>3))),1)
}

resTable
colnames(resTable) <- c("\\pi", "Bias", "Emp. SE", "Est. SE", "Raghu df 95% CI","Z 95% CI")
xtable::xtable(resTable, digits=c(0,2,3,3,3,1,1))