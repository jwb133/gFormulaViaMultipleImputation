#------------------------------------------------------------------------------------------------------------
# R code accompanying "G-formula for causal inference via multiple imputation" 
# Author: Emily Granger
# Date last modified: 22/08/2024

#This script includes analysis code for the illustrative example of the use of the G-formula via MI approach.
#The analysis investigates the effects of DNase and hypertonic saline used in combination compared to the use
#of DNase alone on lung function, in people with cystic fibrosis. 

#We use data from the UK CF Registry. Data are available following application to the UK CF Registry Research 
#Committee. 
#https://www.cysticfibrosis.org.uk/the-work-we-do/uk-cf-registry/apply-for-data-from-the-uk-cf-registry.   

#Code is provided for the following analyses: 

# 1a analysis on a complete cases using the gFormulaMI package
# 1b analysis on a complete cases using the gfoRmula package
# 2a analysis with missing data handled using multiple imputation (via the gFormulaMI package)
# 2b analysis with missing data handled using complete case analysis (via the gfoRmula package)
#-----------------------------------------------------------------------------------------------------------

#Load libraries
library(tidyr)
library(dplyr)
library(mice)
library(gFormulaMI)
library(data.table)
library(gfoRmula)

#Load data
rm(list=ls())
load(".../analysis_dat.RData")

#-----------------------------------------------------------------------------------------------------------
# 0. Description of the data
#-----------------------------------------------------------------------------------------------------------

#analysis_dat.RData is a workspace that includes four dataframes: 
# (1) dat_wide (data on all included individuals, wide format)
# (2) dat_long (data on all included individuals, long format)
# (3) dat_cc_wide (data on complete cases only, wide format)
# (4) dat_cc_long (data on complete cases only, long format)

#dat_long contains 23795 observations on 18 variables. There are 4759 individuals and each individual has 
#5 rows of data, one for each year of follow-up. The variables included in this dataframe are:

#Patient ID: "patient_id"
#Follow-up year: "fu_year" (takes values 1-5)

#Indicator for treatment group: "trt_grp" (factor: DNase, Hypertonic saline, Both, Neither)

#Outcome variable: "gli_percpredfev" 

#Confounding variables:
#Sex: "dmg_sex"
#CFTR Genotype: "genotype_class" (factor: high, low or none_assigned)
#Ethnicity: "ethnicity" (binary variable: white or non-white)
#Baseline age: "Bage" (age when fu_year==1)
#Rate of decline in FEV1%: "slopes" (estimate of the rate of decline in FEV1% prior to study entry)
#Past FEV1%: "lag_fev" (at time k, this is equal to FEV1% measured at time k-1)
#Baseline FEV: "Bfev" (FEV1% measured when fu_year==1)
#P.aeruginosa: "s05culturespeciespseudoaeruginos" (binary: yes/no)
#Staphylococcus aureus: "s05culturespeciesstaph" (binary: yes/no)
#NTMs: "ntm_combined" (binary: yes/no)
#IV hospital admissions: "hospIV" (binary: yes/no for any hospital admissions during the past year)
#BMI: "lag_bmi_zscore" (at time k, this is equal to the BMI z-score measured at time k-1)
#Pancreatic insufficiency: "PI" (binary: yes/no)
#Number of days on IV antibiotics: "IVdays_cat" (Factor: 1 (0 days), 2 (1-14 days), 3(15-28 days), 4(28+ days) )

#dat_wide contains 4759 observations (one for each individual) on 51 variables.

#Time-constant confounders: "dmg_sex"; "genotype_class"; "ethnicity"; "Bage"; "slopes"; "Bfev"
#Time-varying confounders:"trt_grp.k"; "gli_percpredfev.k"; "s05culturespeciespseudoaeruginos.k"; "s05culturespeciesstaph.k"; 
#                       "ntm_combined.k"; "hospIV.k"; "lag_bmi_zscore.k"; "PI.k"; "IVdays_cat.k"

#Each of the time-varying confounders have five variables (k=1,2,3,4,5) where the kth variable represents the confounder
#measured at time k. 

#-------------------------------------------------------------------------------
# 1. Analysis on complete cases
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1a. Using gFormulaMI package

#Remove id and imputation number from wide dataset
dat_cc_wide = subset(dat_cc_wide, select = -c(.id, .imp))

#Start timer
start.time<-proc.time()

#Impute synthetic datasets under four treatment regimes
set.seed(04111990)
imps <- gFormulaImpute(data=dat_cc_wide, M=200, trtVars=c("trt_grp.1",
                       "trt_grp.2","trt_grp.3","trt_grp.4","trt_grp.5"),
                        trtRegimes=list(c(1,1,1,1,1),c(2,2,2,2,2),c(3,3,3,3,3),
                        c(4,4,4,4,4)))

#Fit outcome model
fits1 <- with(imps, lm(gli_percpredfev.1~factor(regime)))
fits2 <- with(imps, lm(gli_percpredfev.2~factor(regime)))
fits3 <- with(imps, lm(gli_percpredfev.3~factor(regime)))
fits4 <- with(imps, lm(gli_percpredfev.4~factor(regime)))
fits5 <- with(imps, lm(gli_percpredfev.5~factor(regime)))

#Use synthetic imputation combination rules
results1<-syntheticPool(fits1)
results2<-syntheticPool(fits2)
results3<-syntheticPool(fits3)
results4<-syntheticPool(fits4)
results5<-syntheticPool(fits5)

#Display results
results1; results2; results3; results4; results5

#End timer
proc.time() - start.time

#Estimate Monte Carlo Standard Errors
sqrt(results1[,4]/200)
sqrt(results2[,4]/200)
sqrt(results3[,4]/200)
sqrt(results4[,4]/200)
sqrt(results5[,4]/200)

#-------------------------------------------------------------------------------
# 1b. Using gfoRmula package

#------------------------------------------------------------
#gfoRmula package required time variable to start at 0
dat_cc_long$fu_year<-dat_cc_long$fu_year-1
setDT(dat_cc_long)

#----------------------
#Analysis

#Start timer
start.time<-proc.time()

#Define inputs for gfoRmula
id<-"id"
time_name<-"fu_year"
covnames<-c("trt_grp","s05culturespeciespseudoaeruginos", "s05culturespeciesstaph",
            "hospIV", "IVdays_cat", "PI", "ntm_combined",
            "lag_fev", "lag_bmi_zscore")
basecovs<-c("genotype_class", "ethnicity", "Bfev", "Bage", "dmg_sex","slopes")
outcome_name<-"gli_percpredfev"
covtypes<-c("categorical","binary","binary","binary",
            "categorical","binary","binary","normal", "normal")
histories<-c(lagged)
histvars<-list(c("trt_grp","s05culturespeciespseudoaeruginos", "s05culturespeciesstaph",
                 "hospIV",  "IVdays_cat", "PI", "ntm_combined",
                 "lag_fev", "lag_bmi_zscore"))

covparams<-list(covmodels=c(
  
  trt_grp~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  s05culturespeciespseudoaeruginos~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  s05culturespeciesstaph~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  hospIV~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  IVdays_cat~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  PI~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  ntm_combined~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  lag_fev~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  lag_bmi_zscore~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined
  
))
#Define interventions
intvars <- list("trt_grp", "trt_grp", "trt_grp", "trt_grp")
int_descript <- c("DNase only", "HS only", "No treatment", "HS & DNase")
nsimul <- 100000
int_times <- list(list(0:4), list=(0:4), list=(0:4), list=(0:4))
interventions <- list(list(c(static, rep("1",5))),list(c(static, rep("2",5))),list(c(static, rep("3",5))),list(c(static, rep("4",5))))

#Outcome model
ymodel<-gli_percpredfev~as.factor(trt_grp)+as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
  as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
  s05culturespeciespseudoaeruginos+lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
  s05culturespeciesstaph+lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
  hospIV+lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
  as.factor(IVdays_cat)+as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
  lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
  lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
  PI+lag1_PI+lag2_PI+lag3_PI+lag4_PI+
  ntm_combined+lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year

#Use gformula_continuous_eof
gform <- gformula_continuous_eof(obs_data = dat_cc_long,
                                 id = id,
                                 time_name = time_name,
                                 basecovs = basecovs,
                                 covnames = covnames,
                                 outcome_name = outcome_name,
                                 covtypes = covtypes,
                                 covparams = covparams, 
                                 ymodel = ymodel,
                                 intvars = intvars,
                                 interventions = interventions,
                                 int_descript = int_descript,
                                 histories = histories, histvars = histvars,
                                 nsimul = nsimul, 
                                 sim_data_b = TRUE,
                                 seed = 04111990)

#Calculate effect estimates
EY45<-mean(gform$sim_data$'HS & DNase'$Ey[gform$sim_data$'HS & DNase'$fu_year==4], na.rm=TRUE)
EY15<-mean(gform$sim_data$'DNase only'$Ey[gform$sim_data$'DNase only'$fu_year==4], na.rm=TRUE)
CE43.5<-EY45-EY15

EY44<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==4], na.rm=TRUE)
EY14<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==4], na.rm=TRUE)
CE43.4<-EY44-EY14

EY43<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==3], na.rm=TRUE)
EY13<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==3], na.rm=TRUE)
CE43.3<-EY43-EY13

EY42<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==2], na.rm=TRUE)
EY12<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==2], na.rm=TRUE)
CE43.2<-EY42-EY12

EY41<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==1], na.rm=TRUE)
EY11<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==1], na.rm=TRUE)
CE43.1<-EY41-EY11

#Save estimated outcomes at each time point for MCSE estimates
Y.DNHS.5 <- gform$sim_data$'HS & DNase'$Ey[gform$sim_data$'HS & DNase'$fu_year==4]
Y.DNHS.4 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==4]
Y.DNHS.3 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==3]
Y.DNHS.2 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==2]
Y.DNHS.1 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==1]

Y.DN.5 <- gform$sim_data$'DNase only'$Ey[gform$sim_data$'DNase only'$fu_year==4]
Y.DN.4 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==4]
Y.DN.3 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==3]
Y.DN.2 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==2]
Y.DN.1 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==1]

#bootstrap.R runs the above analysis 1000 times in a loop to obtain 
#bootstrap estimates of the SE
source(".../bootstrap.R")

#End timer
proc.time()-start.time

#bootstrap.R results in 5 vectors: ce43.1, ce43.2, ce43.3, ce43.4, ce43.5
#each vector contains 1000 bootstrap estimates of the treatment effects for years 1-5
c(CE43.1, CE43.1-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.1),CE43.1+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.1))
c(CE43.2, CE43.2-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.2),CE43.2+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.2))
c(CE43.3, CE43.3-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.3),CE43.3+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.3))
c(CE43.4, CE43.4-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.4),CE43.4+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.4))
c(CE43.5, CE43.5-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.5),CE43.5+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.5))

#calculate Monte-Carlo SE of difference in means assuming they are two independent samples
sqrt(var(Y.DNHS.5)/nsimul + var(Y.DN.5)/nsimul)
sqrt(var(Y.DNHS.4)/nsimul + var(Y.DN.4)/nsimul)
sqrt(var(Y.DNHS.3)/nsimul + var(Y.DN.3)/nsimul)
sqrt(var(Y.DNHS.2)/nsimul + var(Y.DN.2)/nsimul)
sqrt(var(Y.DNHS.1)/nsimul + var(Y.DN.1)/nsimul)

#-------------------------------------------------------------------------------
# 2. Analysis on data with missingness 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2a. Multiple imputation approach using gFormulaMI

#Start timer
start.time<-proc.time()

#Impute missing data
impdata<-mice(dat_wide, defaultMethod = c("norm.predict", "logreg", "polyreg", 
                                          "polr"), m=200, seed=24042023)

#Impute synthetic datasets under four treatment regimes
set.seed(04111990)
imps <- gFormulaImpute(data=impdata, M=200, trtVars=c("trt_grp.1","trt_grp.2",
                                            "trt_grp.3","trt_grp.4","trt_grp.5"),
                                            trtRegimes=list(c(1,1,1,1,1),
                                            c(2,2,2,2,2),c(3,3,3,3,3), c(4,4,4,4,4)))

#Fit outcome model
fits1 <- with(imps, lm(gli_percpredfev.1~factor(regime)))
fits2 <- with(imps, lm(gli_percpredfev.2~factor(regime)))
fits3 <- with(imps, lm(gli_percpredfev.3~factor(regime)))
fits4 <- with(imps, lm(gli_percpredfev.4~factor(regime)))
fits5 <- with(imps, lm(gli_percpredfev.5~factor(regime)))

#Use synthetic imputation combination rules
results1<-syntheticPool(fits1)
results2<-syntheticPool(fits2)
results3<-syntheticPool(fits3)
results4<-syntheticPool(fits4)
results5<-syntheticPool(fits5)

#Display results 
results1; results2; results3; results4; results5

#End timer
proc.time() - start.time

#Estimate Monte Carlo Standard Errors
sqrt(results1[,4]/200)
sqrt(results2[,4]/200)
sqrt(results3[,4]/200)
sqrt(results4[,4]/200)
sqrt(results5[,4]/200)

#-------------------------------------------------------------------------------
# 2b. Complete cases analysis using gfoRmula


#------------------------------------------------------------
#gfoRmula package required time variable to start at 0
dat_long$fu_year<-dat_cc_long$fu_year-1
setDT(dat_long)

#----------------------
#Analysis

#Start timer
start.time<-proc.time()

#Define inputs for gfoRmula
id<-"id"
time_name<-"fu_year"
covnames<-c("trt_grp","s05culturespeciespseudoaeruginos", "s05culturespeciesstaph",
            "hospIV", "IVdays_cat", "PI", "ntm_combined",
            "lag_fev", "lag_bmi_zscore")
basecovs<-c("genotype_class", "ethnicity", "Bfev", "Bage", "dmg_sex","slopes")
outcome_name<-"gli_percpredfev"
covtypes<-c("categorical","binary","binary","binary",
            "categorical","binary","binary","normal", "normal")
histories<-c(lagged)
histvars<-list(c("trt_grp","s05culturespeciespseudoaeruginos", "s05culturespeciesstaph",
                 "hospIV",  "IVdays_cat", "PI", "ntm_combined",
                 "lag_fev", "lag_bmi_zscore"))

covparams<-list(covmodels=c(
  
  trt_grp~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  s05culturespeciespseudoaeruginos~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  s05culturespeciesstaph~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  hospIV~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  IVdays_cat~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  PI~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  ntm_combined~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  lag_fev~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag_bmi_zscore+lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined,
  
  lag_bmi_zscore~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag1_lag_fev+lag2_lag_fev+lag3_lag_fev+
    lag1_lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined
  
))
#Define interventions
intvars <- list("trt_grp", "trt_grp", "trt_grp", "trt_grp")
int_descript <- c("DNase only", "HS only", "No treatment", "HS & DNase")
nsimul <- 100000
int_times <- list(list(0:4), list=(0:4), list=(0:4), list=(0:4))
interventions <- list(list(c(static, rep("1",5))),list(c(static, rep("2",5))),list(c(static, rep("3",5))),list(c(static, rep("4",5))))

#Outcome model
ymodel<-gli_percpredfev~as.factor(trt_grp)+as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
  as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
  s05culturespeciespseudoaeruginos+lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
  s05culturespeciesstaph+lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
  hospIV+lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
  as.factor(IVdays_cat)+as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
  lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
  lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
  PI+lag1_PI+lag2_PI+lag3_PI+lag4_PI+
  ntm_combined+lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year

#Use gformula_continuous_eof
gform <- gformula_continuous_eof(obs_data = dat_cc_long,
                                 id = id,
                                 time_name = time_name,
                                 basecovs = basecovs,
                                 covnames = covnames,
                                 outcome_name = outcome_name,
                                 covtypes = covtypes,
                                 covparams = covparams, 
                                 ymodel = ymodel,
                                 intvars = intvars,
                                 interventions = interventions,
                                 int_descript = int_descript,
                                 histories = histories, histvars = histvars,
                                 nsimul = nsimul, 
                                 sim_data_b = TRUE,
                                 seed = 04111990)

#Calculate effect estimates
EY45<-mean(gform$sim_data$'HS & DNase'$Ey[gform$sim_data$'HS & DNase'$fu_year==4], na.rm=TRUE)
EY15<-mean(gform$sim_data$'DNase only'$Ey[gform$sim_data$'DNase only'$fu_year==4], na.rm=TRUE)
CE43.5<-EY45-EY15

EY44<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==4], na.rm=TRUE)
EY14<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==4], na.rm=TRUE)
CE43.4<-EY44-EY14

EY43<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==3], na.rm=TRUE)
EY13<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==3], na.rm=TRUE)
CE43.3<-EY43-EY13

EY42<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==2], na.rm=TRUE)
EY12<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==2], na.rm=TRUE)
CE43.2<-EY42-EY12

EY41<-mean(gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==1], na.rm=TRUE)
EY11<-mean(gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==1], na.rm=TRUE)
CE43.1<-EY41-EY11

#Save estimated outcomes at each time point for MCSE estimates
Y.DNHS.5 <- gform$sim_data$'HS & DNase'$Ey[gform$sim_data$'HS & DNase'$fu_year==4]
Y.DNHS.4 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==4]
Y.DNHS.3 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==3]
Y.DNHS.2 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==2]
Y.DNHS.1 <- gform$sim_data$'HS & DNase'$lag_fev[gform$sim_data$'HS & DNase'$fu_year==1]

Y.DN.5 <- gform$sim_data$'DNase only'$Ey[gform$sim_data$'DNase only'$fu_year==4]
Y.DN.4 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==4]
Y.DN.3 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==3]
Y.DN.2 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==2]
Y.DN.1 <- gform$sim_data$'DNase only'$lag_fev[gform$sim_data$'DNase only'$fu_year==1]

#bootstrap.R runs the above analysis 1000 times in a loop to obtain 
#bootstrap estimates of the SE
source(".../bootstrap.R")

#End timer
proc.time()-start.time

#bootstrap.R results in 5 vectors: ce43.1, ce43.2, ce43.3, ce43.4, ce43.5
#each vector contains 1000 bootstrap estimates of the treatment effects for years 1-5
c(CE43.1, CE43.1-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.1),CE43.1+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.1))
c(CE43.2, CE43.2-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.2),CE43.2+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.2))
c(CE43.3, CE43.3-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.3),CE43.3+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.3))
c(CE43.4, CE43.4-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.4),CE43.4+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.4))
c(CE43.5, CE43.5-qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.5),CE43.5+qnorm(p=.05/2, lower.tail=FALSE)*sd(ce43.5))

#calculate Monte-Carlo SE of difference in means assuming they are two independent samples
sqrt(var(Y.DNHS.5)/nsimul + var(Y.DN.5)/nsimul)
sqrt(var(Y.DNHS.4)/nsimul + var(Y.DN.4)/nsimul)
sqrt(var(Y.DNHS.3)/nsimul + var(Y.DN.3)/nsimul)
sqrt(var(Y.DNHS.2)/nsimul + var(Y.DN.2)/nsimul)
sqrt(var(Y.DNHS.1)/nsimul + var(Y.DN.1)/nsimul)
