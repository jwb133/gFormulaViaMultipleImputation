#------------------------------------------------------------------------------------------------------------
# R code accompanying "G-formula for causal inference via multiple imputation" 
# Author: Emily Granger
# Date last modified: 09/11/2023

#This script includes analysis code for the illustrative example of the use of the G-formula via MI approach.
#The analysis investigates the effects of DNase and hypertonic saline used in combination compared to the use
#of DNase alone on lung function, in people with cystic fibrosis. 

#We use data from the UK CF Registry. Data are available following application to the UK CF Registry Research 
#Committee. 
#https://www.cysticfibrosis.org.uk/the-work-we-do/uk-cf-registry/apply-for-data-from-the-uk-cf-registry.   

#Code is provided for the following analyses: 

# 1a analysis on a complete dataset using the gFormulaMI package
# 1b analysis on a complete dataset using the gfoRmula package
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

#analysis_dat.RData is a workspace that includes two dataframes: dat_wide and dat_long

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
# 1. Analysis on complete data (using an imputed dataset)
#-------------------------------------------------------------------------------
#Describe data that's prepared in MI-gFormulaImpute-200 (workspace should include
# two datasets, one in long format and one in wide)

#Impute data using mice
impdata<-mice(dat_wide, defaultMethod = c("norm.predict", "logreg", "polyreg", 
                                          "polr"), m=1, seed=24042023)

#Extract imputed data set 1 (wide format for gFormulaMI package)
imputed_data<-complete(impdata, action="long", include=FALSE)
imputed_dat1_wide<-imputed_data[imputed_data$.imp==1,]

#Convert imputed data set 1 to long format (for gfoRmula package)
imputed_dat1_long<-reshape(data = imputed_dat1_wide, direction = "long", 
                      varying = 9:53, timevar = "fu_year"
)

#-------------------------------------------------------------------------------
# 1a. Using gFormulaMI package

#Remove id and imputation number from wide dataset
imputed_dat1_wide = subset(imputed_dat1_wide, select = -c(.id, .imp))

#Start timer
start.time<-proc.time()

#Impute synthetic datasets under four treatment regimes
set.seed(04111990)
imps <- gFormulaImpute(data=imputed_dat1_wide, M=200, trtVars=c("trt_grp.1",
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
#Formatting required to prepare dataset for use with gfoRmula

#Convert variables to either factor or numeric
imputed_dat1_long<-imputed_dat1_long %>%
  mutate(
    trt_grp = as.factor(trt_grp),
    dmg_sex = as.numeric(as.character(dmg_sex)),
    genotype_class = as.factor(genotype_class),
    ethnicity = as.numeric(ethnicity),
    Bfev = as.numeric(Bfev),
    Bage = as.numeric(Bage),
    slopes = as.numeric(slopes),
    gli_percpredfev = as.numeric(gli_percpredfev),
    lag_bmi_zscore = as.numeric(lag_bmi_zscore),
    IVdays_cat = as.factor(IVdays_cat),
    hospIV = as.numeric(as.character(hospIV)),
    PI = as.numeric(as.character(PI)),
    s05culturespeciespseudoaeruginos = as.numeric(as.character(s05culturespeciespseudoaeruginos)),
    s05culturespeciesstaph = as.numeric(as.character(s05culturespeciesstaph)),
    ntm_combined = as.numeric(as.character(ntm_combined))
  )

#Order data
imputed_dat1_long<-imputed_dat1_long[order(imputed_dat1_long$.id,imputed_dat1_long$fu_year),]

#Time variable must start from 0 and increase in increments of 1
imputed_dat1_long$fu_year<-imputed_dat1_long$fu_year-1

#Create lag_fev
imputed_dat1_long<- imputed_dat1_long %>% group_by(.id) %>% mutate(lag_fev = dplyr::lag(gli_percpredfev, n = 1, default = NA))
imputed_dat1_long$lag_fev[imputed_dat1_long$fu_year==1]<-imputed_dat1_long$Bfev[imputed_dat1_long$fu_year==1]

#Dataset must be set as data.table
setDT(imputed_dat1_long)

#----------------------
#Analysis

#Start timer
start.time<-proc.time()

#Define inputs for gfoRmula
id<-".id"
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

#Covariate models
covparams<-list(covmodels=c(
  
  trt_grp~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  s05culturespeciespseudoaeruginos~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  s05culturespeciesstaph~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  hospIV~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  IVdays_cat~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  PI~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  ntm_combined~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  lag_fev~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  lag_bmi_zscore~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year
  
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
gform <- gformula_continuous_eof(obs_data = imputed_dat1_long,
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
                                 nsamples=10,
                                 boot_diag=TRUE,
                                 seed = 04111990)

#Display results
gform

#End timer
proc.time()-start.time

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
# 2b. Complete cases using gfoRmula


#------------------------------------------------------------
#Formatting required to prepare dataset for use with gfoRmula
dat_long$fu_year<-dat_long$fu_year-1
setDT(dat_long)

#----------------------
#Analysis

#Start timer
start.time<-proc.time()

#Define inputs for gfoRmula
id<-"patient_id"
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

#Covariate models
covparams<-list(covmodels=c(
  
  trt_grp~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  s05culturespeciespseudoaeruginos~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  s05culturespeciesstaph~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  hospIV~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  IVdays_cat~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  PI~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  ntm_combined~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  lag_fev~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag_bmi_zscore+lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year,
  
  lag_bmi_zscore~as.factor(lag1_trt_grp)+as.factor(lag2_trt_grp)+as.factor(lag3_trt_grp)+as.factor(lag4_trt_grp)+
    as.factor(genotype_class)+as.factor(ethnicity)+Bfev+Bage+dmg_sex+slopes+
    lag1_s05culturespeciespseudoaeruginos+lag2_s05culturespeciespseudoaeruginos+lag3_s05culturespeciespseudoaeruginos+lag4_s05culturespeciespseudoaeruginos+
    lag1_s05culturespeciesstaph+lag2_s05culturespeciesstaph+lag3_s05culturespeciesstaph+lag4_s05culturespeciesstaph+
    lag1_hospIV+lag2_hospIV+lag3_hospIV+lag4_hospIV+
    as.factor(lag1_IVdays_cat)+as.factor(lag2_IVdays_cat)+as.factor(lag3_IVdays_cat)+as.factor(lag4_IVdays_cat)+
    lag_fev+lag2_lag_fev+lag3_lag_fev+lag4_lag_fev+
    lag2_lag_bmi_zscore+lag3_lag_bmi_zscore+lag4_lag_bmi_zscore+
    lag1_PI+lag2_PI+lag3_PI+lag4_PI+
    lag1_ntm_combined+lag2_ntm_combined+lag3_ntm_combined+lag4_ntm_combined+fu_year
  
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
gform <- gformula_continuous_eof(obs_data = dat_long,
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
                                 nsamples=10,
                                 boot_diag=TRUE,
                                 seed = 04111990)

#Display results
gform

#End timer
proc.time()-start.time
