# 2TF-GWAS
It is a script for two steps fitting-GWAS of quantitative traits in R in order to find high-acurated genotype-phenotype associations. First step: quality control of canine visceral leishmaniasis traits using epidemiological traits as covariates in a linear regression model. Second step: residual normalization from the model selected in first step.

setwd("/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test")

install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)

source("boxcox.R")
pheno<-read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/PHENO_DEF2.txt", header=T, sep="\t", stringsAsFactors = F)
#View(pheno)

#Testing traits with enviromental traits and MSL as covariate to look if MSL is associated to traits

#Modelo logístico ajustado para covariáveis ambientais e MSL

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ ORIGIN,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXorigin_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ SEX,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXsex_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ AGE,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXage_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ COLLAR,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXcollar_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ VACCINE,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXvaccine_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ TREATMENT,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXtreatment_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ MSL,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXMSL_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
model<-glm(CLINIC ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXalltraits_residual"] <- model$residuals

barplot(table(pheno$CLINIC))
boxplot(table(pheno$CLINIC))
hist(table(pheno$CLINIC))
rug(pheno$IGA)
summary(pheno$IGA)
model<-glm(CLINIC ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"CLINICXalltraitsMSL_residual"] <- model$residuals
pheno$CLINIC.t<-boxcox.ols(pheno$CLINIC)
hist(pheno$CLINIC.t)
rug(pheno$CLINIC.t)
summary(pheno$CLINIC.t)
model<-lm(CLINIC.t ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
rug(model$residuals)
pheno[names(model$residuals),"CLINIC.tXalltraitsMSL_residual"] <- model$residuals
View(pheno)
View(model$residuals)



#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXorigin_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXsex_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXage_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXcollar_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXvaccine_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXtreatment_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXMSL_residual"] <- model$residuals

barplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXallcovars_residual"] <- model$residuals
pheno$STAGING.t<-boxcox.ols(pheno$STAGING)
hist(pheno$STAGING.t)
rug(pheno$STAGING.t)
summary(pheno$STAGING.t)
model<-lm(STAGING.t ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
rug(model$residuals)
pheno[names(model$residuals),"STAGING.tXalltraits_no_MSL_residual"] <- model$residuals

barplot(table(pheno$STAGING))
boxplot(table(pheno$STAGING))
hist(table(pheno$STAGING))
model<-lm(STAGING ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"STAGINGXallcovarsMSL_residual"] <- model$residuals
pheno$STAGING.t<-boxcox.ols(pheno$STAGING)
hist(pheno$STAGING.t)
rug(pheno$STAGING.t)
summary(pheno$STAGING.t)
model<-lm(STAGING.t ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
rug(model$residuals)
pheno[names(model$residuals),"STAGING.tXalltraits_residual"] <- model$residuals



#Modelo linear ajustado para covariáveis ambientais e MSL

hist(table(pheno$LOG.PL.LFND))
barplot(table(pheno$LOG.PL.LFND))
boxplot(table(pheno$LOG.PL.LFND))
LOG.PL.LFND_out_rm <- pheno$LOG.PL.LFND[!pheno$LOG.PL.LFND %in% boxplot.stats(pheno$LOG.PL.LFND)$out]
length(pheno$LOG.PL.LFND)
length(LOG.PL.LFND_out_rm)
length(pheno$LOG.PL.LFND) - length(LOG.PL.LFND_out_rm)
boxplot(LOG.PL.LFND_out_rm)
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXorigin_residual"] <- model$residuals
pheno$LOG.PL.LFND.t<-boxcox.ols(pheno$LOG.PL.LFND)
hist(pheno$LOG.PL.LFND.t)
rug(pheno$LOG.PL.LFND.t)
summary(pheno$LOG.PL.LFND.t)
model<-lm(LOG.PL.LFND.t ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
rug(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFND.tXorigin_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXsex_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXorigin_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXcollar_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXvaccine_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXtreatment_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXMSL_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXallcovars_residual"] <- model$residuals

barplot(table(pheno$LOG.PL.LFND))
hist(table(pheno$LOG.PL.LFND))
model<-lm(LOG.PL.LFND ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LOG.PL.LFNDXallcovarsMSL_residual"] <- model$residuals

#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXorigin_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXsex_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXage_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXcolar_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXvaccine_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXtreatment_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXMSL_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXallcovars_residual"] <- model$residuals

barplot(table(pheno$IgG.UE))
hist(table(pheno$IgG.UE))
model<-lm(IgG.UE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.UEXallcovarsMSL_residual"] <- model$residuals

#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXorigin_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXsex_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXage_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXcollar_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXvaccine_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXtreatment_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXMSL_residual"] <- model$residuals

barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXallcovars_residual"] <- model$residuals


barplot(table(pheno$IgA))
hist(table(pheno$IgA))
model<-lm(IgA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgAXallcovarsMSL_residual"] <- model$residuals

#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXorigin_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXsex_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXage_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXcollar_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXvaccine_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXtreatment_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXMSL_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXallcovars_residual"] <- model$residuals

barplot(table(pheno$IgM))
hist(table(pheno$IgM))
model<-lm(IgM ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgMXallcovarsMSL_residual"] <- model$residuals

#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXorigin_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXsex_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXage_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXcollar_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXvaccine_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXtreatment_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXMSL_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXallcovars_residual"] <- model$residuals

barplot(table(pheno$IgE))
hist(table(pheno$IgE))
model<-lm(IgE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgEXallcovarsMSL_residual"] <- model$residuals

#Modelo linear ajustado para covariáveis ambientais e MSL

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXorigin_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXsex_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.UE ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXage_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXcollar_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXvaccine_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXtreatment_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXMSL_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXallcovars_residual"] <- model$residuals

barplot(table(pheno$IgG.SALIVA))
hist(table(pheno$IgG.SALIVA))
model<-lm(IgG.SALIVA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IgG.SALIVAXallcovarsMSL_residual"] <- model$residuals

#Modelo linear INDURATION ajustado para covariáveis ambientais

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXorigin_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXsex_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXage_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXcollar_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXvaccine_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXtreatment_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXMSL_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXallcovars_residual"] <- model$residuals

barplot(table(pheno$IDRM.VALUE..mm2.))
hist(table(pheno$IDRM.VALUE..mm2.))
model<-lm(IDRM.VALUE..mm2. ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"INDURATIONXallcovarsMSL_residual"] <- model$residuals

#Modelo logístico para LST ajustado para covariáveis ambientais

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ ORIGIN,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXorigin_residual"] <- model$residuals


barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ SEX,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)pheno[names(model$residuals),"LSTXsex_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ AGE,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXage_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ COLLAR,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXcollar_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ VACCINE,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXvaccine_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ TREATMENT,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXtreatment_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ MSL,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXMSL_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXallcovars_residual"] <- model$residuals

barplot(table(pheno$IDRM.STATUS))
hist(table(pheno$IDRM.STATUS))
model<-glm(IDRM.STATUS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno,family = binomial(link = "logit"))
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"LSTXallcovarsMSL_residual"] <- model$residuals


#Modelo linear IFN.NN ajustado para covariáveis ambientais
barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXorigin_residual"] <- model$residuals


barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXsex_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXage_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXcollar_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXvaccine_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXtreatment_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXMSL_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXallcovars_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$IFN.NN))
model<-lm(IFN.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IFNXallcovarsMSL_residual"] <- model$residuals


#Modelo linear TNF.NN ajustado para covariáveis ambientais

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXorigin_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXsex_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXage_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXcollar_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXvaccine_residual"] <- model$residuals

barplot(table(pheno$IFN.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXtreatment_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXMSL_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXallcovars_residual"] <- model$residuals

barplot(table(pheno$TNF.NN))
hist(table(pheno$TNF.NN))
model<-lm(TNF.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TNFAXallcovarsMSL_residual"] <- model$residuals

#Modelo linear IL10.NN ajustado para covariáveis ambientais

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xorigin_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xsex_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xage_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xcollar_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xvaccine_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xtreatment_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10XMSL_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10Xallcovars_residual"] <- model$residuals

barplot(table(pheno$IL10.NN))
hist(table(pheno$IL10.NN))
model<-lm(IL10.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL10XallcovarsMSL_residual"] <- model$residuals


#Modelo linear IL4.NN ajustado para covariáveis ambientais

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xorigin_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xsex_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xage_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xcollar_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xvaccine_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xtreatment_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4XMSL_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4Xallcovars_residual"] <- model$residuals

barplot(table(pheno$IL4.NN))
hist(table(pheno$IL4.NN))
model<-lm(IL4.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"IL4XallcovarsMSL_residual"] <- model$residuals


#Modelo linear TGFB ajustado para covariáveis ambientais

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXorigin_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUE_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXage_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUE_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXvaccine_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXtreatment_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXMSL_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXallcovars_residual"] <- model$residuals

barplot(table(pheno$TGF.VALUE))
hist(table(pheno$TGF.VALUE))
model<-lm(TGF.VALUE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TGF.VALUEXallcovarsMSL_residual"] <- model$residuals

#Modelo linear NO ajustado para covariáveis ambientais

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXorigin_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXsex_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXage_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXcollar_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXvaccine_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXtreatment_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXMSL_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXallcovars_residual"] <- model$residuals

barplot(table(pheno$NO.NN))
hist(table(pheno$NO.NN))
model<-lm(NO.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"NOXallcovarsMSL_residual"] <- model$residuals

#Modelo linear SOD ajustado para covariáveis ambientais

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXorigin_residual"] <- model$residuals


barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXsex_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXage_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXcollar_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXvaccine_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXtreatment_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXMSL_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXallcovars_residual"] <- model$residuals

barplot(table(pheno$SOD.NN))
hist(table(pheno$SOD.NN))
model<-lm(SOD.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"SODXallcovars_residual"] <- model$residuals

#Modelo linear TOS ajustado para covariáveis ambientais

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXorigin_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXsex_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXage_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXcollar_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXvaccine_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXtreatment_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXMSL_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXallcovars_residual"] <- model$residuals

barplot(table(pheno$TOS))
hist(table(pheno$TOS))
model<-lm(TOS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TOCXallcovarsMSL_residual"] <- model$residuals


#Modelo linear TAS ajustado para covariáveis ambientais

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXorigin_residual"] <- model$residuals


barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXsex_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXage_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXcollar_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXvaccine_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXorigin_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXMSL_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXallcovars_residual"] <- model$residuals

barplot(table(pheno$TAS))
hist(table(pheno$TAS))
model<-lm(TAS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"TACXallcovarsMSL_residual"] <- model$residuals


#Criando um subset MDA para normalização
MDA <- select(pheno,IID,FID,MDA)
#View(MDA)

#Remoção dos outliers
hist(MDA$MDA)
boxplot(MDA$MDA)
length(MDA$MDA)
MDA$MDA[MDA$MDA %in% boxplot.stats(MDA$MDA)$out]
length(MDA$MDA[MDA$MDA %in% boxplot.stats(MDA$MDA)$out])
MDA_out_rm <- subset(MDA, MDA < 83.13)
nrow(MDA_out_rm)
#View(MDA_out_rm)
MDA_out_rm$out_rm <- MDA_out_rm$MDA
#View(MDA_out_rm)
MDA_out_rm <- within(MDA_out_rm, rm(MDA))
#View(MDAXallcovars_no_MSL_res_out_rm)
length(MDA$MDA)
length(MDA_out_rm$out_rm)
length(MDA$MDA) - length(MDA_out_rm$out_rm)
boxplot(MDA_out_rm$out_rm)
hist(MDA_out_rm$out_rm)
rug(MDA_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
MDA <- full_join(MDA,MDA_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(MDA)
#View(MDAXallcovars_no_MSL_res)

#Transformação boxcox
MDA$MDA.t<-boxcox.ols(MDA$MDA)
boxplot(MDA$MDA.t)
hist(MDA$MDA.t)
rug(MDA$MDA.t)
summary(MDA$MDA.t)

#Transformação logarítmica
MDA$LOG <- log10(MDA$MDA)
boxplot(MDA$LOG)
hist(MDA$LOG)
rug(MDA$LOG)
summary(MDA$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$MDA)
hist(pheno$MDA)
model<-lm(MDA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(MDA$MDA)
hist(model$residuals)
rug(models$residuals)
MDA[names(model$residuals),"MDAXallcovars_no_MSL_RES"] <- model$residuals
#View(MDA)

#Remoção dos outliers
hist(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES)
MDA$MDAXallcovars_no_MSL_RES[MDA$MDAXallcovars_no_MSL_RES %in% boxplot.stats(MDA$MDAXallcovars_no_MSL_RES)$out]
length(MDA$MDAXallcovars_no_MSL_RES[MDA$MDAXallcovars_no_MSL_RES %in% boxplot.stats(MDA$MDAXallcovars_no_MSL_RES)$out])
MDAXalcovars_no_MSL_out_rm <- subset(MDA, MDAXallcovars_no_MSL_RES < 57.41991)
nrow(MDAXalcovars_no_MSL_out_rm)
MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm <- MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES
MDAXalcovars_no_MSL_out_rm <- within(MDAXalcovars_no_MSL_out_rm, rm(MDAXallcovars_no_MSL_RES))
length(MDA$MDAXallcovars_no_MSL_RES)
length(MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
length(MDA$MDAXallcovars_no_MSL_RES - MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
boxplot(MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
hist(MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
rug(MDAXalcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
MDA <- full_join(MDA,MDAXalcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(MDA)
View(MDA)


#Transformação boxcox
MDA$MDAXallcovars_no_MSL_RES.t<-boxcox.ols(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES.t)
hist(MDA$MDAXallcovars_no_MSL_RES.t)
rug(MDA$MDAXallcovars_no_MSL_RES.t)
summary(MDA$MDAXallcovars_no_MSL_RES.t)

#Transformação logarítmica
MDA$MDAXallcovars_no_MSL_RES.LOG <- log10(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES.LOG)
hist(MDA$MDAXallcovars_no_MSL_RES.LOG)
rug(MDA$MDAXallcovars_no_MSL_RES.LOG)
summary(MDA$MDAXallcovars_no_MSL_RES.LOG)

View(MDA)

write.table(MDA,"MDA_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


model<-lm(MDA ~ ORIGIN,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXorigin_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ SEX,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXsex_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ AGE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXage_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ COLLAR,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXcollar_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ VACCINE,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXvaccine_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXtreatment_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXMSL_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXallcovars_residual"] <- model$residuals

barplot(table(pheno$MDA))
hist(table(pheno$MDA))
model<-lm(MDA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT + MSL,data = pheno)
summary(model)
hist(model$residuals)
pheno[names(model$residuals),"MDAXallcovarsMSL_residual"] <- model$residuals



#Normalização do resíduo do trait + Kinship + covariáveis

clinic_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/CLINIC.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(clinic_res$Residual)
boxplot(clinic_res$Residual)
clinic_res$CLINIC_res_out_rm <- clinic_res$Residual[!clinic_res$Residual %in% boxplot.stats(clinic_res$Residual)$out]
length(clinic_res$Residual)
length(clinic_res$CLINIC_res_out_rm)
length(clinic_res$Residual) - length(clinic_res$CLINIC_res_out_rm)
boxplot(clinic_res$CLINIC_res_out_rm)
hist(clinic_res$CLINIC_res_out_rm)
rug(clinic_res$CLINIC_res_out_rm)

#Transformação boxcox
clinic_res$Residual.t<-boxcox.ols(clinic_res$Residual)
boxplot(clinic_res$Residual.t)
hist(clinic_res$Residual.t)
rug(clinic_res$Residual.t)
summary(clinic_res$Residual.t)

#Transformação logarítmica
clinic_res$LOG <- log10(clinic_res$Residual)
boxplot(clinic_res$LOG)
hist(clinic_res$LOG)
rug(clinic_res$LOG)
summary(clinic_res$LOG)

#View(clinic_res)

write.table(clinic_res,"CLINIC_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

clinicXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/CLINIC.allcovars_no_MSL_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(clinicXallcovars_no_MSL_res$Residual)
boxplot(clinicXallcovars_no_MSL_res$Residual)
clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm <- clinicXallcovars_no_MSL_res$Residual[!clinicXallcovars_no_MSL_res$Residual %in% boxplot.stats(clinicXallcovars_no_MSL_res$Residual)$out]
length(clinicXallcovars_no_MSL_res$Residual)
length(clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm)
length(clinicXallcovars_no_MSL_res$Residual) - length(clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm)
boxplot(clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm)
hist(clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm)
rug(clinicXallcovars_no_MSL_res$CLINIC_Xallcovars_no_MSL_res_out_rm)

#Transformação boxcox
clinicXallcovars_no_MSL_res$Residual.t<-boxcox.ols(clinicXallcovars_no_MSL_res$Residual)
boxplot(clinicXallcovars_no_MSL_res$Residual.t)
hist(clinicXallcovars_no_MSL_res$Residual.t)
rug(clinicXallcovars_no_MSL_res$Residual.t)
summary(clinicXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
clinicXallcovars_no_MSL_res$LOG <- log10(clinicXallcovars_no_MSL_res$Residual)
boxplot(clinicXallcovars_no_MSL_res$LOG)
hist(clinicXallcovars_no_MSL_res$LOG)
rug(clinicXallcovars_no_MSL_res$LOG)
summary(clinicXallcovars_no_MSL_res$LOG)

#View(clinicXallcovars_no_MSL_res)

write.table(clinic_res,"CLINIC_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


clinicXallcovars_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/CLINIC.all_covariates_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(clinicXallcovars_res$Residual)
boxplot(clinicXallcovars_res$Residual)
clinicXallcovars_res$CLINICXallcovars_res_out_rm <- clinicXallcovars_res$Residual[!clinicXallcovars_res$Residual %in% boxplot.stats(clinicXallcovars_res$Residual)$out]
length(clinicXallcovars_res$Residual)
length(clinicXallcovars_res$CLINICXallcovars_res_out_rm)
length(clinicXallcovars_res$Residual) - length(clinicXallcovars_res$CLINICXallcovars_res_out_rm)
boxplot(clinicXallcovars_res$CLINICXallcovars_res_out_rm)
hist(clinicXallcovars_res$CLINICXallcovars_res_out_rm)
rug(clinicXallcovars_res$CLINICXallcovars_res_out_rm)

#Transformação boxcox
clinicXallcovars_res$Residual.t<-boxcox.ols(clinicXallcovars_res$Residual)
boxplot(clinicXallcovars_res$Residual.t)
hist(clinicXallcovars_res$Residual.t)
rug(clinicXallcovars_res$Residual.t)
summary(clinicXallcovars_res$Residual.t)

#Transformação logarítmica
clinicXallcovars_res$LOG <- log10(clinicXallcovars_res$Residual)
boxplot(clinicXallcovars_res$LOG)
hist(clinicXallcovars_res$LOG)
rug(clinicXallcovars_res$LOG)
summary(clinicXallcovars_res$LOG)

#View(clinicXallcovars_res)

write.table(clinic_res,"CLINIC_ALLCOVARS_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


stagingXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/STAGING.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(stagingXallcovars_no_MSL_res$Residual)
boxplot(stagingXallcovars_no_MSL_res$Residual)
stagingXallcovars_no_MSL_res$Residual[stagingXallcovars_no_MSL_res$Residual %in% boxplot.stats(stagingXallcovars_no_MSL_res$Residual)$out]
length(stagingXallcovars_no_MSL_res$Residual[stagingXallcovars_no_MSL_res$Residual %in% boxplot.stats(stagingXallcovars_no_MSL_res$Residual)$out])
StagingXallcovars_no_MSL_res_out_rm <- subset(stagingXallcovars_no_MSL_res, Residual > -1.483464)
nrow(StagingXallcovars_no_MSL_res_out_rm)
#View(StagingXallcovars_no_MSL_res_out_rm)
StagingXallcovars_no_MSL_res_out_rm$res_out_rm <- StagingXallcovars_no_MSL_res_out_rm$Residual
#View(StagingXallcovars_no_MSL_res_out_rm)
StagingXallcovars_no_MSL_res_out_rm <- within(StagingXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(StagingXallcovars_no_MSL_res_out_rm)
length(stagingXallcovars_no_MSL_res$Residual)
length(StagingXallcovars_no_MSL_res_out_rm$res_out_rm)
length(stagingXallcovars_no_MSL_res$Residual) - length(StagingXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(StagingXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(StagingXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(StagingXallcovars_no_MSL_res_out_rm$res_out_rm)

stagingXallcovars_no_MSL_res <- full_join(stagingXallcovars_no_MSL_res,StagingXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(stagingXallcovars_no_MSL_res)
#View(stagingXallcovars_no_MSL_res)

#Transformação boxcox
stagingXallcovars_no_MSL_res$Residual.t<-boxcox.ols(stagingXallcovars_no_MSL_res$Residual)
boxplot(stagingXallcovars_no_MSL_res$Residual.t)
hist(stagingXallcovars_no_MSL_res$Residual.t)
rug(stagingXallcovars_no_MSL_res$Residual.t)
summary(stagingXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
stagingXallcovars_no_MSL_res$LOG <- log10(stagingXallcovars_no_MSL_res$Residual)
boxplot(stagingXallcovars_no_MSL_res$LOG)
hist(stagingXallcovars_no_MSL_res$LOG)
rug(stagingXallcovars_no_MSL_res$LOG)
summary(stagingXallcovars_no_MSL_res$LOG)

#View(stagingXallcovars_no_MSL_res)

write.table(stagingXallcovars_no_MSL_res,"STAGING_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

pl.lfndXorigin_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/LOG.PL.LFND.ORIGIN.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(pl.lfndXorigin_res$Residual)
boxplot(pl.lfndXorigin_res$Residual)
length(pl.lfndXorigin_res$Residual)
pl.lfndXorigin_res$Residual[pl.lfndXorigin_res$Residual %in% boxplot.stats(pl.lfndXorigin_res$Residual)$out]
length(pl.lfndXorigin_res$Residual[pl.lfndXorigin_res$Residual %in% boxplot.stats(pl.lfndXorigin_res$Residual)$out])
pl.lfndXorigin_res_out_rm <- subset(pl.lfndXorigin_res, Residual > -1.653004 & Residual < 4.751120)
nrow(pl.lfndXorigin_res_out_rm)
#View(pl.lfndXorigin_res_out_rm)
pl.lfndXorigin_res_out_rm$res_out_rm <- pl.lfndXorigin_res_out_rm$Residual
#View(pl.lfndXorigin_res_out_rm)
pl.lfndXorigin_res_out_rm <- within(pl.lfndXorigin_res_out_rm, rm(Residual))
#View(pl.lfndXorigin_res_out_rm)
length(pl.lfndXorigin_res$Residual)
length(pl.lfndXorigin_res_out_rm$res_out_rm)
length(pl.lfndXorigin_res$Residual) - length(pl.lfndXorigin_res_out_rm$res_out_rm)
boxplot(pl.lfndXorigin_res_out_rm$res_out_rm)
hist(pl.lfndXorigin_res_out_rm$res_out_rm)
rug(pl.lfndXorigin_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
pl.lfndXorigin_res <- full_join(pl.lfndXorigin_res,pl.lfndXorigin_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(pl.lfndXorigin_res)
#View(pl.lfndXorigin_res)

#Transformação boxcox
pl.lfndXorigin_res$Residual.t<-boxcox.ols(pl.lfndXorigin_res$Residual)
boxplot(pl.lfndXorigin_res$Residual.t)
hist(pl.lfndXorigin_res$Residual.t)
rug(pl.lfndXorigin_res$Residual.t)
summary(pl.lfndXorigin_res$Residual.t)

#Transformação logarítmica
pl.lfndXorigin_res$LOG <- log10(pl.lfndXorigin_res$Residual)
boxplot(pl.lfndXorigin_res$LOG)
hist(pl.lfndXorigin_res$LOG)
rug(pl.lfndXorigin_res$LOG)
summary(pl.lfndXorigin_res$LOG)

#View(pl.lfndXorigin_res)

write.table(pl.lfndXorigin_res,"PL.LFND_ORIGIN_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


pl.bloodXMSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/PL.BLOOD.MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(pl.bloodXMSL_res$Residual)
boxplot(pl.bloodXMSL_res$Residual)
length(pl.bloodXMSL_res$Residual)
pl.bloodXMSL_res$Residual[pl.bloodXMSL_res$Residual %in% boxplot.stats(pl.bloodXMSL_res$Residual)$out]
length(pl.bloodXMSL_res$Residual[pl.bloodXMSL_res$Residual %in% boxplot.stats(pl.bloodXMSL_res$Residual)$out])
#order(pl.bloodXMSL_res$Residual[pl.bloodXMSL_res$Residual %in% boxplot.stats(pl.bloodXMSL_res$Residual)$out], decreasing = TRUE, na.last = TRUE)
pl.bloodXMSL_res_out_rm <- subset(pl.bloodXMSL_res, Residual > -30.85968 & Residual < 37.32495)
nrow(pl.bloodXMSL_res_out_rm)
#View(pl.bloodXMSL_res_out_rm)
pl.bloodXMSL_res_out_rm$res_out_rm <- pl.bloodXMSL_res_out_rm$Residual
#View(pl.bloodXMSL_res_out_rm)
pl.bloodXMSL_res_out_rm <- within(pl.bloodXMSL_res_out_rm, rm(Residual))
#View(pl.bloodXMSL_res_out_rm)
length(pl.bloodXMSL_res$Residual)
length(pl.bloodXMSL_res_out_rm$res_out_rm)
length(pl.bloodXMSL_res$Residual) - length(pl.bloodXMSL_res_out_rm$res_out_rm)
boxplot(pl.bloodXMSL_res_out_rm$res_out_rm)
hist(pl.bloodXMSL_res_out_rm$res_out_rm)
rug(pl.bloodXMSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
pl.bloodXMSL_res <- full_join(pl.bloodXMSL_res,pl.bloodXMSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(pl.bloodXMSL_res)
#View(pl.bloodXMSL_res)

#Transformação boxcox
pl.bloodXMSL_res$Residual.t<-boxcox.ols(pl.bloodXMSL_res$Residual)
boxplot(pl.bloodXMSL_res$Residual.t)
hist(pl.bloodXMSL_res$Residual.t)
rug(pl.bloodXMSL_res$Residual.t)
summary(pl.bloodXMSL_res$Residual.t)

#Transformação logarítmica
pl.bloodXMSL_res$LOG <- log10(pl.bloodXMSL_res$Residual)
boxplot(pl.bloodXMSL_res$LOG)
hist(pl.bloodXMSL_res$LOG)
rug(pl.bloodXMSL_res$LOG)
summary(pl.bloodXMSL_res$LOG)

#View(pl.bloodXMSL_res)

write.table(pl.bloodXMSL_res,"PL.BLOOD_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

IgGXallcovars_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IgG.UE.all_covariates.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IgGXallcovars_res$Residual)
boxplot(IgGXallcovars_res$Residual)
length(IgGXallcovars_res$Residual)
IgGXallcovars_res$Residual[IgGXallcovars_res$Residual %in% boxplot.stats(IgGXallcovars_res$Residual)$out]
length(IgGXallcovars_res$Residual[IgGXallcovars_res$Residual %in% boxplot.stats(IgGXallcovars_res$Residual)$out])
IgGXallcovars_res_out_rm <- subset(IgGXallcovars_res, Residual > -5.897590 & Residual < 4.922174)
nrow(IgGXallcovars_res_out_rm)
#View(IgGXallcovars_res_out_rm)
IgGXallcovars_res_out_rm$res_out_rm <- IgGXallcovars_res_out_rm$Residual
#View(IgGXallcovars_res_out_rm)
IgGXallcovars_res_out_rm <- within(IgGXallcovars_res_out_rm, rm(Residual))
#View(IgGXallcovars_res_out_rm)
length(IgGXallcovars_res$Residual)
length(IgGXallcovars_res_out_rm$res_out_rm)
length(IgGXallcovars_res$Residual) - length(IgGXallcovars_res_out_rm$res_out_rm)
boxplot(IgGXallcovars_res_out_rm$res_out_rm)
hist(IgGXallcovars_res_out_rm$res_out_rm)
rug(IgGXallcovars_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgGXallcovars_res <- full_join(IgGXallcovars_res,IgGXallcovars_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgGXallcovars_res)
#View(IgGXallcovars_res)

#Transformação boxcox
IgGXallcovars_res$Residual.t<-boxcox.ols(IgGXallcovars_res$Residual)
boxplot(IgGXallcovars_res$Residual.t)
hist(IgGXallcovars_res$Residual.t)
rug(IgGXallcovars_res$Residual.t)
summary(IgGXallcovars_res$Residual.t)

#Transformação logarítmica
IgGXallcovars_res$LOG <- log10(IgGXallcovars_res$Residual)
boxplot(IgGXallcovars_res$LOG)
hist(IgGXallcovars_res$LOG)
rug(IgGXallcovars_res$LOG)
summary(IgGXallcovars_res$LOG)

#View(IgGXallcovars_res)

write.table(IgGXallcovars_res,"IgG_ALLCOVARS_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


IgAXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IgA.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IgAXallcovars_no_MSL_res$Residual)
boxplot(IgAXallcovars_no_MSL_res$Residual)
length(IgAXallcovars_no_MSL_res$Residual)
IgAXallcovars_no_MSL_res$Residual[IgAXallcovars_no_MSL_res$Residual %in% boxplot.stats(IgAXallcovars_no_MSL_res$Residual)$out]
length(IgAXallcovars_no_MSL_res$Residual[IgAXallcovars_no_MSL_res$Residual %in% boxplot.stats(IgAXallcovars_no_MSL_res$Residual)$out])
IgAXallcovars_no_MSL_res_out_rm <- subset(IgAXallcovars_no_MSL_res, Residual > -13.81384 & Residual < 10.74568)
nrow(IgAXallcovars_no_MSL_res_out_rm)
#View(IgAXallcovars_no_MSL_res_out_rm)
IgAXallcovars_no_MSL_res_out_rm$res_out_rm <- IgAXallcovars_no_MSL_res_out_rm$Residual
#View(IgAXallcovars_no_MSL_res_out_rm)
IgAXallcovars_no_MSL_res_out_rm <- within(IgAXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(IgAXallcovars_no_MSL_res_out_rm)
length(IgAXallcovars_no_MSL_res$Residual)
length(IgAXallcovars_no_MSL_res_out_rm$res_out_rm)
length(IgAXallcovars_no_MSL_res$Residual) - length(IgAXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(IgAXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(IgAXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(IgAXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgAXallcovars_no_MSL_res <- full_join(IgAXallcovars_no_MSL_res,IgAXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgAXallcovars_no_MSL_res)
#View(IgAXallcovars_no_MSL_res)

#Transformação boxcox
IgAXallcovars_no_MSL_res$Residual.t<-boxcox.ols(IgAXallcovars_no_MSL_res$Residual)
boxplot(IgAXallcovars_no_MSL_res$Residual.t)
hist(IgAXallcovars_no_MSL_res$Residual.t)
rug(IgAXallcovars_no_MSL_res$Residual.t)
summary(IgAXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
IgAXallcovars_no_MSL_res$LOG <- log10(IgAXallcovars_no_MSL_res$Residual)
boxplot(IgAXallcovars_no_MSL_res$LOG)
hist(IgAXallcovars_no_MSL_res$LOG)
rug(IgAXallcovars_no_MSL_res$LOG)
summary(IgAXallcovars_no_MSL_res$LOG)

#View(IgAXallcovars_no_MSL_res)

write.table(IgAXallcovars_no_MSL_res,"IgA_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


IgEXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IgE.allcovars_no_MS.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IgEXallcovars_no_MSL_res$Residual)
boxplot(IgEXallcovars_no_MSL_res$Residual)
length(IgEXallcovars_no_MSL_res$Residual)
IgEXallcovars_no_MSL_res$Residual[IgEXallcovars_no_MSL_res$Residual %in% boxplot.stats(IgEXallcovars_no_MSL_res$Residual)$out]
length(IgEXallcovars_no_MSL_res$Residual[IgEXallcovars_no_MSL_res$Residual %in% boxplot.stats(IgEXallcovars_no_MSL_res$Residual)$out])
IgEXallcovars_no_MSL_res_out_rm <- subset(IgEXallcovars_no_MSL_res, Residual < 54.17923)
nrow(IgEXallcovars_no_MSL_res_out_rm)
#View(IgEXallcovars_no_MSL_res_out_rm)
IgEXallcovars_no_MSL_res_out_rm$res_out_rm <- IgEXallcovars_no_MSL_res_out_rm$Residual
#View(IgEXallcovars_no_MSL_res_out_rm)
IgEXallcovars_no_MSL_res_out_rm <- within(IgEXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(IgEXallcovars_no_MSL_res_out_rm)
length(IgEXallcovars_no_MSL_res$Residual)
length(IgEXallcovars_no_MSL_res_out_rm$res_out_rm)
length(IgEXallcovars_no_MSL_res$Residual) - length(IgEXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(IgEXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(IgEXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(IgEXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgEXallcovars_no_MSL_res <- full_join(IgEXallcovars_no_MSL_res,IgEXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgEXallcovars_no_MSL_res)
#View(IgEXallcovars_no_MSL_res)

#Transformação boxcox
IgEXallcovars_no_MSL_res$Residual.t<-boxcox.ols(IgEXallcovars_no_MSL_res$Residual)
boxplot(IgEXallcovars_no_MSL_res$Residual.t)
hist(IgEXallcovars_no_MSL_res$Residual.t)
rug(IgEXallcovars_no_MSL_res$Residual.t)
summary(IgEXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
IgEXallcovars_no_MSL_res$LOG <- log10(IgEXallcovars_no_MSL_res$Residual)
boxplot(IgEXallcovars_no_MSL_res$LOG)
hist(IgEXallcovars_no_MSL_res$LOG)
rug(IgEXallcovars_no_MSL_res$LOG)
summary(IgEXallcovars_no_MSL_res$LOG)

#View(IgEXallcovars_no_MSL_res)

write.table(IgEXallcovars_no_MSL_res,"IgE_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


IgMXcollar_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IgM.COLLAR.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IgMXcollar_res$Residual)
boxplot(IgMXcollar_res$Residual)
length(IgMXcollar_res$Residual)
IgMXcollar_res$Residual[IgMXcollar_res$Residual %in% boxplot.stats(IgMXcollar_res$Residual)$out]
length(IgMXcollar_res$Residual[IgMXcollar_res$Residual %in% boxplot.stats(IgMXcollar_res$Residual)$out])
IgMXcollar_res_out_rm <- subset(IgMXcollar_res, Residual < 85.11396)
nrow(IgMXcollar_res_out_rm)
#View(IgMXcollar_res_out_rm)
IgMXcollar_res_out_rm$res_out_rm <- IgMXcollar_res_out_rm$Residual
#View(IgMXcollar_res_out_rm)
IgMXcollar_res_out_rm <- within(IgMXcollar_res_out_rm, rm(Residual))
#View(IgMXcollar_res_out_rm)
length(IgMXcollar_res$Residual)
length(IgMXcollar_res_out_rm$res_out_rm)
length(IgMXcollar_res$Residual) - length(IgMXcollar_res_out_rm$res_out_rm)
boxplot(IgMXcollar_res_out_rm$res_out_rm)
hist(IgMXcollar_res_out_rm$res_out_rm)
rug(IgMXcollar_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgMXcollar_res <- full_join(IgMXcollar_res,IgMXcollar_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgMXcollar_res)
#View(IgMXcollar_res)

#Transformação boxcox
IgMXcollar_res$Residual.t<-boxcox.ols(IgMXcollar_res$Residual)
boxplot(IgMXcollar_res$Residual.t)
hist(IgMXcollar_res$Residual.t)
rug(IgMXcollar_res$Residual.t)
summary(IgMXcollar_res$Residual.t)

#Transformação logarítmica
IgMXcollar_res$LOG <- log10(IgMXcollar_res$Residual)
boxplot(IgMXcollar_res$LOG)
hist(IgMXcollar_res$LOG)
rug(IgMXcollar_res$LOG)
summary(IgMXcollar_res$LOG)

#View(IgMXcollar_res)

write.table(IgMXcollar_res,"IgM_COLLAR_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


IgG.salivaXtreatment_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IgG.SALIVA.TREATMENT.res_col.txt", header=T, sep="\t", stringsAsFactors = F) 

#Remoção dos outliers
hist(IgG.salivaXtreatment_res$Residual)
boxplot(IgG.salivaXtreatment_res$Residual)
length(IgG.salivaXtreatment_res$Residual)
IgG.salivaXtreatment_res$Residual[IgG.salivaXtreatment_res$Residual %in% boxplot.stats(IgG.salivaXtreatment_res$Residual)$out]
length(IgG.salivaXtreatment_res$Residual[IgG.salivaXtreatment_res$Residual %in% boxplot.stats(IgG.salivaXtreatment_res$Residual)$out])
IgG.salivaXtreatment_res_out_rm <- subset(IgG.salivaXtreatment_res, Residual < 39.64364)
nrow(IgG.salivaXtreatment_res_out_rm)
#View(IgG.salivaXtreatment_res_out_rm)
IgG.salivaXtreatment_res_out_rm$res_out_rm <- IgG.salivaXtreatment_res_out_rm$Residual
#View(IgG.salivaXtreatment_res_out_rm)
IgG.salivaXtreatment_res_out_rm <- within(IgG.salivaXtreatment_res_out_rm, rm(Residual))
#View(IgG.salivaXtreatment_res_out_rm)
length(IgG.salivaXtreatment_res$Residual)
length(IgG.salivaXtreatment_res_out_rm$res_out_rm)
length(IgG.salivaXtreatment_res$Residual) - length(IgG.salivaXtreatment_res_out_rm$res_out_rm)
boxplot(IgG.salivaXtreatment_res_out_rm$res_out_rm)
hist(IgG.salivaXtreatment_res_out_rm$res_out_rm)
rug(IgG.salivaXtreatment_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgG.salivaXtreatment_res <- full_join(IgG.salivaXtreatment_res,IgG.salivaXtreatment_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgG.salivaXtreatment_res)
#View(IgG.salivaXtreatment_res)

#Transformação boxcox
IgG.salivaXtreatment_res$Residual.t<-boxcox.ols(IgG.salivaXtreatment_res$Residual)
boxplot(IgG.salivaXtreatment_res$Residual.t)
hist(IgG.salivaXtreatment_res$Residual.t)
rug(IgG.salivaXtreatment_res$Residual.t)
summary(IgG.salivaXtreatment_res$Residual.t)

#Transformação logarítmica
IgG.salivaXtreatment_res$LOG <- log10(IgG.salivaXtreatment_res$Residual)
boxplot(IgG.salivaXtreatment_res$LOG)
hist(IgG.salivaXtreatment_res$LOG)
rug(IgG.salivaXtreatment_res$LOG)
summary(IgG.salivaXtreatment_res$LOG)

#View(IgG.salivaXtreatment_res)

write.table(IgG.salivaXtreatment_res,"IgG.SALIVA_TREATMENT_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


LSTXorigin_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/LST.ORIGIN.res_col.txt", header=T, sep="\t", stringsAsFactors = F) 

#Remoção dos outliers
hist(LSTXorigin_res$Residual)
boxplot(LSTXorigin_res$Residual)
length(LSTXorigin_res$Residual)
LSTXorigin_res$Residual[LSTXorigin_res$Residual %in% boxplot.stats(LSTXorigin_res$Residual)$out]
length(LSTXorigin_res$Residual[LSTXorigin_res$Residual %in% boxplot.stats(LSTXorigin_res$Residual)$out])
LSTXorigin_res_out_rm <- subset(LSTXorigin_res, Residual < 0.653855)
nrow(LSTXorigin_res_out_rm)
#View(LSTXorigin_res_out_rm)
LSTXorigin_res_out_rm$res_out_rm <- LSTXorigin_res_out_rm$Residual
#View(LSTXorigin_res_out_rm)
LSTXorigin_res_out_rm <- within(LSTXorigin_res_out_rm, rm(Residual))
#View(LSTXorigin_res_out_rm)
length(LSTXorigin_res$Residual)
length(LSTXorigin_res_out_rm$res_out_rm)
length(LSTXorigin_res$Residual) - length(LSTXorigin_res_out_rm$res_out_rm)
boxplot(LSTXorigin_res_out_rm$res_out_rm)
hist(LSTXorigin_res_out_rm$res_out_rm)
rug(LSTXorigin_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
LSTXorigin_res <- full_join(LSTXorigin_res,LSTXorigin_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(LSTXorigin_res)
#View(LSTXorigin_res)

#Transformação boxcox
LSTXorigin_res$Residual.t<-boxcox.ols(LSTXorigin_res$Residual)
boxplot(LSTXorigin_res$Residual.t)
hist(LSTXorigin_res$Residual.t)
rug(LSTXorigin_res$Residual.t)
summary(LSTXorigin_res$Residual.t)

#Transformação logarítmica
LSTXorigin_res$LOG <- log10(LSTXorigin_res$Residual)
boxplot(LSTXorigin_res$LOG)
hist(LSTXorigin_res$LOG)
rug(LSTXorigin_res$LOG)
summary(LSTXorigin_res$LOG)

#View(LSTXorigin_res)

write.table(LSTXorigin_res,"LST_ORIGIN_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


indurationXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/INDURATION.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F) 

#Remoção dos outliers
hist(indurationXallcovars_no_MSL_res$Residual)
boxplot(indurationXallcovars_no_MSL_res$Residual)
length(indurationXallcovars_no_MSL_res$Residual)
indurationXallcovars_no_MSL_res$Residual[indurationXallcovars_no_MSL_res$Residual %in% boxplot.stats(indurationXallcovars_no_MSL_res$Residual)$out]
length(indurationXallcovars_no_MSL_res$Residual[indurationXallcovars_no_MSL_res$Residual %in% boxplot.stats(indurationXallcovars_no_MSL_res$Residual)$out])
indurationXallcovars_no_MSL_res_out_rm <- subset(indurationXallcovars_no_MSL_res, Residual < 72.50773)
nrow(indurationXallcovars_no_MSL_res_out_rm)
#View(indurationXallcovars_no_MSL_res_out_rm)
indurationXallcovars_no_MSL_res_out_rm$res_out_rm <- indurationXallcovars_no_MSL_res_out_rm$Residual
#View(indurationXallcovars_no_MSL_res_out_rm)
indurationXallcovars_no_MSL_res_out_rm <- within(indurationXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(indurationXallcovars_no_MSL_res_out_rm)
length(indurationXallcovars_no_MSL_res$Residual)
length(indurationXallcovars_no_MSL_res_out_rm$res_out_rm)
length(indurationXallcovars_no_MSL_res$Residual) - length(indurationXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(indurationXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(indurationXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(indurationXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
indurationXallcovars_no_MSL_res <- full_join(indurationXallcovars_no_MSL_res,indurationXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(indurationXallcovars_no_MSL_res)
#View(indurationXallcovars_no_MSL_res)

#Transformação boxcox
indurationXallcovars_no_MSL_res$Residual.t<-boxcox.ols(indurationXallcovars_no_MSL_res$Residual)
boxplot(indurationXallcovars_no_MSL_res$Residual.t)
hist(indurationXallcovars_no_MSL_res$Residual.t)
rug(indurationXallcovars_no_MSL_res$Residual.t)
summary(indurationXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
indurationXallcovars_no_MSL_res$LOG <- log10(indurationXallcovars_no_MSL_res$Residual)
boxplot(indurationXallcovars_no_MSL_res$LOG)
hist(indurationXallcovars_no_MSL_res$LOG)
rug(indurationXallcovars_no_MSL_res$LOG)
summary(indurationXallcovars_no_MSL_res$LOG)

#View(indurationXallcovars_no_MSL_res)

write.table(indurationXallcovars_no_MSL_res,"INDURATION_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

PIXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/PI.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F) 

#Remoção dos outliers
hist(PIXallcovars_no_MSL_res$Residual)
boxplot(PIXallcovars_no_MSL_res$Residual)
length(PIXallcovars_no_MSL_res$Residual)
PIXallcovars_no_MSL_res$Residual[PIXallcovars_no_MSL_res$Residual %in% boxplot.stats(PIXallcovars_no_MSL_res$Residual)$out]
length(PIXallcovars_no_MSL_res$Residual[PIXallcovars_no_MSL_res$Residual %in% boxplot.stats(PIXallcovars_no_MSL_res$Residual)$out])
PIXallcovars_no_MSL_res_out_rm <- subset(PIXallcovars_no_MSL_res, Residual < 84.09875)
nrow(PIXallcovars_no_MSL_res_out_rm)
#View(PIXallcovars_no_MSL_res_out_rm)
PIXallcovars_no_MSL_res_out_rm$res_out_rm <- PIXallcovars_no_MSL_res_out_rm$Residual
#View(PIXallcovars_no_MSL_res_out_rm)
PIXallcovars_no_MSL_res_out_rm <- within(PIXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(PIXallcovars_no_MSL_res_out_rm)
length(PIXallcovars_no_MSL_res$Residual)
length(PIXallcovars_no_MSL_res_out_rm$res_out_rm)
length(PIXallcovars_no_MSL_res$Residual) - length(PIXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(PIXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(PIXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(PIXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
PIXallcovars_no_MSL_res <- full_join(PIXallcovars_no_MSL_res,PIXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(PIXallcovars_no_MSL_res)
#View(PIXallcovars_no_MSL_res)

#Transformação boxcox
PIXallcovars_no_MSL_res$Residual.t<-boxcox.ols(PIXallcovars_no_MSL_res$Residual)
boxplot(PIXallcovars_no_MSL_res$Residual.t)
hist(PIXallcovars_no_MSL_res$Residual.t)
rug(PIXallcovars_no_MSL_res$Residual.t)
summary(PIXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
PIXallcovars_no_MSL_res$LOG <- log10(PIXallcovars_no_MSL_res$Residual)
boxplot(PIXallcovars_no_MSL_res$LOG)
hist(PIXallcovars_no_MSL_res$LOG)
rug(PIXallcovars_no_MSL_res$LOG)
summary(PIXallcovars_no_MSL_res$LOG)

#View(PIXallcovars_no_MSL_res)

write.table(PIXallcovars_no_MSL_res,"PI_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


PI.STATUSXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/PI.STATUS.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
# hist(PI.STATUSXallcovars_no_MSL_res$Residual)
# boxplot(PI.STATUSXallcovars_no_MSL_res$Residual)
# length(PI.STATUSXallcovars_no_MSL_res$Residual)
# PI.STATUSXallcovars_no_MSL_res$Residual[PI.STATUSXallcovars_no_MSL_res$Residual %in% boxplot.stats(PI.STATUSXallcovars_no_MSL_res$Residual)$out]
# length(PI.STATUSXallcovars_no_MSL_res$Residual[PI.STATUSXallcovars_no_MSL_res$Residual %in% boxplot.stats(PI.STATUSXallcovars_no_MSL_res$Residual)$out])
# PI.STATUSXallcovars_no_MSL_res_out_rm <- PI.STATUSXallcovars_no_MSL_res$Residual
# nrow(PI.STATUSXallcovars_no_MSL_res_out_rm)
# #View(PI.STATUSXallcovars_no_MSL_res_out_rm)
# PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm <- PI.STATUSXallcovars_no_MSL_res_out_rm$Residual
# #View(PI.STATUSXallcovars_no_MSL_res_out_rm)
# PI.STATUSXallcovars_no_MSL_res_out_rm <- within(PI.STATUSXallcovars_no_MSL_res_out_rm, rm(Residual))
# #View(PI.STATUSXallcovars_no_MSL_res_out_rm)
# length(PI.STATUSXallcovars_no_MSL_res$Residual)
# length(PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm)
# length(PI.STATUSXallcovars_no_MSL_res$Residual) - length(PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm)
# boxplot(PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm)
# hist(PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm)
# rug(PI.STATUSXallcovars_no_MSL_res_out_rm$res_out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# PI.STATUSXallcovars_no_MSL_res <- full_join(PI.STATUSXallcovars_no_MSL_res,PI.STATUSXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(PI.STATUSXallcovars_no_MSL_res)
# View(PI.STATUSXallcovars_no_MSL_res)


#Transformação boxcox
PI.STATUSXallcovars_no_MSL_res$Residual.t<-boxcox.ols(PI.STATUSXallcovars_no_MSL_res$Residual)
boxplot(PI.STATUSXallcovars_no_MSL_res$Residual.t)
hist(PI.STATUSXallcovars_no_MSL_res$Residual.t)
rug(PI.STATUSXallcovars_no_MSL_res$Residual.t)
summary(PI.STATUSXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
PI.STATUSXallcovars_no_MSL_res$LOG <- log10(PI.STATUSXallcovars_no_MSL_res$Residual)
boxplot(PI.STATUSXallcovars_no_MSL_res$LOG)
hist(PI.STATUSXallcovars_no_MSL_res$LOG)
rug(PI.STATUSXallcovars_no_MSL_res$LOG)
summary(PI.STATUSXallcovars_no_MSL_res$LOG)

#View(PI.STATUSXallcovars_no_MSL_res)

write.table(PI.STATUSXallcovars_no_MSL_res,"PI.STATUS_ALLCOVARS_no_MSL_RESIDUAL_PHENOm.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

IFNGXvaccine_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IFNG.VACCINE.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IFNGXvaccine_res$Residual)
boxplot(IFNGXvaccine_res$Residual)
length(IFNGXvaccine_res$Residual)
IFNGXvaccine_res$Residual[IFNGXvaccine_res$Residual %in% boxplot.stats(IFNGXvaccine_res$Residual)$out]
length(IFNGXvaccine_res$Residual[IFNGXvaccine_res$Residual %in% boxplot.stats(IFNGXvaccine_res$Residual)$out])
IFNGXvaccine_res_out_rm <- subset(IFNGXvaccine_res, Residual > -544.9409 & Residual < 388.2073)
nrow(IFNGXvaccine_res_out_rm)
#View(IFNGXvaccine_res_out_rm)
IFNGXvaccine_res_out_rm$res_out_rm <- IFNGXvaccine_res_out_rm$Residual
#View(IFNGXvaccine_res_out_rm)
IFNGXvaccine_res_out_rm <- within(IFNGXvaccine_res_out_rm, rm(Residual))
#View(IFNGXvaccine_res_out_rm)
length(IFNGXvaccine_res$Residual)
length(IFNGXvaccine_res_out_rm$res_out_rm)
length(IFNGXvaccine_res$Residual) - length(IFNGXvaccine_res_out_rm$res_out_rm)
boxplot(IFNGXvaccine_res_out_rm$res_out_rm)
hist(IFNGXvaccine_res_out_rm$res_out_rm)
rug(IFNGXvaccine_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IFNGXvaccine_res <- full_join(IFNGXvaccine_res,IFNGXvaccine_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IFNGXvaccine_res)
#View(IFNGXvaccine_res)

#Transformação boxcox
IFNGXvaccine_res$Residual.t<-boxcox.ols(IFNGXvaccine_res$Residual)
boxplot(IFNGXvaccine_res$Residual.t)
hist(IFNGXvaccine_res$Residual.t)
rug(IFNGXvaccine_res$Residual.t)
summary(IFNGXvaccine_res$Residual.t)

#Transformação logarítmica
IFNGXvaccine_res$LOG <- log10(IFNGXvaccine_res$Residual)
boxplot(IFNGXvaccine_res$LOG)
hist(IFNGXvaccine_res$LOG)
rug(IFNGXvaccine_res$LOG)
summary(IFNGXvaccine_res$LOG)

#View(IFNGXvaccine_res)

write.table(IFNGXvaccine_res,"IFNG_VACCINE_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


TNFAXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/TNFA.cov.allcovars_no_MSL.indi.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(TNFAXallcovars_no_MSL_res$Residual)
boxplot(TNFAXallcovars_no_MSL_res$Residual)
length(TNFAXallcovars_no_MSL_res$Residual)
TNFAXallcovars_no_MSL_res$Residual[TNFAXallcovars_no_MSL_res$Residual %in% boxplot.stats(TNFAXallcovars_no_MSL_res$Residual)$out]
length(TNFAXallcovars_no_MSL_res$Residual[TNFAXallcovars_no_MSL_res$Residual %in% boxplot.stats(TNFAXallcovars_no_MSL_res$Residual)$out])
TNFAXallcovars_no_MSL_res_out_rm <- subset(TNFAXallcovars_no_MSL_res, Residual > -2.922735 & Residual < 3.038026)
nrow(TNFAXallcovars_no_MSL_res_out_rm)
#View(TNFAXallcovars_no_MSL_res_out_rm)
TNFAXallcovars_no_MSL_res_out_rm$res_out_rm <- TNFAXallcovars_no_MSL_res_out_rm$Residual
#View(TNFAXallcovars_no_MSL_res_out_rm)
TNFAXallcovars_no_MSL_res_out_rm <- within(TNFAXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(TNFAXallcovars_no_MSL_res_out_rm)
length(TNFAXallcovars_no_MSL_res$Residual)
length(TNFAXallcovars_no_MSL_res_out_rm$res_out_rm)
length(TNFAXallcovars_no_MSL_res$Residual) - length(TNFAXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(TNFAXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(TNFAXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(TNFAXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TNFAXallcovars_no_MSL_res <- full_join(TNFAXallcovars_no_MSL_res,TNFAXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TNFAXallcovars_no_MSL_res)
#View(TNFAXallcovars_no_MSL_res)

#Transformação boxcox
TNFAXallcovars_no_MSL_res$Residual.t<-boxcox.ols(TNFAXallcovars_no_MSL_res$Residual)
boxplot(TNFAXallcovars_no_MSL_res$Residual.t)
hist(TNFAXallcovars_no_MSL_res$Residual.t)
rug(TNFAXallcovars_no_MSL_res$Residual.t)
summary(TNFAXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
TNFAXallcovars_no_MSL_res$LOG <- log10(TNFAXallcovars_no_MSL_res$Residual)
boxplot(TNFAXallcovars_no_MSL_res$LOG)
hist(TNFAXallcovars_no_MSL_res$LOG)
rug(TNFAXallcovars_no_MSL_res$LOG)
summary(TNFAXallcovars_no_MSL_res$LOG)

#View(TNFAXallcovars_no_MSL_res)

write.table(TNFAXallcovars_no_MSL_res,"TNFA_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

IL10Xallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IL10.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IL10Xallcovars_no_MSL_res$Residual)
boxplot(IL10Xallcovars_no_MSL_res$Residual)
length(IL10Xallcovars_no_MSL_res$Residual)
IL10Xallcovars_no_MSL_res$Residual[IL10Xallcovars_no_MSL_res$Residual %in% boxplot.stats(IL10Xallcovars_no_MSL_res$Residual)$out]
length(IL10Xallcovars_no_MSL_res$Residual[IL10Xallcovars_no_MSL_res$Residual %in% boxplot.stats(IL10Xallcovars_no_MSL_res$Residual)$out])
IL10Xallcovars_no_MSL_res_out_rm <- subset(IL10Xallcovars_no_MSL_res, Residual > -17.51029 & Residual < 11.50499)
nrow(IL10Xallcovars_no_MSL_res_out_rm)
#View(IL10Xallcovars_no_MSL_res_out_rm)
IL10Xallcovars_no_MSL_res_out_rm$res_out_rm <- IL10Xallcovars_no_MSL_res_out_rm$Residual
#View(IL10Xallcovars_no_MSL_res_out_rm)
IL10Xallcovars_no_MSL_res_out_rm <- within(IL10Xallcovars_no_MSL_res_out_rm, rm(Residual))
#View(IL10Xallcovars_no_MSL_res_out_rm)
length(IL10Xallcovars_no_MSL_res$Residual)
length(IL10Xallcovars_no_MSL_res_out_rm$res_out_rm)
length(IL10Xallcovars_no_MSL_res$Residual) - length(IL10Xallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(IL10Xallcovars_no_MSL_res_out_rm$res_out_rm)
hist(IL10Xallcovars_no_MSL_res_out_rm$res_out_rm)
rug(IL10Xallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL10Xallcovars_no_MSL_res <- full_join(IL10Xallcovars_no_MSL_res,IL10Xallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL10Xallcovars_no_MSL_res)
#View(IL10Xallcovars_no_MSL_res)

#Transformação boxcox
IL10Xallcovars_no_MSL_res$Residual.t<-boxcox.ols(IL10Xallcovars_no_MSL_res$Residual)
boxplot(IL10Xallcovars_no_MSL_res$Residual.t)
hist(IL10Xallcovars_no_MSL_res$Residual.t)
rug(IL10Xallcovars_no_MSL_res$Residual.t)
summary(IL10Xallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
IL10Xallcovars_no_MSL_res$LOG <- log10(IL10Xallcovars_no_MSL_res$Residual)
boxplot(IL10Xallcovars_no_MSL_res$LOG)
hist(IL10Xallcovars_no_MSL_res$LOG)
rug(IL10Xallcovars_no_MSL_res$LOG)
summary(IL10Xallcovars_no_MSL_res$LOG)

#View(IL10Xallcovars_no_MSL_res)

write.table(IL10Xallcovars_no_MSL_res,"IL10_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


IL4Xallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/IL4.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(IL4Xallcovars_no_MSL_res$Residual)
boxplot(IL4Xallcovars_no_MSL_res$Residual)
length(IL4Xallcovars_no_MSL_res$Residual)
IL4Xallcovars_no_MSL_res$Residual[IL4Xallcovars_no_MSL_res$Residual %in% boxplot.stats(IL4Xallcovars_no_MSL_res$Residual)$out]
length(IL4Xallcovars_no_MSL_res$Residual[IL4Xallcovars_no_MSL_res$Residual %in% boxplot.stats(IL4Xallcovars_no_MSL_res$Residual)$out])
IL4Xallcovars_no_MSL_res_out_rm <- subset(IL4Xallcovars_no_MSL_res, Residual < 192.4286)
nrow(IL4Xallcovars_no_MSL_res_out_rm)
#View(IL4Xallcovars_no_MSL_res_out_rm)
IL4Xallcovars_no_MSL_res_out_rm$res_out_rm <- IL4Xallcovars_no_MSL_res_out_rm$Residual
#View(IL4Xallcovars_no_MSL_res_out_rm)
IL4Xallcovars_no_MSL_res_out_rm <- within(IL4Xallcovars_no_MSL_res_out_rm, rm(Residual))
#View(IL4Xallcovars_no_MSL_res_out_rm)
length(IL4Xallcovars_no_MSL_res$Residual)
length(IL4Xallcovars_no_MSL_res_out_rm$res_out_rm)
length(IL4Xallcovars_no_MSL_res$Residual) - length(IL4Xallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(IL4Xallcovars_no_MSL_res_out_rm$res_out_rm)
hist(IL4Xallcovars_no_MSL_res_out_rm$res_out_rm)
rug(IL4Xallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL4Xallcovars_no_MSL_res <- full_join(IL4Xallcovars_no_MSL_res,IL4Xallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL4Xallcovars_no_MSL_res)
#View(IL4Xallcovars_no_MSL_res)

#Transformação boxcox
IL4Xallcovars_no_MSL_res$Residual.t<-boxcox.ols(IL4Xallcovars_no_MSL_res$Residual)
boxplot(IL4Xallcovars_no_MSL_res$Residual.t)
hist(IL4Xallcovars_no_MSL_res$Residual.t)
rug(IL4Xallcovars_no_MSL_res$Residual.t)
summary(IL4Xallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
IL4Xallcovars_no_MSL_res$LOG <- log10(IL4Xallcovars_no_MSL_res$Residual)
boxplot(IL4Xallcovars_no_MSL_res$LOG)
hist(IL4Xallcovars_no_MSL_res$LOG)
rug(IL4Xallcovars_no_MSL_res$LOG)
summary(IL4Xallcovars_no_MSL_res$LOG)

#View(IL4Xallcovars_no_MSL_res)

write.table(IL4Xallcovars_no_MSL_res,"IL4_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

TGFBXallcovars_no_MSL_res <- read.table(file ="/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/TGFB.allcovars_no_MSL.res_col.txt", header=T, sep="\t", dec =".")

#Remoção dos outliers
hist(TGFBXallcovars_no_MSL_res$Residual)
boxplot(TGFBXallcovars_no_MSL_res$Residual)
length(TGFBXallcovars_no_MSL_res$Residual)
TGFBXallcovars_no_MSL_res$Residual[TGFBXallcovars_no_MSL_res$Residual %in% boxplot.stats(TGFBXallcovars_no_MSL_res$Residual)$out]
length(TGFBXallcovars_no_MSL_res$Residual[TGFBXallcovars_no_MSL_res$Residual %in% boxplot.stats(TGFBXallcovars_no_MSL_res$Residual)$out])
TGFBXallcovars_no_MSL_res_out_rm <- subset(TGFBXallcovars_no_MSL_res, Residual > -85.95065 & Residual < 82.37950)
nrow(TGFBXallcovars_no_MSL_res_out_rm)
#View(TGFBXallcovars_no_MSL_res_out_rm)
TGFBXallcovars_no_MSL_res_out_rm$res_out_rm <- TGFBXallcovars_no_MSL_res_out_rm$Residual
#View(TGFBXallcovars_no_MSL_res_out_rm)
TGFBXallcovars_no_MSL_res_out_rm <- within(TGFBXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(TGFBXallcovars_no_MSL_res_out_rm)
length(TGFBXallcovars_no_MSL_res$Residual)
length(TGFBXallcovars_no_MSL_res_out_rm$res_out_rm)
length(TGFBXallcovars_no_MSL_res$Residual) - length(TGFBXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(TGFBXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(TGFBXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(TGFBXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TGFBXallcovars_no_MSL_res <- full_join(TGFBXallcovars_no_MSL_res,TGFBXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TGFBXallcovars_no_MSL_res)
#View(TGFBXallcovars_no_MSL_res)

#Transformação boxcox
TGFBXallcovars_no_MSL_res$Residual.t<-boxcox.ols(TGFBXallcovars_no_MSL_res$Residual)
boxplot(TGFBXallcovars_no_MSL_res$Residual.t)
hist(TGFBXallcovars_no_MSL_res$Residual.t)
rug(TGFBXallcovars_no_MSL_res$Residual.t)
summary(TGFBXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
TGFBXallcovars_no_MSL_res$LOG <- log10(TGFBXallcovars_no_MSL_res$Residual)
boxplot(TGFBXallcovars_no_MSL_res$LOG)
hist(TGFBXallcovars_no_MSL_res$LOG)
rug(TGFBXallcovars_no_MSL_res$LOG)
summary(TGFBXallcovars_no_MSL_res$LOG)

#View(TGFBXallcovars_no_MSL_res)

write.table(TGFBXallcovars_no_MSL_res,"TGFB_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

NO_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/NO.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(NO_res$Residual)
boxplot(NO_res$Residual)
length(NO_res$Residual)
NO_res$Residual[NO_res$Residual %in% boxplot.stats(NO_res$Residual)$out]
length(NO_res$Residual[NO_res$Residual %in% boxplot.stats(NO_res$Residual)$out])
NO_res_out_rm <- subset(NO_res, Residual > -2.165642 & Residual < 1.811385)
nrow(NO_res_out_rm)
#View(NO_res_out_rm)
NO_res_out_rm$res_out_rm <- NO_res_out_rm$Residual
#View(NO_res_out_rm)
NO_res_out_rm <- within(NO_res_out_rm, rm(Residual))
#View(NO_res_out_rm)
length(NO_res$Residual)
length(NO_res_out_rm$res_out_rm)
length(NO_res$Residual) - length(NO_res_out_rm$res_out_rm)
boxplot(NO_res_out_rm$res_out_rm)
hist(NO_res_out_rm$res_out_rm)
rug(NO_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
NO_res <- full_join(NO_res,NO_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(NO_res)
#View(NO_res)

# #Transformação boxcox
# NO_res$Residual.t <- boxcox.ols(NO_res$Residual)
# boxplot(NO_res$Residual.t)
# hist(NO_res$Residual.t)
# rug(NO_res$Residual.t)
# summary(NO_res$Residual.t)

#Transformação logarítmica
NO_res$LOG <- log10(NO_res$Residual)
boxplot(NO_res$LOG)
hist(NO_res$LOG)
rug(NO_res$LOG)
summary(NO_res$LOG)

#View(NO_res)

write.table(NO_res,"NO_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

NOXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/NO.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(NOXallcovars_no_MSL_res$Residual)
boxplot(NOXallcovars_no_MSL_res$Residual)
length(NOXallcovars_no_MSL_res$Residual)
NOXallcovars_no_MSL_res$Residual[NOXallcovars_no_MSL_res$Residual %in% boxplot.stats(NOXallcovars_no_MSL_res$Residual)$out]
length(NOXallcovars_no_MSL_res$Residual[NOXallcovars_no_MSL_res$Residual %in% boxplot.stats(NOXallcovars_no_MSL_res$Residual)$out])
NOXallcovars_no_MSL_res_out_rm <- subset(NOXallcovars_no_MSL_res, Residual > -3.040702 & Residual < 2.787357)
nrow(NOXallcovars_no_MSL_res_out_rm)
#View(NOXallcovars_no_MSL_res_out_rm)
NOXallcovars_no_MSL_res_out_rm$res_out_rm <- NOXallcovars_no_MSL_res_out_rm$Residual
#View(NOXallcovars_no_MSL_res_out_rm)
NOXallcovars_no_MSL_res_out_rm <- within(NOXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(NOXallcovars_no_MSL_res_out_rm)
length(NOXallcovars_no_MSL_res$Residual)
length(NOXallcovars_no_MSL_res_out_rm$res_out_rm)
length(NOXallcovars_no_MSL_res$Residual) - length(NOXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(NOXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(NOXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(NOXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
NOXallcovars_no_MSL_res <- full_join(NOXallcovars_no_MSL_res,NOXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(NOXallcovars_no_MSL_res)
#View(NOXallcovars_no_MSL_res)

#Transformação boxcox
NOXallcovars_no_MSL_res$Residual.t<-boxcox.ols(NOXallcovars_no_MSL_res$Residual)
boxplot(NOXallcovars_no_MSL_res$Residual.t)
hist(NOXallcovars_no_MSL_res$Residual.t)
rug(NOXallcovars_no_MSL_res$Residual.t)
summary(NOXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
NOXallcovars_no_MSL_res$LOG <- log10(NOXallcovars_no_MSL_res$Residual)
boxplot(NOXallcovars_no_MSL_res$LOG)
hist(NOXallcovars_no_MSL_res$LOG)
rug(NOXallcovars_no_MSL_res$LOG)
summary(NOXallcovars_no_MSL_res$LOG)

#View(NOXallcovars_no_MSL_res)

write.table(NOXallcovars_no_MSL_res,"NO_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


SOD_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/SOD.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(SOD_res$Residual)
boxplot(SOD_res$Residual)
length(SOD_res$Residual)
SOD_res$Residual[SOD_res$Residual %in% boxplot.stats(SOD_res$Residual)$out]
length(SOD_res$Residual[SOD_res$Residual %in% boxplot.stats(SOD_res$Residual)$out])
SOD_res_out_rm <- subset(SOD_res, Residual < 0.493914)
nrow(SOD_res_out_rm)
#View(SOD_res_out_rm)
SOD_res_out_rm$res_out_rm <- SOD_res_out_rm$Residual
#View(SOD_res_out_rm)
SOD_res_out_rm <- within(SOD_res_out_rm, rm(Residual))
#View(SOD_res_out_rm)
length(SOD_res$Residual)
length(SOD_res_out_rm$res_out_rm)
length(SOD_res$Residual) - length(SOD_res_out_rm$res_out_rm)
boxplot(SOD_res_out_rm$res_out_rm)
hist(SOD_res_out_rm$res_out_rm)
rug(SOD_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
SOD_res <- full_join(SOD_res,SOD_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(SOD_res)
#View(SOD_res)

#Transformação boxcox
SOD_res$Residual.t<-boxcox.ols(SOD_res$Residual)
boxplot(SOD_res$Residual.t)
hist(SOD_res$Residual.t)
rug(SOD_res$Residual.t)
summary(SOD_res$Residual.t)

#Transformação logarítmica
SOD_res$LOG <- log10(SOD_res$Residual)
boxplot(SOD_res$LOG)
hist(SOD_res$LOG)
rug(SOD_res$LOG)
summary(SOD_res$LOG)

#View(SOD_res)

write.table(SOD_res,"SOD_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

SODXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/SOD.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(SODXallcovars_no_MSL_res$Residual)
boxplot(SODXallcovars_no_MSL_res$Residual)
length(SODXallcovars_no_MSL_res$Residual)
SODXallcovars_no_MSL_res$Residual[SODXallcovars_no_MSL_res$Residual %in% boxplot.stats(SODXallcovars_no_MSL_res$Residual)$out]
length(SODXallcovars_no_MSL_res$Residual[SODXallcovars_no_MSL_res$Residual %in% boxplot.stats(SODXallcovars_no_MSL_res$Residual)$out])
SODXallcovars_no_MSL_res_out_rm <- subset(SODXallcovars_no_MSL_res, Residual > -0.503816 & Residual < 0.661213)
nrow(SODXallcovars_no_MSL_res_out_rm)
#View(SODXallcovars_no_MSL_res_out_rm)
SODXallcovars_no_MSL_res_out_rm$res_out_rm <- SODXallcovars_no_MSL_res_out_rm$Residual
#View(SODXallcovars_no_MSL_res_out_rm)
SODXallcovars_no_MSL_res_out_rm <- within(SODXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(SODXallcovars_no_MSL_res_out_rm)
length(SODXallcovars_no_MSL_res$Residual)
length(SODXallcovars_no_MSL_res_out_rm$res_out_rm)
length(SODXallcovars_no_MSL_res$Residual) - length(SODXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(SODXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(SODXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(SODXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
SODXallcovars_no_MSL_res <- full_join(SODXallcovars_no_MSL_res,SODXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(SODXallcovars_no_MSL_res)
#View(SODXallcovars_no_MSL_res)

#Transformação boxcox
SODXallcovars_no_MSL_res$Residual.t<-boxcox.ols(SODXallcovars_no_MSL_res$Residual)
boxplot(SODXallcovars_no_MSL_res$Residual.t)
hist(SODXallcovars_no_MSL_res$Residual.t)
rug(SODXallcovars_no_MSL_res$Residual.t)
summary(SODXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
SODXallcovars_no_MSL_res$LOG <- log10(SODXallcovars_no_MSL_res$Residual)
boxplot(SODXallcovars_no_MSL_res$LOG)
hist(SODXallcovars_no_MSL_res$LOG)
rug(SODXallcovars_no_MSL_res$LOG)
summary(SODXallcovars_no_MSL_res$LOG)

#View(SODXallcovars_no_MSL_res)

write.table(SODXallcovars_no_MSL_res,"SOD_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


TOC_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/TOS.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(TOC_res$Residual)
boxplot(TOC_res$Residual)
length(TOC_res$Residual)
TOC_res$Residual[TOC_res$Residual %in% boxplot.stats(TOC_res$Residual)$out]
length(TOC_res$Residual[TOC_res$Residual %in% boxplot.stats(TOC_res$Residual)$out])
TOC_res_out_rm <- subset(TOC_res, Residual > -92.47320)
nrow(TOC_res_out_rm)
#View(TOC_res_out_rm)
TOC_res_out_rm$res_out_rm <- TOC_res_out_rm$Residual
#View(TOC_res_out_rm)
TOC_res_out_rm <- within(TOC_res_out_rm, rm(Residual))
#View(TOC_res_out_rm)
length(TOC_res$Residual)
length(TOC_res_out_rm$res_out_rm)
length(TOC_res$Residual) - length(TOC_res_out_rm$res_out_rm)
boxplot(TOC_res_out_rm$res_out_rm)
hist(TOC_res_out_rm$res_out_rm)
rug(TOC_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TOC_res <- full_join(TOC_res,TOC_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TOC_res)
#View(TOC_res)

#Transformação boxcox
TOC_res$Residual.t<-boxcox.ols(TOC_res$Residual)
boxplot(TOC_res$Residual.t)
hist(TOC_res$Residual.t)
rug(TOC_res$Residual.t)
summary(TOC_res$Residual.t)

#Transformação logarítmica
TOC_res$LOG <- log10(TOC_res$Residual)
boxplot(TOC_res$LOG)
hist(TOC_res$LOG)
rug(TOC_res$LOG)
summary(TOC_res$LOG)

#View(TOC_res)

write.table(TOC_res,"TOC_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

TOCXcollar_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/TOS.COLLAR.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(TOCXcollar_res$Residual)
boxplot(TOCXcollar_res$Residual)
length(TOCXcollar_res$Residual)
TOCXcollar_res$Residual[TOCXcollar_res$Residual %in% boxplot.stats(TOCXcollar_res$Residual)$out]
length(TOCXcollar_res$Residual[TOCXcollar_res$Residual %in% boxplot.stats(TOCXcollar_res$Residual)$out])
TOCXcollar_res_out_rm <- subset(TOCXcollar_res, Residual > -92.10292)
nrow(TOCXcollar_res_out_rm)
#View(TOCXcollar_res_out_rm)
TOCXcollar_res_out_rm$res_out_rm <- TOCXcollar_res_out_rm$Residual
#View(TOCXcollar_res_out_rm)
TOCXcollar_res_out_rm <- within(TOCXcollar_res_out_rm, rm(Residual))
#View(TOCXcollar_res_out_rm)
length(TOCXcollar_res$Residual)
length(TOCXcollar_res_out_rm$res_out_rm)
length(TOCXcollar_res$Residual) - length(TOCXcollar_res_out_rm$res_out_rm)
boxplot(TOCXcollar_res_out_rm$res_out_rm)
hist(TOCXcollar_res_out_rm$res_out_rm)
rug(TOCXcollar_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TOCXcollar_res <- full_join(TOCXcollar_res,TOCXcollar_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TOCXcollar_res)
#View(TOCXcollar_res)

#Transformação boxcox
TOCXcollar_res$Residual.t<-boxcox.ols(TOCXcollar_res$Residual)
boxplot(TOCXcollar_res$Residual.t)
hist(TOCXcollar_res$Residual.t)
rug(TOCXcollar_res$Residual.t)
summary(TOCXcollar_res$Residual.t)

#Transformação logarítmica
TOCXcollar_res$LOG <- log10(TOCXcollar_res$Residual)
boxplot(TOCXcollar_res$LOG)
hist(TOCXcollar_res$LOG)
rug(TOCXcollar_res$LOG)
summary(TOCXcollar_res$LOG)

#View(TOCXcollar_res)

write.table(TOCXcollar_res,"TOC_COLLAR_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

TACXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/TAS.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(TACXallcovars_no_MSL_res$Residual)
boxplot(TACXallcovars_no_MSL_res$Residual)
length(TACXallcovars_no_MSL_res$Residual)
TACXallcovars_no_MSL_res$Residual[TACXallcovars_no_MSL_res$Residual %in% boxplot.stats(TACXallcovars_no_MSL_res$Residual)$out]
length(TACXallcovars_no_MSL_res$Residual[TACXallcovars_no_MSL_res$Residual %in% boxplot.stats(TACXallcovars_no_MSL_res$Residual)$out])
TACXallcovars_no_MSL_res_out_rm <- subset(TACXallcovars_no_MSL_res, Residual > -0.105914 & Residual < 0.106002)
nrow(TACXallcovars_no_MSL_res_out_rm)
#View(TACXallcovars_no_MSL_res_out_rm)
TACXallcovars_no_MSL_res_out_rm$res_out_rm <- TACXallcovars_no_MSL_res_out_rm$Residual
#View(TACXallcovars_no_MSL_res_out_rm)
TACXallcovars_no_MSL_res_out_rm <- within(TACXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(TACXallcovars_no_MSL_res_out_rm)
length(TACXallcovars_no_MSL_res$Residual)
length(TACXallcovars_no_MSL_res_out_rm$res_out_rm)
length(TACXallcovars_no_MSL_res$Residual) - length(TACXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(TACXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(TACXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(TACXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TACXallcovars_no_MSL_res <- full_join(TACXallcovars_no_MSL_res,TACXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TACXallcovars_no_MSL_res)
#View(TACXallcovars_no_MSL_res)

#Transformação boxcox
TACXallcovars_no_MSL_res$Residual.t<-boxcox.ols(TACXallcovars_no_MSL_res$Residual)
boxplot(TACXallcovars_no_MSL_res$Residual.t)
hist(TACXallcovars_no_MSL_res$Residual.t)
rug(TACXallcovars_no_MSL_res$Residual.t)
summary(TACXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
TACXallcovars_no_MSL_res$LOG <- log10(TACXallcovars_no_MSL_res$Residual)
boxplot(TACXallcovars_no_MSL_res$LOG)
hist(TACXallcovars_no_MSL_res$LOG)
rug(TACXallcovars_no_MSL_res$LOG)
summary(TACXallcovars_no_MSL_res$LOG)

#View(TACXallcovars_no_MSL_res)

write.table(TACXallcovars_no_MSL_res,"TAC_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

MDAXallcovars_no_MSL_res <- read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/Analise_Projeto/GWAS_H-W_test/Heritability_H-W_FILT/HERITABILITY_ALLTRAITS_23-3-2024/Modelos_selecionados_para_GWAS/MDA.allcovars_no_MSL.res_col.txt", header=T, sep="\t", stringsAsFactors = F)

#Remoção dos outliers
hist(MDAXallcovars_no_MSL_res$Residual)
boxplot(MDAXallcovars_no_MSL_res$Residual)
length(MDAXallcovars_no_MSL_res$Residual)
MDAXallcovars_no_MSL_res$Residual[MDAXallcovars_no_MSL_res$Residual %in% boxplot.stats(MDAXallcovars_no_MSL_res$Residual)$out]
length(MDAXallcovars_no_MSL_res$Residual[MDAXallcovars_no_MSL_res$Residual %in% boxplot.stats(MDAXallcovars_no_MSL_res$Residual)$out])
MDAXallcovars_no_MSL_res_out_rm <- subset(MDAXallcovars_no_MSL_res, Residual < 59.31210)
nrow(MDAXallcovars_no_MSL_res_out_rm)
#View(MDAXallcovars_no_MSL_res_out_rm)
MDAXallcovars_no_MSL_res_out_rm$res_out_rm <- MDAXallcovars_no_MSL_res_out_rm$Residual
#View(MDAXallcovars_no_MSL_res_out_rm)
MDAXallcovars_no_MSL_res_out_rm <- within(MDAXallcovars_no_MSL_res_out_rm, rm(Residual))
#View(MDAXallcovars_no_MSL_res_out_rm)
length(MDAXallcovars_no_MSL_res$Residual)
length(MDAXallcovars_no_MSL_res_out_rm$res_out_rm)
length(MDAXallcovars_no_MSL_res$Residual) - length(MDAXallcovars_no_MSL_res_out_rm$res_out_rm)
boxplot(MDAXallcovars_no_MSL_res_out_rm$res_out_rm)
hist(MDAXallcovars_no_MSL_res_out_rm$res_out_rm)
rug(MDAXallcovars_no_MSL_res_out_rm$res_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
MDAXallcovars_no_MSL_res <- full_join(MDAXallcovars_no_MSL_res,MDAXallcovars_no_MSL_res_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(MDAXallcovars_no_MSL_res)
#View(MDAXallcovars_no_MSL_res)

#Transformação boxcox
MDAXallcovars_no_MSL_res$Residual.t<-boxcox.ols(MDAXallcovars_no_MSL_res$Residual)
boxplot(MDAXallcovars_no_MSL_res$Residual.t)
hist(MDAXallcovars_no_MSL_res$Residual.t)
rug(MDAXallcovars_no_MSL_res$Residual.t)
summary(MDAXallcovars_no_MSL_res$Residual.t)

#Transformação logarítmica
MDAXallcovars_no_MSL_res$LOG <- log10(MDAXallcovars_no_MSL_res$Residual)
boxplot(MDAXallcovars_no_MSL_res$LOG)
hist(MDAXallcovars_no_MSL_res$LOG)
rug(MDAXallcovars_no_MSL_res$LOG)
summary(MDAXallcovars_no_MSL_res$LOG)

#View(MDAXallcovars_no_MSL_res)

write.table(MDAXallcovars_no_MSL_res,"MDA_ALLCOVARS_no_MSL_RESIDUAL_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

############################################################################################################################################################################
#################### Ajustando traits com covariáveis e normalizando apenas trait puro e os modelos escolhidos ############################################################# 
############################################################################################################################################################################

#Criando um subset STAGING para normalização
STAGING <- select(pheno,IID,FID,STAGING)
#View(STAGING)

#Remoção dos outliers
hist(STAGING$STAGING)
boxplot(STAGING$STAGING)
length(STAGING$STAGING)
STAGING$STAGING[STAGING$STAGING %in% boxplot.stats(STAGING$STAGING)$out]
length(STAGING$STAGING[STAGING$STAGING %in% boxplot.stats(STAGING$STAGING)$out])
STAGING_out_rm <- subset(STAGING, STAGING > 1)
nrow(STAGING_out_rm)
#View(STAGING_out_rm)
STAGING_out_rm$out_rm <- STAGING_out_rm$STAGING
#View(STAGING_out_rm)
STAGING_out_rm <- within(STAGING_out_rm, rm(STAGING))
#View(STAGINGXallcovars_no_MSL_res_out_rm)
length(STAGING$STAGING)
length(STAGING_out_rm$out_rm)
length(STAGING$STAGING) - length(STAGING_out_rm$out_rm)
boxplot(STAGING_out_rm$out_rm)
hist(STAGING_out_rm$out_rm)
rug(STAGING_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
STAGING <- full_join(STAGING,STAGING_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(STAGING)
#View(STAGINGXallcovars_no_MSL_res)

#Transformação boxcox
STAGING$STAGING.t<-boxcox.ols(STAGING$STAGING)
boxplot(STAGING$STAGING.t)
hist(STAGING$STAGING.t)
rug(STAGING$STAGING.t)
summary(STAGING$STAGING.t)

#Transformação logarítmica
STAGING$LOG <- log10(STAGING$STAGING)
boxplot(STAGING$LOG)
hist(STAGING$LOG)
rug(STAGING$LOG)
summary(STAGING$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$STAGING)
hist(pheno$STAGING)
model<-lm(STAGING ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(STAGING$STAGING)
hist(model$residuals)
rug(model$residuals)
STAGING[names(model$residuals),"STAGINGXallcovars_no_MSL_RES"] <- model$residuals
#View(STAGING)

#Remoção dos outliers
hist(STAGING$STAGINGXallcovars_no_MSL_RES)
boxplot(STAGING$STAGINGXallcovars_no_MSL_RES)
STAGING$STAGINGXallcovars_no_MSL_RES[STAGING$STAGINGXallcovars_no_MSL_RES %in% boxplot.stats(STAGING$STAGINGXallcovars_no_MSL_RES)$out]
length(STAGING$STAGINGXallcovars_no_MSL_RES[STAGING$STAGINGXallcovars_no_MSL_RES %in% boxplot.stats(STAGING$STAGINGXallcovars_no_MSL_RES)$out])
STAGINGXalcovars_no_MSL_out_rm <- subset(STAGING, STAGINGXallcovars_no_MSL_RES < 2.181968 & STAGINGXallcovars_no_MSL_RES >= -2.087437)
nrow(STAGINGXalcovars_no_MSL_out_rm)
STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm <- STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES
STAGINGXalcovars_no_MSL_out_rm <- within(STAGINGXalcovars_no_MSL_out_rm, rm(STAGINGXallcovars_no_MSL_RES))
length(STAGING$STAGINGXallcovars_no_MSL_RES)
length(STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm)
length(STAGING$STAGINGXallcovars_no_MSL_RES - STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm)
boxplot(STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm)
hist(STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm)
rug(STAGINGXalcovars_no_MSL_out_rm$STAGINGXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
STAGING <- full_join(STAGING,STAGINGXalcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(STAGING)
#View(STAGING)


#Transformação boxcox
STAGING$STAGINGXallcovars_no_MSL_RES.t<-boxcox.ols(STAGING$STAGINGXallcovars_no_MSL_RES)
boxplot(STAGING$STAGINGXallcovars_no_MSL_RES.t)
hist(STAGING$STAGINGXallcovars_no_MSL_RES.t)
rug(STAGING$STAGINGXallcovars_no_MSL_RES.t)
summary(STAGING$STAGINGXallcovars_no_MSL_RES.t)

#Transformação logarítmica
STAGING$STAGINGXallcovars_no_MSL_RES.LOG <- log10(STAGING$STAGINGXallcovars_no_MSL_RES)
boxplot(STAGING$STAGINGXallcovars_no_MSL_RES.LOG)
hist(STAGING$STAGINGXallcovars_no_MSL_RES.LOG)
rug(STAGING$STAGINGXallcovars_no_MSL_RES.LOG)
summary(STAGING$STAGINGXallcovars_no_MSL_RES.LOG)

View(STAGING)

write.table(STAGING,"STAGING_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)




#Criando um subset LOG.PL.LFND para normalização
LOG.PL.LFND <- select(pheno,IID,FID,LOG.PL.LFND)
#View(LOG.PL.LFND)

#Remoção dos outliers
hist(LOG.PL.LFND$LOG.PL.LFND)
boxplot(LOG.PL.LFND$LOG.PL.LFND)
length(LOG.PL.LFND$LOG.PL.LFND)
LOG.PL.LFND$LOG.PL.LFND[LOG.PL.LFND$LOG.PL.LFND %in% boxplot.stats(LOG.PL.LFND$LOG.PL.LFND)$out]
length(LOG.PL.LFND$LOG.PL.LFND[LOG.PL.LFND$LOG.PL.LFND %in% boxplot.stats(LOG.PL.LFND$LOG.PL.LFND)$out])
LOG.PL.LFND_out_rm <- subset(LOG.PL.LFND, LOG.PL.LFND < 7.55)
nrow(LOG.PL.LFND_out_rm)
#View(LOG.PL.LFND_out_rm)
LOG.PL.LFND_out_rm$out_rm <- LOG.PL.LFND_out_rm$LOG.PL.LFND
#View(LOG.PL.LFND_out_rm)
LOG.PL.LFND_out_rm <- within(LOG.PL.LFND_out_rm, rm(LOG.PL.LFND))
#View(LOG.PL.LFNDXorigin_res_out_rm)
length(LOG.PL.LFND$LOG.PL.LFND)
length(LOG.PL.LFND_out_rm$out_rm)
length(LOG.PL.LFND$LOG.PL.LFND) - length(LOG.PL.LFND_out_rm$out_rm)
boxplot(LOG.PL.LFND_out_rm$out_rm)
hist(LOG.PL.LFND_out_rm$out_rm)
rug(LOG.PL.LFND_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
LOG.PL.LFND <- full_join(LOG.PL.LFND,LOG.PL.LFND_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(LOG.PL.LFND)
#View(LOG.PL.LFNDXorigin_res)

#Transformação boxcox
LOG.PL.LFND$LOG.PL.LFND.t<-boxcox.ols(LOG.PL.LFND$LOG.PL.LFND)
boxplot(LOG.PL.LFND$LOG.PL.LFND.t)
hist(LOG.PL.LFND$LOG.PL.LFND.t)
rug(LOG.PL.LFND$LOG.PL.LFND.t)
summary(LOG.PL.LFND$LOG.PL.LFND.t)

#Transformação logarítmica
LOG.PL.LFND$LOG <- log10(LOG.PL.LFND$LOG.PL.LFND)
boxplot(LOG.PL.LFND$LOG)
hist(LOG.PL.LFND$LOG)
rug(LOG.PL.LFND$LOG)
summary(LOG.PL.LFND$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$LOG.PL.LFND)
hist(pheno$LOG.PL.LFND)
model<-lm(LOG.PL.LFND ~ ORIGIN,data = pheno)
summary(model)
length(model$residuals)
length(LOG.PL.LFND$LOG.PL.LFND)
hist(model$residuals)
rug(model$residuals)
LOG.PL.LFND[names(model$residuals),"LOG.PL.LFNDXorigin_RES"] <- model$residuals
#View(LOG.PL.LFND)

#Remoção dos outliers
# hist(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)
# boxplot(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)
# LOG.PL.LFND$LOG.PL.LFNDXorigin_RES[LOG.PL.LFND$LOG.PL.LFNDXorigin_RES %in% boxplot.stats(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)$out]
# length(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES[LOG.PL.LFND$LOG.PL.LFNDXorigin_RES %in% boxplot.stats(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)$out])
# LOG.PL.LFNDXalcovars_no_MSL_out_rm <- subset(LOG.PL.LFND, LOG.PL.LFNDXorigin_RES < 57.41991)
# nrow(LOG.PL.LFNDXalcovars_no_MSL_out_rm)
# LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm <- LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES
# LOG.PL.LFNDXalcovars_no_MSL_out_rm <- within(LOG.PL.LFNDXalcovars_no_MSL_out_rm, rm(LOG.PL.LFNDXorigin_RES))
# length(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)
# length(LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm)
# length(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES - LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm)
# boxplot(LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm)
# hist(LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm)
# rug(LOG.PL.LFNDXalcovars_no_MSL_out_rm$LOG.PL.LFNDXorigin_RES_out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# LOG.PL.LFND <- full_join(LOG.PL.LFND,LOG.PL.LFNDXalcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(LOG.PL.LFND)
# View(LOG.PL.LFND)


#Transformação boxcox
LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.t<-boxcox.ols(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)
boxplot(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.t)
hist(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.t)
rug(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.t)
summary(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.t)

#Transformação logarítmica
LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.LOG <- log10(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES)
boxplot(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.LOG)
hist(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.LOG)
rug(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.LOG)
summary(LOG.PL.LFND$LOG.PL.LFNDXorigin_RES.LOG)

View(LOG.PL.LFND)

write.table(LOG.PL.LFND,"LOG.PL.LFND_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset PL.BLOOD para normalização
PL.BLOOD <- select(pheno,IID,FID,PL.BLOOD)
#View(PL.BLOOD)

#Remoção dos outliers
hist(PL.BLOOD$PL.BLOOD)
boxplot(PL.BLOOD$PL.BLOOD)
length(PL.BLOOD$PL.BLOOD)
PL.BLOOD$PL.BLOOD[PL.BLOOD$PL.BLOOD %in% boxplot.stats(PL.BLOOD$PL.BLOOD)$out]
length(PL.BLOOD$PL.BLOOD[PL.BLOOD$PL.BLOOD %in% boxplot.stats(PL.BLOOD$PL.BLOOD)$out])
PL.BLOOD_out_rm <- subset(PL.BLOOD, PL.BLOOD < 78.5)
nrow(PL.BLOOD_out_rm)
#View(PL.BLOOD_out_rm)
PL.BLOOD_out_rm$out_rm <- PL.BLOOD_out_rm$PL.BLOOD
#View(PL.BLOOD_out_rm)
PL.BLOOD_out_rm <- within(PL.BLOOD_out_rm, rm(PL.BLOOD))
#View(PL.BLOODXcollar_res_out_rm)
length(PL.BLOOD$PL.BLOOD)
length(PL.BLOOD_out_rm$out_rm)
length(PL.BLOOD$PL.BLOOD) - length(PL.BLOOD_out_rm$out_rm)
boxplot(PL.BLOOD_out_rm$out_rm)
hist(PL.BLOOD_out_rm$out_rm)
rug(PL.BLOOD_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
PL.BLOOD <- full_join(PL.BLOOD,PL.BLOOD_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(PL.BLOOD)
#View(PL.BLOODXcollar_res)

#Transformação boxcox
PL.BLOOD$PL.BLOOD.t<-boxcox.ols(PL.BLOOD$PL.BLOOD)
boxplot(PL.BLOOD$PL.BLOOD.t)
hist(PL.BLOOD$PL.BLOOD.t)
rug(PL.BLOOD$PL.BLOOD.t)
summary(PL.BLOOD$PL.BLOOD.t)

#Transformação logarítmica
PL.BLOOD$LOG <- log10(PL.BLOOD$PL.BLOOD)
boxplot(PL.BLOOD$LOG)
hist(PL.BLOOD$LOG)
rug(PL.BLOOD$LOG)
summary(PL.BLOOD$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$PL.BLOOD)
hist(pheno$PL.BLOOD)
model<-lm(PL.BLOOD ~ COLLAR,data = pheno)
summary(model)
length(model$residuals)
length(PL.BLOOD$PL.BLOOD)
hist(model$residuals)
rug(model$residuals)
PL.BLOOD[names(model$residuals),"PL.BLOODXcollar_RES"] <- model$residuals
#View(PL.BLOOD)

#Remoção dos outliers
hist(PL.BLOOD$PL.BLOODXcollar_RES)
boxplot(PL.BLOOD$PL.BLOODXcollar_RES)
PL.BLOOD$PL.BLOODXcollar_RES[PL.BLOOD$PL.BLOODXcollar_RES %in% boxplot.stats(PL.BLOOD$PL.BLOODXcollar_RES)$out]
length(PL.BLOOD$PL.BLOODXcollar_RES[PL.BLOOD$PL.BLOODXcollar_RES %in% boxplot.stats(PL.BLOOD$PL.BLOODXcollar_RES)$out])
PL.BLOODXcollar_out_rm <- subset(PL.BLOOD, PL.BLOODXcollar_RES < 49.21500)
nrow(PL.BLOODXcollar_out_rm)
PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm <- PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES
PL.BLOODXcollar_out_rm <- within(PL.BLOODXcollar_out_rm, rm(PL.BLOODXcollar_RES))
length(PL.BLOOD$PL.BLOODXcollar_RES)
length(PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm)
length(PL.BLOOD$PL.BLOODXcollar_RES - PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm)
boxplot(PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm)
hist(PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm)
rug(PL.BLOODXcollar_out_rm$PL.BLOODXcollar_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
PL.BLOOD <- full_join(PL.BLOOD,PL.BLOODXcollar_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(PL.BLOOD)
View(PL.BLOOD)


#Transformação boxcox
PL.BLOOD$PL.BLOODXcollar_RES.t<-boxcox.ols(PL.BLOOD$PL.BLOODXcollar_RES)
boxplot(PL.BLOOD$PL.BLOODXcollar_RES.t)
hist(PL.BLOOD$PL.BLOODXcollar_RES.t)
rug(PL.BLOOD$PL.BLOODXcollar_RES.t)
summary(PL.BLOOD$PL.BLOODXcollar_RES.t)

#Transformação logarítmica
PL.BLOOD$PL.BLOODXcollar_RES.LOG <- log10(PL.BLOOD$PL.BLOODXcollar_RES)
boxplot(PL.BLOOD$PL.BLOODXcollar_RES.LOG)
hist(PL.BLOOD$PL.BLOODXcollar_RES.LOG)
rug(PL.BLOOD$PL.BLOODXcollar_RES.LOG)
summary(PL.BLOOD$PL.BLOODXcollar_RES.LOG)

View(PL.BLOOD)

write.table(PL.BLOOD,"PL.BLOOD_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IgA para normalização
IgA <- select(pheno,IID,FID,IgA)
#View(IgA)

#Remoção dos outliers
hist(IgA$IgA)
boxplot(IgA$IgA)
length(IgA$IgA)
IgA$IgA[IgA$IgA %in% boxplot.stats(IgA$IgA)$out]
length(IgA$IgA[IgA$IgA %in% boxplot.stats(IgA$IgA)$out])
IgA_out_rm <- subset(IgA, IgA < 51.48)
nrow(IgA_out_rm)
#View(IgA_out_rm)
IgA_out_rm$out_rm <- IgA_out_rm$IgA
#View(IgA_out_rm)
IgA_out_rm <- within(IgA_out_rm, rm(IgA))
#View(IgAXallcovars_no_MSL_res_out_rm)
length(IgA$IgA)
length(IgA_out_rm$out_rm)
length(IgA$IgA) - length(IgA_out_rm$out_rm)
boxplot(IgA_out_rm$out_rm)
hist(IgA_out_rm$out_rm)
rug(IgA_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgA <- full_join(IgA,IgA_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgA)
#View(IgAXallcovars_no_MSL_res)

#Transformação boxcox
IgA$IgA.t<-boxcox.ols(IgA$IgA)
boxplot(IgA$IgA.t)
hist(IgA$IgA.t)
rug(IgA$IgA.t)
summary(IgA$IgA.t)

#Transformação logarítmica
IgA$LOG <- log10(IgA$IgA)
boxplot(IgA$LOG)
hist(IgA$LOG)
rug(IgA$LOG)
summary(IgA$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IgA)
hist(pheno$IgA)
model<-lm(IgA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IgA$IgA)
hist(model$residuals)
rug(model$residuals)
IgA[names(model$residuals),"IgAXallcovars_no_MSL_RES"] <- model$residuals
#View(IgA)

#Remoção dos outliers
hist(IgA$IgAXallcovars_no_MSL_RES)
boxplot(IgA$IgAXallcovars_no_MSL_RES)
IgA$IgAXallcovars_no_MSL_RES[IgA$IgAXallcovars_no_MSL_RES %in% boxplot.stats(IgA$IgAXallcovars_no_MSL_RES)$out]
length(IgA$IgAXallcovars_no_MSL_RES[IgA$IgAXallcovars_no_MSL_RES %in% boxplot.stats(IgA$IgAXallcovars_no_MSL_RES)$out])
IgAXallcovars_no_MSL_out_rm <- subset(IgA, IgAXallcovars_no_MSL_RES < 50.22969)
nrow(IgAXallcovars_no_MSL_out_rm)
IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm <- IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES
IgAXallcovars_no_MSL_out_rm <- within(IgAXallcovars_no_MSL_out_rm, rm(IgAXallcovars_no_MSL_RES))
length(IgA$IgAXallcovars_no_MSL_RES)
length(IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm)
length(IgA$IgAXallcovars_no_MSL_RES - IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm)
boxplot(IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm)
hist(IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm)
rug(IgAXallcovars_no_MSL_out_rm$IgAXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgA <- full_join(IgA,IgAXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgA)
View(IgA)


#Transformação boxcox
IgA$IgAXallcovars_no_MSL_RES.t<-boxcox.ols(IgA$IgAXallcovars_no_MSL_RES)
boxplot(IgA$IgAXallcovars_no_MSL_RES.t)
hist(IgA$IgAXallcovars_no_MSL_RES.t)
rug(IgA$IgAXallcovars_no_MSL_RES.t)
summary(IgA$IgAXallcovars_no_MSL_RES.t)

#Transformação logarítmica
IgA$IgAXallcovars_no_MSL_RES.LOG <- log10(IgA$IgAXallcovars_no_MSL_RES)
boxplot(IgA$IgAXallcovars_no_MSL_RES.LOG)
hist(IgA$IgAXallcovars_no_MSL_RES.LOG)
rug(IgA$IgAXallcovars_no_MSL_RES.LOG)
summary(IgA$IgAXallcovars_no_MSL_RES.LOG)

View(IgA)

write.table(IgA,"IgA_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IgG.UE para normalização
IgG <- select(pheno,IID,FID,IgG.UE)
#View(IgG.UE)

#Remoção dos outliers
hist(IgG$IgG.UE)
boxplot(IgG$IgG.UE)
# length(IgG$IgG.UE)
# IgG$IgG.UE[IgG$IgG.UE %in% boxplot.stats(IgG$IgG.UE)$out]
# length(IgG$IgG.UE[IgG$IgG.UE %in% boxplot.stats(IgG$IgG.UE)$out])
# IgG_out_rm <- subset(IgG, IgG.UE < 83.13)
# nrow(IgG_out_rm)
# #View(IgG_out_rm)
# IgG_out_rm$out_rm <- IgG_out_rm$IgG.UE
# #View(IgG_out_rm)
# IgG_out_rm <- within(IgG_out_rm, rm(IgG.UE))
# #View(IgG_out_rm)
# length(IgG$IgG.UE)
# length(IgG_out_rm$out_rm)
# length(IgG$IgG.UE) - length(IgG_out_rm$out_rm)
# boxplot(IgG_out_rm$out_rm)
# hist(IgG_out_rm$out_rm)
# rug(IgG_out_rm$out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# IgG.UE <- full_join(IgG,IgG_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(IgG)
# #View(IgG)

#Transformação boxcox
IgG$IgG.UE.t<-boxcox.ols(IgG$IgG.UE)
boxplot(IgG$IgG.UE.t)
hist(IgG$IgG.UE.t)
rug(IgG$IgG.UE.t)
summary(IgG$IgG.UE.t)

#Transformação logarítmica
IgG$LOG <- log10(IgG$IgG.UE)
boxplot(IgG$LOG)
hist(IgG$LOG)
rug(IgG$LOG)
summary(IgG$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IgG.UE)
hist(pheno$IgG.UE)
model<-lm(IgG.UE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IgG$IgG.UE)
hist(model$residuals)
rug(model$residuals)
IgG[names(model$residuals),"IgG.UEXallcovariates_RES"] <- model$residuals
#View(IgG)

#Remoção dos outliers
hist(IgG$IgG.UEXallcovariates_RES)
boxplot(IgG$IgG.UEXallcovariates_RES)
IgG$IgG.UEXallcovariates_RES[IgG$IgG.UEXallcovariates_RES %in% boxplot.stats(IgG$IgG.UEXallcovariates_RES)$out]
length(IgG$IgG.UEXallcovariates_RES[IgG$IgG.UEXallcovariates_RES %in% boxplot.stats(IgG$IgG.UEXallcovariates_RES)$out])
IgGXallcovariates_out_rm <- subset(IgG, IgG.UEXallcovariates_RES < 107.846)
nrow(IgGXallcovariates_out_rm)
IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm <- IgGXallcovariates_out_rm$IgG.UEXallcovariates_RES
IgGXallcovariates_out_rm <- within(IgGXallcovariates_out_rm, rm(IgG.UEXallcovariates_RES))
length(IgG$IgG.UEXallcovariates_RES)
length(IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm)
length(IgG$IgG.UEXallcovariates_RES - IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm)
boxplot(IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm)
hist(IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm)
rug(IgGXallcovariates_out_rm$IgGXallcovariates_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgG <- full_join(IgG,IgGXallcovariates_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgG)
View(IgG)


#Transformação boxcox
IgG$IgG.UEXallcovariates_RES.t<-boxcox.ols(IgG$IgG.UEXallcovariates_RES)
boxplot(IgG$IgG.UEXallcovariates_RES.t)
hist(IgG$IgG.UEXallcovariates_RES.t)
rug(IgG$IgG.UEXallcovariates_RES.t)
summary(IgG$IgG.UEXallcovariates_RES.t)

#Transformação logarítmica
IgG$IgG.UEXallcovariates_RES.LOG <- log10(IgG$IgG.UEXallcovariates_RES)
boxplot(IgG$IgG.UEXallcovariates_RES.LOG)
hist(IgG$IgG.UEXallcovariates_RES.LOG)
rug(IgG$IgG.UEXallcovariates_RES.LOG)
summary(IgG$IgG.UEXallcovariates_RES.LOG)

View(IgG)

write.table(IgG,"IgG.UE_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

#Criando um subset IgE para normalização
IgE <- select(pheno,IID,FID,IgE)
#View(IgE)

#Remoção dos outliers
hist(IgE$IgE)
boxplot(IgE$IgE)
length(IgE$IgE)
IgE$IgE[IgE$IgE %in% boxplot.stats(IgE$IgE)$out]
length(IgE$IgE[IgE$IgE %in% boxplot.stats(IgE$IgE)$out])
IgE_out_rm <- subset(IgE, IgE < 92.67)
nrow(IgE_out_rm)
#View(IgE_out_rm)
IgE_out_rm$out_rm <- IgE_out_rm$IgE
#View(IgE_out_rm)
IgE_out_rm <- within(IgE_out_rm, rm(IgE))
#View(IgEXallcovars_no_MSL_res_out_rm)
length(IgE$IgE)
length(IgE_out_rm$out_rm)
length(IgE$IgE) - length(IgE_out_rm$out_rm)
boxplot(IgE_out_rm$out_rm)
hist(IgE_out_rm$out_rm)
rug(IgE_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgE <- full_join(IgE,IgE_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgE)
#View(IgEXallcovars_no_MSL_res)

#Transformação boxcox
IgE$IgE.t<-boxcox.ols(IgE$IgE)
boxplot(IgE$IgE.t)
hist(IgE$IgE.t)
rug(IgE$IgE.t)
summary(IgE$IgE.t)

#Transformação logarítmica
IgE$LOG <- log10(IgE$IgE)
boxplot(IgE$LOG)
hist(IgE$LOG)
rug(IgE$LOG)
summary(IgE$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IgE)
hist(pheno$IgE)
model<-lm(IgE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IgE$IgE)
hist(model$residuals)
rug(model$residuals)
IgE[names(model$residuals),"IgEXallcovars_no_MSL_RES"] <- model$residuals
#View(IgE)

#Remoção dos outliers
hist(IgE$IgEXallcovars_no_MSL_RES)
boxplot(IgE$IgEXallcovars_no_MSL_RES)
IgE$IgEXallcovars_no_MSL_RES[IgE$IgEXallcovars_no_MSL_RES %in% boxplot.stats(IgE$IgEXallcovars_no_MSL_RES)$out]
length(IgE$IgEXallcovars_no_MSL_RES[IgE$IgEXallcovars_no_MSL_RES %in% boxplot.stats(IgE$IgEXallcovars_no_MSL_RES)$out])
IgEXallcovars_no_MSL_out_rm <- subset(IgE, IgEXallcovars_no_MSL_RES < 66.11224)
nrow(IgEXallcovars_no_MSL_out_rm)
IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm <- IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES
IgEXallcovars_no_MSL_out_rm <- within(IgEXallcovars_no_MSL_out_rm, rm(IgEXallcovars_no_MSL_RES))
length(IgE$IgEXallcovars_no_MSL_RES)
length(IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm)
length(IgE$IgEXallcovars_no_MSL_RES - IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm)
boxplot(IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm)
hist(IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm)
rug(IgEXallcovars_no_MSL_out_rm$IgEXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgE <- full_join(IgE,IgEXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgE)
View(IgE)


#Transformação boxcox
IgE$IgEXallcovars_no_MSL_RES.t<-boxcox.ols(IgE$IgEXallcovars_no_MSL_RES)
boxplot(IgE$IgEXallcovars_no_MSL_RES.t)
hist(IgE$IgEXallcovars_no_MSL_RES.t)
rug(IgE$IgEXallcovars_no_MSL_RES.t)
summary(IgE$IgEXallcovars_no_MSL_RES.t)

#Transformação logarítmica
IgE$IgEXallcovars_no_MSL_RES.LOG <- log10(IgE$IgEXallcovars_no_MSL_RES)
boxplot(IgE$IgEXallcovars_no_MSL_RES.LOG)
hist(IgE$IgEXallcovars_no_MSL_RES.LOG)
rug(IgE$IgEXallcovars_no_MSL_RES.LOG)
summary(IgE$IgEXallcovars_no_MSL_RES.LOG)

View(IgE)

write.table(IgE,"IgE_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IgM para normalização
IgM <- select(pheno,IID,FID,IgM)
#View(IgM)

#Remoção dos outliers
hist(IgM$IgM)
boxplot(IgM$IgM)
length(IgM$IgM)
IgM$IgM[IgM$IgM %in% boxplot.stats(IgM$IgM)$out]
length(IgM$IgM[IgM$IgM %in% boxplot.stats(IgM$IgM)$out])
IgM_out_rm <- subset(IgM, IgM < 162.50)
nrow(IgM_out_rm)
#View(IgM_out_rm)
IgM_out_rm$out_rm <- IgM_out_rm$IgM
#View(IgM_out_rm)
IgM_out_rm <- within(IgM_out_rm, rm(IgM))
#View(IgMXcollar_res_out_rm)
length(IgM$IgM)
length(IgM_out_rm$out_rm)
length(IgM$IgM) - length(IgM_out_rm$out_rm)
boxplot(IgM_out_rm$out_rm)
hist(IgM_out_rm$out_rm)
rug(IgM_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgM <- full_join(IgM,IgM_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgM)
#View(IgMXcollar_res)

#Transformação boxcox
IgM$IgM.t<-boxcox.ols(IgM$IgM)
boxplot(IgM$IgM.t)
hist(IgM$IgM.t)
rug(IgM$IgM.t)
summary(IgM$IgM.t)

#Transformação logarítmica
IgM$LOG <- log10(IgM$IgM)
boxplot(IgM$LOG)
hist(IgM$LOG)
rug(IgM$LOG)
summary(IgM$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IgM)
hist(pheno$IgM)
model<-lm(IgM ~ COLLAR,data = pheno)
summary(model)
length(model$residuals)
length(IgM$IgM)
hist(model$residuals)
rug(model$residuals)
IgM[names(model$residuals),"IgMXcollar_RES"] <- model$residuals
#View(IgM)

#Remoção dos outliers
hist(IgM$IgMXcollar_RES)
boxplot(IgM$IgMXcollar_RES)
IgM$IgMXcollar_RES[IgM$IgMXcollar_RES %in% boxplot.stats(IgM$IgMXcollar_RES)$out]
length(IgM$IgMXcollar_RES[IgM$IgMXcollar_RES %in% boxplot.stats(IgM$IgMXcollar_RES)$out])
IgMXcollar_out_rm <- subset(IgM, IgMXcollar_RES < 103.1916)
nrow(IgMXcollar_out_rm)
IgMXcollar_out_rm$IgMXcollar_RES_out_rm <- IgMXcollar_out_rm$IgMXcollar_RES
IgMXcollar_out_rm <- within(IgMXcollar_out_rm, rm(IgMXcollar_RES))
length(IgM$IgMXcollar_RES)
length(IgMXcollar_out_rm$IgMXcollar_RES_out_rm)
length(IgM$IgMXcollar_RES - IgMXcollar_out_rm$IgMXcollar_RES_out_rm)
boxplot(IgMXcollar_out_rm$IgMXcollar_RES_out_rm)
hist(IgMXcollar_out_rm$IgMXcollar_RES_out_rm)
rug(IgMXcollar_out_rm$IgMXcollar_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgM <- full_join(IgM,IgMXcollar_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgM)
View(IgM)


#Transformação boxcox
IgM$IgMXcollar_RES.t<-boxcox.ols(IgM$IgMXcollar_RES)
boxplot(IgM$IgMXcollar_RES.t)
hist(IgM$IgMXcollar_RES.t)
rug(IgM$IgMXcollar_RES.t)
summary(IgM$IgMXcollar_RES.t)

#Transformação logarítmica
IgM$IgMXcollar_RES.LOG <- log10(IgM$IgMXcollar_RES)
boxplot(IgM$IgMXcollar_RES.LOG)
hist(IgM$IgMXcollar_RES.LOG)
rug(IgM$IgMXcollar_RES.LOG)
summary(IgM$IgMXcollar_RES.LOG)

View(IgM)

write.table(IgM,"IgM_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IgG.SALIVA para normalização
IgG.SALIVA <- select(pheno,IID,FID,IgG.SALIVA)
#View(IgG.SALIVA)

#Remoção dos outliers
hist(IgG.SALIVA$IgG.SALIVA)
boxplot(IgG.SALIVA$IgG.SALIVA)
length(IgG.SALIVA$IgG.SALIVA)
IgG.SALIVA$IgG.SALIVA[IgG.SALIVA$IgG.SALIVA %in% boxplot.stats(IgG.SALIVA$IgG.SALIVA)$out]
length(IgG.SALIVA$IgG.SALIVA[IgG.SALIVA$IgG.SALIVA %in% boxplot.stats(IgG.SALIVA$IgG.SALIVA)$out])
IgG.SALIVA_out_rm <- subset(IgG.SALIVA, IgG.SALIVA < 109.23)
nrow(IgG.SALIVA_out_rm)
#View(IgG.SALIVA_out_rm)
IgG.SALIVA_out_rm$out_rm <- IgG.SALIVA_out_rm$IgG.SALIVA
#View(IgG.SALIVA_out_rm)
IgG.SALIVA_out_rm <- within(IgG.SALIVA_out_rm, rm(IgG.SALIVA))
#View(IgG.SALIVAXtreatment_res_out_rm)
length(IgG.SALIVA$IgG.SALIVA)
length(IgG.SALIVA_out_rm$out_rm)
length(IgG.SALIVA$IgG.SALIVA) - length(IgG.SALIVA_out_rm$out_rm)
boxplot(IgG.SALIVA_out_rm$out_rm)
hist(IgG.SALIVA_out_rm$out_rm)
rug(IgG.SALIVA_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgG.SALIVA <- full_join(IgG.SALIVA,IgG.SALIVA_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgG.SALIVA)
#View(IgG.SALIVAXtreatment_res)

#Transformação boxcox
IgG.SALIVA$IgG.SALIVA.t<-boxcox.ols(IgG.SALIVA$IgG.SALIVA)
boxplot(IgG.SALIVA$IgG.SALIVA.t)
hist(IgG.SALIVA$IgG.SALIVA.t)
rug(IgG.SALIVA$IgG.SALIVA.t)
summary(IgG.SALIVA$IgG.SALIVA.t)

#Transformação logarítmica
IgG.SALIVA$LOG <- log10(IgG.SALIVA$IgG.SALIVA)
boxplot(IgG.SALIVA$LOG)
hist(IgG.SALIVA$LOG)
rug(IgG.SALIVA$LOG)
summary(IgG.SALIVA$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IgG.SALIVA)
hist(pheno$IgG.SALIVA)
model<-lm(IgG.SALIVA ~ TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IgG.SALIVA$IgG.SALIVA)
hist(model$residuals)
rug(model$residuals)
IgG.SALIVA[names(model$residuals),"IgG.SALIVAXtreatment_RES"] <- model$residuals
#View(IgG.SALIVA)

#Remoção dos outliers
hist(IgG.SALIVA$IgG.SALIVAXtreatment_RES)
boxplot(IgG.SALIVA$IgG.SALIVAXtreatment_RES)
IgG.SALIVA$IgG.SALIVAXtreatment_RES[IgG.SALIVA$IgG.SALIVAXtreatment_RES %in% boxplot.stats(IgG.SALIVA$IgG.SALIVAXtreatment_RES)$out]
length(IgG.SALIVA$IgG.SALIVAXtreatment_RES[IgG.SALIVA$IgG.SALIVAXtreatment_RES %in% boxplot.stats(IgG.SALIVA$IgG.SALIVAXtreatment_RES)$out])
IgG.SALIVAXtreatment_out_rm <- subset(IgG.SALIVA, IgG.SALIVAXtreatment_RES < 90.20486)
nrow(IgG.SALIVAXtreatment_out_rm)
IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm <- IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES
IgG.SALIVAXtreatment_out_rm <- within(IgG.SALIVAXtreatment_out_rm, rm(IgG.SALIVAXtreatment_RES))
length(IgG.SALIVA$IgG.SALIVAXtreatment_RES)
length(IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm)
length(IgG.SALIVA$IgG.SALIVAXtreatment_RES - IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm)
boxplot(IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm)
hist(IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm)
rug(IgG.SALIVAXtreatment_out_rm$IgG.SALIVAXtreatment_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IgG.SALIVA <- full_join(IgG.SALIVA,IgG.SALIVAXtreatment_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IgG.SALIVA)
View(IgG.SALIVA)


#Transformação boxcox
IgG.SALIVA$IgG.SALIVAXtreatment_RES.t<-boxcox.ols(IgG.SALIVA$IgG.SALIVAXtreatment_RES)
boxplot(IgG.SALIVA$IgG.SALIVAXtreatment_RES.t)
hist(IgG.SALIVA$IgG.SALIVAXtreatment_RES.t)
rug(IgG.SALIVA$IgG.SALIVAXtreatment_RES.t)
summary(IgG.SALIVA$IgG.SALIVAXtreatment_RES.t)

#Transformação logarítmica
IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG <- log10(IgG.SALIVA$IgG.SALIVAXtreatment_RES)
boxplot(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)
hist(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)
rug(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)
summary(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)

IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG[IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG %in% boxplot.stats(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)$out]
length(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG[IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG %in% boxplot.stats(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)$out])

# length(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)
# IgG.SALIVA(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG[IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG <= -1.0863598]) <- "NA"
# nrow(IgG.SALIVA$IgG.SALIVAXtreatment_RES.LOG)

View(IgG.SALIVA)


write.table(IgG.SALIVA,"IgG.SALIVA_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset INDURATION para normalização
INDURATION <- select(pheno,IID,FID,IDRM.VALUE..mm2.)
#View(INDURATION)

# #Remoção dos outliers
# hist(INDURATION$IDRM.VALUE..mm2.)
# boxplot(INDURATION$IDRM.VALUE..mm2.)
# length(INDURATION$IDRM.VALUE..mm2.)
# INDURATION$IDRM.VALUE..mm2.[INDURATION$IDRM.VALUE..mm2. %in% boxplot.stats(INDURATION$IDRM.VALUE..mm2.)$out]
# length(INDURATION$IDRM.VALUE..mm2.[INDURATION$IDRM.VALUE..mm2. %in% boxplot.stats(INDURATION$IDRM.VALUE..mm2.)$out])
# INDURATION_out_rm <- subset(INDURATION, IDRM.VALUE..mm2. < 25)
# nrow(INDURATION_out_rm)
# #View(INDURATION_out_rm)
# INDURATION_out_rm$out_rm <- INDURATION_out_rm$IDRM.VALUE..mm2.
# #View(INDURATION_out_rm)
# INDURATION_out_rm <- within(INDURATION_out_rm, rm(IDRM.VALUE..mm2.))
# #View(INDURATIONXallcovars_no_MSL_res_out_rm)
# length(INDURATION$IDRM.VALUE..mm2.)
# length(INDURATION_out_rm$out_rm)
# length(INDURATION$IDRM.VALUE..mm2.) - length(INDURATION_out_rm$out_rm)
# boxplot(INDURATION_out_rm$out_rm)
# hist(INDURATION_out_rm$out_rm)
# rug(INDURATION_out_rm$out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# INDURATION <- full_join(INDURATION,INDURATION_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(INDURATION)
# #View(INDURATIONXallcovars_no_MSL_res)

#Transformação boxcox
INDURATION$IDRM.VALUE..mm2..t<-boxcox.ols(INDURATION$IDRM.VALUE..mm2.)
boxplot(INDURATION$IDRM.VALUE..mm2..t)
hist(INDURATION$IDRM.VALUE..mm2..t)
rug(INDURATION$IDRM.VALUE..mm2..t)
summary(INDURATION$IDRM.VALUE..mm2..t)

#Transformação logarítmica
INDURATION$LOG <- log10(INDURATION$IDRM.VALUE..mm2.)
boxplot(INDURATION$LOG)
hist(INDURATION$LOG)
rug(INDURATION$LOG)
summary(INDURATION$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IDRM.VALUE..mm2.)
hist(pheno$IDRM.VALUE..mm2.)
model<-lm(IDRM.VALUE..mm2. ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(INDURATION$IDRM.VALUE..mm2.)
hist(model$residuals)
rug(model$residuals)
INDURATION[names(model$residuals),"INDURATIONXallcovars_no_MSL_RES"] <- model$residuals
#View(INDURATION)

#Remoção dos outliers
hist(INDURATION$INDURATIONXallcovars_no_MSL_RES)
boxplot(INDURATION$INDURATIONXallcovars_no_MSL_RES)
INDURATION$INDURATIONXallcovars_no_MSL_RES[INDURATION$INDURATIONXallcovars_no_MSL_RES %in% boxplot.stats(INDURATION$INDURATIONXallcovars_no_MSL_RES)$out]
length(INDURATION$INDURATIONXallcovars_no_MSL_RES[INDURATION$INDURATIONXallcovars_no_MSL_RES %in% boxplot.stats(INDURATION$INDURATIONXallcovars_no_MSL_RES)$out])
INDURATIONXallcovars_no_MSL_out_rm <- subset(INDURATION, INDURATIONXallcovars_no_MSL_RES < 75.85043)
nrow(INDURATIONXallcovars_no_MSL_out_rm)
INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm <- INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES
INDURATIONXallcovars_no_MSL_out_rm <- within(INDURATIONXallcovars_no_MSL_out_rm, rm(INDURATIONXallcovars_no_MSL_RES))
length(INDURATION$INDURATIONXallcovars_no_MSL_RES)
length(INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm)
length(INDURATION$INDURATIONXallcovars_no_MSL_RES - INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm)
boxplot(INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm)
hist(INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm)
rug(INDURATIONXallcovars_no_MSL_out_rm$INDURATIONXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
INDURATION <- full_join(INDURATION,INDURATIONXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(INDURATION)
View(INDURATION)

#Transformação boxcox
INDURATION$INDURATIONXallcovars_no_MSL_RES.t<-boxcox.ols(INDURATION$INDURATIONXallcovars_no_MSL_RES)
boxplot(INDURATION$INDURATIONXallcovars_no_MSL_RES.t)
hist(INDURATION$INDURATIONXallcovars_no_MSL_RES.t)
rug(INDURATION$INDURATIONXallcovars_no_MSL_RES.t)
summary(INDURATION$INDURATIONXallcovars_no_MSL_RES.t)

#Transformação logarítmica
INDURATION$INDURATIONXallcovars_no_MSL_RES.LOG <- log10(INDURATION$INDURATIONXallcovars_no_MSL_RES)
boxplot(INDURATION$INDURATIONXallcovars_no_MSL_RES.LOG)
hist(INDURATION$INDURATIONXallcovars_no_MSL_RES.LOG)
rug(INDURATION$INDURATIONXallcovars_no_MSL_RES.LOG)
summary(INDURATION$INDURATIONXallcovars_no_MSL_RES.LOG)

View(INDURATION)

write.table(INDURATION,"INDURATION_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)



#Criando um subset PI para normalização
PI <- select(pheno,IID,FID,IP.NN)
#View(PI)

#Remoção dos outliers
hist(PI$IP.NN)
boxplot(PI$IP.NN)
length(PI$IP.NN)
PI$IP.NN[PI$IP.NN %in% boxplot.stats(PI$IP.NN)$out]
length(PI$IP.NN[PI$IP.NN %in% boxplot.stats(PI$IP.NN)$out])
PI_out_rm <- subset(PI, IP.NN < 148.25)
nrow(PI_out_rm)
#View(PI_out_rm)
PI_out_rm$out_rm <- PI_out_rm$IP.NN
#View(PI_out_rm)
PI_out_rm <- within(PI_out_rm, rm(IP.NN))
#View(PIXallcovars_no_MSL_res_out_rm)
length(PI$IP.NN)
length(PI_out_rm$out_rm)
length(PI$IP.NN) - length(PI_out_rm$out_rm)
boxplot(PI_out_rm$out_rm)
hist(PI_out_rm$out_rm)
rug(PI_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
PI <- full_join(PI,PI_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(PI)
#View(PIXallcovars_no_MSL_res)

# #Transformação boxcox
# PI$IP.NN.t<-boxcox.ols(PI$IP.NN)
# boxplot(PI$IP.NN.t)
# hist(PI$IP.NN.t)
# rug(PI$IP.NN.t)
# summary(PI$IP.NN.t)

#Transformação logarítmica
PI$LOG <- log10(PI$IP.NN)
boxplot(PI$LOG)
hist(PI$LOG)
rug(PI$LOG)
summary(PI$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IP.NN)
hist(pheno$IP.NN)
model<-lm(IP.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(PI$IP.NN)
hist(model$residuals)
rug(model$residuals)
PI[names(model$residuals),"PIXallcovars_no_MSL_RES"] <- model$residuals
View(PI)

#Remoção dos outliers
hist(PI$PIXallcovars_no_MSL_RES)
boxplot(PI$PIXallcovars_no_MSL_RES)
PI$PIXallcovars_no_MSL_RES[PI$PIXallcovars_no_MSL_RES %in% boxplot.stats(PI$PIXallcovars_no_MSL_RES)$out]
length(PI$PIXallcovars_no_MSL_RES[PI$PIXallcovars_no_MSL_RES %in% boxplot.stats(PI$PIXallcovars_no_MSL_RES)$out])
PIXallcovars_no_MSL_out_rm <- subset(PI, PIXallcovars_no_MSL_RES < 98.34209)
nrow(PIXallcovars_no_MSL_out_rm)
PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm <- PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES
PIXallcovars_no_MSL_out_rm <- within(PIXallcovars_no_MSL_out_rm, rm(PIXallcovars_no_MSL_RES))
length(PI$PIXallcovars_no_MSL_RES)
length(PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm)
length(PI$PIXallcovars_no_MSL_RES - PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm)
boxplot(PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm)
hist(PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm)
rug(PIXallcovars_no_MSL_out_rm$PIXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
PI <- full_join(PI,PIXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(PI)
View(PI)


#Transformação boxcox
PI$PIXallcovars_no_MSL_RES.t<-boxcox.ols(PI$PIXallcovars_no_MSL_RES)
boxplot(PI$PIXallcovars_no_MSL_RES.t)
hist(PI$PIXallcovars_no_MSL_RES.t)
rug(PI$PIXallcovars_no_MSL_RES.t)
summary(PI$PIXallcovars_no_MSL_RES.t)

#Transformação logarítmica
PI$PIXallcovars_no_MSL_RES_LOG <- log10(PI$PIXallcovars_no_MSL_RES)
boxplot(PI$PIXallcovars_no_MSL_RES_LOG)
hist(PI$PIXallcovars_no_MSL_RES_LOG)
rug(PI$PIXallcovars_no_MSL_RES_LOG)
summary(PI$PIXallcovars_no_MSL_RES_LOG)

View(PI)

write.table(PI,"PI_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IFNG para normalização
IFNG <- select(pheno,IID,FID,IFN.NN)
#View(IFNG)

#Remoção dos outliers
hist(IFNG$IFN.NN)
boxplot(IFNG$IFN.NN)
length(IFNG$IFN.NN)
IFNG$IFN.NN[IFNG$IFN.NN %in% boxplot.stats(IFNG$IFN.NN)$out]
length(IFNG$IFN.NN[IFNG$IFN.NN %in% boxplot.stats(IFNG$IFN.NN)$out])
IFNG_out_rm <- subset(IFNG, IFN.NN < 1186.05)
nrow(IFNG_out_rm)
#View(IFNG_out_rm)
IFNG_out_rm$out_rm <- IFNG_out_rm$IFN.NN
#View(IFNG_out_rm)
IFNG_out_rm <- within(IFNG_out_rm, rm(IFN.NN))
#View(IFNGXvaccine_res_out_rm)
length(IFNG$IFN.NN)
length(IFNG_out_rm$out_rm)
length(IFNG$IFN.NN) - length(IFNG_out_rm$out_rm)
boxplot(IFNG_out_rm$out_rm)
hist(IFNG_out_rm$out_rm)
rug(IFNG_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IFNG <- full_join(IFNG,IFNG_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IFNG)
#View(IFNGXvaccine_res)

# #Transformação boxcox
# IFNG$IFN.NN.t<-boxcox.ols(IFNG$IFN.NN)
# boxplot(IFNG$IFN.NN.t)
# hist(IFNG$IFN.NN.t)
# rug(IFNG$IFN.NN.t)
# summary(IFNG$IFN.NN.t)

#Transformação logarítmica
IFNG$LOG <- log10(IFNG$IFN.NN)
boxplot(IFNG$LOG)
hist(IFNG$LOG)
rug(IFNG$LOG)
summary(IFNG$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IFN.NN)
hist(pheno$IFN.NN)
model<-lm(IFN.NN ~ VACCINE,data = pheno)
summary(model)
length(model$residuals)
length(IFNG$IFN.NN)
hist(model$residuals)
rug(model$residuals)
IFNG[names(model$residuals),"IFNGXvaccine_RES"] <- model$residuals
View(IFNG)

#Remoção dos outliers
hist(IFNG$IFNGXvaccine_RES)
boxplot(IFNG$IFNGXvaccine_RES)
IFNG$IFNGXvaccine_RES[IFNG$IFNGXvaccine_RES %in% boxplot.stats(IFNG$IFNGXvaccine_RES)$out]
length(IFNG$IFNGXvaccine_RES[IFNG$IFNGXvaccine_RES %in% boxplot.stats(IFNG$IFNGXvaccine_RES)$out])
IFNGXvaccine_out_rm <- subset(IFNG, IFNGXvaccine_RES < 253.9045 & IFNGXvaccine_RES > -1380.0429)
nrow(IFNGXvaccine_out_rm)
IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm <- IFNGXvaccine_out_rm$IFNGXvaccine_RES
IFNGXvaccine_out_rm <- within(IFNGXvaccine_out_rm, rm(IFNGXvaccine_RES))
length(IFNG$IFNGXvaccine_RES)
length(IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm)
length(IFNG$IFNGXvaccine_RES - IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm)
boxplot(IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm)
hist(IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm)
rug(IFNGXvaccine_out_rm$IFNGXvaccine_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IFNG <- full_join(IFNG,IFNGXvaccine_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IFNG)
View(IFNG)


#Transformação boxcox
IFNG$IFNGXvaccine_RES.t<-boxcox.ols(IFNG$IFNGXvaccine_RES)
boxplot(IFNG$IFNGXvaccine_RES.t)
hist(IFNG$IFNGXvaccine_RES.t)
rug(IFNG$IFNGXvaccine_RES.t)
summary(IFNG$IFNGXvaccine_RES.t)

#Transformação logarítmica
IFNG$IFNGXvaccine_RES_LOG <- log10(IFNG$IFNGXvaccine_RES)
boxplot(IFNG$IFNGXvaccine_RES_LOG)
hist(IFNG$IFNGXvaccine_RES_LOG)
rug(IFNG$IFNGXvaccine_RES_LOG)
summary(IFNG$IFNGXvaccine_RES_LOG)

View(IFNG)

write.table(IFNG,"IFNG_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset TNFA para normalização
TNFA <- select(pheno,IID,FID,TNF.NN)
#View(TNFA)

#Remoção dos outliers
hist(TNFA$TNF.NN)
boxplot(TNFA$TNF.NN)
length(TNFA$TNF.NN)
TNFA$TNF.NN[TNFA$TNF.NN %in% boxplot.stats(TNFA$TNF.NN)$out]
length(TNFA$TNF.NN[TNFA$TNF.NN %in% boxplot.stats(TNFA$TNF.NN)$out])
TNFA_out_rm <- subset(TNFA, TNF.NN < 70.79)
nrow(TNFA_out_rm)
#View(TNFA_out_rm)
TNFA_out_rm$out_rm <- TNFA_out_rm$TNF.NN
#View(TNFA_out_rm)
TNFA_out_rm <- within(TNFA_out_rm, rm(TNF.NN))
#View(TNFAXallcovars_no_MSL_res_out_rm)
length(TNFA$TNF.NN)
length(TNFA_out_rm$out_rm)
length(TNFA$TNF.NN) - length(TNFA_out_rm$out_rm)
boxplot(TNFA_out_rm$out_rm)
hist(TNFA_out_rm$out_rm)
rug(TNFA_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TNFA <- full_join(TNFA,TNFA_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TNFA)
#View(TNFAXallcovars_no_MSL_res)

#Transformação boxcox
TNFA$TNF.NN.t<-boxcox.ols(TNFA$TNF.NN)
boxplot(TNFA$TNF.NN.t)
hist(TNFA$TNF.NN.t)
rug(TNFA$TNF.NN.t)
summary(TNFA$TNF.NN.t)

#Transformação logarítmica
TNFA$LOG <- log10(TNFA$TNF.NN)
boxplot(TNFA$LOG)
hist(TNFA$LOG)
rug(TNFA$LOG)
summary(TNFA$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$TNF.NN)
hist(pheno$TNF.NN)
model<-lm(TNF.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(TNFA$TNF.NN)
hist(model$residuals)
rug(model$residuals)
TNFA[names(model$residuals),"TNFAXallcovars_no_MSL_RES"] <- model$residuals
View(TNFA)

#Remoção dos outliers
hist(TNFA$TNFAXallcovars_no_MSL_RES)
boxplot(TNFA$TNFAXallcovars_no_MSL_RES)
TNFA$TNFAXallcovars_no_MSL_RES[TNFA$TNFAXallcovars_no_MSL_RES %in% boxplot.stats(TNFA$TNFAXallcovars_no_MSL_RES)$out]
length(TNFA$TNFAXallcovars_no_MSL_RES[TNFA$TNFAXallcovars_no_MSL_RES %in% boxplot.stats(TNFA$TNFAXallcovars_no_MSL_RES)$out])
TNFAXallcovars_no_MSL_out_rm <- subset(TNFA, TNFAXallcovars_no_MSL_RES < 73.13531)
nrow(TNFAXallcovars_no_MSL_out_rm)
TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm <- TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES
TNFAXallcovars_no_MSL_out_rm <- within(TNFAXallcovars_no_MSL_out_rm, rm(TNFAXallcovars_no_MSL_RES))
length(TNFA$TNFAXallcovars_no_MSL_RES)
length(TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm)
length(TNFA$TNFAXallcovars_no_MSL_RES - TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm)
boxplot(TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm)
hist(TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm)
rug(TNFAXallcovars_no_MSL_out_rm$TNFAXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TNFA <- full_join(TNFA,TNFAXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TNFA)
View(TNFA)


#Transformação boxcox
TNFA$TNFAXallcovars_no_MSL_RES.t<-boxcox.ols(TNFA$TNFAXallcovars_no_MSL_RES)
boxplot(TNFA$TNFAXallcovars_no_MSL_RES.t)
hist(TNFA$TNFAXallcovars_no_MSL_RES.t)
rug(TNFA$TNFAXallcovars_no_MSL_RES.t)
summary(TNFA$TNFAXallcovars_no_MSL_RES.t)

#Transformação logarítmica
TNFA$TNFAXallcovars_no_MSL_RES_LOG <- log10(TNFA$TNFAXallcovars_no_MSL_RES)
boxplot(TNFA$TNFAXallcovars_no_MSL_RES_LOG)
hist(TNFA$TNFAXallcovars_no_MSL_RES_LOG)
rug(TNFA$TNFAXallcovars_no_MSL_RES_LOG)
summary(TNFA$TNFAXallcovars_no_MSL_RES_LOG)

View(TNFA)

write.table(TNFA,"TNFA_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset IL10 para normalização
IL10 <- select(pheno,IID,FID,IL10.NN)
#View(IL10)

#Remoção dos outliers
hist(IL10$IL10.NN)
boxplot(IL10$IL10.NN)
length(IL10$IL10.NN)
IL10$IL10.NN[IL10$IL10.NN %in% boxplot.stats(IL10$IL10.NN)$out]
length(IL10$IL10.NN[IL10$IL10.NN %in% boxplot.stats(IL10$IL10.NN)$out])
IL10_out_rm <- subset(IL10, IL10.NN < 72.44)
nrow(IL10_out_rm)
#View(IL10_out_rm)
IL10_out_rm$out_rm <- IL10_out_rm$IL10.NN
#View(IL10_out_rm)
IL10_out_rm <- within(IL10_out_rm, rm(IL10.NN))
#View(IL10Xallcovars_no_MSL_res_out_rm)
length(IL10$IL10.NN)
length(IL10_out_rm$out_rm)
length(IL10$IL10.NN) - length(IL10_out_rm$out_rm)
boxplot(IL10_out_rm$out_rm)
hist(IL10_out_rm$out_rm)
rug(IL10_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL10 <- full_join(IL10,IL10_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL10)
#View(IL10Xallcovars_no_MSL_res)

# #Transformação boxcox
# IL10$IL10.NN.t<-boxcox.ols(IL10$IL10.NN)
# boxplot(IL10$IL10.NN.t)
# hist(IL10$IL10.NN.t)
# rug(IL10$IL10.NN.t)
# summary(IL10$IL10.NN.t)

#Transformação logarítmica
IL10$LOG <- log10(IL10$IL10.NN)
boxplot(IL10$LOG)
hist(IL10$LOG)
rug(IL10$LOG)
summary(IL10$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IL10.NN)
hist(pheno$IL10.NN)
model<-lm(IL10.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IL10$IL10.NN)
hist(model$residuals)
rug(model$residuals)
IL10[names(model$residuals),"IL10Xallcovars_no_MSL_RES"] <- model$residuals
View(IL10)

#Remoção dos outliers
hist(IL10$IL10Xallcovars_no_MSL_RES)
boxplot(IL10$IL10Xallcovars_no_MSL_RES)
IL10$IL10Xallcovars_no_MSL_RES[IL10$IL10Xallcovars_no_MSL_RES %in% boxplot.stats(IL10$IL10Xallcovars_no_MSL_RES)$out]
length(IL10$IL10Xallcovars_no_MSL_RES[IL10$IL10Xallcovars_no_MSL_RES %in% boxplot.stats(IL10$IL10Xallcovars_no_MSL_RES)$out])
IL10Xallcovars_no_MSL_out_rm <- subset(IL10, IL10Xallcovars_no_MSL_RES < 62.47707)
nrow(IL10Xallcovars_no_MSL_out_rm)
IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm <- IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES
IL10Xallcovars_no_MSL_out_rm <- within(IL10Xallcovars_no_MSL_out_rm, rm(IL10Xallcovars_no_MSL_RES))
length(IL10$IL10Xallcovars_no_MSL_RES)
length(IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm)
length(IL10$IL10Xallcovars_no_MSL_RES - IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm)
boxplot(IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm)
hist(IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm)
rug(IL10Xallcovars_no_MSL_out_rm$IL10Xallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL10 <- full_join(IL10,IL10Xallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL10)
View(IL10)


#Transformação boxcox
IL10$IL10Xallcovars_no_MSL_RES.t<-boxcox.ols(IL10$IL10Xallcovars_no_MSL_RES)
boxplot(IL10$IL10Xallcovars_no_MSL_RES.t)
hist(IL10$IL10Xallcovars_no_MSL_RES.t)
rug(IL10$IL10Xallcovars_no_MSL_RES.t)
summary(IL10$IL10Xallcovars_no_MSL_RES.t)

#Transformação logarítmica
IL10$IL10Xallcovars_no_MSL_RES_LOG <- log10(IL10$IL10Xallcovars_no_MSL_RES)
boxplot(IL10$IL10Xallcovars_no_MSL_RES_LOG)
hist(IL10$IL10Xallcovars_no_MSL_RES_LOG)
rug(IL10$IL10Xallcovars_no_MSL_RES_LOG)
summary(IL10$IL10Xallcovars_no_MSL_RES_LOG)

View(IL10)

write.table(IL10,"IL10_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)



#Criando um subset IL4 para normalização
IL4 <- select(pheno,IID,FID,IL4.NN)
#View(IL4)

#Remoção dos outliers
hist(IL4$IL4.NN)
boxplot(IL4$IL4.NN)
length(IL4$IL4.NN)
IL4$IL4.NN[IL4$IL4.NN %in% boxplot.stats(IL4$IL4.NN)$out]
length(IL4$IL4.NN[IL4$IL4.NN %in% boxplot.stats(IL4$IL4.NN)$out])
IL4_out_rm <- subset(IL4, IL4.NN < 391.69)
nrow(IL4_out_rm)
#View(IL4_out_rm)
IL4_out_rm$out_rm <- IL4_out_rm$IL4.NN
#View(IL4_out_rm)
IL4_out_rm <- within(IL4_out_rm, rm(IL4.NN))
#View(IL4Xallcovars_no_MSL_res_out_rm)
length(IL4$IL4.NN)
length(IL4_out_rm$out_rm)
length(IL4$IL4.NN) - length(IL4_out_rm$out_rm)
boxplot(IL4_out_rm$out_rm)
hist(IL4_out_rm$out_rm)
rug(IL4_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL4 <- full_join(IL4,IL4_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL4)
#View(IL4Xallcovars_no_MSL_res)

#Transformação boxcox
IL4$IL4.NN.t<-boxcox.ols(IL4$IL4.NN)
boxplot(IL4$IL4.NN.t)
hist(IL4$IL4.NN.t)
rug(IL4$IL4.NN.t)
summary(IL4$IL4.NN.t)

#Transformação logarítmica
IL4$LOG <- log10(IL4$IL4.NN)
boxplot(IL4$LOG)
hist(IL4$LOG)
rug(IL4$LOG)
summary(IL4$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$IL4.NN)
hist(pheno$IL4.NN)
model<-lm(IL4.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(IL4$IL4.NN)
hist(model$residuals)
rug(model$residuals)
IL4[names(model$residuals),"IL4Xallcovars_no_MSL_RES"] <- model$residuals
View(IL4)

#Remoção dos outliers
hist(IL4$IL4Xallcovars_no_MSL_RES)
boxplot(IL4$IL4Xallcovars_no_MSL_RES)
IL4$IL4Xallcovars_no_MSL_RES[IL4$IL4Xallcovars_no_MSL_RES %in% boxplot.stats(IL4$IL4Xallcovars_no_MSL_RES)$out]
length(IL4$IL4Xallcovars_no_MSL_RES[IL4$IL4Xallcovars_no_MSL_RES %in% boxplot.stats(IL4$IL4Xallcovars_no_MSL_RES)$out])
IL4Xallcovars_no_MSL_out_rm <- subset(IL4, IL4Xallcovars_no_MSL_RES < 210.0491)
nrow(IL4Xallcovars_no_MSL_out_rm)
IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm <- IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES
IL4Xallcovars_no_MSL_out_rm <- within(IL4Xallcovars_no_MSL_out_rm, rm(IL4Xallcovars_no_MSL_RES))
length(IL4$IL4Xallcovars_no_MSL_RES)
length(IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm)
length(IL4$IL4Xallcovars_no_MSL_RES - IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm)
boxplot(IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm)
hist(IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm)
rug(IL4Xallcovars_no_MSL_out_rm$IL4Xallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
IL4 <- full_join(IL4,IL4Xallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(IL4)
View(IL4)


#Transformação boxcox
IL4$IL4Xallcovars_no_MSL_RES.t<-boxcox.ols(IL4$IL4Xallcovars_no_MSL_RES)
boxplot(IL4$IL4Xallcovars_no_MSL_RES.t)
hist(IL4$IL4Xallcovars_no_MSL_RES.t)
rug(IL4$IL4Xallcovars_no_MSL_RES.t)
summary(IL4$IL4Xallcovars_no_MSL_RES.t)

#Transformação logarítmica
IL4$IL4Xallcovars_no_MSL_RES_LOG <- log10(IL4$IL4Xallcovars_no_MSL_RES)
boxplot(IL4$IL4Xallcovars_no_MSL_RES_LOG)
hist(IL4$IL4Xallcovars_no_MSL_RES_LOG)
rug(IL4$IL4Xallcovars_no_MSL_RES_LOG)
summary(IL4$IL4Xallcovars_no_MSL_RES_LOG)

View(IL4)

write.table(IL4,"IL4_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset TGFB para normalização
TGFB <- select(pheno,IID,FID,TGF.VALUE)
#View(TGFB)

#Remoção dos outliers
hist(TGFB$TGF.VALUE)
boxplot(TGFB$TGF.VALUE)
length(TGFB$TGF.VALUE)
TGFB$TGF.VALUE[TGFB$TGF.VALUE %in% boxplot.stats(TGFB$TGF.VALUE)$out]
length(TGFB$TGF.VALUE[TGFB$TGF.VALUE %in% boxplot.stats(TGFB$TGF.VALUE)$out])
TGFB_out_rm <- subset(TGFB, TGF.VALUE < 148.25)
nrow(TGFB_out_rm)
#View(TGFB_out_rm)
TGFB_out_rm$out_rm <- TGFB_out_rm$TGF.VALUE
#View(TGFB_out_rm)
TGFB_out_rm <- within(TGFB_out_rm, rm(TGF.VALUE))
#View(TGFBXallcovars_no_MSL_res_out_rm)
length(TGFB$TGF.VALUE)
length(TGFB_out_rm$out_rm)
length(TGFB$TGF.VALUE) - length(TGFB_out_rm$out_rm)
boxplot(TGFB_out_rm$out_rm)
hist(TGFB_out_rm$out_rm)
rug(TGFB_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TGFB <- full_join(TGFB,TGFB_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TGFB)
#View(TGFBXallcovars_no_MSL_res)

# #Transformação boxcox
# TGFB$TGF.VALUE.t<-boxcox.ols(TGFB$TGF.VALUE)
# boxplot(TGFB$TGF.VALUE.t)
# hist(TGFB$TGF.VALUE.t)
# rug(TGFB$TGF.VALUE.t)
# summary(TGFB$TGF.VALUE.t)

#Transformação logarítmica
TGFB$LOG <- log10(TGFB$TGF.VALUE)
boxplot(TGFB$LOG)
hist(TGFB$LOG)
rug(TGFB$LOG)
summary(TGFB$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$TGF.VALUE)
hist(pheno$TGF.VALUE)
model<-lm(TGF.VALUE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(TGFB$TGF.VALUE)
hist(model$residuals)
rug(model$residuals)
TGFB[names(model$residuals),"TGFBXallcovars_no_MSL_RES"] <- model$residuals
View(TGFB)

#Remoção dos outliers
hist(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES)
TGFB$TGFBXallcovars_no_MSL_RES[TGFB$TGFBXallcovars_no_MSL_RES %in% boxplot.stats(TGFB$TGFBXallcovars_no_MSL_RES)$out]
length(TGFB$TGFBXallcovars_no_MSL_RES[TGFB$TGFBXallcovars_no_MSL_RES %in% boxplot.stats(TGFB$TGFBXallcovars_no_MSL_RES)$out])
TGFBXallcovars_no_MSL_out_rm <- subset(TGFB, TGFBXallcovars_no_MSL_RES < 98.34209)
nrow(TGFBXallcovars_no_MSL_out_rm)
TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm <- TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES
TGFBXallcovars_no_MSL_out_rm <- within(TGFBXallcovars_no_MSL_out_rm, rm(TGFBXallcovars_no_MSL_RES))
length(TGFB$TGFBXallcovars_no_MSL_RES)
length(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
length(TGFB$TGFBXallcovars_no_MSL_RES - TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
boxplot(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
hist(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
rug(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TGFB <- full_join(TGFB,TGFBXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TGFB)
View(TGFB)


#Transformação boxcox
TGFB$TGFBXallcovars_no_MSL_RES.t<-boxcox.ols(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES.t)
hist(TGFB$TGFBXallcovars_no_MSL_RES.t)
rug(TGFB$TGFBXallcovars_no_MSL_RES.t)
summary(TGFB$TGFBXallcovars_no_MSL_RES.t)

#Transformação logarítmica
TGFB$TGFBXallcovars_no_MSL_RES_LOG <- log10(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
hist(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
rug(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
summary(TGFB$TGFBXallcovars_no_MSL_RES_LOG)

View(TGFB)

write.table(TGFB,"TGFB_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)



#Criando um subset TGFB para normalização
TGFB <- select(pheno,IID,FID,TGF.VALUE)
#View(TGFB)

#Remoção dos outliers
hist(TGFB$TGF.VALUE)
boxplot(TGFB$TGF.VALUE)
length(TGFB$TGF.VALUE)
TGFB$TGF.VALUE[TGFB$TGF.VALUE %in% boxplot.stats(TGFB$TGF.VALUE)$out]
length(TGFB$TGF.VALUE[TGFB$TGF.VALUE %in% boxplot.stats(TGFB$TGF.VALUE)$out])
TGFB_out_rm <- subset(TGFB, TGF.VALUE < 110.88 & TGF.VALUE > -87.29)
nrow(TGFB_out_rm)
#View(TGFB_out_rm)
TGFB_out_rm$out_rm <- TGFB_out_rm$TGF.VALUE
#View(TGFB_out_rm)
TGFB_out_rm <- within(TGFB_out_rm, rm(TGF.VALUE))
#View(TGFBXallcovars_no_MSL_res_out_rm)
length(TGFB$TGF.VALUE)
length(TGFB_out_rm$out_rm)
length(TGFB$TGF.VALUE) - length(TGFB_out_rm$out_rm)
boxplot(TGFB_out_rm$out_rm)
hist(TGFB_out_rm$out_rm)
rug(TGFB_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TGFB <- full_join(TGFB,TGFB_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TGFB)
#View(TGFBXallcovars_no_MSL_res)

#Transformação boxcox
TGFB$TGF.VALUE.t<-boxcox.ols(TGFB$TGF.VALUE)
boxplot(TGFB$TGF.VALUE.t)
hist(TGFB$TGF.VALUE.t)
rug(TGFB$TGF.VALUE.t)
summary(TGFB$TGF.VALUE.t)

#Transformação logarítmica
TGFB$LOG <- log10(TGFB$TGF.VALUE)
boxplot(TGFB$LOG)
hist(TGFB$LOG)
rug(TGFB$LOG)
summary(TGFB$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$TGF.VALUE)
hist(pheno$TGF.VALUE)
model<-lm(TGF.VALUE ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(TGFB$TGF.VALUE)
hist(model$residuals)
rug(model$residuals)
TGFB[names(model$residuals),"TGFBXallcovars_no_MSL_RES"] <- model$residuals
View(TGFB)

#Remoção dos outliers
hist(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES)
TGFB$TGFBXallcovars_no_MSL_RES[TGFB$TGFBXallcovars_no_MSL_RES %in% boxplot.stats(TGFB$TGFBXallcovars_no_MSL_RES)$out]
length(TGFB$TGFBXallcovars_no_MSL_RES[TGFB$TGFBXallcovars_no_MSL_RES %in% boxplot.stats(TGFB$TGFBXallcovars_no_MSL_RES)$out])
TGFBXallcovars_no_MSL_out_rm <- subset(TGFB, TGFBXallcovars_no_MSL_RES < 86.26450 & TGFBXallcovars_no_MSL_RES > -81.83218)
nrow(TGFBXallcovars_no_MSL_out_rm)
TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm <- TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES
TGFBXallcovars_no_MSL_out_rm <- within(TGFBXallcovars_no_MSL_out_rm, rm(TGFBXallcovars_no_MSL_RES))
length(TGFB$TGFBXallcovars_no_MSL_RES)
length(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
length(TGFB$TGFBXallcovars_no_MSL_RES - TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
boxplot(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
hist(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
rug(TGFBXallcovars_no_MSL_out_rm$TGFBXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TGFB <- full_join(TGFB,TGFBXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TGFB)
View(TGFB)


#Transformação boxcox
TGFB$TGFBXallcovars_no_MSL_RES.t<-boxcox.ols(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES.t)
hist(TGFB$TGFBXallcovars_no_MSL_RES.t)
rug(TGFB$TGFBXallcovars_no_MSL_RES.t)
summary(TGFB$TGFBXallcovars_no_MSL_RES.t)

#Transformação logarítmica
TGFB$TGFBXallcovars_no_MSL_RES_LOG <- log10(TGFB$TGFBXallcovars_no_MSL_RES)
boxplot(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
hist(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
rug(TGFB$TGFBXallcovars_no_MSL_RES_LOG)
summary(TGFB$TGFBXallcovars_no_MSL_RES_LOG)

View(TGFB)

write.table(TGFB,"TGFB_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset NO para normalização
NO <- select(pheno,IID,FID,NO.NN)
#View(NO)

#Remoção dos outliers
hist(NO$NO.NN)
boxplot(NO$NO.NN)
length(NO$NO.NN)
NO$NO.NN[NO$NO.NN %in% boxplot.stats(NO$NO.NN)$out]
length(NO$NO.NN[NO$NO.NN %in% boxplot.stats(NO$NO.NN)$out])
NO_out_rm <- subset(NO, NO.NN > 0 & NO.NN < 4.48)
nrow(NO_out_rm)
#View(NO_out_rm)
NO_out_rm$out_rm <- NO_out_rm$NO.NN
#View(NO_out_rm)
NO_out_rm <- within(NO_out_rm, rm(NO.NN))
#View(NOXallcovars_no_MSL_res_out_rm)
length(NO$NO.NN)
length(NO_out_rm$out_rm)
length(NO$NO.NN) - length(NO_out_rm$out_rm)
boxplot(NO_out_rm$out_rm)
hist(NO_out_rm$out_rm)
rug(NO_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
NO <- full_join(NO,NO_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(NO)
#View(NOXallcovars_no_MSL_res)

#Transformação boxcox
NO$NO.NN.t<-boxcox.ols(NO$NO.NN)
boxplot(NO$NO.NN.t)
hist(NO$NO.NN.t)
rug(NO$NO.NN.t)
summary(NO$NO.NN.t)

#Transformação logarítmica
NO$LOG <- log10(NO$NO.NN)
boxplot(NO$LOG)
hist(NO$LOG)
rug(NO$LOG)
summary(NO$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$NO.NN)
hist(pheno$NO.NN)
model<-lm(NO.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(NO$NO.NN)
hist(model$residuals)
rug(model$residuals)
NO[names(model$residuals),"NOXallcovars_no_MSL_RES"] <- model$residuals
View(NO)

#Remoção dos outliers
hist(NO$NOXallcovars_no_MSL_RES)
boxplot(NO$NOXallcovars_no_MSL_RES)
NO$NOXallcovars_no_MSL_RES[NO$NOXallcovars_no_MSL_RES %in% boxplot.stats(NO$NOXallcovars_no_MSL_RES)$out]
length(NO$NOXallcovars_no_MSL_RES[NO$NOXallcovars_no_MSL_RES %in% boxplot.stats(NO$NOXallcovars_no_MSL_RES)$out])
NOXallcovars_no_MSL_out_rm <- subset(NO, NOXallcovars_no_MSL_RES < 3.032773)
nrow(NOXallcovars_no_MSL_out_rm)
NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm <- NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES
NOXallcovars_no_MSL_out_rm <- within(NOXallcovars_no_MSL_out_rm, rm(NOXallcovars_no_MSL_RES))
length(NO$NOXallcovars_no_MSL_RES)
length(NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm)
length(NO$NOXallcovars_no_MSL_RES - NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm)
boxplot(NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm)
hist(NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm)
rug(NOXallcovars_no_MSL_out_rm$NOXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
NO <- full_join(NO,NOXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(NO)
View(NO)


#Transformação boxcox
NO$NOXallcovars_no_MSL_RES.t<-boxcox.ols(NO$NOXallcovars_no_MSL_RES)
boxplot(NO$NOXallcovars_no_MSL_RES.t)
hist(NO$NOXallcovars_no_MSL_RES.t)
rug(NO$NOXallcovars_no_MSL_RES.t)
summary(NO$NOXallcovars_no_MSL_RES.t)

#Transformação logarítmica
NO$NOXallcovars_no_MSL_RES_LOG <- log10(NO$NOXallcovars_no_MSL_RES)
boxplot(NO$NOXallcovars_no_MSL_RES_LOG)
hist(NO$NOXallcovars_no_MSL_RES_LOG)
rug(NO$NOXallcovars_no_MSL_RES_LOG)
summary(NO$NOXallcovars_no_MSL_RES_LOG)

View(NO)

write.table(NO,"NO_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset H2O2 para normalização
H2O2 <- select(pheno,IID,FID,H2O2.NN)
#View(H2O2)

#Remoção dos outliers
hist(H2O2$H2O2.NN)
boxplot(H2O2$H2O2.NN)
length(H2O2$H2O2.NN)
H2O2$H2O2.NN[H2O2$H2O2.NN %in% boxplot.stats(H2O2$H2O2.NN)$out]
length(H2O2$H2O2.NN[H2O2$H2O2.NN %in% boxplot.stats(H2O2$H2O2.NN)$out])
H2O2_out_rm <- subset(H2O2, H2O2.NN < 0.65 & H2O2.NN > -0.66)
nrow(H2O2_out_rm)
#View(H2O2_out_rm)
H2O2_out_rm$out_rm <- H2O2_out_rm$H2O2.NN
#View(H2O2_out_rm)
H2O2_out_rm <- within(H2O2_out_rm, rm(H2O2.NN))
#View(H2O2Xtreatment_res_out_rm)
length(H2O2$H2O2.NN)
length(H2O2_out_rm$out_rm)
length(H2O2$H2O2.NN) - length(H2O2_out_rm$out_rm)
boxplot(H2O2_out_rm$out_rm)
hist(H2O2_out_rm$out_rm)
rug(H2O2_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
H2O2 <- full_join(H2O2,H2O2_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(H2O2)
#View(H2O2Xtreatment_res)

#Transformação boxcox
H2O2$H2O2.NN.t<-boxcox.ols(H2O2$H2O2.NN)
boxplot(H2O2$H2O2.NN.t)
hist(H2O2$H2O2.NN.t)
rug(H2O2$H2O2.NN.t)
summary(H2O2$H2O2.NN.t)

#Transformação logarítmica
H2O2$LOG <- log10(H2O2$H2O2.NN)
boxplot(H2O2$LOG)
hist(H2O2$LOG)
rug(H2O2$LOG)
summary(H2O2$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$H2O2.NN)
hist(pheno$H2O2.NN)
model<-lm(H2O2.NN ~ TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(H2O2$H2O2.NN)
hist(model$residuals)
rug(model$residuals)
H2O2[names(model$residuals),"H2O2Xtreatment_RES"] <- model$residuals
View(H2O2)

#Remoção dos outliers
hist(H2O2$H2O2Xtreatment_RES)
boxplot(H2O2$H2O2Xtreatment_RES)
H2O2$H2O2Xtreatment_RES[H2O2$H2O2Xtreatment_RES %in% boxplot.stats(H2O2$H2O2Xtreatment_RES)$out]
length(H2O2$H2O2Xtreatment_RES[H2O2$H2O2Xtreatment_RES %in% boxplot.stats(H2O2$H2O2Xtreatment_RES)$out])
H2O2Xtreatment_out_rm <- subset(H2O2, H2O2Xtreatment_RES < 0.5815833 & H2O2Xtreatment_RES > -0.7284167)
nrow(H2O2Xtreatment_out_rm)
H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm <- H2O2Xtreatment_out_rm$H2O2Xtreatment_RES
H2O2Xtreatment_out_rm <- within(H2O2Xtreatment_out_rm, rm(H2O2Xtreatment_RES))
length(H2O2$H2O2Xtreatment_RES)
length(H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm)
length(H2O2$H2O2Xtreatment_RES - H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm)
boxplot(H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm)
hist(H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm)
rug(H2O2Xtreatment_out_rm$H2O2Xtreatment_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
H2O2 <- full_join(H2O2,H2O2Xtreatment_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(H2O2)
View(H2O2)


#Transformação boxcox
H2O2$H2O2Xtreatment_RES.t<-boxcox.ols(H2O2$H2O2Xtreatment_RES)
boxplot(H2O2$H2O2Xtreatment_RES.t)
hist(H2O2$H2O2Xtreatment_RES.t)
rug(H2O2$H2O2Xtreatment_RES.t)
summary(H2O2$H2O2Xtreatment_RES.t)

#Transformação logarítmica
H2O2$H2O2Xtreatment_RES_LOG <- log10(H2O2$H2O2Xtreatment_RES)
boxplot(H2O2$H2O2Xtreatment_RES_LOG)
hist(H2O2$H2O2Xtreatment_RES_LOG)
rug(H2O2$H2O2Xtreatment_RES_LOG)
summary(H2O2$H2O2Xtreatment_RES_LOG)

View(H2O2)

write.table(H2O2,"H2O2_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset SOD para normalização
SOD <- select(pheno,IID,FID,SOD.NN)
#View(SOD)

# #Remoção dos outliers
# hist(SOD$SOD.NN)
# boxplot(SOD$SOD.NN)
# length(SOD$SOD.NN)
# SOD$SOD.NN[SOD$SOD.NN %in% boxplot.stats(SOD$SOD.NN)$out]
# length(SOD$SOD.NN[SOD$SOD.NN %in% boxplot.stats(SOD$SOD.NN)$out])
# SOD_out_rm <- subset(SOD, SOD.NN < 148.25)
# nrow(SOD_out_rm)
# #View(SOD_out_rm)
# SOD_out_rm$out_rm <- SOD_out_rm$SOD.NN
# #View(SOD_out_rm)
# SOD_out_rm <- within(SOD_out_rm, rm(SOD.NN))
# #View(SODXallcovars_no_MSL_res_out_rm)
# length(SOD$SOD.NN)
# length(SOD_out_rm$out_rm)
# length(SOD$SOD.NN) - length(SOD_out_rm$out_rm)
# boxplot(SOD_out_rm$out_rm)
# hist(SOD_out_rm$out_rm)
# rug(SOD_out_rm$out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# SOD <- full_join(SOD,SOD_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(SOD)
# #View(SODXallcovars_no_MSL_res)

#Transformação boxcox
SOD$SOD.NN.t<-boxcox.ols(SOD$SOD.NN)
boxplot(SOD$SOD.NN.t)
hist(SOD$SOD.NN.t)
rug(SOD$SOD.NN.t)
summary(SOD$SOD.NN.t)

#Transformação logarítmica
SOD$LOG <- log10(SOD$SOD.NN)
boxplot(SOD$LOG)
hist(SOD$LOG)
rug(SOD$LOG)
summary(SOD$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$SOD.NN)
hist(pheno$SOD.NN)
model<-lm(SOD.NN ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(SOD$SOD.NN)
hist(model$residuals)
rug(model$residuals)
SOD[names(model$residuals),"SODXallcovars_no_MSL_RES"] <- model$residuals
View(SOD)

# #Remoção dos outliers
# hist(SOD$SODXallcovars_no_MSL_RES)
# boxplot(SOD$SODXallcovars_no_MSL_RES)
# SOD$SODXallcovars_no_MSL_RES[SOD$SODXallcovars_no_MSL_RES %in% boxplot.stats(SOD$SODXallcovars_no_MSL_RES)$out]
# length(SOD$SODXallcovars_no_MSL_RES[SOD$SODXallcovars_no_MSL_RES %in% boxplot.stats(SOD$SODXallcovars_no_MSL_RES)$out])
# SODXallcovars_no_MSL_out_rm <- subset(SOD, SODXallcovars_no_MSL_RES < 98.34209)
# nrow(SODXallcovars_no_MSL_out_rm)
# SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm <- SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES
# SODXallcovars_no_MSL_out_rm <- within(SODXallcovars_no_MSL_out_rm, rm(SODXallcovars_no_MSL_RES))
# length(SOD$SODXallcovars_no_MSL_RES)
# length(SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm)
# length(SOD$SODXallcovars_no_MSL_RES - SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm)
# boxplot(SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm)
# hist(SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm)
# rug(SODXallcovars_no_MSL_out_rm$SODXallcovars_no_MSL_RES_out_rm)
# #Para juntar o data frame sem os outliers no data frame to resíduo total
# SOD <- full_join(SOD,SODXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
# nrow(SOD)
# View(SOD)


#Transformação boxcox
SOD$SODXallcovars_no_MSL_RES.t<-boxcox.ols(SOD$SODXallcovars_no_MSL_RES)
boxplot(SOD$SODXallcovars_no_MSL_RES.t)
hist(SOD$SODXallcovars_no_MSL_RES.t)
rug(SOD$SODXallcovars_no_MSL_RES.t)
summary(SOD$SODXallcovars_no_MSL_RES.t)

#Transformação logarítmica
SOD$SODXallcovars_no_MSL_RES_LOG <- log10(SOD$SODXallcovars_no_MSL_RES)
boxplot(SOD$SODXallcovars_no_MSL_RES_LOG)
hist(SOD$SODXallcovars_no_MSL_RES_LOG)
rug(SOD$SODXallcovars_no_MSL_RES_LOG)
summary(SOD$SODXallcovars_no_MSL_RES_LOG)

View(SOD)

write.table(SOD,"SOD_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)


#Criando um subset TOC para normalização
TOC <- select(pheno,IID,FID,TOS)
#View(TOC)

#Remoção dos outliers
hist(TOC$TOS)
boxplot(TOC$TOS)
length(TOC$TOS)
TOC$TOS[TOC$TOS %in% boxplot.stats(TOC$TOS)$out]
length(TOC$TOS[TOC$TOS %in% boxplot.stats(TOC$TOS)$out])
TOC_out_rm <- subset(TOC, TOS > 34.04)
nrow(TOC_out_rm)
#View(TOC_out_rm)
TOC_out_rm$out_rm <- TOC_out_rm$TOS
#View(TOC_out_rm)
TOC_out_rm <- within(TOC_out_rm, rm(TOS))
#View(TOCXcollar_res_out_rm)
length(TOC$TOS)
length(TOC_out_rm$out_rm)
length(TOC$TOS) - length(TOC_out_rm$out_rm)
boxplot(TOC_out_rm$out_rm)
hist(TOC_out_rm$out_rm)
rug(TOC_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TOC <- full_join(TOC,TOC_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TOC)
#View(TOCXcollar_res)

#Transformação boxcox
TOC$TOS.t<-boxcox.ols(TOC$TOS)
boxplot(TOC$TOS.t)
hist(TOC$TOS.t)
rug(TOC$TOS.t)
summary(TOC$TOS.t)

#Transformação logarítmica
TOC$LOG <- log10(TOC$TOS)
boxplot(TOC$LOG)
hist(TOC$LOG)
rug(TOC$LOG)
summary(TOC$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$TOS)
hist(pheno$TOS)
model<-lm(TOS ~ COLLAR,data = pheno)
summary(model)
length(model$residuals)
length(TOC$TOS)
hist(model$residuals)
rug(model$residuals)
TOC[names(model$residuals),"TOCXcollar_RES"] <- model$residuals
View(TOC)

#Remoção dos outliers
hist(TOC$TOCXcollar_RES)
boxplot(TOC$TOCXcollar_RES)
TOC$TOCXcollar_RES[TOC$TOCXcollar_RES %in% boxplot.stats(TOC$TOCXcollar_RES)$out]
length(TOC$TOCXcollar_RES[TOC$TOCXcollar_RES %in% boxplot.stats(TOC$TOCXcollar_RES)$out])
TOCXcollar_out_rm <- subset(TOC, TOCXcollar_RES > -113.8620)
nrow(TOCXcollar_out_rm)
TOCXcollar_out_rm$TOCXcollar_RES_out_rm <- TOCXcollar_out_rm$TOCXcollar_RES
TOCXcollar_out_rm <- within(TOCXcollar_out_rm, rm(TOCXcollar_RES))
length(TOC$TOCXcollar_RES)
length(TOCXcollar_out_rm$TOCXcollar_RES_out_rm)
length(TOC$TOCXcollar_RES - TOCXcollar_out_rm$TOCXcollar_RES_out_rm)
boxplot(TOCXcollar_out_rm$TOCXcollar_RES_out_rm)
hist(TOCXcollar_out_rm$TOCXcollar_RES_out_rm)
rug(TOCXcollar_out_rm$TOCXcollar_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TOC <- full_join(TOC,TOCXcollar_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TOC)
View(TOC)


#Transformação boxcox
TOC$TOCXcollar_RES.t<-boxcox.ols(TOC$TOCXcollar_RES)
boxplot(TOC$TOCXcollar_RES.t)
hist(TOC$TOCXcollar_RES.t)
rug(TOC$TOCXcollar_RES.t)
summary(TOC$TOCXcollar_RES.t)

#Transformação logarítmica
TOC$TOCXcollar_RES_LOG <- log10(TOC$TOCXcollar_RES)
boxplot(TOC$TOCXcollar_RES_LOG)
hist(TOC$TOCXcollar_RES_LOG)
rug(TOC$TOCXcollar_RES_LOG)
summary(TOC$TOCXcollar_RES_LOG)

View(TOC)

write.table(TOC,"TOC_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)



#Criando um subset TACpara normalização
TAC<- select(pheno,IID,FID,TAS)
#View(TAC)

#Remoção dos outliers
hist(TAC$TAS)
boxplot(TAC$TAS)
length(TAC$TAS)
TAC$TAS[TAC$TAS %in% boxplot.stats(TAC$TAS)$out]
length(TAC$TAS[TAC$TAS %in% boxplot.stats(TAC$TAS)$out])
TAC_out_rm <- subset(TAC, TAS < 0.76 & TAS > 0.29)
nrow(TAC_out_rm)
#View(TAC_out_rm)
TAC_out_rm$out_rm <- TAC_out_rm$TAS
#View(TAC_out_rm)
TAC_out_rm <- within(TAC_out_rm, rm(TAS))
#View(TACXallcovars_no_MSL_res_out_rm)
length(TAC$TAS)
length(TAC_out_rm$out_rm)
length(TAC$TAS) - length(TAC_out_rm$out_rm)
boxplot(TAC_out_rm$out_rm)
hist(TAC_out_rm$out_rm)
rug(TAC_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TAC<- full_join(TAC,TAC_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TAC)
#View(TACXallcovars_no_MSL_res)

#Transformação boxcox
TAC$TAS.t<-boxcox.ols(TAC$TAS)
boxplot(TAC$TAS.t)
hist(TAC$TAS.t)
rug(TAC$TAS.t)
summary(TAC$TAS.t)

#Transformação logarítmica
TAC$LOG <- log10(TAC$TAS)
boxplot(TAC$LOG)
hist(TAC$LOG)
rug(TAC$LOG)
summary(TAC$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$TAS)
hist(pheno$TAS)
model<-lm(TAS ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(TAC$TAS)
hist(model$residuals)
rug(model$residuals)
TAC[names(model$residuals),"TACXallcovars_no_MSL_RES"] <- model$residuals
View(TAC)

#Remoção dos outliers
hist(TAC$TACXallcovars_no_MSL_RES)
boxplot(TAC$TACXallcovars_no_MSL_RES)
TAC$TACXallcovars_no_MSL_RES[TAC$TACXallcovars_no_MSL_RES %in% boxplot.stats(TAC$TACXallcovars_no_MSL_RES)$out]
length(TAC$TACXallcovars_no_MSL_RES[TAC$TACXallcovars_no_MSL_RES %in% boxplot.stats(TAC$TACXallcovars_no_MSL_RES)$out])
TACXallcovars_no_MSL_out_rm <- subset(TAC, TACXallcovars_no_MSL_RES < 0.1796365 & TACXallcovars_no_MSL_RES > -0.1916669)
nrow(TACXallcovars_no_MSL_out_rm)
TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm <- TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES
TACXallcovars_no_MSL_out_rm <- within(TACXallcovars_no_MSL_out_rm, rm(TACXallcovars_no_MSL_RES))
length(TAC$TACXallcovars_no_MSL_RES)
length(TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm)
length(TAC$TACXallcovars_no_MSL_RES - TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm)
boxplot(TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm)
hist(TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm)
rug(TACXallcovars_no_MSL_out_rm$TACXallcovars_no_MSL_RES_out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
TAC<- full_join(TAC,TACXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(TAC)
View(TAC)


#Transformação boxcox
TAC$TACXallcovars_no_MSL_RES.t<-boxcox.ols(TAC$TACXallcovars_no_MSL_RES)
boxplot(TAC$TACXallcovars_no_MSL_RES.t)
hist(TAC$TACXallcovars_no_MSL_RES.t)
rug(TAC$TACXallcovars_no_MSL_RES.t)
summary(TAC$TACXallcovars_no_MSL_RES.t)

#Transformação logarítmica
TAC$TACXallcovars_no_MSL_RES_LOG <- log10(TAC$TACXallcovars_no_MSL_RES)
boxplot(TAC$TACXallcovars_no_MSL_RES_LOG)
hist(TAC$TACXallcovars_no_MSL_RES_LOG)
rug(TAC$TACXallcovars_no_MSL_RES_LOG)
summary(TAC$TACXallcovars_no_MSL_RES_LOG)

View(TAC)

write.table(TAC,"TAC_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)



#Criando um subset MDA para normalização
MDA<- select(pheno,IID,FID,MDA)
#View(MDA)

#Remoção dos outliers
hist(MDA$MDA)
boxplot(MDA$MDA)
length(MDA$MDA)
MDA$MDA[MDA$MDA %in% boxplot.stats(MDA$MDA)$out]
length(MDA$MDA[MDA$MDA %in% boxplot.stats(MDA$MDA)$out])
MDA_out_rm <- subset(MDA, MDA < 83.13)
nrow(MDA_out_rm)
#View(MDA_out_rm)
MDA_out_rm$out_rm <- MDA_out_rm$MDA
#View(MDA_out_rm)
MDA_out_rm <- within(MDA_out_rm, rm(MDA))
#View(MDAXallcovars_no_MSL_res_out_rm)
length(MDA$MDA)
length(MDA_out_rm$out_rm)
length(MDA$MDA) - length(MDA_out_rm$out_rm)
boxplot(MDA_out_rm$out_rm)
hist(MDA_out_rm$out_rm)
rug(MDA_out_rm$out_rm)
#Para juntar o data frame sem os outliers no data frame to resíduo total
MDA<- full_join(MDA,MDA_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(MDA)
#View(MDAXallcovars_no_MSL_res)

#Transformação boxcox
MDA$MDA.t<-boxcox.ols(MDA$MDA)
boxplot(MDA$MDA.t)
hist(MDA$MDA.t)
rug(MDA$MDA.t)
summary(MDA$MDA.t)

#Transformação logarítmica
MDA$LOG <- log10(MDA$MDA)
boxplot(MDA$LOG)
hist(MDA$LOG)
rug(MDA$LOG)
summary(MDA$LOG)

#Para normalizar o modelo de covariáveis selecionado
boxplot(pheno$MDA)
hist(pheno$MDA)
model<-lm(MDA ~ ORIGIN + SEX + AGE + COLLAR + VACCINE + TREATMENT,data = pheno)
summary(model)
length(model$residuals)
length(MDA$MDA)
hist(model$residuals)
rug(model$residuals)
MDA[names(model$residuals),"MDAXallcovars_no_MSL_RES"] <- model$residuals
View(MDA)

#Remoção dos outliers
hist(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES)
MDA$MDAXallcovars_no_MSL_RES[MDA$MDAXallcovars_no_MSL_RES %in% boxplot.stats(MDA$MDAXallcovars_no_MSL_RES)$out]
length(MDA$MDAXallcovars_no_MSL_RES[MDA$MDAXallcovars_no_MSL_RES %in% boxplot.stats(MDA$MDAXallcovars_no_MSL_RES)$out])
MDAXallcovars_no_MSL_out_rm <- subset(MDA, MDAXallcovars_no_MSL_RES < 57.41991)
nrow(MDAXallcovars_no_MSL_out_rm)
MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm <- MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES
MDAXallcovars_no_MSL_out_rm <- within(MDAXallcovars_no_MSL_out_rm, rm(MDAXallcovars_no_MSL_RES))
length(MDA$MDAXallcovars_no_MSL_RES)
length(MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
length(MDA$MDAXallcovars_no_MSL_RES - MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
boxplot(MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
hist(MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)
rug(MDAXallcovars_no_MSL_out_rm$MDAXallcovars_no_MSL_RES_out_rm)

#Para juntar o data frame sem os outliers no data frame to resíduo total
MDA<- full_join(MDA,MDAXallcovars_no_MSL_out_rm,by = NULL, copy = FALSE, suffix = c(".x", ".y"))
nrow(MDA)
View(MDA)


#Transformação boxcox
MDA$MDAXallcovars_no_MSL_RES.t<-boxcox.ols(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES.t)
hist(MDA$MDAXallcovars_no_MSL_RES.t)
rug(MDA$MDAXallcovars_no_MSL_RES.t)
summary(MDA$MDAXallcovars_no_MSL_RES.t)

#Transformação logarítmica
MDA$MDAXallcovars_no_MSL_RES_LOG <- log10(MDA$MDAXallcovars_no_MSL_RES)
boxplot(MDA$MDAXallcovars_no_MSL_RES_LOG)
hist(MDA$MDAXallcovars_no_MSL_RES_LOG)
rug(MDA$MDAXallcovars_no_MSL_RES_LOG)
summary(MDA$MDAXallcovars_no_MSL_RES_LOG)

View(MDA)

write.table(MDA,"MDA_NORMALIZED_PHENO.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
