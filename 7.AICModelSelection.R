#7 - Model Selection using AIC

rm(list=ls());  options(show.error.locations = TRUE); #remove the stuff you did before working on this script

#load packages
library(MuMIn) #AIC and dredge functions. Citation:
#Bart?n, K. 2017. MuMIn, Multi-model inference. R package version 1.40.0 [online]. 
#Available from www.CRAN.R-project.org/package=MuMIn
library(car) #VIF scores
library(lme4) #linear modeling package
library(readr) #directly read in .csv files
library(tidyverse)

#load in data
df <- read.csv("Models/ModelInputAll_081222.csv")

#Nb - LD ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Nb_LD ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Nb_LD ~ SampSize + TreatMonth + 
          Drainage + TreatYear +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
df$TreatMonth <- as.factor(df$TreatMonth)
LD_global <- glm(Nb_LD ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + EstLarvalAge + SampSize * SampDist +
                   YearSinceTreat * TreatMonth + EstLarvalAge * YearSinceTreat,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)

LD_AICModel <- dredge(LD_global)
View(LD_AICModel)
summary(LD_AICModel)
write.csv(LD_AICModel,file = "Models/LD_dredged_AICTable_081922.csv",append = F,quote = F)
LD_ModelAve <- model.avg(LD_AICModel, subset = delta < 2)
View(LD_ModelAve)
summary(LD_ModelAve)


#Nb - SF ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(SF_Nb ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Nb_LD ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
SF_global <- glm(SF_Nb ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + EstLarvalAge + SampSize * SampDist +
                   YearSinceTreat * TreatMonth + EstLarvalAge * YearSinceTreat,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(SF_global)

SF_AICModel <- dredge(SF_global)
View(SF_AICModel)
write.csv(SF_AICModel,file = "Models/SF_dredged_AICTable_081922.csv",append = F,quote = F)
SF_ModelAve <- model.avg(SF_AICModel, subset = delta < 2)
summary(SF_ModelAve)

#Ns - Chao ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Ns_Chao ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Ns_Chao ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
Ns_global <- glm(Ns_Chao ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + EstLarvalAge + SampSize * SampDist +
                   YearSinceTreat * TreatMonth + EstLarvalAge * YearSinceTreat,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(Ns_global)

Ns_AICModel <- dredge(Ns_global)
View(Ns_AICModel)
write.csv(Ns_AICModel,file = "Models/Ns_dredged_AICTable_081922.csv",append = F,quote = F)
Ns_ModelAve <- model.avg(Ns_AICModel, subset = delta < 2)
summary(Ns_ModelAve)

#Vk ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Vk ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Ns_Chao ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist + EstLarvalAge,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
Vk_global <- glm(Vk ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + EstLarvalAge + SampSize * SampDist +
                   YearSinceTreat * TreatMonth + EstLarvalAge * YearSinceTreat,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(Vk_global)

Vk_AICModel <- dredge(Vk_global)
View(Vk_AICModel)
write.csv(Vk_AICModel,file = "Models/Vk_dredged_AICTable_081922.csv",append = F,quote = F)

Vk_ModelAve <- model.avg(Vk_AICModel, subset = delta < 2)
summary(Vk_ModelAve)
