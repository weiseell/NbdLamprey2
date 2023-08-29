##5 - Model Selection using AIC
## Code adapted from Doug Larson

rm(list=ls());  options(show.error.locations = TRUE); #remove the stuff you did before working on this script

#load packages
library(MuMIn) #AIC and dredge functions. Citation:
#Bart?n, K. 2017. MuMIn, Multi-model inference. R package version 1.40.0 [online]. 
#Available from www.CRAN.R-project.org/package=MuMIn
library(car) #VIF scores
library(lme4) #linear modeling package
library(readr) #directly read in .csv files
library(tidyverse)
library(ggpubr)
#load in data
df <- read.csv("Models/ModelInputAll_120222.csv")

#Nb - LD ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Nb_LD ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Nb_LD ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist,  
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
                   SampDist + SampSize * SampDist +
                   YearSinceTreat * TreatMonth,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)

LD_AICModel <- dredge(LD_global)
View(LD_AICModel)
write.csv(LD_AICModel,file = "Models/LD_dredged_AICTable_030123.csv",append = F,quote = F)
LD_ModelAve <- model.avg(LD_AICModel, subset = delta < 2)
View(LD_ModelAve)
summary(LD_ModelAve)

#summarizing averaged coefficients
LD_ModelAve[["coefArray"]]

coefout1 <- data.frame(estimate = "Nb - LD",
                      varname = c("Drainage","Sample Size","Sampling \nDistance"),
           coefficient = c(0.022009,0.559521,-4.088938),
           stderr = c(0.007874,0.161669,5.64836),
           diffzero = c("Y","Y","N"))

#plotting coefficients
#Nb - SF ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(SF_Nb ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Nb_LD ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
SF_global <- glm(SF_Nb ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + SampSize * SampDist +
                   YearSinceTreat * TreatMonth,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(SF_global)

SF_AICModel <- dredge(SF_global)
View(SF_AICModel)
write.csv(SF_AICModel,file = "Models/SF_dredged_AICTable_030123.csv",append = F,quote = F)
SF_ModelAve <- model.avg(SF_AICModel, subset = delta < 2)
summary(SF_ModelAve)

coefout2 <- data.frame(estimate = "Nb - SF",
                       varname = c("Drainage","Sample Size","Sampling \nDistance"),
                       coefficient = c(0.03191,0.82661,-5.40557),
                       stderr = c(0.01254,0.25233,8.40593),
                       diffzero = c("Y","Y","N"))

#Ns - Chao ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Ns_Chao ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Ns_Chao ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
Ns_global <- glm(Ns_Chao ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + SampSize * SampDist +
                   YearSinceTreat * TreatMonth,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(Ns_global)

Ns_AICModel <- dredge(Ns_global)
View(Ns_AICModel)
write.csv(Ns_AICModel,file = "Models/Ns_dredged_AICTable_030123.csv",append = F,quote = F)
Ns_ModelAve <- model.avg(Ns_AICModel, subset = delta < 2)
summary(Ns_ModelAve)

coefout3 <- data.frame(estimate = "Ns",
                       varname = c("Sampling \nDistance","Sample Size","Sample Size:\nSampling Distance","Years Since \nTFM Treatment","Drainage"),
                       coefficient = c(-11.221045,1.970903,-0.094592,-7.389010,0.002842),
                       stderr = c(17.968862,0.548067,0.106730,14.623939,0.008078),
                       diffzero = c("N","Y","N","N","N"))

#Vk ####
#colinearity check
#Variance Inflation - Run on your fixed effects ONLY to determine if you have multicolinearity in your data.
vif(glm(Vk ~ SampSize + TreatYear +
          TreatMonth + Drainage + SampSites + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#recheck after removing colinear variables
vif(glm(Ns_Chao ~ SampSize + TreatMonth + 
          Drainage + YearSinceTreat +
          SampDist,  
        data = df, 
        family = gaussian 
        (link = "identity"), 
        na.action = na.fail))

#run global model
#!#Run twice - first run with only single fixed effects
#second run with interactions between terms in the best model
Vk_global <- glm(Vk ~ SampSize + TreatMonth + 
                   Drainage + YearSinceTreat +
                   SampDist + SampSize * SampDist +
                   YearSinceTreat * TreatMonth,  
                 data = df, 
                 family = gaussian 
                 (link = "identity"), 
                 na.action = na.fail)
summary(Vk_global)

Vk_AICModel <- dredge(Vk_global)
View(Vk_AICModel)
write.csv(Vk_AICModel,file = "Models/Vk_dredged_AICTable_030123.csv",append = F,quote = F)
Vk_ModelAve <- model.avg(Vk_AICModel, subset = delta < 2)
summary(Vk_ModelAve)

coefout4 <- data.frame(estimate = "Vk",
                       varname = c("Sampling Distance","Sample Size","Years Since \n TFM Treatment","Sample Size:\nSampling Distance"),
                       coefficient = c(13.21744,-0.13621,12.26493,-0.02253),
                       stderr = c(6.45261,0.19147,8.55307,0.03652),
                       diffzero = c("Y","N","Y","N"))

##Plotting coefficiencts #####
coefoutall <- rbind(coefout1,coefout2,coefout3,coefout4)

estnames <- c("Nb - Linkage Disequilibrium",
              "Nb - Sibship Frequency",
              "Extrapolated Ns",
              "Vk")
names(estnames) <- unique(coefoutall$estimate)

tiff(filename = "Figures/CoeffVisualization_GLM_030123.tiff",width = 6.5,height = 8,units = "in",res = 100)
ggplot(coefoutall,aes(x = coefficient,y = varname,color = diffzero))+
  facet_wrap(~estimate, scales = "free_y",labeller = as_labeller(estnames))+
  geom_point()+
  geom_errorbar(aes(xmin = coefficient-stderr,xmax = coefficient+stderr),width = 0.2)+
  geom_vline(xintercept = 0,linetype = "dashed")+
  geom_label(label = coefoutall$coefficient,nudge_y = 0.3)+
  theme_pubclean()+
  theme(legend.position = "none",axis.title = element_blank())
dev.off()
