#load libraries
library(tidyverse)
library(qpcR)
source("Homebrew/multiplot.R")
#read in data
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
df1 <- merge(Nb_Ns,Nc)
###global to single factor models
input <- df1[,c("Pop","LD_Nb","SF_Nb","Ns_Chao","Vk","NcAge1","Drainage","SampSites","SampSize","YearSinceTreat")]
input <- na.omit(input)
factors <- input[,c("LD_Nb","SF_Nb","Ns_Chao","Vk")]
responses <- input[,c("NcAge1","Drainage","SampSites","SampSize","YearSinceTreat")]

###checking correlations between factors
plot(responses)

##examine models with Nc and TrapToMouth distance ####
#creating blank data frame to save all AICc values for comparison
AICc_table <- data.frame(matrix(nrow = 9,ncol = 5))
colnames(AICc_table) <- c("ModelName",colnames(factors))
AICc_table$ModelName <- c("global",colnames(responses),"sampling","environmental","intercept")
pval_table <- data.frame(matrix(nrow = 9,ncol = 5))
colnames(pval_table) <- c("ModelName",colnames(factors))
pval_table$ModelName <- c("global",colnames(responses),"sampling","environmental","intercept")
facnames <- colnames(factors)

##loop to test all potential models
#run global models first
modeltmp <- glm(factors$Ns_Chao~responses$YearSinceTreat+responses$Drainage+responses$NcAge1+responses$SampSites+responses$SampSize)
AICc_table[1,4] <- AICc(modeltmp) 
pval_table[1,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Vk~responses$YearSinceTreat+responses$Drainage+responses$NcAge1+responses$SampSites+responses$SampSize)
AICc_table[1,5] <- AICc(modeltmp)  
pval_table[1,5] <- summary(modeltmp)$coefficients[,4][2]

#run intercept models
modeltmp <- glm(factors$Ns_Chao~1)
AICc_table[9,4] <- AICc(modeltmp) 
pval_table[9,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Vk~1)
AICc_table[9,5] <- AICc(modeltmp)  
pval_table[9,5] <- summary(modeltmp)$coefficients[,4][2]

for (a in 1:length(factors)) {
  f <- factors[,a]
  fname <- facnames[a]
  tmp <- cbind(factors[,a],responses)
  colnames(tmp) <- c(fname,colnames(responses))
  count <- 0
  for (fact1 in colnames(responses)) {
    count <- count+1
    modeltmp <- glm(f~tmp[,fact1])
    AICc_table[1+count,a+1] <- AICc(modeltmp) 
    pval_table[1+count,a+1] <- summary(modeltmp)$coefficients[,4][2]
  }
  #sampling model
  modeltmp <- glm(f~SampSites+SampSize,data = tmp)
  AICc_table[7,a+1] <- AICc(modeltmp) 
  pval_table[7,a+1] <- summary(modeltmp)$coefficients[,4][2]
  #environmental model
  #modeltmp <- glm(f~Drainage+TrapToMouth,data = tmp)
  #AICc_table[8,a+1] <- AICc(modeltmp) 
  #pval_table[8,a+1] <- summary(modeltmp)$coefficients[,4][2]
  #biotic model
  modeltmp <- glm(f~NcAge1+YearSinceTreat,data = tmp)
  AICc_table[8,a+1] <- AICc(modeltmp) 
  pval_table[8,a+1] <- summary(modeltmp)$coefficients[,4][2]
}

#examine AICc for the best model in each factor
apply(AICc_table, 2, min)

#create another table with Akeike weights
Chao_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$Ns_Chao)
Chao_AICc$AkaikeWeights <- akaike.weights(Chao_AICc$AICc_value)$weights
Vk_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$Vk)
Vk_AICc$AkaikeWeights <- akaike.weights(Vk_AICc$AICc_value)$weights

#Chao
modelChao <- glm(factors$Ns_Chao~responses$SampSites)
summary(modelChao)
plotdat <- data.frame(Ns_Chao=factors$Ns_Chao,SampSites=responses$SampSites)
sitesChao_plot <- ggplot(plotdat,aes(x=SampSites,y=Ns_Chao))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "glm", se = F,col = "darkblue")+
  stat_regline_equation(aes(label = ..rr.label..))+
  xlab("Sample Sites")+
  ylab("Ns - Chao")

tiff(filename = "Figures/LMPlots.tiff",height = 8,width = 6,units = "in",res = 400)
multiplot(cols = 2,sampLD_plot,sitesSF_plot,sizeChao_plot,drainLD_plot,sizeSF_plot,sitesChao_plot)
dev.off()

modelChao2 <- glm(factors$Ns_Chao~responses$Drainage)
summary(modelChao2)
plot(factors$Ns_Chao,responses$Drainage,
     pch = 16, 
     xlab = "log(Ns - Chao)",  
     ylab = "Drainage (ha)")
modelChao3 <- glm(factors$Ns_Chao~responses$SampSize+responses$SampSites)
summary(modelChao3)
modelChao4 <- glm(factors$Ns_Chao~responses$Drainage)
summary(modelChao4)
modelChao5 <- glm(factors$Ns_Chao~1)
summary(modelChao5)
modelChao8 <- glm(factors$Ns_Chao~responses$YearSinceTreat)
summary(modelChao8)

#Vk
modelVk <- glm(factors$Vk~1)
summary(modelVk)
modelVk2 <- glm(factors$Vk~responses$YearSinceTreat)
summary(modelVk2)
modelVk3 <- glm(factors$Vk~responses$Drainage)
summary(modelVk3)
modelVk4 <- glm(factors$Vk~responses$NcAge1)
summary(modelVk4)
modelVk5 <- glm(factors$Vk~responses$SampSites)
summary(modelVk5)





