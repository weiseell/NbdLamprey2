#load libraries
library(tidyverse)
library(qpcR)
#read in data
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
df1 <- merge(Nb_Ns,Nc)
###global to single factor models
input <- df1[,c("LD_Nb","SF_Nb","Ns_Chao","Vk","YearSinceTreat","Drainage","TrapToMouthKm","SampDist","SampSize","NcAge1")]
input <- na.omit(input)
factors <- input[,c("LD_Nb","SF_Nb","Ns_Chao","Vk")]
responses <- input[,c("YearSinceTreat","Drainage","TrapToMouthKm","SampDist","SampSize","NcAge1")]
input <- df1 %>% 
  select(LD_Nb,SF_Nb,Ns_Chao,Vk,YearSinceTreat,Drainage,TrapToMouthKm,SampDist,SampSize,NcAge1)
input <- na.omit(input)
factors <- input %>% 
  select(LD_Nb,SF_Nb,Ns_Chao,Vk)
responses <- input %>% 
  select(YearSinceTreat,Drainage,TrapToMouthKm,SampDist,SampSize,NcAge1)

###checking correlations between factors
plot(responses)

#log transform factors
factors$LD_Nb <- log(factors$LD_Nb)
factors$SF_Nb <- log(factors$SF_Nb)
factors$Ns_Chao <- log(factors$Ns_Chao)
factors$Vk <- log(factors$Vk)
##examine models with Nc and TrapToMouth distance ####
#creating blank data frame to save all AICc values for comparison
AICc_table <- data.frame(matrix(nrow = 11,ncol = 5))
colnames(AICc_table) <- c("ModelName",colnames(factors))
AICc_table$ModelName <- c("global",colnames(responses),"sampling","environmental","biotic","intercept")
pval_table <- data.frame(matrix(nrow = 11,ncol = 5))
colnames(pval_table) <- c("ModelName",colnames(factors))
pval_table$ModelName <- c("global",colnames(responses),"sampling","environmental","biotic","intercept")
facnames <- colnames(factors)

##loop to test all potential models
#run global models first
modeltmp <- glm(factors$LD_Nb~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
AICc_table[1,2] <- AICc(modeltmp) 
pval_table[1,2] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$SF_Nb~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
AICc_table[1,3] <- AICc(modeltmp) 
pval_table[1,3] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Ns_Chao~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
AICc_table[1,4] <- AICc(modeltmp) 
pval_table[1,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Vk~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize+responses$NcAge1)
AICc_table[1,5] <- AICc(modeltmp)  
pval_table[1,5] <- summary(modeltmp)$coefficients[,4][2]

#run intercept models
modeltmp <- glm(factors$LD_Nb~1)
AICc_table[11,2] <- AICc(modeltmp) 
pval_table[11,2] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$SF_Nb~1)
AICc_table[11,3] <- AICc(modeltmp) 
pval_table[11,3] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Ns_Chao~1)
AICc_table[11,4] <- AICc(modeltmp) 
pval_table[11,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors$Vk~1)
AICc_table[11,5] <- AICc(modeltmp)  
pval_table[11,5] <- summary(modeltmp)$coefficients[,4][2]

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
  modeltmp <- glm(f~SampDist+SampSize,data = tmp)
  AICc_table[8,a+1] <- AICc(modeltmp) 
  pval_table[8,a+1] <- summary(modeltmp)$coefficients[,4][2]
  #environmental model
  modeltmp <- glm(f~Drainage+TrapToMouthKm,data = tmp)
  AICc_table[9,a+1] <- AICc(modeltmp) 
  pval_table[9,a+1] <- summary(modeltmp)$coefficients[,4][2]
  #biotic model
  modeltmp <- glm(f~NcAge1+YearSinceTreat,data = tmp)
  AICc_table[10,a+1] <- AICc(modeltmp) 
  pval_table[10,a+1] <- summary(modeltmp)$coefficients[,4][2]
}

#examine AICc for the best model in each factor
apply(AICc_table, 2, min)

#create another table with Akeike weights
LD_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$LD_Nb)
LD_AICc <- LD_AICc[-c(7,10),]
LD_AICc$AkaikeWeights <- akaike.weights(LD_AICc$AICc_value)$weights
SF_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$SF_Nb)
SF_AICc <- SF_AICc[-c(7,10),]
SF_AICc$AkaikeWeights <- akaike.weights(SF_AICc$AICc_value)$weights
Chao_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$Ns_Chao)
Chao_AICc$AkaikeWeights <- akaike.weights(Chao_AICc$AICc_value)$weights
Vk_AICc <- data.frame(ModelName = AICc_table$ModelName,AICc_value = AICc_table$Vk)
Vk_AICc$AkaikeWeights <- akaike.weights(Vk_AICc$AICc_value)$weights

##run best model from each group
modelSF <- glm(factors$SF_Nb~responses$SampSize)
summary(modelSF)
modelChao <- glm(factors$Ns_Chao~responses$SampSize)
summary(modelChao)
modelVk <- glm(factors$Vk~responses$YearSinceTreat)
summary(modelVk)
##run models that are similar AIC to best model
#LD
modelLD <- glm(factors$LD_Nb~responses$SampSize)
summary(modelLD)

plotdat <- data.frame(LD_Nb=factors$LD_Nb,SampSize=responses$SampSize)
ggplot(plotdat,aes(x=SampSize,y=LD_Nb))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "glm", se = F,col = "darkblue")+
  xlab("Sample Size")+
  ylab("log(Nb - LD)")


modelLD2 <- glm(factors$LD_Nb~responses$Drainage)
summary(modelLD2)
plotdat <- data.frame(LD_Nb=factors$LD_Nb,Drainage=responses$Drainage)
ggplot(plotdat,aes(x=Drainage,y=LD_Nb))+
  geom_point()+
  theme_bw()+
  xlab("Drainage (ha)")+
  ylab("log(Nb - LD)")

modelLD2 <- glm(factors$LD_Nb~responses$Drainage)
summary(modelLD2)

modelLD3 <- glm(factors$LD_Nb~responses$YearSinceTreat)
summary(modelLD3)
modelLD4 <- glm(factors$LD_Nb~responses$SampDist+responses$SampSize)
summary(modelLD4)
modelLD5 <- glm(factors$LD_Nb~responses$SampDist)
summary(modelLD5)
modelLD6 <- glm(factors$LD_Nb~responses$TrapToMouthKm)
summary(modelLD6)
modelLD7 <- glm(factors$LD_Nb~responses$Drainage+responses$TrapToMouthKm)
summary(modelLD7)
modelLD8 <- glm(factors$LD_Nb~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
summary(modelLD8)

#SF
modelSF <- glm(factors$SF_Nb~responses$SampSize)
summary(modelSF)
plot(factors$SF_Nb,responses$SampSize,
     pch = 16, 
     xlab = "log(Nb - SF)",  
     ylab = "Sample Size")

modelSF2 <- glm(factors$SF_Nb~responses$Drainage)
summary(modelSF2)
plotdat <- data.frame(SF_Nb=factors$SF_Nb,Drainage=responses$Drainage)
ggplot(plotdat,aes(x=Drainage,y=SF_Nb))+
  geom_point()+
  theme_bw()+
  xlab("Drainage (ha)")+
  ylab("log(Nb - SF)")

modelSF2 <- glm(factors$SF_Nb~responses$Drainage)
summary(modelSF2)
modelSF3 <- glm(factors$SF_Nb~responses$YearSinceTreat)
summary(modelSF3)
modelSF4 <- glm(factors$SF_Nb~responses$SampDist+responses$SampSize)
summary(modelSF4)
modelSF5 <- glm(factors$SF_Nb~responses$SampDist)
summary(modelSF5)
modelSF6 <- glm(factors$SF_Nb~responses$TrapToMouthKm)
summary(modelSF6)
modelSF7 <- glm(factors$SF_Nb~responses$Drainage+responses$TrapToMouthKm)
summary(modelSF7)
modelSF8 <- glm(factors$SF_Nb~responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
summary(modelSF8)

#Chao
modelChao <- glm(factors$Ns_Chao~responses$SampSize)
summary(modelChao)
plotdat <- data.frame(Ns_Chao=factors$Ns_Chao,SampSize=responses$SampSize)
ggplot(plotdat,aes(x=SampSize,y=Ns_Chao))+
  geom_point()+
  theme_bw()+
  #geom_smooth(method = "glm", se = F,col = "darkblue")+
  xlab("Sample Size")+
  ylab("log(Ns - Chao)")

modelChao2 <- glm(factors$Ns_Chao~responses$Drainage)
summary(modelChao2)
plot(factors$Ns_Chao,responses$Drainage,
     pch = 16, 
     xlab = "log(Ns - Chao)",  
     ylab = "Drainage (ha)")

modelChao2 <- glm(factors$Ns_Chao~responses$Drainage)
summary(modelChao2)
modelChao3 <- glm(factors$Ns_Chao~responses$YearSinceTreat)
summary(modelChao3)
modelChao4 <- glm(factors$Ns_Chao~responses$SampDist+responses$SampSize)
summary(modelChao4)
modelChao5 <- glm(factors$Ns_Chao~responses$SampDist)
summary(modelChao5)
plotdat <- data.frame(Ns_Chao=factors$Ns_Chao,SampDist=responses$SampDist)
ggplot(plotdat,aes(x=SampDist,y=Ns_Chao))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "glm", se = F,col = "darkblue")+
  xlab("Sampling Distance")+
  ylab("log(Ns - Chao)")

modelChao6 <- glm(factors$Ns_Chao~responses$TrapToMouthKm)
summary(modelChao6)
modelChao7 <- glm(factors$Ns_Chao~responses$NcAge1)
summary(modelChao7)
modelChao8 <- glm(factors$Ns_Chao~responses$Drainage+responses$TrapToMouthKm)
summary(modelChao8)
modelChao9 <- glm(factors$Ns_Chao~responses$NcAge1+responses$YearSinceTreat)
summary(modelChao9)
modelChao10 <- glm(factors$Ns_Chao~responses$NcAge1+responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
summary(modelChao10)

#Vk
modelVk <- glm(factors$Vk~responses$SampSize)
summary(modelVk)
modelVk2 <- glm(factors$Vk~responses$Drainage)
summary(modelVk2)
modelVk3 <- glm(factors$Vk~responses$YearSinceTreat)
summary(modelVk3)
modelVk4 <- glm(factors$Vk~responses$SampDist+responses$SampSize)
summary(modelVk4)
modelVk5 <- glm(factors$Vk~responses$SampDist)
summary(modelVk5)
modelVk6 <- glm(factors$Vk~responses$TrapToMouthKm)
summary(modelVk6)
modelVk7 <- glm(factors$Vk~responses$NcAge1)
summary(modelVk7)
modelVk8 <- glm(factors$Vk~responses$Drainage+responses$TrapToMouthKm)
summary(modelVk8)
modelVk9 <- glm(factors$Vk~responses$NcAge1+responses$YearSinceTreat)
summary(modelVk9)
modelVk10 <- glm(factors$Vk~responses$NcAge1+responses$YearSinceTreat+responses$Drainage+responses$TrapToMouthKm+responses$SampDist+responses$SampSize)
summary(modelVk10)

#separate larval population models ####
###global to single factor models
input2 <- df1[,c("LD_Nb","SF_Nb","Ns_Chao","Vk","YearSinceTreat","Drainage","SampDist","SampSize","LarvalPop")]
input2 <- na.omit(input2)
factors2 <- input2[,c("LD_Nb","SF_Nb","Ns_Chao","Vk")]
responses2 <- input2[,c("YearSinceTreat","Drainage","SampDist","SampSize","LarvalPop")]

#rerun global models with this set
AICc_table2 <- data.frame(matrix(nrow = 9,ncol = 5))
colnames(AICc_table2) <- c("ModelName",colnames(factors2))

AICc_table2$ModelName <- c("global",colnames(responses2),"sampling","biotic","intercept")
pval_table2 <- data.frame(matrix(nrow = 9,ncol = 5))
colnames(pval_table2) <- c("ModelName",colnames(factors2))
pval_table2$ModelName <- c("global",colnames(responses2),"sampling","biotic","intercept")

AICc_table2$ModelName <- c("global",colnames(responses2),"sampling","environmental","biotic")
pval_table2 <- data.frame(matrix(nrow = 9,ncol = 5))
colnames(pval_table2) <- c("ModelName",colnames(factors2))
pval_table2$ModelName <- c("global",colnames(responses2),"sampling","environmental","biotic")

facnames <- colnames(factors2)
modeltmp <- glm(factors2$LD_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,2] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,2] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$SF_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,3] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,3] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$Ns_Chao~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,4] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$Vk~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,5] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,5] <- summary(modeltmp)$coefficients[,4][2]
#run intercept models
modeltmp <- glm(factors2$LD_Nb~1)
AICc_table2[9,2] <- AICc(modeltmp) 
pval_table2[9,2] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$SF_Nb~1)
AICc_table2[9,3] <- AICc(modeltmp) 
pval_table2[9,3] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$Ns_Chao~1)
AICc_table2[9,4] <- AICc(modeltmp) 
pval_table2[9,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$Vk~1)
AICc_table2[9,5] <- AICc(modeltmp)  
pval_table2[9,5] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$LD_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,3] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,3] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$LD_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,4] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,4] <- summary(modeltmp)$coefficients[,4][2]

modeltmp <- glm(factors2$LD_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
AICc_table2[1,5] <- AIC(modeltmp)+((2*5*(5+1))/(nrow(modeltmp$model)-5-1)) 
pval_table2[1,5] <- summary(modeltmp)$coefficients[,4][2]


#single factor and hypotheses models
for (a in 1:length(factors2)) {
  f <- factors2[,a]
  fname <- facnames[a]
  tmp <- cbind(factors2[,a],responses2)
  colnames(tmp) <- c(fname,colnames(responses2))
  count <- 0
  for (fact1 in colnames(responses2)) {
    count <- count+1
    modeltmp <- glm(f~tmp[,fact1])
    AICc_table2[1+count,a+1] <- AIC(modeltmp)+((2*1*(1+1))/(nrow(modeltmp$model)-1-1))
    pval_table2[1+count,a+1] <- summary(modeltmp)$coefficients[,4][2]
  }
  #sampling model
  modeltmp <- glm(f~SampDist+SampSize,data = tmp)
  AICc_table2[7,a+1] <- AIC(modeltmp)+((2*2*(2+1))/(nrow(modeltmp$model)-2-1))
  pval_table2[7,a+1] <- summary(modeltmp)$coefficients[,4][2]

  #biotic model
  modeltmp <- glm(f~LarvalPop+YearSinceTreat,data = tmp)
  AICc_table2[8,a+1] <- AIC(modeltmp)+((2*2*(2+1))/(nrow(modeltmp$model)-2-1))
  pval_table2[8,a+1] <- summary(modeltmp)$coefficients[,4][2]

  #environmental model
  modeltmp <- glm(f~Drainage,data = tmp)
  AICc_table2[8,a+1] <- AIC(modeltmp)+((2*2*(2+1))/(nrow(modeltmp$model)-2-1))
  pval_table2[8,a+1] <- summary(modeltmp)$coefficients[,4][2]
  #biotic model
  modeltmp <- glm(f~LarvalPop+YearSinceTreat,data = tmp)
  AICc_table2[9,a+1] <- AIC(modeltmp)+((2*2*(2+1))/(nrow(modeltmp$model)-2-1))
  pval_table2[9,a+1] <- summary(modeltmp)$coefficients[,4][2]

}

#examine AICc for the best model in each factor
apply(AICc_table2, 2, min)


#Akaike Weights for larval pop models
LD_AICc <- data.frame(ModelName = AICc_table2$ModelName,AICc_value = AICc_table2$LD_Nb)
LD_AICc$AkaikeWeights <- akaike.weights(LD_AICc$AICc_value)$weights
SF_AICc <- data.frame(ModelName = AICc_table2$ModelName,AICc_value = AICc_table2$SF_Nb)
SF_AICc$AkaikeWeights <- akaike.weights(SF_AICc$AICc_value)$weights
Chao_AICc <- data.frame(ModelName = AICc_table2$ModelName,AICc_value = AICc_table2$Ns_Chao)
Chao_AICc$AkaikeWeights <- akaike.weights(Chao_AICc$AICc_value)$weights
Vk_AICc <- data.frame(ModelName = AICc_table2$ModelName,AICc_value = AICc_table2$Vk)
Vk_AICc$AkaikeWeights <- akaike.weights(Vk_AICc$AICc_value)$weights

#Run best model for each factor
modelLD2.1 <- glm(factors2$LD_Nb~responses2$Drainage)
summary(modelLD2.1)
modelSF2.1 <- glm(factors2$SF_Nb~responses2$Drainage)
summary(modelSF2.1)
modelChao2.1 <- glm(factors2$Ns_Chao~responses2$YearSinceTreat)
summary(modelChao2.1)
modelVk2.1 <- glm(factors2$Vk~responses2$SampDist+responses2$SampSize)
summary(modelVk2.1)

#run other models that performed similarly well
modelLD2.1 <- glm(factors2$LD_Nb~responses2$Drainage)
summary(modelLD2.1)
modelLD2.2 <- glm(factors2$LD_Nb~responses2$LarvalPop)
summary(modelLD2.2)
modelLD2.3 <- glm(factors2$LD_Nb~responses2$YearSinceTreat)
summary(modelLD2.3)
modelLD2.4 <- glm(factors2$LD_Nb~responses2$LarvalPop+responses2$YearSinceTreat)
summary(modelLD2.4)
modelLD2.5 <- glm(factors2$LD_Nb~responses2$SampDist)
summary(modelLD2.5)
modelLD2.6 <- glm(factors2$LD_Nb~responses2$SampSize)
summary(modelLD2.6)
modelLD2.7 <- glm(factors2$LD_Nb~responses2$SampDist+responses2$SampSize)
summary(modelLD2.7)
modelLD2.8 <- glm(factors2$LD_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
summary(modelLD2.8)

modelSF2.1 <- glm(factors2$SF_Nb~responses2$Drainage)
summary(modelSF2.1)
modelSF2.2 <- glm(factors2$SF_Nb~responses2$LarvalPop)
summary(modelSF2.2)
modelSF2.3 <- glm(factors2$SF_Nb~responses2$YearSinceTreat)
summary(modelSF2.3)
modelSF2.4 <- glm(factors2$SF_Nb~responses2$LarvalPop+responses2$YearSinceTreat)
summary(modelSF2.4)
modelSF2.5 <- glm(factors2$SF_Nb~responses2$SampDist)
summary(modelSF2.5)
modelSF2.6 <- glm(factors2$SF_Nb~responses2$SampSize)
summary(modelSF2.6)
modelSF2.7 <- glm(factors2$SF_Nb~responses2$SampDist+responses2$SampSize)
summary(modelSF2.7)
modelSF2.8 <- glm(factors2$SF_Nb~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
summary(modelSF2.8)

modelChao2.1 <- glm(factors2$Ns_Chao~responses2$Drainage)
summary(modelChao2.1)
modelChao2.2 <- glm(factors2$Ns_Chao~responses2$LarvalPop)
summary(modelChao2.2)
modelChao2.3 <- glm(factors2$Ns_Chao~responses2$YearSinceTreat)
summary(modelChao2.3)
modelChao2.4 <- glm(factors2$Ns_Chao~responses2$LarvalPop+responses2$YearSinceTreat)
summary(modelChao2.4)
modelChao2.5 <- glm(factors2$Ns_Chao~responses2$SampDist)
summary(modelChao2.5)
modelChao2.6 <- glm(factors2$Ns_Chao~responses2$SampSize)
summary(modelChao2.6)
modelChao2.7 <- glm(factors2$Ns_Chao~responses2$SampDist+responses2$SampSize)
summary(modelChao2.7)
modelChao2.8 <- glm(factors2$Ns_Chao~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
summary(modelChao2.8)

modelVk2.1 <- glm(factors2$Vk~responses2$Drainage)
summary(modelVk2.1)
modelVk2.2 <- glm(factors2$Vk~responses2$LarvalPop)
summary(modelVk2.2)
modelVk2.3 <- glm(factors2$Vk~responses2$YearSinceTreat)
summary(modelVk2.3)
modelVk2.4 <- glm(factors2$Vk~responses2$LarvalPop+responses2$YearSinceTreat)
summary(modelVk2.4)
modelVk2.5 <- glm(factors2$Vk~responses2$SampDist)
summary(modelVk2.5)
modelVk2.6 <- glm(factors2$Vk~responses2$SampSize)
summary(modelVk2.6)
modelVk2.7 <- glm(factors2$Vk~responses2$SampDist+responses2$SampSize)
summary(modelVk2.7)
modelVk2.8 <- glm(factors2$Vk~responses2$YearSinceTreat+responses2$Drainage+responses2$LarvalPop+responses2$SampDist+responses2$SampSize)
summary(modelVk2.8)

modelLD2.2 <- glm(factors2$LD_Nb~responses2$YearSinceTreat)
summary(modelLD2.2)
modelLD2.3 <- glm(factors2$LD_Nb~responses2$SampSize)
summary(modelLD2.3)

modelSF2.2 <- glm(factors2$SF_Nb~responses2$YearSinceTreat)
summary(modelSF2.2)
modelSF2.3 <- glm(factors2$SF_Nb~responses2$SampSize)
summary(modelSF2.3)

modelVk2.2 <- glm(factors2$Vk~responses2$YearSinceTreat)
summary(modelVk2.2)
modelSF2.3 <- glm(factors2$Vk~responses2$SampDist)
summary(modelSF2.3)


