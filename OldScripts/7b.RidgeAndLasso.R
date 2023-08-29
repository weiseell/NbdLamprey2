##Ridge and Lasso Model Evaluation Method

##load data
#library load
library(tidyverse)
library(glmnet)

#data load
env <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
#Nb:Nc ratio calculations and correlations
Nb_Ns1 <- Nb_Ns %>% 
  select(Pop,SampSize,LD_Nb,PwoP_Nb,SF_Nb,Ns_Chao,Ns_Jackknife,kbar,Vk)

df <- Nb_Ns1 %>% 
  full_join(env,by="Pop")

df1 <- df %>% 
  mutate(NbNc_LD=LD_Nb/NcAge1,NbNc_SF=SF_Nb/NcAge1,NbNc_PwoP=PwoP_Nb/NcAge1)

#try without Brule river due to lampricide treatment

#generate models for a single factor
factors <- df1 %>% 
  select(NbNc_PwoP,NbNc_LD,NbNc_SF,Ns_Chao,Ns_Jackknife,kbar,Vk)
responses <- df1 %>% 
  select(YearSinceTreat,Drainage,TrapToMouthKm,SampDist,SampSize,NcAge1,LarvalPop)

##check colinearity of variables
df2 <- cbind(df1$Vk,responses)
df2 <- df2 %>% rename(Vk='df1$Vk')
plot(responses)

#check normality of variables
#!# if not normal, transformations attempted
shapiro.test(log(df2$Vk))
shapiro.test(log(df2$Drainage))
shapiro.test(log(df2$TrapToMouthKm))
shapiro.test(df2$SampDist)
shapiro.test(log(df2$NcAge1))
shapiro.test(log(df2$LarvalPop))
##transform data 
df2$Vk_log <- log(df2$Vk)
df2$Drainage_log <- log(df2$Drainage)
df2$TrapToMouthKm_log <- log(df2$TrapToMouthKm)
df2$NcAge1_log <- log(df2$NcAge1)
df2$LarvalPop_log <- log(df2$LarvalPop)

#Vk Models ####
#subset data to remove locations with missing data
rldata <- df2[which(!(is.na(df2$TrapToMouthKm) | is.na(df2$NcAge1))),]
#subset df2 to just response and factors
rldata <- rldata %>% select(Vk,NcAge1,SampSize,TrapToMouthKm,Drainage,YearSinceTreat)
varmtx <- model.matrix(Vk~.-1, data=rldata)
response <- rldata$Vk

##lasso models
# alpha=1 means lasso regression. 
lasso <- glmnet(scale(varmtx), response,nlambda = 1000,standardize = F)
plot(lasso, xvar = "lambda", label = TRUE)
plot(lasso, xvar = "dev", label = TRUE)
print(lasso)

#coefficients
coef(lasso,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.lasso <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.lasso)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.1se")

plot(lasso, xvar = "lambda", label=T)
lbs_fun(ridge, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

##elastic net model
enet <- glmnet(scale(varmtx), response,alpha = 0.2)
plot(enet, xvar = "lambda", label = TRUE)
plot(enet, xvar = "dev", label = TRUE)
print(enet)

#coefficients
coef(enet,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.net <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.net)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(enet, xvar = "lambda", label=T)
lbs_fun(enet, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

#Nb Models - LD Method####
#try with and without larval population variable
#subset data to remove locations with missing data
#subset df2 to just response and factors
rldata <- df1 %>% select(LD_Nb,SampSize,Drainage,YearSinceTreat,LarvalPop)
rldata <- rldata[which(!(is.na(rldata$LarvalPop))),]

varmtx <- model.matrix(LD_Nb~.-1, data=rldata)
response <- rldata$LD_Nb

##lasso models
# alpha=1 means lasso regression. 
lasso <- glmnet(scale(varmtx), response,nlambda = 1000,standardize = T)
plot(lasso, xvar = "lambda", label = TRUE)
plot(lasso, xvar = "dev", label = TRUE)
print(lasso)

#coefficients
coef(lasso,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.lasso <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.lasso)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(lasso, xvar = "lambda", label=T)
lbs_fun(lasso, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

##elastic net model
enet <- glmnet(scale(varmtx), response,alpha = 0.2)
plot(enet, xvar = "lambda", label = TRUE)
plot(enet, xvar = "dev", label = TRUE)
print(enet)

#coefficients
coef(enet,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.net <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.net)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(enet, xvar = "lambda", label=T)
lbs_fun(enet, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

#Nb Models - PwoP Method####
#try with and without larval population variable
#subset data to remove locations with missing data
#subset df2 to just response and factors
rldata <- df1 %>% select(PwoP_Nb,SampSize,Drainage,YearSinceTreat)
#rldata <- rldata[which(!(is.na(rldata$LarvalPop))),]

varmtx <- model.matrix(PwoP_Nb~.-1, data=rldata)
response <- rldata$PwoP_Nb

##lasso models
# alpha=1 means lasso regression. 
lasso <- glmnet(scale(varmtx), response,nlambda = 1000,standardize = T)
plot(lasso, xvar = "lambda", label = TRUE)
plot(lasso, xvar = "dev", label = TRUE)
print(lasso)

#coefficients
coef(lasso,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.lasso <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.lasso)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(lasso, xvar = "lambda", label=T)
lbs_fun(lasso, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

##elastic net model
enet <- glmnet(scale(varmtx), response,alpha = 0.2)
plot(enet, xvar = "lambda", label = TRUE)
plot(enet, xvar = "dev", label = TRUE)
print(enet)

#coefficients
coef(enet,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.net <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.net)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(enet, xvar = "lambda", label=T)
lbs_fun(enet, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)



#Ns Models - Chao Method####
#try with and without larval population variable
#subset data to remove locations with missing data
#subset df2 to just response and factors
rldata <- df1 %>% select(Ns_Chao,SampSize,Drainage,YearSinceTreat,LarvalPop)
rldata <- rldata[which(!(is.na(rldata$LarvalPop))),]

varmtx <- model.matrix(Ns_Chao~.-1, data=rldata)
response <- rldata$Ns_Chao

##lasso models
# alpha=1 means lasso regression. 
lasso <- glmnet(scale(varmtx), response,nlambda = 1000,standardize = T)
plot(lasso, xvar = "lambda", label = TRUE)
plot(lasso, xvar = "dev", label = TRUE)
print(lasso)

#coefficients
coef(lasso,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.lasso <- cv.glmnet(varmtx, response,type.measure = "mse", nfolds = 12,grouped = F)
plot(cv.lasso)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(lasso, xvar = "lambda", label=T)
lbs_fun(lasso, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

##elastic net model
enet <- glmnet(scale(varmtx), response,alpha = 0.2)
plot(enet, xvar = "lambda", label = TRUE)
plot(enet, xvar = "dev", label = TRUE)
print(enet)

#coefficients
coef(enet,s=0.1)
# Cross validation to find the optimal lambda penalization
cv.net <- cv.glmnet(varmtx, response,grouped = F)
plot(cv.net)
opt.lam = c(cv.lasso$lambda.min, cv.lasso$lambda.1se)
coef(cv.lasso, s = "lambda.min")

plot(enet, xvar = "lambda", label=T)
lbs_fun(enet, offset_x = -2)
abline(v=cv.lasso$lambda.min, col = "red", lty=2)
abline(v=cv.lasso$lambda.1se, col="blue", lty=2)

