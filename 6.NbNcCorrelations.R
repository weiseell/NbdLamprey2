##Nb:Nc Ratio correlations and plots

#Goals:

#library load

#data load
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")

#Nb:Nc ratio calculations and correlations
Nb_Ns1 <- Nb_Ns %>% 
  select(Pop,LD_Nb,PwoP_Nb,SF_Nb,Ns_Chao,Ns_Jackknife,kbar,Vk)

Nc <- Nc %>% 
  select(Pop,ModelNc)

df <- merge(Nb_Ns1,Nc)
df <- df %>% 
  filter(ModelNc != "---")
df$ModelNc <- as.numeric(df$ModelNc)
df$SF_Nb <- as.numeric(df$SF_Nb)
df$PwoP_Nb <- as.numeric(df$PwoP_Nb)
df$Ns_Chao <- as.numeric(df$Ns_Chao)
df$Ns_Jackknife <- as.numeric(df$Ns_Jackknife)
df$kbar <- as.numeric(df$kbar)
df$Vk <- as.numeric(df$Vk)

#checking scatterplots before calculating models
plot(df$ModelNc,df$LD_Nb)
plot(df$ModelNc,df$PwoP_Nb)
plot(df$ModelNc,df$SF_Nb)
plot(df$ModelNc,df$Ns_Chao)
plot(df$ModelNc,df$Ns_Jackknife)
plot(df$ModelNc,df$kbar)
plot(df$ModelNc,df$Vk)



model_LD <- lm(LD_Nb~ModelNc,data = df)
model_SF <- lm(SF_Nb~ModelNc,data = df)
model_PwoP <- lm(PwoP_Nb~ModelNc,data = df)
model_Chao <- lm(Ns_Chao~ModelNc,data = df)
model_Jack <- lm(Ns_Chao~ModelNc,data = df)
model_kbar <- lm(kbar~ModelNc,data = df)
model_Vk <- lm(Vk~ModelNc,data = df)

summary(model_LD)
summary(model_SF)
summary(model_PwoP)
summary(model_Chao)
summary(model_Jack)
summary(model_kbar)
summary(model_Vk)


#checking scatterplots before calculating models
plot(df$ModelNc,df$LD_Nb)
abline(model_LD)
plot(df$ModelNc,df$PwoP_Nb)
abline(model_PwoP)
plot(df$ModelNc,df$SF_Nb)
abline(model_SF)
plot(df$ModelNc,df$Ns_Chao)
abline(model_Chao)
plot(df$ModelNc,df$Ns_Jackknife)
abline(model_Jack)
plot(df$ModelNc,df$kbar)
abline(model_kbar)
plot(df$ModelNc,df$Vk)
abline(model_Vk)







