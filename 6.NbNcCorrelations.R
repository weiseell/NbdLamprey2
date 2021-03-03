##Nb:Nc Ratio correlations and plots

#Goals:

#library load

#data load
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")

#Nb:Nc ratio calculations and correlations
Nb_Ns1 <- Nb_Ns %>% 
  select(Pop,LD_Nb,PwoP_Nb,SF_Nb,Ns_Chao,Ns_Jackknife,kbar,Vk)

Nc <- Nc %>% 
  select(Pop,NcAge1)

df <- merge(Nb_Ns1,Nc)
df <- df[-which(is.na(df$NcAge1)),]

#checking scatterplots before calculating models
plot(df$NcAge1,df$LD_Nb)
plot(df$NcAge1,df$PwoP_Nb)
plot(df$NcAge1,df$SF_Nb)
plot(df$NcAge1,df$Ns_Chao)
plot(df$NcAge1,df$Ns_Jackknife)
plot(df$NcAge1,df$kbar)
plot(df$NcAge1,df$Vk)

#remove Brule due to 2018 lampricide between spawning and collection
df <- df[-5,]

model_LD <- cor.test(df$NcAge1,df$LD_Nb)
model_SF <- cor.test(df$NcAge1,df$SF_Nb)
model_PwoP <- cor.test(df$NcAge1,df$PwoP_Nb)
model_Chao <- cor.test(df$NcAge1,df$Ns_Chao)
model_Jack <- cor.test(df$NcAge1,df$Ns_Jack)
model_kbar <- cor.test(df$NcAge1,df$kbar)
model_Vk <- cor.test(df$NcAge1,df$Vk)

#no model significance

#Spearman Correlations
cor.test(df$NcAge1,df$LD_Nb,method = c("spearman"))
cor.test(df$NcAge1,df$SF_Nb,method = c("spearman"))
cor.test(df$NcAge1,df$PwoP_Nb,method = c("spearman"))
cor.test(df$NcAge1,df$Ns_Chao,method = c("spearman"))
cor.test(df$NcAge1,df$Ns_Jackknife,method = c("spearman"))
cor.test(df$NcAge1,df$kbar,method = c("spearman"))
cor.test(df$NcAge1,df$Vk,method = c("spearman"))
