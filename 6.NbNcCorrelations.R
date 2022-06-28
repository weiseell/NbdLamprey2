##Nb:Nc Ratio correlations and plots

#Goals:

#library load
library(tidyverse)
#data load
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")

#Nb:Nc ratio calculations and correlations
Nb_Ns1 <- Nb_Ns %>% 
  select(Pop,Nb_LD,SF_Nb,Ns_Chao,Ns_Jackknife,kbar,Vk,SampSize) %>% 
  rename(LD_Nb = Nb_LD,genest = Pop) %>% 
  mutate(Pop1 = genest) %>% 
  separate(Pop1,c("Pop","num"),sep = 3) %>%
  select(-num) %>% 
  slice(-c(11,16,21))

Nc <- Nc %>% 
  select(Pop,NcAge1,AdultRemove)

df <- Nb_Ns1 %>% full_join(Nc,by = "Pop")
df$SFNc_ratio <- df$SF_Nb/df$NcAge1
df$LDNc_ratio <- df$LD_Nb/df$NcAge1
df <- df[-which(is.na(df$NcAge1)),]

#checking normality
#if normality assumption is violated, transformations are tested
shapiro.test(log(df$LD_Nb))
shapiro.test(log(df$SF_Nb))
shapiro.test(df$Ns_Chao)
shapiro.test(log(df$Ns_Jackknife))
shapiro.test(df$kbar)
shapiro.test(df$Vk)
shapiro.test(df$NcAge1)


df$LD_Nb_log <- log(df$LD_Nb)
df$SF_Nb_log <- log(df$SF_Nb)
df$Ns_Chao_log <- log(df$Ns_Chao)
df$Ns_Jackknife_log <- log(df$Ns_Jackknife)
df$NcAge1_log <- log(df$NcAge1)

#checking scatterplots before calculating models
plot(df$NcAge1,df$LD_Nb)
plot(df$NcAge1,df$SF_Nb)
plot(df$NcAge1,df$Ns_Chao)
plot(df$NcAge1,df$Ns_Jackknife)

#calculating models with full data set:
cor.test(df$LD_Nb,df$NcAge1)
cor.test(df$SF_Nb,df$NcAge1)
cor.test(df$Ns_Chao,df$NcAge1)
cor.test(df$Ns_Jackknife,df$NcAge1)
cor.test(df$SF_Nb,df$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(df$LD_Nb,df$NcAge1,method = c("spearman"))
cor.test(df$SF_Nb,df$NcAge1,method = c("spearman"))
cor.test(df$Ns_Chao,df$NcAge1,method = c("spearman"))
cor.test(df$Ns_Jackknife,df$NcAge1,method = c("spearman"))


#remove outlier Brule - age 1 birth year was a treatment year - collection is much lower than Nc 2018 estimate
df1 <- df %>% 
  slice(-4)

#checking scatterplots before calculating models
plot(df1$NcAge1,df1$LD_Nb)
plot(df1$NcAge1,df1$SF_Nb)
plot(df1$NcAge1,df1$Ns_Chao)
plot(df1$NcAge1,df1$Ns_Jackknife)
plot(df1$SF_Nb,df1$Ns_Chao)
#calculating models with full data set:
cor.test(df1$LD_Nb,df1$NcAge1)
cor.test(df1$SF_Nb,df1$NcAge1)
cor.test(df1$Ns_Chao,df1$NcAge1)
cor.test(df1$Ns_Jackknife,df1$NcAge1)
cor.test(df1$SF_Nb,df1$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(df1$LD_Nb,df1$NcAge1,method = c("spearman"))
cor.test(df1$SF_Nb,df1$NcAge1,method = c("spearman"))
cor.test(df1$Ns_Chao,df1$NcAge1,method = c("spearman"))
cor.test(df1$Ns_Jackknife,df1$NcAge1,method = c("spearman"))

#remove locations that meet the following criteria:
#1.larval sample size < 50

Nb_Ns2 <- df1[which(!(df1$SampSize < 50)),]

shapiro.test(Nb_Ns2$LD_Nb)
shapiro.test(Nb_Ns2$SF_Nb)
shapiro.test(Nb_Ns2$Ns_Chao)
shapiro.test(Nb_Ns2$Ns_Jackknife)
shapiro.test(Nb_Ns2$NcAge1)

cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1,method = c("spearman"))

## correction for caught lamprey
#checking normality
#if normality assumption is violated, transformations are tested
dfcatch <- df1 %>% 
  select(Pop,LD_Nb,SF_Nb,Ns_Jackknife,Ns_Chao,NcAge1,AdultRemove) %>% 
  mutate(corNc=NcAge1-AdultRemove)


shapiro.test(log(df$LD_Nb))
shapiro.test(log(df$SF_Nb))
shapiro.test(log(df$Ns_Chao))
shapiro.test(log(df$Ns_Jackknife))
shapiro.test(log(df$kbar))
shapiro.test(log(df$Vk))

#checking scatterplots before calculating models
plot(dfcatch$corNc,dfcatch$LD_Nb)
plot(dfcatch$corNc,dfcatch$SF_Nb)
plot(dfcatch$corNc,dfcatch$Ns_Chao)

cor.test(dfcatch$LD_Nb,dfcatch$corNc)
cor.test(dfcatch$SF_Nb,dfcatch$corNc)
cor.test(dfcatch$Ns_Chao,dfcatch$corNc)
cor.test(dfcatch$Ns_Jackknife,dfcatch$corNc)
cor.test(dfcatch$SF_Nb,dfcatch$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(dfcatch$LD_Nb,dfcatch$corNc,method = c("spearman"))
cor.test(dfcatch$SF_Nb,dfcatch$corNc,method = c("spearman"))
cor.test(dfcatch$Ns_Chao,dfcatch$corNc,method = c("spearman"))
cor.test(dfcatch$Ns_Jackknife,dfcatch$corNc,method = c("spearman"))

##figure plotting Nb and Nc correlations
tiff(filename = "Figures/NbNcCorrelationPlot.tiff",height = 3,width = 7,units = "in",res = 200)
par(mfrow = c(1,3),mai=c(.6,.5,.4,.4),mgp=c(1.5,.5,0))
plot(df$NcAge1,df$LD_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Linkage Disequilibrium",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - LD",side = 2,line = 2.5,cex = 0.75)

plot(df$NcAge1,df$SF_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Sibship Frequency",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - SF",side = 2,line = 2.5,cex = 0.75)

plot(df$NcAge1,df$Ns_Chao,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Ns - Chao",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Ns - Chao",side = 2,line = 2.5,cex = 0.75)
dev.off()

tiff(filename = "Figures/NbNcCorrelationPlot_sampcorrection.tiff",height = 3,width = 7,units = "in",res = 200)
par(mfrow = c(1,3),mai=c(.6,.5,.4,.4),mgp=c(1.5,.5,0))
plot(Nb_Ns2$NcAge1,Nb_Ns2$LD_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Linkage Disequilibrium",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - LD",side = 2,line = 2.5,cex = 0.75)

plot(Nb_Ns2$NcAge1,Nb_Ns2$SF_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Sibship Frequency",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - SF",side = 2,line = 2.5,cex = 0.75)

plot(Nb_Ns2$NcAge1,Nb_Ns2$Ns_Chao,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Ns - Chao",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Ns - Chao",side = 2,line = 2.5,cex = 0.75)
dev.off()

tiff(filename = "Figures/NbNcCorrelationPlot_catchcorrection.tiff",height = 5,width = 7,units = "in",res = 200)
par(mfrow = c(1,3),mai=c(.6,.5,.4,.4),mgp=c(1.5,.5,0))
plot(dfcatch$corNc,dfcatch$LD_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Linkage Disequilibrium",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - LD",side = 2,line = 2.5,cex = 0.75)
plot(dfcatch$corNc,dfcatch$SF_Nb,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Nb - Sibship Frequency",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Nb - SF",side = 2,line = 2.5,cex = 0.75)
plot(dfcatch$corNc,dfcatch$Ns_Chao,xlab = "",ylab = "",xaxt = "n",yaxt = "n",main = "Ns - Chao",pch=16,col = "darkgrey")
axis(side = 2,las = 2,mgp = c(3, 0.75, 0))
axis(side = 1,las = 2,mgp = c(3, 0.75, 0))
mtext("Nc", side=1, line=3.5,cex = 0.75)
mtext("Ns - Chao",side = 2,line = 2.5,cex = 0.75)
dev.off()
