##Nb:Nc Ratio correlations and plots

#Goals:

#library load
library(tidyverse)
#data load
Nc <- read.table("Models/env_data_all_locations.txt",header = T,sep = "\t")
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
#Nb:Nc ratio calculations and correlations
Nb_Ns1 <- Nb_Ns %>% 
  select(Pop,LD_Nb,PwoP_Nb,SF_Nb,Ns_Chao,Ns_Jackknife,kbar,Vk)

Nc <- Nc %>% 
  select(Pop,NcAge1,AdultRemove)

df <- merge(Nb_Ns1,Nc)
df <- df[-which(is.na(df$NcAge1)),]

#checking normality
#if normality assumption is violated, transformations are tested
shapiro.test(log(df$LD_Nb))
shapiro.test(log(df$SF_Nb))
shapiro.test(log(df$PwoP_Nb))
shapiro.test(df$Ns_Chao)
shapiro.test(log(df$Ns_Jackknife))
shapiro.test(df$kbar)
shapiro.test(df$Vk)
shapiro.test(df$NcAge1)


df$LD_Nb_log <- log(df$LD_Nb)
df$SF_Nb_log <- log(df$SF_Nb)
df$PwoP_Nb_log <- log(df$PwoP_Nb)
df$Ns_Chao_log <- log(df$Ns_Chao)
df$Ns_Jackknife_log <- log(df$Ns_Jackknife)
df$NcAge1_log <- log(df$NcAge1)

#checking scatterplots before calculating models
plot(df$NcAge1,df$LD_Nb)
plot(df$NcAge1_log,df$PwoP_Nb_log)
plot(df$NcAge1_log,df$SF_Nb_log)
plot(df$NcAge1_log,df$Ns_Chao_log)
plot(df$NcAge1_log,df$Ns_Jackknife_log)

#remove locations that meet the following criteria:
#1.larval sample size < 50

Nb_Ns2 <- Nb_Ns[which(!(Nb_Ns$SampSize < 50)),]
Nb_Ns2 <- merge(Nb_Ns2,Nc)
Nb_Ns2 <- Nb_Ns2[-which(is.na(Nb_Ns2$NcAge1)),]

shapiro.test(Nb_Ns2$LD_Nb)
shapiro.test(Nb_Ns2$SF_Nb)
shapiro.test(Nb_Ns2$PwoP_Nb)
shapiro.test(Nb_Ns2$Ns_Chao)
shapiro.test(Nb_Ns2$Ns_Jackknife)
shapiro.test(Nb_Ns2$NcAge1)

cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1,method = c("spearman"))


## correction for caught lamprey
#checking normality
#if normality assumption is violated, transformations are tested
dfcatch <- df %>% 
  select(Pop,LD_Nb_log,SF_Nb_log,PwoP_Nb_log,Ns_Jackknife_log,Ns_Chao_log,NcAge1,AdultRemove) %>% 
  mutate(corNc=NcAge1-AdultRemove)


shapiro.test(log(df$LD_Nb))
shapiro.test(log(df$SF_Nb))
shapiro.test(log(df$PwoP_Nb))
shapiro.test(log(df$Ns_Chao))
shapiro.test(log(df$Ns_Jackknife))
shapiro.test(log(df$kbar))
shapiro.test(log(df$Vk))
shapiro.test(dfcatch$corNc)

#checking scatterplots before calculating models
plot(dfcatch$corNc,dfcatch$LD_Nb_log)
plot(dfcatch$corNc,dfcatch$PwoP_Nb_log)
plot(dfcatch$corNc,dfcatch$SF_Nb_log)
plot(dfcatch$corNc,dfcatch$Ns_Chao_lo)
plot(dfcatch$corNc,dfcatch$PwoP_Nb_log)


#checking scatterplots before calculating models
plot(df$NcAge1,df$LD_Nb)
plot(df$NcAge1_log,df$PwoP_Nb_log)
plot(df$NcAge1_log,df$SF_Nb_log)
plot(df$NcAge1_log,df$Ns_Chao_log)
plot(df$NcAge1_log,df$Ns_Jackknife_log)

#remove locations that meet the following criteria:
#1.larval sample size < 50

Nb_Ns2 <- Nb_Ns[which(!(Nb_Ns$SampSize < 50)),]
Nb_Ns2 <- merge(Nb_Ns2,Nc)
Nb_Ns2 <- Nb_Ns2[-which(is.na(Nb_Ns2$NcAge1)),]

shapiro.test(Nb_Ns2$LD_Nb)
shapiro.test(Nb_Ns2$SF_Nb)
shapiro.test(Nb_Ns2$PwoP_Nb)
shapiro.test(Nb_Ns2$Ns_Chao)
shapiro.test(Nb_Ns2$Ns_Jackknife)
shapiro.test(Nb_Ns2$NcAge1)

cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1)
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(Nb_Ns2$LD_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$PwoP_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$SF_Nb,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Jackknife,Nb_Ns2$NcAge1,method = c("spearman"))


## correction for caught lamprey
#checking normality
#if normality assumption is violated, transformations are tested
dfcatch <- df %>% 
  select(Pop,LD_Nb_log,SF_Nb_log,PwoP_Nb_log,Ns_Jackknife_log,Ns_Chao_log,NcAge1,AdultRemove) %>% 
  mutate(corNc=NcAge1-AdultRemove)


shapiro.test(log(df$LD_Nb))
shapiro.test(log(df$SF_Nb))
shapiro.test(log(df$PwoP_Nb))
shapiro.test(log(df$Ns_Chao))
shapiro.test(log(df$Ns_Jackknife))
shapiro.test(log(df$kbar))
shapiro.test(log(df$Vk))
shapiro.test(dfcatch$corNc)

dfcatch$corNc_log <- log(dfcatch$corNc)
#checking scatterplots before calculating models
plot(df$NcAge1,df$LD_Nb)
plot(df$NcAge1,df$PwoP_Nb)
plot(df$NcAge1,df$SF_Nb)
plot(df$NcAge1,df$Ns_Chao)
plot(df$NcAge1,df$Ns_Jackknife)
plot(df$NcAge1,df$kbar)
plot(df$NcAge1,df$Vk)

#remove Brule due to 2018 lampricide between spawning and collection
cor.test(dfcatch$LD_Nb_log,dfcatch$corNc)
cor.test(dfcatch$PwoP_Nb_log,dfcatch$corNc)
cor.test(dfcatch$SF_Nb_log,dfcatch$corNc)
cor.test(dfcatch$Ns_Chao_log,dfcatch$corNc)
cor.test(dfcatch$Ns_Jackknife_log,dfcatch$corNc)
cor.test(dfcatch$PwoP_Nb_log,dfcatch$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(dfcatch$LD_Nb_log,dfcatch$corNc_log,method = c("spearman"))
cor.test(dfcatch$PwoP_Nb_log,dfcatch$corNc_log,method = c("spearman"))
cor.test(dfcatch$SF_Nb_log,dfcatch$corNc_log,method = c("spearman"))
cor.test(dfcatch$Ns_Chao_log,dfcatch$corNc_log,method = c("spearman"))
cor.test(dfcatch$Ns_Jackknife_log,dfcatch$corNc_log,method = c("spearman"))

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

tiff(filename = "Figures/NbNcCorrelationPlot_correction.tiff",height = 5,width = 7,units = "in",res = 200)
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
