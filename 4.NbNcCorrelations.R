##Nb:Nc Ratio correlations and plots
#!# This script is still a bit chaos, not sure how helpful it will be to other people

#Goals:

#library load
library(tidyverse)
library(ggpubr)

#load data
df <- read.csv("Models/ModelInputAll_120222.csv")
#remove NAs from data set
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

#checking scatterplots before calculating models
plot(df$NcAge1,df$Nb_LD)
plot(df$NcAge1,df$Nb_SF)
plot(df$NcAge1,df$Ns_Chao)

#calculating models with full data set:
cor.test(df$Nb_LD,df$NcAge1)
cor.test(df$Nb_SF,df$NcAge1)
cor.test(df$Ns_Chao,df$NcAge1)
cor.test(df$Nb_SF,df$Ns_Chao)
#no model significance

#Spearman Correlations
cor.test(df$Nb_LD,df$NcAge1,method = c("spearman"))
cor.test(df$Nb_SF,df$NcAge1,method = c("spearman"))
cor.test(df$Ns_Chao,df$NcAge1,method = c("spearman"))

#remove outlier Brule - age 1 birth year was a treatment year - collection is much lower than Nc 2018 estimate
df1 <- df %>% 
  slice(-4)

#checking scatterplots before calculating models
plot(df1$NcAge1,df1$Nb_LD)
plot(df1$NcAge1,df1$Nb_SF)
plot(df1$NcAge1,df1$Ns_Chao)
plot(df1$Nb_SF,df1$Ns_Chao)
#calculating models with full data set:
cor.test(df1$Nb_LD,df1$NcAge1,alternative = "g")
cor.test(df1$Nb_SF,df1$NcAge1,alternative = "g")
cor.test(df1$Ns_Chao,df1$NcAge1,alternative = "g")
cor.test(df1$Nb_SF,df1$Ns_Chao,alternative = "g")

#make correlation result labels for plot
reslabels <- data.frame(label = c("cor = 0.122", 
                                  "cor = 0.138", 
                                  "cor = 0.192"),
                        estimate = c("Nb_LD","Nb_SF","Ns_Chao"))
#no model significance
corplot1 <- df1 %>% 
  select(Pop,NcAge1,Nb_LD,Nb_SF,Ns_Chao) %>% 
  gather(key = "estimate",value = "value",-Pop,-NcAge1)

plotnames <- c("Nb - Linkage Disequilibrium", "Nb - Sibship Frequency", "Extrapolated Ns")
names(plotnames) <- unique(corplot1$estimate) 

p1 <- ggplot(corplot1, aes(x = value,y = NcAge1))+
  facet_wrap(~estimate,scales = "free_x",labeller = as_labeller(plotnames))+
  geom_point()+
  theme_bw()+
  stat_cor(mapping = aes(x = value, y = NcAge1),alternative = "g",cor.coef.name = "R")+
  ggtitle("All Streams (N = 10)")+
  labs(x = "Nb or Ns Point Estimate",y = "Mark-recapture estimate (Nc)")

#Spearman Correlations
cor.test(df1$Nb_LD,df1$NcAge1,method = c("spearman"))
cor.test(df1$Nb_SF,df1$NcAge1,method = c("spearman"))
cor.test(df1$Ns_Chao,df1$NcAge1,method = c("spearman"))
cor.test(df1$Nb_SF,df1$Ns_Chao,method = c("spearman"))

#remove locations that meet the following criteria:
#1.larval sample size < 50

Nb_Ns2 <- df1[which(!(df1$SampSize < 50)),]

#checking scatterplots before calculating models
plot(Nb_Ns2$NcAge1,Nb_Ns2$Nb_LD)
plot(Nb_Ns2$NcAge1,Nb_Ns2$Nb_SF)
plot(Nb_Ns2$NcAge1,Nb_Ns2$Ns_Chao)
plot(Nb_Ns2$Nb_SF,Nb_Ns2$Ns_Chao)

shapiro.test(Nb_Ns2$LD_Nb)
shapiro.test(Nb_Ns2$SF_Nb)
shapiro.test(Nb_Ns2$Ns_Chao)
shapiro.test(Nb_Ns2$Ns_Jackknife)
shapiro.test(Nb_Ns2$NcAge1)

cor.test(Nb_Ns2$Nb_LD,Nb_Ns2$NcAge1,alternative = "g")
cor.test(Nb_Ns2$Nb_SF,Nb_Ns2$NcAge1,alternative = "g")
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1,alternative = "g")
cor.test(Nb_Ns2$Nb_SF,Nb_Ns2$Ns_Chao,alternative = "g")

#Spearman Correlations
cor.test(Nb_Ns2$Nb_LD,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Nb_SF,Nb_Ns2$NcAge1,method = c("spearman"))
cor.test(Nb_Ns2$Ns_Chao,Nb_Ns2$NcAge1,method = c("spearman"))

corplot2 <- Nb_Ns2 %>% 
  select(Pop,NcAge1,Nb_LD,Nb_SF,Ns_Chao) %>% 
  gather(key = "estimate",value = "value",-Pop,-NcAge1)

#make correlation result labels for plot
reslabels2 <- data.frame(label = c("cor = 0.652", "cor = 0.649", "cor = 0.738"),
                        estimate = c("Nb_LD","Nb_SF","Ns_Chao"))

p2 <- ggplot(corplot2,aes(x = value,y = NcAge1))+
  facet_wrap(~estimate,scales = "free_x",labeller = as_labeller(plotnames))+
  geom_point()+
  geom_smooth(method = lm,mapping = aes(x = value, y = NcAge1))+
  stat_cor(mapping = aes(x = value, y = NcAge1),alternative = "g",cor.coef.name = "R")+
  theme_bw()+
  ggtitle("Streams with Sufficient Samples (n > 50 samples; N = 7)")+
  labs(x = "Nb or Ns Point Estimate",y = "Mark-recapture estimate (Nc)")

#putting the plots together
tiff(filename = "Figures/CorrelationPlot_030123.tiff",width = 12,height = 8,units = "in",res = 100)
ggarrange(p1,p2, 
          labels = c("A", "B"),
          ncol = 1)
dev.off()
#model significance

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
