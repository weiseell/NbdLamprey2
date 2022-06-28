#Alternative Plot ideas

#load data, packages, and source code
#load libraries
library(tidyverse)
library(ggpubr)
library(wesanderson)

#load functions
source("Homebrew/family.plot.R")
source("Homebrew/multiplot.R")
source("Homebrew/Ns_calc.R")
#load Nb data and length data for plot annotation
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
load("Summaries/Alllocs_cohorts.rda")

##Making a subplot for the Middle River with scatter plot, boxplots, family diagram, and Ns figure ####
#Make one for Middle, Muskegon, Bad, Sterling?
tmp <- subset(alllocs_cohorts,alllocs_cohorts$loc == "MIR")
#making scatter plot
scatter <- ggplot(tmp,aes(x = Length,y = Weight,color = as.factor(Class)))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = wes_palette("IsleofDogs1"))+
  ylab("Weight (g)")+
  xlab("Length (mm)")+
  ggtitle("Length and Weight with Mixture Analysis Cohorts")+
  theme(legend.position = "none")

#making boxplots
#make a full-sibling family column
tmp <- tmp %>% mutate(fullsib = paste(MotherID,FatherID,sep = "_"))
#remove families with less than 3 individuals
tmp1 <- tmp[tmp$fullsib %in% names(which(table(tmp$fullsib)>=3)),]
boxplots <- ggplot(tmp1,aes(x=Length,group=fullsib,y=fullsib))+
  geom_boxplot(alpha=0.3,show.legend = F)+
  geom_jitter(aes(color = as.character(Class)),position=position_jitter(0.1))+
  scale_color_manual(values = wes_palette("IsleofDogs1"))+
  theme_bw()+
  #scale_fill_gradient(low="red", high="white",name = "Cluster \n Likelihood")+
  ylab("Full-Sibling Family")+
  xlab("Length (mm)")+
  ggtitle("Length of Larvae sorted by full-sibling family")+
  theme(legend.position = "none")

#writing Ns figure
#mir 2015
#read in file
df <- readLines(paste0("SoftwareOutput/MIR2015.Output.data.BestCluster"))
#separate file into usable data frame
df <- strsplit(df,"\\s+")
df1 <- matrix(unlist(df),byrow = T)
df1 <- df1[df1 != ""]
df1 <- matrix(df1,ncol = 5,byrow = T)
df1 <- as.data.frame(df1)
colnames(df1) <- df1[1,]
df1 <- df1[-1,]

#calculate extrapolated Ns
Ns_tmp <- Ns_calc(df1)

#isolating data to plot
df2 <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps <- cbind(df2$sites,reps)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 1),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 1)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
ns2015 <- ggplot(df2,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao,col = "darkred")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=-1,size = 3)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=1.5,size = 3)+
  #geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=1.5,size = 3)+
  xlab("Number of offspring sampled")+
  ylab("Ns")+
  ylim(0,max(Ns_tmp[[2]]$chao,Ns_tmp[[2]]$jack1)+5)+
  ggtitle("Ns - Middle River 2015")

#mir 2016
#read in file
df <- readLines(paste0("SoftwareOutput/MIR2016.Output.data.BestCluster"))
#separate file into usable data frame
df <- strsplit(df,"\\s+")
df1 <- matrix(unlist(df),byrow = T)
df1 <- df1[df1 != ""]
df1 <- matrix(df1,ncol = 5,byrow = T)
df1 <- as.data.frame(df1)
colnames(df1) <- df1[1,]
df1 <- df1[-1,]

#calculate extrapolated Ns
Ns_tmp <- Ns_calc(df1)

#isolating data to plot
df2 <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps <- cbind(df2$sites,reps)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 1),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 1)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
ns2016 <- ggplot(df2,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao,col = "darkred")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=-1,size = 3)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=1.5,size = 3)+
  #geom_hline(yintercept=Ns_tmp[[2]]$jack1,col = "darkblue")+
  #geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=1.5,size = 3)+
  xlab("Number of offspring sampled")+
  ylab("Ns")+
  ylim(0,max(Ns_tmp[[2]]$chao,Ns_tmp[[2]]$jack1)+5)+
  ggtitle("Ns - Middle River 2016")

#putting the plots together
scatter
boxplots
reconped
ns2015
ns2016

ggarrange(scatter,boxplots, labels = c("A","B"),
          ggarrange(ns2015,ns2016,labels = c("C","D"),ncol = 2),
             nrow=3)




#making family diagram
fam <- readLines(paste0("SoftwareOutput/MIR.Output.data.BestCluster"))
#separate file into usable data frame
fam <- strsplit(fam,"\\s+")
fam <- matrix(unlist(fam),byrow = T)
fam <- fam[fam != ""]
fam <- matrix(fam,ncol = 5,byrow = T)
fam <- as.data.frame(fam)
colnames(fam) <- fam[1,]
fam <- fam[-1,]
#plot pedigree
family.plot(fam,title = "Reconstructed Pedigree")

##Adding all the Ns Accumlations as single lines? ####



























