#Family Diagrams for all locations

#load libraries
library(tidyverse)
#load functions
source("Homebrew/family.plot.R")
source("Homebrew/multiplot.R")
source("Homebrew/Ns_calc.R")
#load Nb data and length data for plot annotation
Nb_Ns <- read.table("Output/genetic.estimates.txt",header = T,sep = "\t")
load("Summaries/Alllocs_cohorts.rda")

#remove Grand River from length data since it was not sequenced
#make vectors of locs and plotnames
locs <- unique(alllocs_cohorts$loc)
plotnames <- c("Bad River",
               "Betsie River",
               "Betsy River",
               "Brule River",
               "Cattaraugus River",
               "Pigeon River",
               "East Au Gres River",
               "Ford River",
               "Manistique River",
               "Manistee River",
               "Middle River",
               "Misery River",
               "Muskegon River",
               "Ocqueoc River",
               "Sterling River",
               "Swan Creek",
               "Tahquamenon River",
               "Two-Hearted River")
names(plotnames) <- locs

#loop to generate pedigree plots
tiff(filename = "Figures/FamilyPlots.tiff",height = 15,width = 12,units = "in",res = 400)
par(mfrow = c(5,4),mai = c(0.3, 0.3, 0.5, 0.5))

for (i in 1:length(locs)) {
  #read in file
  tmp <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  tmp <- strsplit(tmp,"\\s+")
  tmp1 <- matrix(unlist(tmp),byrow = T)
  tmp1 <- tmp1[tmp1 != ""]
  tmp1 <- matrix(tmp1,ncol = 5,byrow = T)
  tmp1 <- as.data.frame(tmp1)
  colnames(tmp1) <- tmp1[1,]
  tmp1 <- tmp1[-1,]
  #plot pedigree
  family.plot(tmp1,title = plotnames[i])
}
dev.off()

#loop to generate boxplots
boxplots <- list()
i <- 1
plotn <- 0
for (i in 1:length(locs)) {
  #read in file
  tmp <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  tmp <- strsplit(tmp,"\\s+")
  tmp1 <- matrix(unlist(tmp),byrow = T)
  tmp1 <- tmp1[tmp1 != ""]
  tmp1 <- matrix(tmp1,ncol = 5,byrow = T)
  tmp1 <- as.data.frame(tmp1)
  colnames(tmp1) <- tmp1[1,]
  tmp1 <- tmp1[-1,]
  
  #combine family data with length data
  tmp1 <- merge(tmp1,alllocs_cohorts)
  #remove families with less than 3 individuals
  #tmp1 <- tmp1[tmp1$ClusterIndex %in% names(which(table(tmp1$ClusterIndex)>=3)),]
  #check to make sure there's more than one family
  #plot pedigree
  if(length(unique(tmp1$ClusterIndex)) > 1){
    plotn <- plotn + 1
    boxplots[[plotn]] <- ggplot(tmp1,aes(x=ClusterIndex,group=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
      geom_boxplot(alpha=0.3,show.legend = F)+
      theme_bw()+
      scale_fill_gradient(low="red", high="white",name = "Cluster \n Likelihood")+
      xlab("Cluster")+
      ylab("Length (mm)")+
      ggtitle(plotnames[i])
  }
}
tiff(filename = "Figures/LengthBoxPlotsNoFamLimit011722.tiff",height = 12,width = 15,units = "in",res = 200)
multiplot(cols = 4,boxplots[[1]],boxplots[[5]],boxplots[[9]],
          boxplots[[13]],boxplots[[2]],boxplots[[6]],
          boxplots[[10]],boxplots[[14]],boxplots[[17]],boxplots[[3]],
          boxplots[[7]],boxplots[[11]],boxplots[[15]],boxplots[[18]],
          boxplots[[4]],boxplots[[8]],boxplots[[12]],
          boxplots[[16]])
dev.off()

#loop to generate Ns accumulation curves
#calculation loop
Nsplots <- list()
for (i in 1:length(locs)) {
  print(i)
  #read in file
  df <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
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
  Nsplots[[i]] <- ggplot(df2,aes(x=sites,y=richness))+
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
    ylab("Number of parent genotypes")+
    ylim(0,max(Ns_tmp[[2]]$chao,Ns_tmp[[2]]$jack1)+5)+
    ggtitle(plotnames[i])
  
  
}

tiff(filename = "Figures/NsAccumPlots.tiff",height = 14,width = 9,units = "in",res = 400)
multiplot(cols = 3,Nsplots[[1]],Nsplots[[4]],Nsplots[[7]],
          Nsplots[[10]],Nsplots[[13]],Nsplots[[16]],
          Nsplots[[2]],Nsplots[[5]],Nsplots[[8]],
          Nsplots[[11]],Nsplots[[14]],Nsplots[[17]],
          Nsplots[[3]],Nsplots[[6]],Nsplots[[9]],
          Nsplots[[12]],Nsplots[[15]],Nsplots[[18]])
dev.off()
