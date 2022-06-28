#Family Diagrams for all locations

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
alllocs_cohorts$ClusterIndex <- factor(alllocs_cohorts$ClusterIndex)
#loop to generate pedigree plots
tiff(filename = "Figures/FamilyPlots.tiff",height = 12,width = 15,units = "in",res = 400)
par(mfrow = c(3,6),mai = c(0.3, 0.3, 0.5, 0.5))

for (i in 1:length(locs)) {
  #read in file
  #tmp <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  #tmp <- strsplit(tmp,"\\s+")
  #tmp1 <- matrix(unlist(tmp),byrow = T)
  #tmp1 <- tmp1[tmp1 != ""]
  #tmp1 <- matrix(tmp1,ncol = 5,byrow = T)
  #tmp1 <- as.data.frame(tmp1)
  #colnames(tmp1) <- tmp1[1,]
  #tmp1 <- tmp1[-1,]
  tmp <- subset(alllocs_cohorts,alllocs_cohorts$loc == locs[i])
  #plot pedigree
  pedigree.plot(tmp,title = plotnames[i])
}

dev.off()

#loop to generate boxplots and scatterplots
boxplots <- list()
scatters <- list()
i <- 1
plotn <- 0
for (i in 1:length(locs)) {
  loc <- locs[i]
  tmp <- subset(alllocs_cohorts,alllocs_cohorts$loc == locs[i])
  #make a full-sibling family column
  tmp <- tmp %>% mutate(fullsib = paste(MotherID,FatherID,sep = "_"))
  #making the scatter plot
  scatters[[i]] <- ggplot(tmp,aes(x = Length,y = Weight,color = as.factor(Class)))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values = wes_palette("IsleofDogs1"))+
    ylab("Weight (g)")+
    xlab("Length (mm)")+
    #ggtitle(plotnames[i])+
    ggtitle("Length and Weight with \n Mixture Analysis Cohorts")+
    theme(legend.position = "none")
  #remove families with less than 3 individuals
  tmp1 <- tmp[tmp$fullsib %in% names(which(table(tmp$fullsib)>=3)),]
  #check to make sure there's more than one family
  #plot pedigree
  if(length(unique(tmp1$fullsib)) > 1){
    plotn <- plotn + 1
    boxplots[[plotn]] <- ggplot(tmp1,aes(x=fullsib,group=fullsib,y=Length))+
      geom_boxplot(alpha=0.3,show.legend = F)+
      geom_jitter(aes(color = as.character(Class)),position=position_jitter(0.1))+
      scale_color_manual(values = wes_palette("IsleofDogs1"))+
      theme_bw()+
      #scale_fill_gradient(low="red", high="white",name = "Cluster \n Likelihood")+
      xlab("Full-Sibling Family")+
      ylab("Length (mm)")+
      ggtitle("Full-sibling Boxplots")+
      #ggtitle(plotnames[i])+
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none")
  }
}

tiff(filename = "Figures/LengthBoxPlots052522.tiff",height = 12,width = 15,units = "in",res = 200)
multiplot(cols = 4,boxplots[[1]],boxplots[[5]],boxplots[[9]],
          boxplots[[13]],boxplots[[2]],boxplots[[6]],
          boxplots[[10]],boxplots[[14]],boxplots[[3]],
          boxplots[[7]],boxplots[[11]],boxplots[[15]],
          boxplots[[4]],boxplots[[8]],boxplots[[12]],
          boxplots[[16]])
dev.off()

#alternate plot method with  single data frame instead of a loop - doesn't quite work right yet
alllocs_cohorts1 <- alllocs_cohorts %>% mutate(fullsib = paste(MotherID,FatherID,sep = "_"))
alllocs_cohorts1 <- alllocs_cohorts1 %>% 
  arrange(loc,as.integer(ClusterIndex))
ggplot(alllocs_cohorts1,aes(x=Length,group=ClusterIndex,y=ClusterIndex))+
  geom_boxplot(alpha=0.3,show.legend = F)+
  geom_jitter(aes(color = as.character(Class)),position=position_jitter(0.1))+
  scale_color_manual(values = wes_palette("Darjeeling2"))+
  facet_wrap(~loc,scales = "free", nrow = 8)+
  theme_bw()+
  xlab("Cluster")+
  ylab("Length (mm)")+
  ggtitle("Boxplots")+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")

#loop to generate Ns accumulation curves
#calculation loop
locs1 <- c("BAD","BEI","BET","BRL",
                  "CAT","CHE","EAG","FOR",
                  "MAI","MAN","MIR2015","MIR2016",
                  "MIS","MUS","OCQ2018","OCQ2019","STE",
                  "SWN","TAQ","TWO2017","TWO2018")
plotnames1 <- c("Bad River",
               "Betsie River",
               "Betsy River",
               "Brule River",
               "Cattaraugus River",
               "Pigeon River",
               "East Au Gres River",
               "Ford River",
               "Manistique River",
               "Manistee River",
               "Middle River - 2015",
               "Middle River - 2016",
               "Misery River",
               "Muskegon River",
               "Ocqueoc River - 2018",
               "Ocqueoc River - 2019",
               "Sterling River",
               "Swan Creek",
               "Tahquamenon River",
               "Two-Hearted River - 2017",
               "Two-Hearted River - 2018")
names(plotnames1) <- locs1
Nsplots <- list()
for (i in 1:length(locs1)) {
  print(i)
  #read in file
  df <- readLines(paste0("SoftwareOutput/",locs1[i],".Output.data.BestCluster"))
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
    ylab("Ns")+
    ylim(0,max(Ns_tmp[[2]]$chao,Ns_tmp[[2]]$jack1)+5)+
    #ggtitle(plotnames1[i])
    ggtitle("Accumulated Ns")
  
  
}

tiff(filename = "Figures/NsAccumPlots.tiff",height = 14,width = 9,units = "in",res = 400)
multiplot(cols = 3,Nsplots[[1]],Nsplots[[4]],Nsplots[[7]],
          Nsplots[[10]],Nsplots[[13]],Nsplots[[16]],Nsplots[[19]],
          Nsplots[[2]],Nsplots[[5]],Nsplots[[8]],
          Nsplots[[11]],Nsplots[[14]],Nsplots[[17]],Nsplots[[20]],
          Nsplots[[3]],Nsplots[[6]],Nsplots[[9]],
          Nsplots[[12]],Nsplots[[15]],Nsplots[[18]],Nsplots[[21]])
dev.off()

##Scatter plots with mixture analysis results


##making per-stream plots from the already generated plots here
pdf(file = "Figures/PerStreamPlots062822.pdf")
badplot <- ggarrange(boxplots[[1]], labels = c("A"),
          ggarrange(scatters[[1]],Nsplots[[1]],labels = c("B","C"),ncol = 2),
          nrow=2)
annotate_figure(badplot,top = text_grob("Bad River",face = "bold",size = 18))
beiplot <- ggarrange(boxplots[[2]], labels = c("A"),
                     ggarrange(scatters[[2]],Nsplots[[2]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(beiplot,top = text_grob("Betsie River",face = "bold",size = 18))
betplot <- ggarrange(boxplots[[3]], labels = c("A"),
                     ggarrange(scatters[[3]],Nsplots[[3]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(betplot,top = text_grob("Betsy River",face = "bold",size = 18))
brlplot <- ggarrange(scatters[[4]],Nsplots[[4]],labels = c("A","B"),ncol = 2)
annotate_figure(brlplot,top = text_grob("Brule River",face = "bold",size = 18))
catplot <- ggarrange(boxplots[[4]], labels = c("A"),
                     ggarrange(scatters[[5]],Nsplots[[5]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(catplot,top = text_grob("Cattaraugus River",face = "bold",size = 18))
cheplot <- ggarrange(boxplots[[5]], labels = c("A"),
                     ggarrange(scatters[[6]],Nsplots[[6]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(cheplot,top = text_grob("Pigeon River",face = "bold",size = 18))
eagplot <- ggarrange(boxplots[[6]], labels = c("A"),
                     ggarrange(scatters[[7]],Nsplots[[7]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(eagplot,top = text_grob("East Au Gres River",face = "bold",size = 18))
forplot <- ggarrange(boxplots[[7]], labels = c("A"),
                     ggarrange(scatters[[8]],Nsplots[[8]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(forplot,top = text_grob("Ford River",face = "bold",size = 18))
maiplot <- ggarrange(boxplots[[8]], labels = c("A"),
                     ggarrange(scatters[[9]],Nsplots[[9]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(maiplot,top = text_grob("Manistique River",face = "bold",size = 18))
manplot <- ggarrange(boxplots[[9]], labels = c("A"),
                     ggarrange(scatters[[10]],Nsplots[[10]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(manplot,top = text_grob("Manistee River",face = "bold",size = 18))
mirplot <- ggarrange(scatters[[11]],boxplots[[10]], labels = c("A","B"),
          ggarrange(Nsplots[[11]],Nsplots[[12]],labels = c("C","D"),ncol = 2),
          nrow=3)
annotate_figure(mirplot,top = text_grob("Middle River",face = "bold",size = 18))
misplot <- ggarrange(boxplots[[11]], labels = c("A"),
                     ggarrange(scatters[[12]],Nsplots[[13]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(misplot,top = text_grob("Misery River",face = "bold",size = 18))
musplot <- ggarrange(scatters[[13]],Nsplots[[14]],labels = c("A","B"),ncol = 2)
annotate_figure(musplot,top = text_grob("Muskegon River",face = "bold",size = 18))
ocqplot <- ggarrange(scatters[[14]],boxplots[[12]], labels = c("A","B"),
                     ggarrange(Nsplots[[15]],Nsplots[[16]],labels = c("C","D"),ncol = 2),
                     nrow=3)
annotate_figure(ocqplot,top = text_grob("Ocqueoc River",face = "bold",size = 18))
steplot <- ggarrange(boxplots[[13]], labels = c("A"),
                     ggarrange(scatters[[15]],Nsplots[[17]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(steplot,top = text_grob("Sterling River",face = "bold",size = 18))
swnplot <- ggarrange(boxplots[[14]], labels = c("A"),
                     ggarrange(scatters[[16]],Nsplots[[18]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(swnplot,top = text_grob("Swan Creek",face = "bold",size = 18))
taqplot <- ggarrange(boxplots[[15]], labels = c("A"),
                     ggarrange(scatters[[17]],Nsplots[[19]],labels = c("B","C"),ncol = 2),
                     nrow=2)
annotate_figure(taqplot,top = text_grob("Tahquamenon River",face = "bold",size = 18))
twoplot <- ggarrange(scatters[[18]],boxplots[[16]], labels = c("A","B"),
                     ggarrange(Nsplots[[20]],Nsplots[[21]],labels = c("C","D"),ncol = 2),
                     nrow=3)
annotate_figure(twoplot,top = text_grob("Two-Hearted River",face = "bold",size = 18))
dev.off()
