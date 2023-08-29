#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/pedigree.plot.R")
##reading in data
load("AgingModels/MclustLW_results.rda")

#load in data
classes_all <- classes_all %>% 
  rename(OffspringID=Indiv)
#note - BestConfig files were reformatted to be tab delimited and 
#special characters in the file were removed prior to load
#identify locations with multiple inferred cohorts
classes_all %>% 
  group_by(loc) %>% 
  summarise(nclust=length(unique(Class)),ss=n(),max_len = max(Length))
locs <- unique(classes_all$loc)
best_config <- data.frame(matrix(ncol=5,nrow = 0))
#read in pedigree data for locations with multiple inferred cohorts
for (i in 1:length(locs)) {
  print(i)
  df <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  df <- strsplit(df,"\\s+")
  df1 <- matrix(unlist(df),byrow = T)
  df1 <- df1[df1 != ""]
  df1 <- matrix(df1,ncol = 5,byrow = T)
  df1 <- as.data.frame(df1)
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]
  
  best_config <- rbind(best_config,df1)
}
#adding ages to model results
df <- merge(classes_all,best_config)
df_Classsum <- df %>% 
  group_by(loc,Class) %>% 
  summarise(ss=n(),max_len = max(Length))

write.table(df,file = "AgingModels/lw_cohorts.txt",append = F,quote = F,sep = "\t",row.names = F)

#loop to generate boxplots for all locations
boxplots_fullsib <- list()
i <- 1
plotn <- 0
df$fullsib <- paste(df$MotherID,df$FatherID,sep = "_")
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
for (i in 1:length(locs)) {
  loc <- locs[i]
  tmp <- subset(df,df$loc == locs[i])
  #check to make sure there's more than one family
  #tmp <- tmp[tmp$fullsib %in% names(which(table(tmp$fullsib)>=3)),]
  #plot pedigree
  if(length(unique(tmp$fullsib)) > 1){
    plotn <- plotn + 1
    boxplots_fullsib[[plotn]] <- ggplot(tmp,aes(x=fullsib,group=fullsib,y=Length))+
      geom_boxplot(alpha=0.3,show.legend = F)+
      geom_jitter(aes(color = as.character(Class)),position=position_jitter(0.1))+
      scale_color_manual(values = wes_palette("IsleofDogs1"))+
      theme_bw()+
      #scale_fill_gradient(low="red", high="white",name = "Cluster \n Likelihood")+
      xlab("Full-Sibling Group")+
      ylab("Length (mm)")+
      ggtitle(plotnames[i])+
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none")
  }
}
pdf(file = "Figures/BoxplotsFullSib062022.pdf")
boxplots_fullsib[[1]]
boxplots_fullsib[[2]]
boxplots_fullsib[[3]]
boxplots_fullsib[[4]]
boxplots_fullsib[[5]]
boxplots_fullsib[[6]]
boxplots_fullsib[[7]]
boxplots_fullsib[[8]]
boxplots_fullsib[[9]]
boxplots_fullsib[[10]]
boxplots_fullsib[[11]]
boxplots_fullsib[[12]]
boxplots_fullsib[[13]]
boxplots_fullsib[[14]]
boxplots_fullsib[[15]]
boxplots_fullsib[[16]]
boxplots_fullsib[[17]]
boxplots_fullsib[[18]]
dev.off()

##Comparing Length assignments and Colony family structure for each location
#BAD - One cohort
bad.1 <- subset(df,df$loc == "BAD" & df$Class == 1)
bad.2 <- subset(df,df$loc == "BAD" & df$Class == 2)
bad <- rbind(bad.1,bad.2)
table(bad$fullsib)
table(bad.2$fullsib%in%bad.1$fullsib)
table(bad.1$fullsib%in%bad.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(bad,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bad,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bad$EstAge <- "BAD_age1"
#BEI - one, length range is very narrow
bei.1 <- subset(df,df$loc == "BEI" & df$Class == 1)
bei.2 <- subset(df,df$loc == "BEI" & df$Class == 2)
bei <- rbind(bei.1,bei.2)
table(bei$fullsib)
table(bei.2$fullsib%in%bei.1$fullsib)
table(bei.1$fullsib%in%bei.2$fullsib)
#!#Overlap of clusters between length cohorts, but some separation as well?????
ggplot(bei,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bei,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bei$EstAge <- "BEI_age1"

#BET - one, length range is very narrow
bet.1 <- subset(df,df$loc == "BET" & df$Class == 1)
bet.2 <- subset(df,df$loc == "BET" & df$Class == 2)
bet.3 <- subset(df,df$loc == "BET" & df$Class == 3)
bet <- rbind(bet.1,bet.2,bet.3)
table(bet$fullsib)
table(bet.2$fullsib%in%bet.1$fullsib)
table(bet.1$fullsib%in%bet.3$fullsib)
table(bet.2$fullsib%in%bet.3$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(bet,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bet,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bet$EstAge <- "BET_age1"

#BRL - one cohort
brl.1 <- subset(df,df$loc == "BRL" & df$Class == 1)
brl.2 <- subset(df,df$loc == "BRL" & df$Class == 2)
brl <- rbind(brl.1,brl.2)
table(brl$fullsib)
table(brl.2$fullsib%in%brl.1$fullsib)
table(brl.1$fullsib%in%brl.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(brl,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(brl,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
brl$EstAge <- "BRL_age1"
#CAT - one cohort
cat.1 <- subset(df,df$loc == "CAT" & df$Class == 1)
cat.2 <- subset(df,df$loc == "CAT" & df$Class == 2)
cat <- rbind(cat.1,cat.2)
table(cat$fullsib)
table(cat.2$fullsib%in%cat.1$fullsib)
table(cat.1$fullsib%in%cat.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(cat,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(cat,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
cat$EstAge <- "CAT_age1"

#CHE - one cohort
che.1 <- subset(df,df$loc == "CHE" & df$Class == 1)
che.2 <- subset(df,df$loc == "CHE" & df$Class == 2)
che <- rbind(che.1,che.2)
table(che$fullsib)
table(che.2$fullsib%in%che.1$fullsib)
table(che.1$fullsib%in%che.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(che,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(che,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
che$EstAge <- "CHE_age1"
#EAG - one cohort
eag.1 <- subset(df,df$loc == "EAG" & df$Class == 1)
eag.2 <- subset(df,df$loc == "EAG" & df$Class == 2)
eag.3 <- subset(df,df$loc == "EAG" & df$Class == 3)
eag <- rbind(eag.1,eag.2,eag.3)
table(eag$fullsib)
table(eag.2$fullsib%in%eag.1$fullsib)
table(eag.1$fullsib%in%eag.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(eag,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(eag,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
eag$EstAge <- "EAG_age1"

#FORD - likely one cluster?
ford.1 <- subset(df,df$loc == "FOR" & df$Class == 1)
ford.2 <- subset(df,df$loc == "FOR" & df$Class == 2)
ford <- rbind(ford.1,ford.2)
table(ford$fullsib)
table(ford.2$fullsib%in%ford.1$fullsib)
table(ford.1$fullsib%in%ford.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(ford,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(ford,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
ford$EstAge <- "FOR_age1"
#MAI - one cluster
mai.1 <- subset(df,df$loc == "MAI" & df$Class == 1)
mai.2 <- subset(df,df$loc == "MAI" & df$Class == 2)
mai <- rbind(mai.1,mai.2)
table(mai$fullsib)
table(mai.2$fullsib%in%mai.1$fullsib)
table(mai.1$fullsib%in%mai.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(mai,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mai,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
mai$EstAge <- "MAI_age1"

#MAN -  lots of overlap
#one cohort, length range is very narrow
man.1 <- subset(df,df$loc == "MAN" & df$Class == 1)
man.2 <- subset(df,df$loc == "MAN" & df$Class == 2)
man <- rbind(man.1,man.2)
table(man$fullsib)
table(man.2$fullsib%in%man.1$fullsib)
table(man.1$fullsib%in%man.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(man,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(man,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
man$EstAge <- "MAN_age1"

#MIR - 3 cohorts, vast majority in age 1 cohort
mir.1 <- subset(df,df$loc == "MIR" & df$Class == 1)
mir.2 <- subset(df,df$loc == "MIR" & df$Class == 2)
mir.3 <- subset(df,df$loc == "MIR" & df$Class == 3)
mir.4 <- subset(df,df$loc == "MIR" & df$Class == 4)
mir <- rbind(mir.1,mir.2,mir.3,mir.4)
table(mir$fullsib)
table(mir.2$fullsib%in%mir.1$fullsib)
table(mir.3$fullsib%in%mir.1$fullsib)
table(mir.4$fullsib%in%mir.1$fullsib)
table(mir.2$fullsib%in%mir.3$fullsib)
table(mir.2$fullsib%in%mir.4$fullsib)
table(mir.3$fullsib%in%mir.4$fullsib)
table(mir.3$fullsib)
table(mir.4$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(mir,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mir,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mir2015_clusts <- unique(mir$fullsib[which(mir$Length > 50)])
mir2015 <- mir[which(mir$fullsib %in% mir2015_clusts),]
mir2015$EstAge <- "MIR_age2"
mir2016 <- mir[which(!(mir$fullsib %in% mir2015_clusts)),]
mir2016$EstAge <- "MIR_age1"

#MIS - one cohort
mis.1 <- subset(df,df$loc == "MIS" & df$Class == 1)
mis.2 <- subset(df,df$loc == "MIS" & df$Class == 2)
mis.3 <- subset(df,df$loc == "MIS" & df$Class == 3)
mis <- rbind(mis.1,mis.2,mis.3)
table(mis$fullsib)
table(mis.2$fullsib%in%mis.1$fullsib)
table(mis.3$fullsib%in%mis.1$fullsib)
table(mis.2$fullsib%in%mis.3$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(mis,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mis,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mis$EstAge <- "MIS_age1"
#MUS
mus.1 <- subset(df,df$loc == "MUS" & df$Class == 1)
mus.2 <- subset(df,df$loc == "MUS" & df$Class == 2)
mus <- rbind(mus.1,mus.2)
table(mus$fullsib)
table(mus.2$fullsib%in%mus.1$fullsib)
table(mus.1$fullsib%in%mus.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(mus,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mus,aes(x=fullsib,y=Weight,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mus$EstAge <- "MUS_age2"
#OCQ - definitely two clusters, maybe 3?
ocq.1 <- subset(df,df$loc == "OCQ" & df$Class == 1)
ocq.2 <- subset(df,df$loc == "OCQ" & df$Class == 2)
ocq.3 <- subset(df,df$loc == "OCQ" & df$Class == 3)
ocq.4 <- subset(df,df$loc == "OCQ" & df$Class == 4)
ocq <- rbind(ocq.1,ocq.2,ocq.3,ocq.4)
table(ocq$fullsib)
table(ocq.2$fullsib%in%ocq.1$fullsib)
table(ocq.3$fullsib%in%ocq.1$fullsib)
table(ocq.4$fullsib%in%ocq.1$fullsib)
table(ocq.2$fullsib%in%ocq.3$fullsib)
table(ocq.2$fullsib%in%ocq.4$fullsib)
table(ocq.3$fullsib%in%ocq.4$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(ocq,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(ocq,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

ocq2018 <- ocq[which(ocq$Length >= 38 | ocq$fullsib == 36),]
ocq2018$EstAge <- "OCQ_age1"
ocq2019 <- ocq[which(ocq$Length < 38 & ocq$fullsib != 36),]
ocq2019$EstAge <- "OCQ_age0"

#STE - one cohort from LBAM, no follow-up needed
ste <- subset(df,df$loc == "STE")
ste$EstAge <- "STE_age1"
#SWN - one cohort
swn.1 <- subset(df,df$loc == "SWN" & df$Class == 1)
swn.2 <- subset(df,df$loc == "SWN" & df$Class == 2)
swn <- rbind(swn.1,swn.2)
table(swn$fullsib)
table(swn.2$fullsib%in%swn.1$fullsib)
table(swn.1$fullsib%in%swn.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(swn,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(swn,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
swn$EstAge <- "SWN_age1"

#TAQ - one cohort
taq.1 <- subset(df,df$loc == "TAQ" & df$Class == 1)
taq.2 <- subset(df,df$loc == "TAQ" & df$Class == 2)
taq <- rbind(taq.1,taq.2)
table(taq$fullsib)
table(taq.2$fullsib%in%taq.1$fullsib)
table(taq.1$fullsib%in%taq.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(taq,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(taq,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
taq$EstAge <- "TAQ_age1"
#TWO - two cohorts
two.1 <- subset(df,df$loc == "TWO" & df$Class == 1)
two.2 <- subset(df,df$loc == "TWO" & df$Class == 2)
two <- rbind(two.1,two.2)
table(two$fullsib)
table(two.2$fullsib%in%two.1$fullsib)
table(two.1$fullsib%in%two.2$fullsib)
#!#Overlap of clusters between length cohorts
ggplot(two,aes(x=fullsib,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(two,aes(x=fullsib,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

two2017 <- two[which(two$Length >= 50 & two$fullsib != 2),]
two2017$EstAge <- "TWO_age2"
two2018 <- two[which(two$Length < 50 | two$fullsib == 2),]
two2018$EstAge <- "TWO_age1"

#merge together all cohorts
alllocs_cohorts <- rbind(bad,bei,bet,brl,cat,
      che,eag,ford,mai,man,
      mir2015,mir2016,mis,mus,
      ocq2018,ocq2019,ste,swn,taq,
      two2017,two2018)
save(alllocs_cohorts,file = "Summaries/Alllocs_cohorts.rda")  

##Make figures combining all locations
alllocs_cohorts$ClusterIndex <- factor(alllocs_cohorts$ClusterIndex)
##make cohort colony files for all locations with multiple cohorts
#mir
#load genotypes
mir_snps <- read.table("Input/MIR_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
mir_colsnps <- merge(Colony_out[["MIR"]],mir_snps)
#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
mir_colinput <- mir_colsnps %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(mir_colsnps))) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
mir_colformat <- vcf_colony(mir_colinput)
#making the markers file for COLONY
SNPs <- colnames(mir_colformat)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
mir2015_col <- mir_colformat[which(mir_colformat$id %in% mir2015$OffspringID),]
mir2016_col <- mir_colformat[which(mir_colformat$id %in% mir2016$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = mir2015_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/MIR2015_colony2.dat"))
colonydat_create(moms = NA,dads = NA,kids = mir2016_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/MIR2016_colony2.dat"))

#ocq
#load genotypes
ocq_snps <- read.table("Input/OCQ_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
ocq_colsnps <- merge(Colony_out[["OCQ"]],ocq_snps)
#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
ocq_colinput <- ocq_colsnps %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(ocq_colsnps))) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
ocq_colformat <- vcf_colony(ocq_colinput)
#making the markers file for COLONY
SNPs <- colnames(ocq_colformat)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
ocq2018_col <- ocq_colformat[which(ocq_colformat$id %in% ocq2018$OffspringID),]
ocq2019_col <- ocq_colformat[which(ocq_colformat$id %in% ocq2019$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = ocq2018_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/OCQ2018_colony2.dat"))
colonydat_create(moms = NA,dads = NA,kids = ocq2019_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/OCQ2019_colony2.dat"))

#two
two_snps <- read.table("Input/TWO_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
two_colsnps <- merge(Colony_out[["TWO"]],two_snps)
#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
two_colinput <- two_colsnps %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(two_colsnps))) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
two_colformat <- vcf_colony(two_colinput)
#making the markers file for COLONY
SNPs <- colnames(two_colformat)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
two2018_col <- two_colformat[which(two_colformat$id %in% two2018$OffspringID),]
two2017_col <- two_colformat[which(two_colformat$id %in% two2017$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = two2018_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/two2018_colony2.dat"))
colonydat_create(moms = NA,dads = NA,kids = two2017_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/two2017_colony2.dat"))








