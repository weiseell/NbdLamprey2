#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/pedigree.plot.R")
##reading in data
load("AgingModels/MclustLW_results.rda")
load("SNPsets/COLONY_SNPsets_062321.rda")
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

##Comparing Length assignments and Colony family structure for each location
#BAD - One cohort
bad.1 <- subset(df,df$loc == "BAD" & df$Class == 1)
bad.2 <- subset(df,df$loc == "BAD" & df$Class == 2)
bad <- rbind(bad.1,bad.2)
table(bad$ClusterIndex)
table(bad.2$ClusterIndex%in%bad.1$ClusterIndex)
table(bad.1$ClusterIndex%in%bad.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(bad,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bad,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bad$EstAge <- "BAD_age1"
#BEI - maybe 2? Hard to tell
bei.1 <- subset(df,df$loc == "BEI" & df$Class == 1)
bei.2 <- subset(df,df$loc == "BEI" & df$Class == 2)
bei <- rbind(bei.1,bei.2)
table(bei$ClusterIndex)
table(bei.2$ClusterIndex%in%bei.1$ClusterIndex)
table(bei.1$ClusterIndex%in%bei.2$ClusterIndex)
#!#Overlap of clusters between length cohorts, but some separation as well?????
ggplot(bei,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bei,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bei$EstAge <- "BEI_age1"

#BET - maybe 2? Hard to tell
bet.1 <- subset(df,df$loc == "BET" & df$Class == 1)
bet.2 <- subset(df,df$loc == "BET" & df$Class == 2)
bet <- rbind(bet.1,bet.2)
table(bet$ClusterIndex)
table(bet.2$ClusterIndex%in%bet.1$ClusterIndex)
table(bet.1$ClusterIndex%in%bet.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(bet,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(bet,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
bet$EstAge <- "BET_age1"

#BRL - one cohort
brl.1 <- subset(df,df$loc == "BRL" & df$Class == 1)
brl.2 <- subset(df,df$loc == "BRL" & df$Class == 2)
brl <- rbind(brl.1,brl.2)
table(brl$ClusterIndex)
table(brl.2$ClusterIndex%in%brl.1$ClusterIndex)
table(brl.1$ClusterIndex%in%brl.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(brl,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(brl,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
brl$EstAge <- "BRL_age1"
#CAT - one cohort
cat.1 <- subset(df,df$loc == "CAT" & df$Class == 1)
cat.2 <- subset(df,df$loc == "CAT" & df$Class == 2)
cat <- rbind(cat.1,cat.2)
table(cat$ClusterIndex)
table(cat.2$ClusterIndex%in%cat.1$ClusterIndex)
table(cat.1$ClusterIndex%in%cat.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(cat,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(cat,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
cat$EstAge <- "CAT_age1"

#CHE - one cohort
che.1 <- subset(df,df$loc == "CHE" & df$Class == 1)
che.2 <- subset(df,df$loc == "CHE" & df$Class == 2)
che <- rbind(che.1,che.2)
table(che$ClusterIndex)
table(che.2$ClusterIndex%in%che.1$ClusterIndex)
table(che.1$ClusterIndex%in%che.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(che,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(che,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
che$EstAge <- "CHE_age1"
#EAG - one cohort
eag.1 <- subset(df,df$loc == "EAG" & df$Class == 1)
eag.2 <- subset(df,df$loc == "EAG" & df$Class == 2)
eag.3 <- subset(df,df$loc == "EAG" & df$Class == 3)
eag <- rbind(eag.1,eag.2,eag.3)
table(eag$ClusterIndex)
table(eag.2$ClusterIndex%in%eag.1$ClusterIndex)
table(eag.1$ClusterIndex%in%eag.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(eag,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(eag,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
eag$EstAge <- "EAG_age1"

#FORD - likely one cluster?
ford.1 <- subset(df,df$loc == "FOR" & df$Class == 1)
ford.2 <- subset(df,df$loc == "FOR" & df$Class == 2)
ford <- rbind(ford.1,ford.2)
table(ford$ClusterIndex)
table(ford.2$ClusterIndex%in%ford.1$ClusterIndex)
table(ford.1$ClusterIndex%in%ford.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(ford,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(ford,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
ford$EstAge <- "FOR_age1"
#MAI - one cluster
mai.1 <- subset(df,df$loc == "MAI" & df$Class == 1)
mai.2 <- subset(df,df$loc == "MAI" & df$Class == 2)
mai <- rbind(mai.1,mai.2)
table(mai$ClusterIndex)
table(mai.2$ClusterIndex%in%mai.1$ClusterIndex)
table(mai.1$ClusterIndex%in%mai.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(mai,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mai,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
mai$EstAge <- "MAI_age1"

#MAN - looks like it could be two clusters, lots of overlap
man.1 <- subset(df,df$loc == "MAN" & df$Class == 1)
man.2 <- subset(df,df$loc == "MAN" & df$Class == 2)
man <- rbind(man.1,man.2)
table(man$ClusterIndex)
table(man.2$ClusterIndex%in%man.1$ClusterIndex)
table(man.1$ClusterIndex%in%man.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(man,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(man,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
man$EstAge <- "MAN_age1"

#MIR - 3 cohorts, vast majority in age 1 cohort
mir.1 <- subset(df,df$loc == "MIR" & df$Class == 1)
mir.2 <- subset(df,df$loc == "MIR" & df$Class == 2)
mir.3 <- subset(df,df$loc == "MIR" & df$Class == 3)
mir.4 <- subset(df,df$loc == "MIR" & df$Class == 4)
mir <- rbind(mir.1,mir.2,mir.3,mir.4)
table(mir$ClusterIndex)
table(mir.2$ClusterIndex%in%mir.1$ClusterIndex)
table(mir.3$ClusterIndex%in%mir.1$ClusterIndex)
table(mir.4$ClusterIndex%in%mir.1$ClusterIndex)
table(mir.2$ClusterIndex%in%mir.3$ClusterIndex)
table(mir.2$ClusterIndex%in%mir.4$ClusterIndex)
table(mir.3$ClusterIndex%in%mir.4$ClusterIndex)
table(mir.3$ClusterIndex)
table(mir.4$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(mir,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mir,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mir2015_clusts <- unique(mir$ClusterIndex[which(mir$Length > 50)])
mir2015 <- mir[which(mir$ClusterIndex %in% mir2015_clusts),]
mir2015$EstAge <- "MIR_age2"
mir2016 <- mir[which(!(mir$ClusterIndex %in% mir2015_clusts)),]
mir2016$EstAge <- "MIR_age1"

#MIS - one cohort
mis.1 <- subset(df,df$loc == "MIS" & df$Class == 1)
mis.2 <- subset(df,df$loc == "MIS" & df$Class == 2)
mis.3 <- subset(df,df$loc == "MIS" & df$Class == 3)
mis <- rbind(mis.1,mis.2,mis.3)
table(mis$ClusterIndex)
table(mis.2$ClusterIndex%in%mis.1$ClusterIndex)
table(mis.3$ClusterIndex%in%mis.1$ClusterIndex)
table(mis.2$ClusterIndex%in%mis.3$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(mis,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mis,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mis$EstAge <- "MIS_age1"
#MUS
mus.1 <- subset(df,df$loc == "MUS" & df$Class == 1)
mus.2 <- subset(df,df$loc == "MUS" & df$Class == 2)
mus <- rbind(mus.1,mus.2)
table(mus$ClusterIndex)
table(mus.2$ClusterIndex%in%mus.1$ClusterIndex)
table(mus.1$ClusterIndex%in%mus.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(mus,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(mus,aes(x=ClusterIndex,y=Weight,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

mus$EstAge <- "MUS_age2"
#OCQ - definitely two clusters, maybe 3?
ocq.1 <- subset(df,df$loc == "OCQ" & df$Class == 1)
ocq.2 <- subset(df,df$loc == "OCQ" & df$Class == 2)
ocq.3 <- subset(df,df$loc == "OCQ" & df$Class == 3)
ocq.4 <- subset(df,df$loc == "OCQ" & df$Class == 4)
ocq <- rbind(ocq.1,ocq.2,ocq.3,ocq.4)
table(ocq$ClusterIndex)
table(ocq.2$ClusterIndex%in%ocq.1$ClusterIndex)
table(ocq.3$ClusterIndex%in%ocq.1$ClusterIndex)
table(ocq.4$ClusterIndex%in%ocq.1$ClusterIndex)
table(ocq.2$ClusterIndex%in%ocq.3$ClusterIndex)
table(ocq.2$ClusterIndex%in%ocq.4$ClusterIndex)
table(ocq.3$ClusterIndex%in%ocq.4$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(ocq,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(ocq,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

ocq2018 <- ocq[which(ocq$Length >= 38 | ocq$ClusterIndex == 36),]
ocq2018$EstAge <- "OCQ_age1"
ocq2019 <- ocq[which(ocq$Length < 38 & ocq$ClusterIndex != 36),]
ocq2019$EstAge <- "OCQ_age0"

#STE - one cohort from LBAM, no follow-up needed
ste <- subset(df,df$loc == "STE")
ste$EstAge <- "STE_age1"
#SWN - one cohort
swn.1 <- subset(df,df$loc == "SWN" & df$Class == 1)
swn.2 <- subset(df,df$loc == "SWN" & df$Class == 2)
swn <- rbind(swn.1,swn.2)
table(swn$ClusterIndex)
table(swn.2$ClusterIndex%in%swn.1$ClusterIndex)
table(swn.1$ClusterIndex%in%swn.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(swn,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(swn,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
swn$EstAge <- "SWN_age1"

#TAQ - one cohort
taq.1 <- subset(df,df$loc == "TAQ" & df$Class == 1)
taq.2 <- subset(df,df$loc == "TAQ" & df$Class == 2)
taq <- rbind(taq.1,taq.2)
table(taq$ClusterIndex)
table(taq.2$ClusterIndex%in%taq.1$ClusterIndex)
table(taq.1$ClusterIndex%in%taq.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(taq,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(taq,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")
taq$EstAge <- "TAQ_age1"
#TWO - two cohorts
two.1 <- subset(df,df$loc == "TWO" & df$Class == 1)
two.2 <- subset(df,df$loc == "TWO" & df$Class == 2)
two <- rbind(two.1,two.2)
table(two$ClusterIndex)
table(two.2$ClusterIndex%in%two.1$ClusterIndex)
table(two.1$ClusterIndex%in%two.2$ClusterIndex)
#!#Overlap of clusters between length cohorts
ggplot(two,aes(x=ClusterIndex,y=Length,color=Class))+
  geom_point()+theme_bw()

ggplot(two,aes(x=ClusterIndex,y=Length,fill=as.numeric(Probability)))+
  geom_boxplot()+theme_bw()+
  scale_fill_gradient(low = "red",high = "white")

two2017 <- two[which(two$Length >= 50 & two$ClusterIndex != 2),]
two2017$EstAge <- "TWO_age2"
two2018 <- two[which(two$Length < 50 | two$ClusterIndex == 2),]
two2018$EstAge <- "TWO_age1"

#merge together all cohorts
alllocs_cohorts <- rbind(bad,bei,bet,brl,cat,
      che,eag,ford,mai,man,
      mir2015,mir2016,mis,mus,
      ocq2018,ocq2019,ste,swn,taq,
      two2017,two2018)
save(alllocs_cohorts,file = "Summaries/Alllocs_cohorts.rda")
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








