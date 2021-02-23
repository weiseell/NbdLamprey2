#NbdLamprey - Script 4: SNP Sets
#Main Objectives:
#1. #1. Make a SNP data set for LD method and pedigree methods of Nb
#2. Make GENEPOP and Colony files for all populations
#load libraries:
library(tidyverse)

#load homebrew functions:
source("Homebrew/LD_filter.R")
source("Homebrew/vcf_genepop.R")
source("Homebrew/genepop_create.R")
#load input data:
load("Genotypes/all_PM_gts_8X.rda")
load("Summaries/SNP_summaries_targets.rda")
load("SNPsets/COLONY_SNPset.rda")
load("SNPsets/NeEst_SNPset.rda")
#Goal 1####
#selcting loci for LD method
#select one SNP per tag
geno_targets <- SNPsumm %>% filter(target != "NonTarget")
LD_select <- LD_filter(geno_targets)
#a few summary histograms to look at the SNP set
hist(LD_select$MAF)
hist(LD_select$pGT)

#generate a match file for chromosomes and SNPs 
#ensures that SNPs on the same chromosome are not compared
chrom_file <- LD_select %>% 
  select(ID) %>% 
  mutate(ID1=ID) %>% 
  separate(ID1,into = c("CHROM","POS"),sep = "-") %>% 
  rename(SNP_name=ID) %>% 
  select(CHROM,SNP_name)
write.table(chrom_file,file = "SNPsets/NeEst_chromfile.txt")
#save LD SNP set
save(LD_select,file = "SNPsets/NeEst_SNPset.rda")
#Goal 2####
#merging SNP sets with the corresponding genotype calls
gt_LD <- merge(LD_select,comb_gt8X)
#making an NeEstimator file
#formatting for input
gp_input <- gt_LD %>%
  select(ID,everything()) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "indiv",value = "gt",-ID) %>% 
  rename(SNP=ID)
#changing formatting from .vcf to genepop
gp_format <- vcf_genepop(gp_input)
gp_format <- gp_format %>% 
  mutate(indiv1=indiv) %>% 
  separate(indiv1,into = c("spp","pop","num"),sep = "_") %>% 
  select(-spp,-num) %>% 
  gather(key = "SNP",value = "gt",-indiv,-pop)

#writing Colony files for pops with multiple cohorts
mirSNP <- comb_gt8X[,which(colnames(comb_gt8X)%in%mir.1$OffspringID)]
mirSNP <- cbind(comb_gt8X[,1:2],mirSNP)
twoSNP <- comb_gt8X[,which(colnames(comb_gt8X)%in%two.1$OffspringID)]
twoSNP <- cbind(comb_gt8X[,1:2],twoSNP)
mirSNP <- merge(Col_select,mirSNP)
twoSNP <- merge(Col_select,twoSNP)
colSNPs <- merge(mirSNP,twoSNP)

col_input <- colSNPs %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(colSNPs))) %>% 
  select(-CHROM:-pGT) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
col_format <- vcf_colony(col_input)
#making the markers file for COLONY
SNPs <- colnames(col_format)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)

pops <- c("MIR","TWO")
i <- 1
for (i in 1:length(pops)) {
  tmp_col <- col_format[grep(pattern = pops[i],x = col_format$id),]
  colonydat_create(moms = NA,dads = NA,kids = tmp_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                   ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                   run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                   prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles/",pops[i],"_021820_colony2.dat"))
  
}


#making the NeEstimator files
mirSNP <- comb_gt8X[,which(colnames(comb_gt8X)%in%mir.1$OffspringID)]
mirSNP <- cbind(comb_gt8X[,1:2],mirSNP)
twoSNP <- comb_gt8X[,which(colnames(comb_gt8X)%in%two.1$OffspringID)]
twoSNP <- cbind(comb_gt8X[,1:2],twoSNP)
mirSNP <- merge(LD_select,mirSNP)
twoSNP <- merge(LD_select,twoSNP)
LDSNPs <- merge(mirSNP,twoSNP)

#formatting for input
gp_input <- LDSNPs %>%
  select(ID,everything()) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "indiv",value = "gt",-ID) %>% 
  rename(SNP=ID)
#changing formatting from .vcf to genepop
gp_format <- vcf_genepop(gp_input)
gp_format <- gp_format %>% 
  mutate(indiv1=indiv) %>% 
  separate(indiv1,into = c("spp","pop","num"),sep = "_") %>% 
  select(-spp,-num) %>% 
  gather(key = "SNP",value = "gt",-indiv,-pop)

genepop_create(gp_format,output_file = "SNPsets/Neestimator_2019_age1.txt",title = "Test file for NeEstimator")





