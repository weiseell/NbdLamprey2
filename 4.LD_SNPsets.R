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

#making the NeEstimator file
genepop_create(gp_format,output_file = "SNPsets/Neestimator_2019.txt",title = "Test file for NeEstimator")





