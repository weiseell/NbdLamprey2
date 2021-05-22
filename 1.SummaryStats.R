#NbdLamprey - Script 1
##Main Objective: Calculate the following summary statistics:
#1. Calculate percent coverage and heterozygosity for all SNPs
#2. Match SNPs to the Rapture Panel and calculate stats for Rapture success

#load libraries:
library(tidyverse)

#load homebrew functions:
source("Homebrew/match_tags.R")

#load input data:
load("Genotypes/all_PM_gts_8X.rda")
load("Rapture_all_genomes.rda")
#Goal 1####
#calculating %coverage, heterozygosity for all loci
#merging with MAF data
#should have a large table with all genotypes and three summary criteria
#allelecounts are per population
AC_summ <- AC_comb %>% 
  group_by(CHROM,POS) %>% 
  summarize(sum_total=sum(N_CHR),
          sum_MinAC=sum(MinAC),
          MAF=sum_MinAC/sum_total) %>% 
  mutate(ID = paste0(CHROM,"-",POS))

#Checking MAF file: if there are MAF > 0.5, fix so they're < 0.5
#check that all MAF values are below 0.5 and saving as an individual file
max(AC_summ$MAF)

MAF <- data.frame(ID = AC_summ$ID, MAF = AC_summ$MAF, stringsAsFactors = F)

#change NAs to empty data
comb_gt8X[is.na(comb_gt8X)] <- "./."
comb_gt8X$ID <- paste0(comb_gt8X$CHROM,"-",comb_gt8X$POS)
#remove empty rows
comb_gt8X <- comb_gt8X[which(!(rowSums(comb_gt8X=="./.")==2075)),]

#calculate heterozygosity
n_indiv <- ncol(comb_gt8X)-1
stats <- data.frame(ID = comb_gt8X$ID,stringsAsFactors = F)
het_counts <- rowSums(comb_gt8X == "0/1")
stats$het <- het_counts/n_indiv
#calculate percent coverage
gt_missing <- rowSums(comb_gt8X == "./.")
stats$pGT <- (1-(gt_missing/n_indiv))
#merge with MAF
stats <- merge(stats,MAF)


#Goal 2####
#matching genotypes to rapture panel

#separating tag IDs to make location ranges on reference genome scaffolds for all tags
rapture1 <- link_file1 %>% 
  rename(min=start2,max=end2) %>% 
  mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
  select(ID, CHROM, min, max)
SNPs <- comb_gt8X %>% 
  select(CHROM,POS)

target <- match_tags(SNPs = SNPs,tags = rapture1)

#merging target results with summary statistics
rapture_class <- target[[1]]
rapture_class <- rapture_class %>% 
  mutate(ID = paste(CHROM,POS,sep = "-")) %>% 
  select(ID,everything())
SNPsumm <- merge(rapture_class,stats)

SNPsumm1 <- subset(SNPsumm,SNPsumm$target != "NonTarget")
mean(SNPsumm1$pGT)
#save two main outputs from this script
#summary stats, target classifications, and on-target SNPs
save(SNPsumm,file = "Summaries/SNP_summaries_targets.rda")
 






