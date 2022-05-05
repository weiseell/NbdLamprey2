
library(tidyverse)
source("Homebrew/match_tags.R")
source("Homebrew/ld_filter.R")
source("Homebrew/genepop_create.R")
source("Homebrew/vcf_genepop.R")
source("Homebrew/match_tags.R")
load("Summaries/Alllocs_cohorts.rda")
load("Rapture_all_genomes.rda")

rapture1 <- link_file1 %>% 
  rename(min=start2,max=end2) %>% 
  mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
  select(ID, CHROM, min, max)

#generating NeEstimator files for all populations
pops <- unique(alllocs_cohorts$loc)
NeEst_out <- list()
for (i in 1:length(pops)) {
  print(i)
  pop <- pops[i]
  gts <- read.table(paste0("Input/",pop,"_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
  gts$ID <- paste0(gts$CHROM,"-",gts$POS)
  
  #calculate het and pGT for SNP selection
  #calculate heterozygosity
  n_indiv <- ncol(gts)-3
  stats <- data.frame(ID = gts$ID,CHROM = gts$CHROM,POS = gts$POS,stringsAsFactors = F)
  het_counts <- rowSums(gts == "0/1")
  stats$het <- het_counts/n_indiv
  #calculate percent coverage
  gt_missing <- rowSums(gts == "./.")
  stats$pGT <- (1-(gt_missing/n_indiv))
  
  #select SNPs
  SNPstmp <- stats %>% select(CHROM,POS)
  targets <- match_tags(SNPs = SNPstmp,tags = rapture1)
  SNPsummtmp <- merge(stats,targets[[1]])
  geno_targets <- SNPsummtmp %>% filter(target != "NonTarget")
  LD_select <- LD_filter(geno_targets)
  #save LD SNP set
  NeEst_out[[pop]] <- LD_select
  
  #generate a match file for chromosomes and SNPs 
  #ensures that SNPs on the same chromosome are not compared
  #chrom_file <- LD_select %>% 
  #  select(ID) %>% 
  #  mutate(ID1=ID) %>% 
  #  separate(ID1,into = c("CHROM","POS"),sep = "-") %>% 
  #  rename(SNP_name=ID) %>% 
  #  select(CHROM,SNP_name)
  #write.table(chrom_file,file = paste0("SNPsets/",pop,"NeEst_chromfile.txt"))

  #merging SNP sets with the corresponding genotype calls
  #gt_LD <- merge(LD_select,gts)
  #making an NeEstimator file
  #formatting for input
  #gp_input <- gt_LD %>%
  #  select(ID,everything()) %>% 
  #  select(-CHROM:-target) %>% 
  #  gather(key = "indiv",value = "gt",-ID) %>% 
  #  rename(SNP=ID)
  #changing formatting from .vcf to genepop
  #gp_format <- vcf_genepop(gp_input)
  #tmp <- alllocs_cohorts %>% 
  #  select(OffspringID,EstAge) %>% 
  #  rename(indiv=OffspringID,pop=EstAge)
  #gp_format <- merge(gp_format,tmp)
  
  #gp_format <- gp_format %>% 
  #  gather(key = "SNP",value = "gt",-indiv,-pop)
  
  #genepop_create(gp_format,output_file = paste0("SNPsets/",pop,"Neestimator_2019.txt"),title = paste(pop,"file for NeEstimator"))
}

save(LD_select,file = paste0("SNPsets/",pop,"_NeEst_SNPset.rda"))




