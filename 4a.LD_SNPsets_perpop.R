


#generating NeEstimator files for all single-cohort populations
pops <- c("BAD","BEI","BET","BRL","CAT","CHE","EAG","FOR","MAI","MIS","MUS","OCQ","STE","SWN","TAQ","TWO")
pops <- c("MAN","MIR")
for (i in 1:length(pops)) {
  print(i)
  pop <- pops[i]
  gts <- read.table(paste0("Input/",pop,"_GTs.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
  gts$ID <- paste0(gts$CHROM,"-",gts$POS)
  MAF <- read.table(paste0("Input/",pop,"_allele.frq"),header = T,sep = "\t",stringsAsFactors = F)
  MAF$ID <- paste(MAF$CHROM,MAF$POS,sep = "-")
  #calculate het and pGT for SNP selection
  #calculate heterozygosity
  n_indiv <- ncol(gts)-3
  stats <- data.frame(ID = gts$ID,stringsAsFactors = F)
  het_counts <- rowSums(gts == "0/1")
  stats$het <- het_counts/n_indiv
  #calculate percent coverage
  gt_missing <- rowSums(gts == "./.")
  stats$pGT <- (1-(gt_missing/n_indiv))
  #merge with MAF
  stats <- merge(stats,MAF)
  SNPsumm <- match_tags(stats,rapture1)
  
  #select SNPs
  geno_targets <- SNPsumm[[1]] %>% filter(target != "NonTarget")
  LD_select <- LD_filter(geno_targets)
  #generate a match file for chromosomes and SNPs 
  #ensures that SNPs on the same chromosome are not compared
  chrom_file <- LD_select %>% 
    select(ID) %>% 
    mutate(ID1=ID) %>% 
    separate(ID1,into = c("CHROM","POS"),sep = "-") %>% 
    rename(SNP_name=ID) %>% 
    select(CHROM,SNP_name)
  write.table(chrom_file,file = paste0("SNPsets/",pop,"NeEst_chromfile.txt"))
  #save LD SNP set
  save(LD_select,file = paste0("SNPsets/",pop,"NeEst_SNPset.rda"))
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
  #gp_format <- gp_format %>% 
  #  mutate(indiv1=indiv) %>% 
  #  separate(indiv1,into = c("spp","pop","num"),sep = "_") %>% 
  #  select(-spp,-num) %>% 
  #  gather(key = "SNP",value = "gt",-indiv,-pop)
  
  #genepop_create(gp_format,output_file = paste0("SNPsets/",pop,"Neestimator_2019_age1.txt"),title = paste(pop,"file for NeEstimator"))
}






