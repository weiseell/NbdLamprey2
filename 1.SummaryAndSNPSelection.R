##Script 1: Summary Statistics and SNP set selection
##Goals:
##1. Generate heterozygostiy and percent genotyped for each SNP
##2. Generate COLONY SNP set and input file for each stream location

#load libraries

#load custom functions
source("Homebrew/COLONY_filter.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/marker_create.R")
source("Homebrew/colonydat_create.R")

#list of populations
pops <- c("BAD","BEI","BET","BRL","CAT","CHE","EAG","FOR","MAI","MAN","MIR","MIS","MUS","OCQ","STE","SWN","TAQ","TWO")

#blank list to store which SNPs are in each COLONY set
Colony_out <- list()

#per population loop
for (i in 1:length(pops)) {
  print(i)
  #read in input files
  pop <- pops[i]
  gts <- read.table(paste0("Input/",pop,"_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
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
  
  prefilter <- stats %>%
    rename(MAF=MinAF) %>% 
    filter(MAF > 0.05 & pGT > 0.8)
  
  #run colony filter function
  Colony_out[[pop]] <- COLONY_filter(prefilter,window = 1000000,pGT_min = 0.8,MAF_min = 0.05)
  
  #Goal 2####
  #merging SNP sets with the corresponding genotype calls
  gt_Col <- merge(Colony_out[[pop]],gts)
  
  #making a COLONY input file
  #changing gt_Col into a three col data frame
  #only ID, Individual, and gt
  col_input <- gt_Col %>% 
    #making generic loci names for COLONY here
    mutate(locus = paste0("L",1:nrow(gt_Col))) %>% 
    select(-CHROM:-pGT) %>% 
    gather(key = "id",value = "gt",-locus)
  #converting gts from .vcf format to COLONY format
  col_format <- vcf_colony(col_input)
  #making the markers file for COLONY
  SNPs <- colnames(col_format)
  SNPs <- SNPs[-1]
  markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
  
  #make colony file with all locations
  tmp_col <- col_format[grep(pattern = pops[i],x = col_format$id),]
  colonydat_create(moms = NA,dads = NA,kids = tmp_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                   ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                   run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                   prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_011222/",pops[i],"_colony2.dat"))
}

#Save COLONY SNP set file
#!# Note: check for length outliers and remake COLONY files without outliers to run for Nb
save(Colony_out,file = "SNPsets/COLONY_SNPsets_011222.rda")




