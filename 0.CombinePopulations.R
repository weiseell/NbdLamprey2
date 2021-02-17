##load in all GT and depth files and merge across populations
BAD_gt <- read.table("Genotypes/BAD_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BEI_gt <- read.table("Genotypes/BEI_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BET_gt <- read.table("Genotypes/BET_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BRL_gt <- read.table("Genotypes/BRL_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
CAT_gt <- read.table("Genotypes/CAT_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
CHE_gt <- read.table("Genotypes/CHE_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
EAG_gt <- read.table("Genotypes/EAG_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
FOR_gt <- read.table("Genotypes/FOR_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MAN_gt <- read.table("Genotypes/MAN_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MAI_gt <- read.table("Genotypes/MAI_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MIR_gt <- read.table("Genotypes/MIR_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MIS_gt <- read.table("Genotypes/MIS_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MSC_gt <- read.table("Genotypes/MSC_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MUS_gt <- read.table("Genotypes/MUS_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
OCQ_gt <- read.table("Genotypes/OCQ_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
STE_gt <- read.table("Genotypes/STE_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
SWN_gt <- read.table("Genotypes/SWN_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
TAQ_gt <- read.table("Genotypes/TAQ_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
TWO_gt <- read.table("Genotypes/TWO_GTs.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)

BAD_depth <- read.table("Genotypes/BAD_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
BEI_depth <- read.table("Genotypes/BEI_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
BET_depth <- read.table("Genotypes/BET_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
BRL_depth <- read.table("Genotypes/BRL_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
CAT_depth <- read.table("Genotypes/CAT_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
CHE_depth <- read.table("Genotypes/CHE_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
EAG_depth <- read.table("Genotypes/EAG_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
FOR_depth <- read.table("Genotypes/FOR_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MAN_depth <- read.table("Genotypes/MAN_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MAI_depth <- read.table("Genotypes/MAI_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MIR_depth <- read.table("Genotypes/MIR_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MIS_depth <- read.table("Genotypes/MIS_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MSC_depth <- read.table("Genotypes/MSC_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
MUS_depth <- read.table("Genotypes/MUS_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
OCQ_depth <- read.table("Genotypes/OCQ_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
STE_depth <- read.table("Genotypes/STE_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
SWN_depth <- read.table("Genotypes/SWN_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
TAQ_depth <- read.table("Genotypes/TAQ_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)
TWO_depth <- read.table("Genotypes/TWO_depth.gdepth",header = T,sep = "\t",stringsAsFactors = F)

BAD_MAF <- read.table("Genotypes/BAD_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
BAD_MAF$pop <- "BAD"
BEI_MAF <- read.table("Genotypes/BEI_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
BET_MAF <- read.table("Genotypes/BET_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
BRL_MAF <- read.table("Genotypes/BRL_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
CAT_MAF <- read.table("Genotypes/CAT_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
CHE_MAF <- read.table("Genotypes/CHE_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
EAG_MAF <- read.table("Genotypes/EAG_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
FOR_MAF <- read.table("Genotypes/FOR_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MAN_MAF <- read.table("Genotypes/MAN_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MAI_MAF <- read.table("Genotypes/MAI_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MIR_MAF <- read.table("Genotypes/MIR_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MIS_MAF <- read.table("Genotypes/MIS_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MSC_MAF <- read.table("Genotypes/MSC_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
MUS_MAF <- read.table("Genotypes/MUS_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
OCQ_MAF <- read.table("Genotypes/OCQ_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
STE_MAF <- read.table("Genotypes/STE_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
SWN_MAF <- read.table("Genotypes/SWN_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
TAQ_MAF <- read.table("Genotypes/TAQ_allele.frq",header = T,sep = "\t",stringsAsFactors = F)
TWO_MAF <- read.table("Genotypes/TWO_allele.frq",header = T,sep = "\t",stringsAsFactors = F)

BAD_AC <- read.table("Genotypes/BAD_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
BEI_AC <- read.table("Genotypes/BEI_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
BET_AC <- read.table("Genotypes/BET_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
BRL_AC <- read.table("Genotypes/BRL_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
CAT_AC <- read.table("Genotypes/CAT_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
CHE_AC <- read.table("Genotypes/CHE_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
EAG_AC <- read.table("Genotypes/EAG_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
FOR_AC <- read.table("Genotypes/FOR_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MAN_AC <- read.table("Genotypes/MAN_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MAI_AC <- read.table("Genotypes/MAI_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MIR_AC <- read.table("Genotypes/MIR_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MIS_AC <- read.table("Genotypes/MIS_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MSC_AC <- read.table("Genotypes/MSC_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
MUS_AC <- read.table("Genotypes/MUS_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
OCQ_AC <- read.table("Genotypes/OCQ_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
STE_AC <- read.table("Genotypes/STE_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
SWN_AC <- read.table("Genotypes/SWN_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
TAQ_AC <- read.table("Genotypes/TAQ_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)
TWO_AC <- read.table("Genotypes/TWO_allele.frq.count",header = T,sep = "\t",stringsAsFactors = F)

BAD_MAF$pop <- "BAD"
BAD_AC$pop <- "BAD"
BEI_MAF$pop <- "BEI"
BEI_AC$pop <- "BEI"
BET_MAF$pop <- "BET"
BET_AC$pop <- "BET"
BRL_MAF$pop <- "BRL"
BRL_AC$pop <- "BRL"
CAT_MAF$pop <- "CAT"
CAT_AC$pop <- "CAT"
CHE_MAF$pop <- "CHE"
CHE_AC$pop <- "CHE"
EAG_MAF$pop <- "EAG"
EAG_AC$pop <- "EAG"
FOR_MAF$pop <- "FOR"
FOR_AC$pop <- "FOR"
MAN_MAF$pop <- "MAN"
MAN_AC$pop <- "MAN"
MAI_MAF$pop <- "MAI"
MAI_AC$pop <- "MAI"
MIR_MAF$pop <- "MIR"
MIR_AC$pop <- "MIR"
MIS_MAF$pop <- "MIS"
MIS_AC$pop <- "MIS"
MSC_MAF$pop <- "MSC"
MSC_AC$pop <- "MSC"
MUS_MAF$pop <- "MUS"
MUS_AC$pop <- "MUS"
OCQ_MAF$pop <- "OCQ"
OCQ_AC$pop <- "OCQ"
STE_MAF$pop <- "STE"
STE_AC$pop <- "STE"
SWN_MAF$pop <- "SWN"
SWN_AC$pop <- "SWN"
TAQ_MAF$pop <- "TAQ"
TAQ_AC$pop <- "TAQ"
TWO_MAF$pop <- "TWO"
TWO_AC$pop <- "TWO"

BAD_idepth <- read.table("Genotypes/BAD_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
BEI_idepth <- read.table("Genotypes/BEI_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
BET_idepth <- read.table("Genotypes/BET_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
BRL_idepth <- read.table("Genotypes/BRL_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
CAT_idepth <- read.table("Genotypes/CAT_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
CHE_idepth <- read.table("Genotypes/CHE_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
EAG_idepth <- read.table("Genotypes/EAG_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
FOR_idepth <- read.table("Genotypes/FOR_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MAN_idepth <- read.table("Genotypes/MAN_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MAI_idepth <- read.table("Genotypes/MAI_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MIR_idepth <- read.table("Genotypes/MIR_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MIS_idepth <- read.table("Genotypes/MIS_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MSC_idepth <- read.table("Genotypes/MSC_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
MUS_idepth <- read.table("Genotypes/MUS_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
OCQ_idepth <- read.table("Genotypes/OCQ_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
STE_idepth <- read.table("Genotypes/STE_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
SWN_idepth <- read.table("Genotypes/SWN_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
TAQ_idepth <- read.table("Genotypes/TAQ_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)
TWO_idepth <- read.table("Genotypes/TWO_indiv_depth.idepth",header = T,sep = "\t",stringsAsFactors = F)

BAD_idepth$pop <- "BAD"
BEI_idepth$pop <- "BEI"
BET_idepth$pop <- "BET"
BRL_idepth$pop <- "BRL"
CAT_idepth$pop <- "CAT"
CHE_idepth$pop <- "CHE"
EAG_idepth$pop <- "EAG"
FOR_idepth$pop <- "FOR"
MAN_idepth$pop <- "MAN"
MAI_idepth$pop <- "MAI"
MIR_idepth$pop <- "MIR"
MIS_idepth$pop <- "MIS"
MSC_idepth$pop <- "MSC"
MUS_idepth$pop <- "MUS"
OCQ_idepth$pop <- "OCQ"
STE_idepth$pop <- "STE"
SWN_idepth$pop <- "SWN"
TAQ_idepth$pop <- "TAQ"
TWO_idepth$pop <- "TWO"

BAD_gt_8X <- read.table("Genotypes/BAD_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BEI_gt_8X <- read.table("Genotypes/BEI_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BET_gt_8X <- read.table("Genotypes/BET_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
BRL_gt_8X <- read.table("Genotypes/BRL_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
CAT_gt_8X <- read.table("Genotypes/CAT_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
CHE_gt_8X <- read.table("Genotypes/CHE_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
EAG_gt_8X <- read.table("Genotypes/EAG_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
FOR_gt_8X <- read.table("Genotypes/FOR_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MAN_gt_8X <- read.table("Genotypes/MAN_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MAI_gt_8X <- read.table("Genotypes/MAI_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MIR_gt_8X <- read.table("Genotypes/MIR_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MIS_gt_8X <- read.table("Genotypes/MIS_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MSC_gt_8X <- read.table("Genotypes/MSC_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
MUS_gt_8X <- read.table("Genotypes/MUS_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
OCQ_gt_8X <- read.table("Genotypes/OCQ_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
STE_gt_8X <- read.table("Genotypes/STE_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
SWN_gt_8X <- read.table("Genotypes/SWN_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
TAQ_gt_8X <- read.table("Genotypes/TAQ_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
TWO_gt_8X <- read.table("Genotypes/TWO_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)

##merge population files (native lamprey not included)
comb_gts <- BAD_gt %>% 
  full_join(BEI_gt, by = c("CHROM","POS")) %>% 
  full_join(BET_gt, by = c("CHROM","POS")) %>% 
  full_join(BRL_gt, by = c("CHROM","POS")) %>% 
  full_join(CAT_gt, by = c("CHROM","POS")) %>% 
  full_join(CHE_gt, by = c("CHROM","POS")) %>% 
  full_join(EAG_gt, by = c("CHROM","POS")) %>% 
  full_join(FOR_gt, by = c("CHROM","POS")) %>% 
  full_join(MAN_gt, by = c("CHROM","POS")) %>% 
  full_join(MAI_gt, by = c("CHROM","POS")) %>% 
  full_join(MIR_gt, by = c("CHROM","POS")) %>% 
  full_join(MIS_gt, by = c("CHROM","POS")) %>% 
  full_join(MUS_gt, by = c("CHROM","POS")) %>% 
  full_join(OCQ_gt, by = c("CHROM","POS")) %>% 
  full_join(STE_gt, by = c("CHROM","POS")) %>% 
  full_join(SWN_gt, by = c("CHROM","POS")) %>% 
  full_join(TAQ_gt, by = c("CHROM","POS")) %>% 
  full_join(TWO_gt, by = c("CHROM","POS"))

comb_depth <- BAD_depth %>% 
  full_join(BEI_depth, by = c("CHROM","POS")) %>% 
  full_join(BET_depth, by = c("CHROM","POS")) %>% 
  full_join(BRL_depth, by = c("CHROM","POS")) %>% 
  full_join(CAT_depth, by = c("CHROM","POS")) %>% 
  full_join(CHE_depth, by = c("CHROM","POS")) %>% 
  full_join(EAG_depth, by = c("CHROM","POS")) %>% 
  full_join(FOR_depth, by = c("CHROM","POS")) %>% 
  full_join(MAN_depth, by = c("CHROM","POS")) %>% 
  full_join(MAI_depth, by = c("CHROM","POS")) %>% 
  full_join(MIR_depth, by = c("CHROM","POS")) %>% 
  full_join(MIS_depth, by = c("CHROM","POS")) %>% 
  full_join(MUS_depth, by = c("CHROM","POS")) %>% 
  full_join(OCQ_depth, by = c("CHROM","POS")) %>% 
  full_join(STE_depth, by = c("CHROM","POS")) %>% 
  full_join(SWN_depth, by = c("CHROM","POS")) %>% 
  full_join(TAQ_depth, by = c("CHROM","POS")) %>% 
  full_join(TWO_depth, by = c("CHROM","POS"))

comb_gt8X <- BAD_gt_8X %>% 
  full_join(BEI_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(BET_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(BRL_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(CAT_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(CHE_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(EAG_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(FOR_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(MAN_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(MAI_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(MIR_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(MIS_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(MUS_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(OCQ_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(STE_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(SWN_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(TAQ_gt_8X, by = c("CHROM","POS")) %>% 
  full_join(TWO_gt_8X, by = c("CHROM","POS"))

comb_gt8X <- comb_gt8X %>% 
  full_join(BET_gt_8X, by = c("CHROM","POS"))
#rbind AC and MAF files
AC_comb <- rbind(BAD_AC,BEI_AC,BET_AC,BRL_AC,CAT_AC,CHE_AC,EAG_AC,FOR_AC,MAI_AC,MAN_AC,MIR_AC,MIS_AC,MUS_AC,OCQ_AC,STE_AC,SWN_AC,TAQ_AC,TWO_AC)
MAF_comb <- rbind(BAD_MAF,BEI_MAF,BET_MAF,BRL_MAF,CAT_MAF,CHE_MAF,EAG_MAF,FOR_MAF,MAI_MAF,MAN_MAF,MIR_MAF,MIS_MAF,MUS_MAF,OCQ_MAF,STE_MAF,SWN_MAF,TAQ_MAF,TWO_MAF)
idepth_comb <- rbind(BAD_idepth,BEI_idepth,BET_idepth,BRL_idepth,CAT_idepth,CHE_idepth,EAG_idepth,FOR_idepth,MAI_idepth,MAN_idepth,MIR_idepth,MIS_idepth,MSC_idepth,MUS_idepth,OCQ_idepth,STE_idepth,SWN_idepth,TAQ_idepth,TWO_idepth)

comb_gt8X <- comb_gt8X %>% 
  mutate(ID=paste0(CHROM,"-",POS))

save(comb_gts,file="Genotypes/all_PM_gts.rda")
save(comb_depth,file="Genotypes/all_PM_depth.rda")
save(comb_gt8X,file="Genotypes/all_PM_gts_8X.rda")
match_file <- read.table("Summaries/chromosome_names.txt",header = F,sep = "\t",na.strings = T)
colnames(match_file) <- c("CHROM","scaf2")

link_file <- read.table("Summaries/rapture_seqs_liftover.txt",header = F,sep = "\t",na.strings = T)
link_file <- link_file[-grep(pattern = "Fail",x = link_file$V4),]
colnames(link_file) <- c("scaf1","start1","end1","arrow","scaf2","start2","end2","mapratio")
link_file1 <- merge(link_file,match_file)
save(link_file1,file="Rapture_all_genomes.rda")

#make a bed file for the rapture sequences
rapture2 <- rapture1 %>% 
  select(CHROM,min,max)

write.table(rapture2,"rapture_seqs.bed",quote = F,append = F,row.names = F)
