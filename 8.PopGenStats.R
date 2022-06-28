#Generating comparative statistics for all stream populations
#expected and observed heterozygosity
#inbreeding coefficient
#Fis calculations
#Fst calculations
#Hierarchal AMOVA among and between lakes*****
#IBD and relationship to Fst
#prep for structure analysis, run in R?

#load libraries
library("adegenet")
library("hierfstat")
library("apex")
library("mmod")
library("ape")
library("pegas")
library("tidyverse")

#load data
load("Summaries/Alllocs_cohorts.rda")
load("SNPsets/LDcomb_SNPs.rda")
load("Input/all_PM_gts_8X.rda")

#load functions?
source("Homebrew/vcf_colony.R")
##Data prep####
#random selection of one full-sibling from each family 
alllocs_cohorts$fullsib <- paste(alllocs_cohorts$MotherID,alllocs_cohorts$FatherID,sep = "_")
sib_select <- alllocs_cohorts %>% 
  group_by(loc,fullsib) %>% 
  slice_sample()

#count sample size for each pop and remove locations with less than 10 individuals
sib_select %>% group_by(loc) %>% count()
sib_select <- subset(sib_select,!(sib_select$loc == "BAD" | sib_select$loc == "CHE" | sib_select$loc == "SWN"))

#merging reduced data set with NeEstimate combined SNP set
selectedSNPs <- merge(LD_select,comb_gt8X)
#selected 500 SNPs with most complete data to simplify calculations
selectedSNPs1 <- selectedSNPs %>% 
  arrange(desc(pGT)) %>% 
  slice_head(n = 1500)
hist(selectedSNPs1$pGT)
hist(selectedSNPs1$het)
table(selectedSNPs1$CHROM)
#merge selected SNPs with selected individuals
selectedSNPs2 <- selectedSNPs1 %>% 
  select(-het:-target) %>% 
  select(-CHROM:-POS) %>% 
  gather(key = "indiv",value = "SNP",-ID) %>% 
  spread(key = "ID",value = "SNP") %>% 
  rename(OffspringID = indiv)
gts_comb <- merge(sib_select,selectedSNPs2)

#converting into geneid and hierfstat objects for stat calculations
SNPs <- gts_comb[,-1:-11]
indiv <- gts_comb[,1]
pop <- gts_comb[,5]
tmpSNPname <- paste("SNP",seq(1,ncol(SNPs)),sep = "_")
gene1 <- df2genind(SNPs, ind.names = indiv, 
                   pop = pop, 
                   loc.names = tmpSNPname,
                   sep = "/",
                   NA.char = "_/_")
gene2 <-genind2loci(gene1)
gene3 <- as.genclone(gene1)
hierfgts <- genind2hierfstat(gene1)
##Per-SNP statistics####
#expected and observed heterozygosity for each SNP and each populations
statsout <- basic.stats(hierfgts)

#means per population
colMeans(statsout[["Fis"]],na.rm = T)
colMeans(statsout[["Ho"]],na.rm = T)
colMeans(statsout[["Hs"]],na.rm = T)

#number of alleles - they're SNPs they should all be two
nAll(gene1)

#allele freq PCA
freqs <- pop.freq(hierfgts)

#!# no interpretable separation
x <- indpca(hierfgts) 
plot(x, cex = 0.7)

#HWE test
#!# slow
hwe <- hw.test(gene1, B = 1000)
##Per-population statistics####
##Fis for each population
#confidence intervals for per population Fis (across 1500 SNP set)
boot.ppfis(hierfgts)

##Population Comparison Analysis####
##between population Fst
fsts <- genet.dist(gene1, method = "WC84")
fsts <- pairwise.WCfst(hierfgts)
fstuncert <- boot.ppfst(hierfgts)
##sort by lake and run an AMOVA
#making a lake column for data
gts_comb$lake <- NA
gts_comb[which(gts_comb$loc == "CAT"),]$lake <- "Erie"
gts_comb[which(gts_comb$loc == "BMR" | 
                 gts_comb$loc == "EAG" | 
                 gts_comb$loc == "OCQ"),]$lake <- "Huron"
gts_comb[which(gts_comb$loc == "BEI" | 
                 gts_comb$loc == "FOR" | 
                 gts_comb$loc == "MAI" | 
                 gts_comb$loc == "MAN" | 
                 gts_comb$loc == "MUS"),]$lake <- "Michigan"
gts_comb[which(gts_comb$loc == "STE"),]$lake <- "Ontario"
gts_comb[which(gts_comb$loc == "BET" | 
                 gts_comb$loc == "BRL" | 
                 gts_comb$loc == "MIR" | 
                 gts_comb$loc == "MIS" | 
                 gts_comb$loc == "TAQ" | 
                 gts_comb$loc == "TWO"),]$lake <- "Superior"
lake <- gts_comb$lake
#test for differentiation based on lake
loci <- hierfgts[,-1]

#test if population is a sig factor on genetic differentiation
#!#significant
popdifftest <- test.g(loci, level = pop,nperm = 1000) 

#test if lake is a sig factor on genetic differentiation
#!#significant
difftest <- test.between(loci, test.lev = pop, rand.unit = lake, nperm = 1000) 

##IBD testing????
distgenEUCL <- dist(gene1, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
gene2 <- genind2genpop(gene1,pop = pop)
distgen <- dist.genpop(gene2,method = 2)
dist(gene1$other)























