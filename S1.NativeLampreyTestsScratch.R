library(tidyverse)
library(adegenet)
library(ape)

nl <- read.table("Genotypes/MSC_GT_8X.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
load("Genotypes/all_PM_gts_8X.rda")
load("Rapture_all_genomes.rda")
load("Summaries/SNP_summaries_targets.rda")
source("Homebrew/match_tags.R")
source("Homebrew/vcf_colony.R")
rapture1 <- link_file1 %>% 
  rename(min=start2,max=end2) %>% 
  mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
  select(ID, CHROM, min, max)

colSums(nl == "./.")
nl <- nl %>% select(-LA_06)
gt_missing <- rowSums(nl == "./.")
nl$pGT <- (1-(gt_missing/26))

nl_SNPs <- nl %>% filter(pGT > 0.95) %>% select(CHROM,POS)
target <- match_tags(SNPs = nl_SNPs,tags = rapture1)
target1 <- target[[2]]
nl_targets <- merge(target1,nl)
SNPsumm <- SNPsumm %>% 
  filter(pGT > 0.95)

test_gts <- comb_gt8X[,sample(3:ncol(comb_gt8X),30,replace = F)]
test_gts <- cbind(comb_gt8X[,1:2],test_gts)
comb_targets <- test_gts %>% 
  mutate(ID=paste(CHROM,POS,sep = "-"))
nl_targets <- nl_targets %>% 
  mutate(ID=paste(CHROM,POS,sep = "-"))
all_SNPs <- comb_targets %>% 
  full_join(nl_targets, by = c("ID","CHROM","POS"))

all_SNPs[is.na(all_SNPs)] <- "0/0"

col1 <- all_SNPs %>% 
  select(ID,everything()) %>% 
  select(-CHROM:-POS) %>% 
  select(-target,-pGT)
row.names(col1) <- col1$ID
SNPs <- rownames(col1)
indiv <- colnames(col1)
indiv <- indiv[-1]
col2 <- col1 %>% select(-ID)

loci <- data.frame(ID=SNPs,stringsAsFactors = F)
loci <- loci %>% 
  separate(ID, into = c("scaffold","pos"),sep = "-")
loc <- data.frame(ID=indiv,stringsAsFactors = F)
locs <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(spp)
PM_subpops <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
#transpose matrix
col2 <- data.frame(t(col2),stringsAsFactors = F)
rownames(col2) <- indiv
colnames(col2) <- SNPs

#check amount of missing data per individual
col2 <- data.frame(lapply(col2, gsub, pattern = "./.", replacement = NA, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/0", replacement = 0, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/1", replacement = 1, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "1/1", replacement = 2, fixed = TRUE))
save(col2,file = "Genotypes/gentype.for.PCA.test.rda")
#convert to genlight object for PCA
snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$spp)
indNames(snp) <- indiv
save(snp,file = "genlight_pca_test.RData")
#running PCA analysis
#getting colors
locs$col <- NA
for(i in 1:length(locs$spp)){
  if(locs$spp[i] == "LA"){locs$col[i] <- "blue"}
  if(locs$spp[i] == "IF"){locs$col[i] <- "dark green"}
  if(locs$spp[i] == "PM"){locs$col[i] <- "purple"}
}

#run pca
pca <- glPca(snp,nf = 10)
save(pca,file = "pca_combined.RData")
#plot pca
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
plot(pca$scores[,1], pca$scores[,3],
     col=locs$col,cex=0.5)
#make a tree




