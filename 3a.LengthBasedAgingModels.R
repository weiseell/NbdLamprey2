#Length-Based Age Models
#1. EM models
#2. Bayes Cluster determining models
#3. Bayes models
#4. Plots
#Written by: Ellie Weise
#Last edited: 06/25/20

#loading libraries
library(tidyverse)
library(bayesmix)
library(bmixture)

#loading homebrew functions

#load in length and weight data
df <- read.table("AgingModels/lw_2019_locations_081220.txt",header = T,sep = "\t",stringsAsFactors = F)

#manipulating data frame
df1 <- df %>% 
  mutate(ID2 = ID_indiv) %>% 
  rename(ID=ID_indiv) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  select(ID,loc,everything()) %>% 
  select(-spp:-num) %>% 
  mutate(samp = paste(loc,Year_collect,sep = "_"))
samples <- unique(df1$samp)

##2. Bayes Cluster determining models ####
RMclust <- vector(mode = "list", length = length(samples))
BDclust <- vector(mode = "list", length = length(samples))
names(RMclust) <- samples
names(BDclust) <- samples
i <- 1
for(i in 1:length(samples)){
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  
  #running BayesMix models with RM criteria
  model_tmp <- BMMmodel(len, k = 10, 
                        priors = list(kind = "independence",
                                      parameter = "priorsUncertain"))
  control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                         burn.in = 10000, n.iter = 500000)
  z <- JAGSrun(len, model = model_tmp, control = control)
  #calculating probabilities
  eta.thresh <- 0.035
  post.k <- rowSums(z$results[,grep("eta", colnames(z$results))] > eta.thresh)
  post.k <- table(post.k)/sum(table(post.k))
  
  RMclust[[samp_i]] <- post.k
  
  #BD-MCMC with bmixture
  bmixt.model <- bmixt(len, k = "unknown", iter = 500000, burnin = 100000, k_max = 4)
  #calculating probabilities
  bmix_vals <- data.frame(kval = bmixt.model$all_k, weight = bmixt.model$all_weights,stringsAsFactors = F)
  props <- bmix_vals %>% group_by(kval) %>% summarise(prop = sum(weight)/sum(bmix_vals$weight))
  BDclust[[samp_i]] <- props
}


##3. Bayes models for all locations ####
i <- 1
samples <- sort(samples)
bestk <- c(2,1,1,2,2,2,2,1,2,1,1,3,2,2,1,1,2,1,2)
Bmodels <- vector(mode = "list", length = length(samples))
names(Bmodels) <- samples

for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  
  if(bestk[i] > 1) {
    #bayesmix model
    model_tmp <- BMMmodel(len,k = bestk[i], 
                          priors = list(kind = "independence",
                                        parameter = "priorsUncertain"))
    control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                           burn.in = 10000, n.iter = 50000)
    
    z <- JAGSrun(len, model = model_tmp, control = control)
    Bmodels[[samp_i]] <- z
  }
}
#sorting for models with multiple cohorts
Bmodels[["BAD_2019"]] <- Sort(Bmodels[["BAD_2019"]],by = "mu")
Bmodels[["BRL_2019"]] <- Sort(Bmodels[["BRL_2019"]],by = "mu")
Bmodels[["CAT_2019"]] <- Sort(Bmodels[["CAT_2019"]],by = "mu")
Bmodels[["CHE_2019"]] <- Sort(Bmodels[["CHE_2019"]],by = "mu")
Bmodels[["EAG_2019"]] <- Sort(Bmodels[["EAG_2019"]],by = "mu")
Bmodels[["MAI_2019"]] <- Sort(Bmodels[["MAI_2019"]],by = "mu")
Bmodels[["MIR_2019"]] <- Sort(Bmodels[["MIR_2019"]],by = "mu")
Bmodels[["MIS_2019"]] <- Sort(Bmodels[["MIS_2019"]],by = "mu")
Bmodels[["MUS_2019"]] <- Sort(Bmodels[["MUS_2019"]],by = "mu")
Bmodels[["MIS_2019"]] <- Sort(Bmodels[["MIS_2019"]],by = "mu")
Bmodels[["SWN_2019"]] <- Sort(Bmodels[["SWN_2019"]],by = "mu")
Bmodels[["TWO_2019"]] <- Sort(Bmodels[["TWO_2019"]],by = "mu")

i <- 1
all_locs_Bayes <- data.frame(matrix(ncol=8,nrow = 0))
for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  zSort <- Bmodels[[samp_i]]
  
  if(bestk[i] > 1) {    
    #Extracting probabilities and assigning classes
    probs <- BMMposteriori(zSort, plot = F)
    probs <- as.data.frame(cbind(probs$data,t(probs$post)),stringsAsFactors = F)
    probs$clust <- apply(probs[,-1], 1, FUN = function(x){which(x == max(x))})
    probs$clust <- paste0("clust",probs$clust)
    
    #assignments from probabilities
    probs <- probs %>% arrange(V1)
    j <- 1
    tmp$clust <- NA
    for (j in 1:length(probs$V1)) {
      tmp[which(tmp$Length == probs$V1[j]),]$clust <- probs[j,]$clust
    }
    if(length(unique(tmp$clust))==1){
      tmp$clust <- "clust1"
    }
  }
  if(bestk[i]==1){
    tmp$clust <- "clust1"
  }
  all_locs_Bayes <- rbind(all_locs_Bayes,tmp)
}

#saving individual assignments
write.table(all_locs_Bayes,file = "AgingModels/lw_Bayes_assignments.txt",sep = "\t",row.names = F,col.names = T,quote = F)

##4. Plotting all models ####
#labels for plot
plotnames <- c("Bad River",
               "Betsie River",
               "Betsy River",
               "Brule River",
               "Cattaraugus River",
               "Pigeon River",
               "East Au Gres River",
               "Ford River",
               "Grand River",
               "Manistique River",
               "Manistee River",
               "Middle River",
               "Misery River",
               "Muskegon River",
               "Ocqueoc River",
               "Sterling River",
               "Swan Creek",
               "Tahquamenon River",
               "Two-Hearted River")
names(plotnames) <- sort(unique(all_locs_Bayes$samp))

all_locs_Bayes <- read.table("AgingModels/lw_Bayes_assignments.txt",header = T,sep = "\t")

tiff(filename = "Figures/LengthHistograms.tiff",height = 9,width = 12,units = "in",res = 200)
ggplot(all_locs_Bayes, aes(x=Length, fill = clust)) +
  facet_wrap(~samp, scales = "free_y",labeller = labeller(samp=plotnames)) +
  geom_histogram(aes(fill = factor(clust)),bins = 50) +
  scale_fill_manual(values = c("#000000","#cccccc","#969696","#636363"),
                    guide = F)+
  labs(x="Length (mm)",y="counts")+
  theme_bw(base_size = 8)
dev.off()
