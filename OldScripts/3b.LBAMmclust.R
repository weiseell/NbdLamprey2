#Mclust mixture analysis
#Last edited 12/13/21

#Load packages
library(mclust)
library(tidyverse)
library(wesanderson)

#Load Dataset
df <- read.table("AgingModels/lw_2019_locations_081220.txt",header = T,sep = "\t")

#adding a loc column
df1 <- df %>% 
  mutate(ID1=ID_indiv) %>% 
  separate(ID1, into = c("spp","loc","num")) %>% 
  select(-spp,-num) %>% 
  select(loc,ID_indiv,Year_collect,Length,Weight)

#remove grand river due to low sample size (not sequenced)
df1 <- subset(df1,df1$loc != "GRN")

locs <- sort(unique(df1$loc))
#loop to run mclust models for all locations to add to current cluster determination models
i <- 4
BICvals <- vector(mode = "list",length = length(locs))
names(BICvals) <- locs
for (i in 1:length(locs)) {
  #select location
  tmp <- subset(df1,df1$loc == locs[i])
  
  #order by length for mclust and select only length and weight data
  tmp <- tmp[order(tmp$Length),]
  tmpmodel <- cbind(tmp$Length,tmp$Weight)
  
  #run mclust with VVV model
  mBIC <- mclustBIC(tmpmodel, G=1:5, modelName= c("VVV"))
  BICvals[[locs[i]]] <- mBIC
  
}

#view all BIC vals for all streams
BICvals

#add mclust results to the other two methods to get bestk for all locations
#make a named vector with all bestk for individual assignments
bestK <- c(2,2,3,2,3,2,3,2,2,2,4,3,2,4,1,2,2,2)
names(bestK) <- locs
models <- vector(mode = "list",length = length(locs))
classes <- vector(mode = "list",length = length(locs))
names(models) <- locs
names(classes) <- locs
#loop for individual assignment
for (i in 1:length(locs)) {
  #select location
  tmp <- subset(df1,df1$loc == locs[i])
  
  #order by length for mclust and select only length and weight data
  tmp <- tmp[order(tmp$Length),]
  tmpmodel <- cbind(tmp$Length,tmp$Weight)
  
  #grab maximum BIC value for individual assignment of pop
  massign <- Mclust(tmpmodel,G = bestK[i], modelNames = c("VVV"))
  tmpmodelout <- data.frame(Indiv = tmp$ID_indiv,
                            Length = tmpmodel[,1],
                            Weight = tmpmodel[,2],
                            Class = massign$classification,
                            loc = locs[i])
  
  #bind classifications to other locations for output
  models[[locs[i]]] <- massign
  classes[[locs[i]]] <- tmpmodelout
}

#bind all classifications together
classes_all <- rbind(classes[["BAD"]],classes[["BEI"]],classes[["BET"]],
                     classes[["BRL"]],classes[["CAT"]],classes[["CHE"]],
                     classes[["EAG"]],classes[["FOR"]],classes[["MAI"]],
                     classes[["MAN"]],classes[["MIR"]],classes[["MIS"]],
                     classes[["MUS"]],classes[["OCQ"]],classes[["STE"]],
                     classes[["SWN"]],classes[["TAQ"]],classes[["TWO"]])

#save individual classifications for downstream analyses
save(classes_all,file = "AgingModels/MclustLW_results.rda")
plotnames <- c("Bad River",
               "Betsie River",
               "Betsy River",
               "Brule River",
               "Cattaraugus River",
               "Pigeon River",
               "East Au Gres River",
               "Ford River",
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
names(plotnames) <- sort(unique(classes_all$loc))
klabels <- data.frame(label = paste0("K = ",bestK), 
                      loc = sort(unique(classes_all$loc)))

tiff(filename = "Figures/LengthHistogramsmclustVVV.tiff",height = 9,width = 12,units = "in",res = 200)
ggplot(classes_all, aes(x=Length))+
  facet_wrap(~loc, ncol = 6,scales = "free_y",labeller = labeller(loc=plotnames))+
  geom_histogram(aes(fill = factor(Class)),bins = 50)+
  geom_text(data = klabels,x = 150, y = 0,aes(label = label))+
  scale_fill_manual(values = wes_palette("IsleofDogs1"),
                    guide = "none")+
  labs(x="Length (mm)",y="counts")+
  ggtitle("Cohort Assignments using VVV models from mclust")+
  theme_bw(base_size = 8)
dev.off()





