#Nb and Ns estimates
#1. Calculate Nb - PwoP method and both Ns estimates using reconstructed pedigree
#2. Extract Nb - LD method from NeEstimator tabular output
#3. Extract Nb - SF method from Colony Ne output

#load libraries
library(tidyverse)

#load functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_uncert.R")
source("Homebrew/Ns_calc.R")

#!#functions require the following dependancies: vegan


##1. Calculate Nb - PwoP method and both Ns estimates using reconstructed pedigree
#!# PwoP estimate not used in publication anymore, kept in to calculate kbar and Vk

pops <- c("BAD","BEI","BET","BRL",
          "CAT","CHE","EAG","FOR",
          "MAI","MAN","MIR2016",
          "MIS","MUS","OCQ","STE",
          "SWN","TAQ","TWO")
i <- 1
Nb_PwoP <- data.frame(matrix(nrow = length(pops),ncol = 7))
colnames(Nb_PwoP) <- c("Pop","SampSize","PwoP_Nb","kbar","Vk", "PwoP_LCI", "PwoP_HCI")
Ns_all <- data.frame(matrix(nrow = length(pops),ncol = 6))
colnames(Ns_all) <- c("Pop","Ns" ,"Ns_Chao","Chao_uncert","Ns_Jackknife","Jackknife_uncert")

#calculation loop
for (i in 1:length(pops)) {
  print(i)
  #read in file
  df <- readLines(paste0("SoftwareOutput/",pops[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  df <- strsplit(df,"\\s+")
  df1 <- matrix(unlist(df),byrow = T)
  df1 <- df1[df1 != ""]
  df1 <- matrix(df1,ncol = 5,byrow = T)
  df1 <- as.data.frame(df1)
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]
  
  #read in configarchive for PwoP uncertainty
  ca_tmp <- readLines(paste0("SoftwareOutput/",pops[i],".Output.data.ConfigArchive"))
  
  #calculate Nb_PwoP
  PwoP_tmp <- PwoP(df1)
  uncert <- PwoP_uncert(ca = ca_tmp,bc = df1)
  names(uncert) <- c("CI_Lower","CI_Upper")
  tmp <- c(pops[i],length(df1$OffspringID),PwoP_tmp,uncert)
  Nb_PwoP[i,] <- tmp
  
  #calculate extrapolated Ns
  Ns_tmp <- Ns_calc(df1)
  Ns_all[i,] <- c(pops[i],Ns_tmp[[3]],Ns_tmp[[2]]$chao,Ns_tmp[[2]]$chao.se,Ns_tmp[[2]]$jack1,Ns_tmp[[2]]$jack1.se)
}


##2. Extract Nb - LD method from NeEstimator tabular output
Nb_LD <- read.table("SoftwareOutput/All_Neestimator_2019_age1LDxLD.txt",header = T)
Nb_LD1 <- Nb_LD %>% 
  filter(CritValue == "0.05") %>% 
  select(Pop,SampSize,Nb_LD,Jack_LCI,Jack_HCI) %>% 
  rename(LD_LCI=Jack_LCI,LD_HCI=Jack_HCI)

#3. Extract Nb - SF method from Colony Ne output
i <- 1
Nb_SF <- data.frame(matrix(nrow = length(pops),ncol = 4))
colnames(Nb_SF) <- c("SF_Nb","SF_LCI","SF_HCI","Pop")
for (i in 1:length(pops)) {
  df <- readLines(paste0("SoftwareOutput/",pops[i],".Output.data.Ne"))
  Ne <- grep(pattern = "Ne",x = df,value = T)
  Ci <- grep(pattern = "CI95",x = df, value = T)
  
  tmp <- data.frame(Nb=gsub(" ","",Ne[1]),CILow=gsub(" ","",Ci[1]),CIHigh=gsub(" ","",Ci[2]))
  tmp1 <- tmp %>% 
    gather(key = "stat",value = "value") %>% 
    separate(value,into = c("tmp","value"),sep = "=") %>% 
    select(-tmp) %>% 
    spread(key = stat,value = value) %>% 
    select(Nb,CILow,CIHigh)
  tmp1$Pop <- pops[i]
  Nb_SF[i,] <- tmp1
}


##combining all estimates into a table, along with Vk and kbar
Nb_Ns <- Nb_LD1 %>% 
  full_join(Nb_PwoP,by = "Pop") %>% 
  full_join(Nb_SF,by="Pop") %>% 
  full_join(Ns_all,by="Pop")

write.table(Nb_Ns,file = "Output/genetic.estimates.txt",append = F,quote = F,sep = "\t")
