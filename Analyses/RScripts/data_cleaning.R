###Data cleaning file to organize datafiles of 
##the 10 oaks+ individuals

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)

#####################
#     Analysis      #
#####################
setwd("C:/Users/eschumacher/Documents/GitHub/Ten_Oaks/Data_Files/Adegent_Files")
allsp_genind <- list.files(path = "Hoban_genind", pattern = "total.gen")

for(gen in 1:length(allsp_genind)){
  
  sp_gen <- read.genepop(paste0("Hoban_genind/", allsp_genind[[gen]]), ncode = 3)
  
  genind2genalex(sp_gen, paste0(gsub("_.*","",allsp_genind[[gen]]), "_genalex.csv"))
  
}

