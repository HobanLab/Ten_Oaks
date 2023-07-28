###################
#     Library     #
###################

library(diveRsity)
library(adegenet)
library(poppr)

#####################
#     Analyses      #
#####################
setwd("C:/Users/eschumacher/Documents/GitHub/Ten_Oaks/Data_Files")

arp2gen("Adegent_Files/QUTO_allpop.arp")

#read in genepop file 
QUTO_allpop_genind <- read.genepop("Adegent_Files/QUTO_allpop.gen",
                                   ncode = 3)

##reduce genind file for individuals with greater than 25% missing data 
QUTO_allpop_gen <- missingno(QUTO_allpop_genind, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##run clone check
#convert genind object to a genclone object 
QUTO_genclone <- as.genclone(QUTO_allpop_gen)

#identify multi-locus genotypes (non-clones)
QUTO_mlg <- mlg.id(QUTO_genclone)

#create clone index 
QUTO_clone_index <- which(sapply(QUTO_mlg, function(x) length(x)>1))

#list clones
QUTO_clone_id <- list()

#create a list of the clone individuals 
for(clones in 1:length(QUTO_clone_index)) QUTO_clone_id[[clones]] <- QUTO_mlg[[QUTO_clone_index[[clones]]]]

#then create genind objects without clones 
QUTO_genind_nocl <- clonecorrect(QUTO_allpop_gen)

#name pops 
levels(QUTO_genind_nocl@pop) <- c("Garden","SR_G1","SR_G2",
                                  "SR_G3", "SR_G5",
                                  "SR_G6","C_PB", "C_GC", "C_MO",
                                  "C_TR", "SCr", "SCl_E1","SCl_E2",
                                  "SCl_E3", "SCl_LTG", "SCl_B",
                                  "SCl_SSG4","SCl_HC")

#analyses
genind2genalex(QUTO_genind_nocl, "CSV_Files/QUTO_nomd_nocl.csv")

##load in by island genotype file 
arp2gen("Data_Files/Adegent_Files/QUTO_byisland_genind.arp")




