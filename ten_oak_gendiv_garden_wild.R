##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(stringr)
library(tidyr)
library(hierfstat)
library(poppr)
library(Demerelate)
library(rworldmap)
library(data.table)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(plotrix)
library(ggpmisc)
library(factoextra)
library(GISTools)
library(raster)
library(rgdal)
library(sp)
library(ggpubr)
library(RColorBrewer)

#####################################
############ Directories ############
#####################################
##data file input paths
tenoak_genind_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\genind_files\\garden_wild"
tenoak_df_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\Oak_Score_Df"

ten_oaks_shared <- "G:\\Shared drives\\Emily_Schumacher\\ten_oaks_gen\\"

species_names <- c("quac", "quaj", "quhi","qupa")

#############################################
############ Conversion Code ################
#############################################
##all pops conversion 
setwd(tenoak_genind_path)

##write in genind files 
allpop_arp_list <- list.files(pattern = ".arp$")

##convert

for(o in 1:length(allpop_arp_list)){
  
  arp2gen(allpop_arp_list[[o]])
  
}

#######################################
############ load in files ############
#######################################
setwd(tenoak_genind_path)

##write in genind files 
ten_oaks_genind_list <- list.files(pattern = ".gen$")

##
tenoak_genind_list <- list()

for(o in 1:length(ten_oaks_genind_list)){
  
  tenoak_genind_list[[o]] <- read.genepop(ten_oaks_genind_list[[o]], ncode = 3)
  
}

##read in data frame 
tenoak_df_list <- list.files(path = tenoak_df_path, 
                             pattern = "_garden_wild_df.csv$")

##load in data frame list
tenoak_dfs <- list()

##loop to load in these datafiles 
for(df in 1:length(tenoak_df_list)){
  
  tenoak_dfs[[df]] <- read.csv(paste0(tenoak_df_path,"\\",tenoak_df_list[[df]]))
  
  
}

##name individuals in genind file 
for(n in 1:length(tenoak_dfs)){
  
  rownames(tenoak_genind_list[[n]]@tab) <- tenoak_dfs[[n]][,1]
  
}

#####remove individuals based on missing data
##list no missing 
tenoaks_nomd_list <- list()

for(i in 1:length(tenoak_genind_list)){
  
  tenoaks_nomd_list[[i]] <- missingno(tenoak_genind_list[[i]], type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) 
  
}
###QUAC - 7 genotypes missing a lot of data 
###QUAJ - no individuals more than 25% 
###QUHI - 3 genotypes with more than 25% missing data 
###QUPA - 1 genotype removed 

####reduce dfs 
##list reduced data frame 
tenoak_garden_wild_nomd_df <- list()

for(i in 1:length(tenoaks_nomd_list)){
  
  tenoak_garden_wild_nomd_df[[i]] <-  tenoak_dfs[[i]][tenoak_dfs[[i]][,1] %in% rownames(tenoaks_nomd_list[[i]]@tab),]
  
  write.csv(tenoak_garden_wild_nomd_df[[i]], paste0(tenoak_df_path, "\\",species_names[[i]], "_nomd_garden_wild_df.csv"))
  
}

#############################################
######### Clone Correct Code ################
#############################################
clone_test <- list()
tenoak_list <- list()
tenoak_clone_index <- list()
tenoak_genind_nocl_list <- list()
popr_nocl <- list()

for(i in 1:length(tenoaks_nomd_list)){
    ##try clone analyses for QUAC
    clone_test[[i]] <- as.genclone(tenoaks_nomd_list[[i]])
    tenoak_list[[i]] <- mlg.id(tenoaks_nomd_list[[i]])
    #Function to pull out individual indices where clone length greater than 1
    tenoak_clone_index[[i]] <- which(sapply(tenoak_list[[i]],function(x) length(x)>1))

  ##remove clones from test doc
  popr_nocl[[i]] <- clonecorrect(clone_test[[i]])
  #genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")
  #Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
  tenoak_genind_nocl_list[[i]] <- genclone2genind(popr_nocl[[i]])

}

################################################
############# Relatedness Code #################
################################################
##list df files 
tenoak_nomd_df_list <- list.files(path = tenoak_df_path, pattern = "_nomd_garden_wild_df.csv$")

##write in data files 
tenoak_rel <- list()

for(df in 1:length(tenoak_nomd_df_list)){
  
  tenoak_rel[[df]] <- read.csv(paste0(tenoak_df_path, "\\", tenoak_nomd_df_list[[df]]))
  
}

####do relatedness code 
##make list of relatedness docs 
tenoak_rel_list <- list()
##half siblings list
tenoak_halfsib_list <- list()

##run relatedness in a loop 
for(o in 1:length(tenoak_rel)){
  
  tenoak_rel_list[[o]] <- Demerelate(tenoak_rel[[o]], object = T, value = "loiselle")
  tenoak_halfsib_list[[o]] <- names(which(unlist(tenoak_rel_list[[o]]$Empirical_Relatedness) > 0.25))
  
}
length(tenoak_halfsib_list[[1]])
length(tenoak_halfsib_list[[2]])
length(tenoak_halfsib_list[[3]])
length(tenoak_halfsib_list[[4]])


################################################
########### Genetic Diversity ##################
################################################
##list out allelic richness 
tenoak_allrich <- list()
tenoak_hexp <- list()

##loop over allelic richness
for(i in 1:length(tenoak_genind_nocl_list)) {
  
    tenoak_allrich[[i]] <- colSums(allelic.richness(tenoak_genind_nocl_list[[i]])$Ar)/length(tenoak_genind_nocl_list[[i]]@loc.n.all)	#remember to divide by number of loci!!
}

##create a table
tenoak_allrich_matrix <- rbind(tenoak_allrich[[1]], tenoak_allrich[[2]], tenoak_allrich[[3]], tenoak_allrich[[4]])

##name rows and columns 
rownames(tenoak_allrich_matrix) <- species_names
colnames(tenoak_allrich_matrix) <- c("Garden", "Wild")

##switch orientation 
tenoak_allrich_df <- rbind(tenoak_allrich_matrix[,1], tenoak_allrich_matrix[,2])

##barplot 
pdf("G:\\Shared drives\\Emily_Schumacher\\ten_oaks_gen\\tenoak_allrich_barplot.pdf")
barplot(tenoak_allrich_df, beside = T, col = c("darkseagreen1","forestgreen"), ylim = c(0,25), 
        ylab = "Allelic Richness", xlab = "Species", names.arg = c("Quercus acerifolia", "Quercus ajoensis", 
                                                                   "Quercus hinckleyi", "Quercus pacifica"))
legend('top', col = c("darkseagreen1","forestgreen"), legend = c("Garden", "Wild"), pch = 15)
abline(h = 0)
dev.off()
