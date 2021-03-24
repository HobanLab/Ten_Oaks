##########################
######## Libraries #######
##########################
library(adegenet)
library(hierfstat)

#####################################
############ Directories ############
#####################################
##data file input paths
allpop_genind_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\genind_files\\allpop\\"
tenoak_df_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\Oak_Score_Df"

##species list 
oak_species_list <- c("QUAC","QUAJ","QUHI","QUPA")

#############################################
############ Conversion Code ################
#############################################
##all pops conversion 
setwd(allpop_genind_path)

##write in genind files 
allpop_arp_list <- list.files(pattern = "_wild_pop.arp$")

##convert

for(o in 1:length(allpop_arp_list)){
  
  arp2gen(allpop_arp_list[[o]])
  
}

#######################################
############ load in files ############
#######################################
setwd(allpop_genind_path)

##write in genind files 
ten_oaks_genind_wild_list <- list.files(pattern = "_wild_pop.gen$")

##conversion list
tenoak_wild_genind_list <- list()

##read in data frame 
tenoak_wild_df_list <- list.files(path = tenoak_df_path, 
                                  pattern = "_wild_pops.csv$")

##no md 

##storage list for data frames 
tenoak_wild_dfs <- list()

##loop to load in genind files
for(o in 1:length(ten_oaks_genind_wild_list)){
  
  ##load in genind files
  tenoak_wild_genind_list[[o]] <- read.genepop(ten_oaks_genind_wild_list[[o]], ncode = 3)
  ##load in data frames
  tenoak_wild_dfs[[o]] <- read.csv(paste0(tenoak_df_path,"\\",tenoak_wild_df_list[[o]]))
  ##name individuals
  rownames(tenoak_wild_genind_list[[o]]@tab) <- tenoak_wild_dfs[[o]][,1]
  ##names populations
  levels(tenoak_wild_genind_list[[o]]@pop) <- unique(tenoak_wild_dfs[[o]][,2])
  
  
}

#########################################
########### Fst Calculations ############
#########################################

##write a loop to do conversions
no_md_list <- list()
hierfstat_list <- list()
pwfst_list <- list()

for(o in 1:length(tenoak_wild_genind_list)){
  
  ##reduce by too much missing data 
  no_md_list[[o]] <- missingno(tenoak_wild_genind_list[[o]], type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)
  
  ##rename populations 
  hierfstat_list[[o]] <- genind2hierfstat(tenoak_wild_genind_list[[o]])
  
  ##calculate pwfst
  pwfst_list[[o]] <- pairwise.neifst(hierfstat_list[[o]])
  
  ##write out pwfst tables 
  write.csv(pwfst_list[[o]],paste0("G:\\Shared drives\\Emily_Schumacher\\ten_oaks_gen\\",oak_species_list[[o]],"_pwfst.csv"))
  
}


