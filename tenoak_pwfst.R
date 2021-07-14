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

###read in three oaks 

###prep data table 
threeoak_df <- matrix(nrow = 3, ncol = 3)
colnames(threeoak_df) <- c("mean","min","max")
rownames(threeoak_df) <- c("QUGE","QUOG","QUBO")
#-----------------------#
#		PREP			#
#-----------------------#

library("adegenet"); library("poppr"); library(hierfstat); library("Demerelate")

#Modification of existing function from adegenet
repool_new<- function(genind_obj,vect_pops){
  genind_obj_sep<-seppop(genind_obj)
  genind_obj_merge<-genind_obj_sep[[vect_pops[1]]]
  for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
  genind_obj_merge
}

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/")
#setwd(paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))	#setwd(paste("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
#setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii")

pop_sizes<-read.csv("3_oak_pop_sizes.csv",header=F)
species_names<-c("Qgeorgiana","Qoglethorpensis","Qboyntonii")
file_names<-c("Qg_wild_w_NC.gen","Qo_wild_w_gos.gen", "Qb_wild_w_ALL.gen")
file_names2<-c("Qg_wild_rel_w_NC.csv", "Qg_wild_rel_w_NC_ESTSSRs.csv", "Qg_wild_rel_w_NC_gSSRs.csv", "Qo_wild_rel_w_gos.csv", "Qb_wild_rel_w_ALL.csv")

min_pop_size_keep<-c(2,5,10)

for (min_p in min_pop_size_keep){
  
  par(mfrow=c(2,3))	#for IBD plots
  
  for (sp in 1:length(species_names)){
    this_species<-species_names[sp]
    
    setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/")
    #------------------------------------------------------------------------------#
    #				SECTION ONE
    #------------------------------------------------------------------------------#
    
    ##############
    #  CLONES	 #
    ##############
    
    #import the data- note individuals with no data at all are dropped
    ade_test<-read.genepop(paste("genetic_data/",file_names[sp],sep=""),ncode=3)
    pop_names<-levels(ade_test@pop)
    
    #First put into poppr format
    popr_test <- as.genclone(ade_test)
    strata(popr_test) <- other(popr_test)$population_hierarchy[-1]
    list_a<-mlg.id(popr_test)
    #Function to pull out individual indices where clone length greater than 1
    clone_index<-which(sapply(list_a,function(x) length(x)>1))
    
    #This removes clones and then saves as new file for Genealex if desired
    popr_nocl<-clonecorrect(popr_test,strata=~Pop)
    #genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")
    
    #Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
    GI_nocl<-genclone2genind(popr_nocl); 	GP_nocl<-genind2genpop(GI_nocl)
    
    ###################
    #  BASIC STATS	  #
    ###################
    
    #narrow down to populations with >10 individuals- call these GP_sub, GI_sub for subset
    pop_keep<- which(as.vector(table(GI_nocl@pop)>=min_p))
    GI_sub<-repool_new(GI_nocl,pop_keep);	GP_sub<-GP_nocl[pop_keep,]
    samp_size<-table(GI_nocl@pop)[pop_keep]
  
    ###################
    #  PAIRWISE FST	  #
    ###################
    
    ##add code to convert 
    GI_fst <- genind2hierfstat(GI_sub)
    sm_fst_mat <- pairwise.neifst(GI_fst)
    rownames(sm_fst_mat)<-pop_names[pop_keep];	colnames(sm_fst_mat)<-pop_names[pop_keep]
    sm_fst_mat[sm_fst_mat==0]<-NA
    threeoak_df[sp,1]<- mean(apply(sm_fst_mat,2,mean,na.rm=T))
    threeoak_df[sp,2] <- min(apply(sm_fst_mat,2,min,na.rm=T))
    threeoak_df[sp,3] <- max(apply(sm_fst_mat,2,max,na.rm=T))
  }
}

    
###havardi

#library(adegenet); library(hierfstat)
#library(parallel);	library(doParallel) #will load foreach
#source("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/Fa_sample_funcs.R")
#source("/home/user/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/Fa_sample_funcs.R")
#colMax <- function(data) sapply(data, max, na.rm = TRUE)

#Get working directory
#which_comp<-getwd()
#if (grepl("shoban",which_comp)) prefix_wd<-"C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/"
#if (grepl("shoban.DE",which_comp)) prefix_wd<-"C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/"
#f (grepl("home",which_comp)) prefix_wd<-"/home/user/Dropbox/Projects/IN_PROGRESS/"
#setwd(prefix_wd)


#--------------------------------------------------------------------------------------------------
#	COMPARE EX SITU AND IN SITU NOW, AND THEN RESAMPLE WILD FOR MINIMUM SAMPLING/ OPTIMAL	#
#-------------------------------------------------------------------------------------------------

#####################################################
#													#
#	CALCULATE GEN DIV CAPTURE (% gen div in 		#
#	in all gardens and in each garden				#
#	AND MAKE .CSVs and R files for plots later		#
#													#
#####################################################

setwd("C:\\Users\\eschumacher\\Documents\\Qhavardii_ex_situ-main")
gard_names<-read.csv("Key_to_Garden_POP.txt")[,2]

wild_files<-c("","_E","_W")
reg_names<-c("all","E","W")

#This will run over a loop of "Include all alleles (n_to_drop=0)" and "Include only alleles present in more than two copies (n_to_drop=2)"
for (n_to_drop in c(0,2)){
  if (n_to_drop==2) n_drop_file<-""
  if (n_to_drop==0) n_drop_file<-"_dr_0"
  
  #####################################################
  #													#
  #	Analysis by region- Total, East, West			#
  #	will output ex_vs_in_situ results csv by reg	#
  #													#
  #####################################################
  
  wild_results<-matrix(nrow=3,ncol=9+1)
  alleles_existing_by_sp<-matrix(nrow=3,ncol=9)
  #The loop below over 'reg' or 'regions' will go through three sets: analysis of all sample, only east samples and only west samples
  #The files are 'all wild pops plus all garden samples', 'all wild pops plus garden samples sourced from E', 'all wild pops plus garden samples sourced from W' 
  #The "wild_sets" will subset the wild reference to only the E and W for the reg=2 and reg=3 analysis
  #1-19 is east, 20-35 is west
  #The "garden_sets" idntify the populations that are garden samples
  #All three "by_garden" files have all wild pop'ns- will subset to east and west below
  #However the "_E" only has seedlings from East and "_W" only has seedlings from West
  #Note only four gardens have West seedlings
  
  wild_sets<-list(list(1,2:19,20:35),list(1,2:19),list(20,21:35))
  garden_sets<-(list(36:43,36:43,36:39))


  
#####################################################
#													#
#	Genetic distances- FST, clusters				#
#													#
#####################################################

  spp_fst_df <- matrix(nrow = 1, ncol = 3)
  
#FST
reg<-1
Spp_tot_genind<-read.genepop(paste0("QH_total_garden_for_FST",wild_files[reg],".gen"),ncode=3)

##now calculate fst 
spp_fst <- genind2hierfstat(Spp_tot_genind)

spp_fst_df[,1] <- mean(sapply(pairwise.neifst(spp_fst),mean,na.rm = T),na.rm = T)
spp_fst_df[,2] <- min(sapply(pairwise.neifst(spp_fst),min,na.rm = T),na.rm = T)
spp_fst_df[,3] <- max(sapply(pairwise.neifst(spp_fst),max,na.rm = T),na.rm = T)


#m<-pairwise.fst(Spp_tot_genind)
#makeSymm <- function(m) {
#  m[upper.tri(m)] <- t(m)[upper.tri(m)]
#  return(m)
 # }
  #a<-matrix(ncol=9,nrow=9)
 # a[lower.tri(a)]<-m
  #a<-makeSymm(a)
 # rowMeans(a,na.rm=T);	mean(a[1,3:9],na.rm=T); mean(a[2,3:9],na.rm=T)	#FST East 0.0066, West 0.0159
  #t.test(a[1,3:9],a[2,3:9],pair=T)	#p=0.049

}
