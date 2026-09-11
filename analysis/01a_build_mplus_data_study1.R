### 01a_build_mplus_data_study1.R — builds the Mplus .dat files for Study 1.
###
### Run before the Study 1 random-slope models in Mplus. See README for order.
###

## Paths come from config.R at the repository root. Edit that file, not this one.
if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
source(file.path(repo_directory, "R/exclusions.R"))
# Ceci Westbrook, 5/2025
# building models for Mplus using Michael Hallquist's MplusAutomation package

# Load required libraries
library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(fmri.pipeline)
library(data.table)
library(parallel)
library(doParallel)
library(lme4)
library(emmeans)
library(sciplot)
library(MplusAutomation)

######################################################################################################
######################################################################################################
#########################################   STUDY 1   ################################################
######################################################################################################
######################################################################################################

# set root directory - change for your needs
rootdir <- rootdir1
repo_directory <- repo_directory

# load some source code that we'll need to use
setwd(rootdir)
source(file.path(repo_directory, "R/get_trial_data.R"))

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
load('MMclock_clock_Aug2023.Rdata')
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)

# load in the hippocampus data, filter to within 5s of RT
load(file.path(rootdir1, "HC_clock_Aug2023.Rdata"))
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
# for the 0-to-4 analyses
# hc <- hc %>% filter(evt_time >= 0 & evt_time < 4)

# Compress data from 12 bins to 2 by averaging across anterior 6 bins and posterior 6 bins to create
# AH and PH averages
hc <- hc %>% group_by(id,run,trial,run_trial,evt_time,HC_region) %>%
  summarize(HCwithin = mean(decon_mean,na.rm=TRUE)) %>% 
  ungroup() # 12 -> 2

# not creating HCbetween as Mplus will estimate within- and between-subject variance

# Merge vmPFC and HC data; remove unscaled within-S variable (decon1)
Q <- inner_join(vmPFC,hc,by=c("id","run","run_trial","evt_time"))

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')

# select only vars of interest and merge into MRI data--can add others in if needed
# not scaling anything for Mplus
behav <- df %>% select(id,run,trial,run_trial,rt_csv,iti_ideal,iti_prev,rt_lag)
Q <- inner_join(behav, Q, by = c("id", "run", "trial", "run_trial")) %>% arrange("id","run","trial","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;

# summarize over evt_time and atlas_value and pivot to create new columns with data averaged across network
# if adding other behav vars, include them in group_by here
Q1 <- Q %>% group_by(id,run,trial,run_trial,network,HC_region,rt_csv,rt_lag) %>% 
  summarize(vmPFC_decon = mean(vmPFC_decon,na.rm=TRUE),
            HCwithin = mean(HCwithin,na.rm=TRUE)) %>% 
  ungroup()
# pivoting to create new vmPFC and HC variables for each network
Q1 <- Q1 %>% group_by(id,run,run_trial) %>% 
  pivot_wider(values_from=c(vmPFC_decon,HCwithin),names_from = 'network') %>% 
  ungroup()

# between-subject-only variables need to be grand-mean scaled.
# within-subject-only variables need to be scaled by ID.
# variables that will be included at both levels should not be scaled because Mplus will do it as 
# part of latent decomposition.
# so it turns out none of our vars should get scaled, except for age and sex!

# add in age and sex variables
demo <- read.table(file=file.path(rootdir, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q1 <- inner_join(Q1,demo,by=c('id'))
Q1$female <- relevel(as.factor(Q1$female),ref='0')
Q1$sex <- ifelse(Q1$female==0,1,
                ifelse(Q1$female==1,0,NA))
Q1$age <- scale(Q1$age)

Q1_AH <- Q1 %>% filter(HC_region == "AH")
Q1_PH <- Q1 %>% filter(HC_region == "PH")

setwd(file.path(rootdir,"Mplus"))

# saving out the data for ease of later use
save(Q1_AH,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_AHforMplus.Rdata"))
save(Q1_PH,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_PHforMplus.Rdata"))

# Now, using Michael Hallquist's MplusAutomation, will prep the dataframes to run in Mplus.
# Download these dataframes and run the models to extract random slopes for vPFC-HC connectivity
# per participant, which will be used as between-subjects variables in the MSEMs.
prepareMplusData(df = Q1_AH, filename = "Study1_HC_vmPFC_clock_AH_forMplus_taa.dat", overwrite = TRUE)
prepareMplusData(df = Q1_PH, filename = "Study1_HC_vmPFC_clock_PH_forMplus_taa.dat", overwrite = TRUE)

# Alright! Was able to run the random slopes models. Now need to extract random slopes and put them
# back into the dataframe. Can use Michael's MplusAutomation for that.

list.files(getwd()) # let's see what we got in here

# let's make a function
getslopes <- function(HC_region, network,fileheader){
  # get filename
  filename <- paste0(fileheader,HC_region,"_",network,".out")
  renamevar <- paste0(HC_region,"_",network,"_rs")
  
  print(renamevar)
  
  #read in model
  wm_rs <- readModels(filename)
  
  ## save the savedata from the mplus output to an object
  slopes <- wm_rs$savedata %>% 
    dplyr::select(ID, RS.Mean) %>% 
    group_by(ID) %>% slice_head() %>% ungroup() %>%
    rename(id = ID, !! sym(renamevar) := RS.Mean)

  return(slopes)
}

# for loading datasets that have already been prepared:
# all timepoints:
#load("Study1_HC_vmPFC_clock_AHforMplus.Rdata")
#load("Study1_HC_vmPFC_clock_PHforMplus.Rdata")

#prepareMplusData(df = Q1_AH, filename = "mmclock_HC_vmPFC_clock_AH_forMplus_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)
#prepareMplusData(df = Q1_PH, filename = "mmclock_HC_vmPFC_clock_PH_forMplus_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)

# ok we ran six models: 3 networks (DMN, CTR, LIM) x 2 HC regions (AH, PH).
fileheader <- "Study1_get_HC_vmPFC_randslopes_"
for(netname in c("CTR","DMN","LIM")){
  Q1_AH <- merge.data.frame(Q1_AH,getslopes("AH",netname,fileheader),by="id")
  Q1_PH <- merge.data.frame(Q1_PH,getslopes("PH",netname,fileheader),by="id")
}

# now saving these again
save(Q1_AH,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_AHforMplus_withrs.Rdata"))
save(Q1_PH,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_PHforMplus_withrs.Rdata"))

# create one omnibus file for the combined model
Q <- Q1_AH %>% select(id, run, trial, run_trial, rt_csv, rt_lag, age, sex, AH_CTR_rs, AH_DMN_rs, AH_LIM_rs)
Q2 <- Q1_PH %>% select(id, run, trial, run_trial, rt_csv, rt_lag, age, sex, PH_CTR_rs, PH_DMN_rs, PH_LIM_rs)
test <- merge.data.frame(Q,Q2,by=c("id","run","trial","run_trial","rt_csv","rt_lag","age","sex"))

prepareMplusData(df = test, filename = "Study1_HC_vmPFC_clock_OMNIBUS_forMplus_taa.dat", overwrite = TRUE)

# And removing an outlier participant
Q1_AH2 <- Q1_AH[Q1_AH$id != 11162,] %>% select(!c(AH_CTR_rs,AH_DMN_rs,AH_LIM_rs))
Q1_PH2 <- Q1_PH[Q1_PH$id != 11162,] %>% select(!c(PH_CTR_rs,PH_DMN_rs,PH_LIM_rs))

save(Q1_AH2,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_sens_AHforMplus.Rdata"))
save(Q1_PH2,file=file.path(rootdir,"Mplus/Study1_HC_vmPFC_clock_sens_PHforMplus.Rdata"))

prepareMplusData(df = Q1_AH2, filename = "Study1_HC_vmPFC_clock_sens_AH_forMplus_taa.dat", overwrite = TRUE)
prepareMplusData(df = Q1_PH2, filename = "Study1_HC_vmPFC_clock_sens_PH_forMplus_taa.dat", overwrite = TRUE)

#load("Study1_HC_vmPFC_clock_sens_AHforMplus.Rdata")
#load("Study1_HC_vmPFC_clock_sens_AHforMplus.Rdata")

fileheader <- "Study1_get_HC_vmPFC_sens_randslopes_"
for(netname in c("CTR","DMN","LIM")){
  Q1_AH2 <- merge.data.frame(Q1_AH2,getslopes("AH",netname,fileheader),by="id")
  Q1_PH2 <- merge.data.frame(Q1_PH2,getslopes("PH",netname,fileheader),by="id")
}

# create one omnibus file for the combined model
Q <- Q1_AH2 %>% select(id, run, trial, run_trial, rt_csv, rt_lag, age, sex, AH_CTR_rs, AH_DMN_rs, AH_LIM_rs)
Q2 <- Q1_PH2 %>% select(id, run, trial, run_trial, rt_csv, rt_lag, age, sex, PH_CTR_rs, PH_DMN_rs, PH_LIM_rs)
test <- merge.data.frame(Q,Q2,by=c("id","run","trial","run_trial","rt_csv","rt_lag","age","sex"))
prepareMplusData(df = test, filename = "Study1_HC_vmPFC_clock_OMNIBUS_sens_forMplus_taa.dat", overwrite = TRUE)
