### 02_study1_models.R — Study 1 (MMClock) multilevel models.
###
### Study 2 lives in 03_study2_models.R, which is the corrected and gated
### version of the same analysis; the Study 2 half of the original combined
### script is archived rather than published, because it still contained the
### protocol-assignment code that misaligned trials.
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
# Ceci Westbrook, 6/2025
# code for vPFC-HC sex differences paper
# building models for analysis using Michael Hallquist's mixed_by function
# Based on code by Andrew Papale

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
library(ggplot2)
library(dplyr)
library(purrr)

# set root directory - change for your needs
rootdir <- rootdir1
repo_directory <- repo_directory

# load some source code that we'll need to use
setwd(rootdir)
source(file.path(repo_directory, "R/get_trial_data.R"))

# load mixed_by function for analyses
source(file.path(repo_directory, "R/mixed_by.R"))
source(file.path(repo_directory, "R/plotting.R"))   # supersedes plot_mixed_by_vmPFC_HC.R
source(file.path(repo_directory, "R/demographics.R"))
source(file.path(repo_directory, "analysis/04_sensitivity_extract.R"))

######################################################################################################
######################################################################################################
#########################################   STUDY 1   ################################################
######################################################################################################
######################################################################################################

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

# Compress data from 12 bins to 2 by averaging across anterior 6 bins and posterior 6 bins to create
# AH and PH averages
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>%
  summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% 
  ungroup() # 12 -> 2

# Create a new scaled within-subjects variable (HCwithin) and a between-subjects
# variable averaged per subject and run (HCbetween)
hc <- hc %>% group_by(id,run) %>%
  mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>%
  ungroup()

# Merge vmPFC and HC data; remove unscaled within-S variable (decon1)
Q <- inner_join(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')

# select and scale variables of interest
df <- df %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev))

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,run,trial,run_trial,rt_lag_sc,trial_neg_inv_sc,iti_ideal,iti_prev,iti_lag_sc,rt_csv,v_entropy_wi,v_max_wi)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# add in age and sex variables
demo <- prep_demo_study1(file.path(rootdir,"fmri/data/MMC_demog.csv"))
Q$id <- as.character(Q$id)
Q <- inner_join(Q,demo,by=c('id'))

#save(Q,file=file.path(rootdir,'mmclock_HC_vmPFC_clock.Rdata'))

s1 <- demo_summary(distinct(as.data.frame(Q), id, .keep_all = TRUE), "Study 1")
print_demo_summary(s1)
#saveRDS(s1, file.path(rootdir, "fmri/mixed_by_output/table1_study1.rds"))

########################
##### Set up models ####
########################

setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(basemodel_formula)
basemodel_formula <- NULL

# determine the base model:
# start with sex and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~sex*HCwithin + (1 | id/run)) 

# sex, hippocampal activity and control vars (trial, iti)
basemodel_formula[[2]] = formula(~sex*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[3]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[4]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[5]] = formula(~sex*HCwithin + iti_lag_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[6]] = formula(~sex*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[7]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + (1 | id/run))
basemodel_formula[[8]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run))

for(i in 1:length(basemodel_formula)){
  assign(paste0("df",i),basemodel_formula[[i]])
  print(get(paste0("df",i)))
  df0 <- get(paste0("df",i))
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0, split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  assign(paste0("ddf",i),ddf)
}

plot_aic_comparison <- function(model_list) {

  # Extract AIC values into a tidy dataframe
  all_df <- imap_dfr(model_list, function(model, name) {
    data.frame(
      Index = seq_len(nrow(model$fit_df)),
      AIC = model$fit_df$AIC,
      Model = name
    )
  })
  
  # Plot
  ggplot(all_df, aes(x = Index, y = AIC, color = Model)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(title = "AIC Values by Model",
         x = "Model Index",
         y = "AIC",
         color = "Model")
}

# can compare models 1, 2, 4, 8; 1, 3, 4, 8; 1, 5, 6, 8; 1, 5, 7, 8; 1, 3, 7, 8; 1, 5, 7, 8 (not all models are nested)

models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf4 = ddf4)
models <- list(ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf4 = ddf4, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf5 = ddf5)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf7 = ddf7)
plot_aic_comparison(models)

# No evidence that including trial or avg hippocampus activity changes model fit. Including ITI worsens it.
# therefore, for the sake of having more control variables accounted for, will include trial num and hippocampal activity.

# Now that we have determined a base model, we can do hypothesis testing:
# rt_lag_sc for RTswings, rt_lag_sc controlled for age, vmax, entropy

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of sex
decode_formula[[1]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding age
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by sex interaction
decode_formula[[3]] = formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~sex*HCwithin + race_white*HCwithin + hispanic_yn*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~sex*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~sex*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~sex*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with vmax
decode_formula[[8]] = formula(~sex*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~sex*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with entropy
decode_formula[[10]] = formula(~sex*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1  | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-sexmodels',i,'.Rdata'))
  
}

#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-sexmodels',j,'.Rdata'))

  setwd(file.path(rootdir,'fmri/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='censored-vmPFC-HC-network-clock',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
}

#########################################################################################################
#########################################################################################################

# Now, we can do the same for age!
setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(basemodel_formula)
basemodel_formula <- NULL

# determine the base model:
# start with age and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~age*HCwithin + (1 | id/run)) 
# determine the base model:
# start with sex and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~sex*HCwithin + (1 | id/run)) 

# sex, hippocampal activity and control vars (trial, iti)
basemodel_formula[[2]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[3]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[4]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[5]] = formula(~age*HCwithin + iti_lag_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[6]] = formula(~age*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[7]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + (1 | id/run))
basemodel_formula[[8]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run))

for(i in 1:length(basemodel_formula)){
  assign(paste0("df",i),basemodel_formula[[i]])
  print(get(paste0("df",i)))
  df0 <- get(paste0("df",i))
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0, split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  )
  assign(paste0("ddf",i),ddf)
}

# can compare models 1, 2, 4, 8; 1, 3, 4, 8; 1, 5, 6, 8; 1, 5, 7, 8; 1, 3, 7, 8; 1, 5, 7, 8 (not all models are nested)

models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf4 = ddf4)
models <- list(ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf4 = ddf4, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf5 = ddf5)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf7 = ddf7)
plot_aic_comparison(models)

# Once again: no evidence that including trial or avg hippocampus activity changes model fit. Including ITI worsens it.
# therefore, for the sake of having more control variables accounted for, will include trial num and hippocampal activity.

# Now that we have determined a base model, we can do hypothesis testing:
# rt_lag_sc for RTswings, rt_lag_sc controlled for age, vmax, entropy

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of age
decode_formula[[1]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding sex
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by sex interaction
decode_formula[[3]] = formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~age*HCwithin + race_white*HCwithin + hispanic_yn*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~age*HCwithin + level_edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~age*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~age*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with vmax
decode_formula[[8]] = formula(~age*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with entropy
decode_formula[[10]] = formula(~age*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~age:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-agemodels-',i,'.Rdata'))
  
}

#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-agemodels-',j,'.Rdata'))
  
  setwd(file.path(rootdir,'fmri/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='censored-',toalign='clock',
                         toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
  
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
}

#############################
#############################
#############################
#########  behavior  ########
#############################
#############################
#############################

# using mixed-by: load in fMRI and MEG data and split by dataset
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,v_max_wi,v_entropy_wi,trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;
demo <- prep_demo_study1(file.path(rootdir, "fmri/data/MMC_demog.csv"))
Qfmri <- inner_join(Q2,demo,by=c('id'))
Qfmri$age <- scale(Qfmri$age)

df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_meg')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,v_max_wi,v_entropy_wi, trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;

Qmeg <- inner_join(Q2,demo,by=c('id'))
#Qmeg$female <- relevel(as.factor(Qmeg$female),ref='0')
Qmeg$age <- scale(Qmeg$age)

Qmeg <- Qmeg %>% mutate(dataset = 'MEG')
Qfmri <- Qfmri %>% mutate(dataset = 'fMRI')
Q3 <- rbind(Qmeg,Qfmri)
#Q3$race <- Q3 %>% select("AmerIndianAlaskan":"White") %>% rowSums()

Q3 <- Q3 %>% filter(rt_csv < 4 & rt_csv > 0.2)
Q3 <- Q3 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

modelfits <- read.csv(file.path(rootdir,'trial_data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv'))
modelfits$id <- as.character(modelfits$id)
Q3 <- inner_join(Q3,modelfits[,c("id","beta")],by="id")

# Save Q3 so the figure scripts can run without rebuilding it. Figure 2 and the
# Study 1 panels of Figure S6 both need it, and it is otherwise session-only.
if (!dir.exists(mbo1)) dir.create(mbo1, recursive = TRUE)
saveRDS(Q3, file.path(mbo1, paste0(curr_date, "_Study1_Q3.rds")))


# look at total earnings by sex
earnings_fmri <- Q3[Q3$dataset=="fMRI",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_MEG <- Q3[Q3$dataset=="MEG",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_fmri$dataset <- "fMRI"
earnings_MEG$dataset <- "MEG"
earnings_all <- rbind(earnings_fmri,earnings_MEG)

# not significant
t.test(earnings_fmri[earnings_fmri$sex==0,"total_earnings"],earnings_fmri[earnings_fmri$sex==1,"total_earnings"])
t.test(earnings_MEG[earnings_MEG$sex==0,"total_earnings"],earnings_MEG[earnings_MEG$sex==1,"total_earnings"])

# plot:
tab <- earnings_all %>% group_by(sex,dataset) %>% summarise(TotalEarnings = mean(total_earnings),se = se(total_earnings))
tab$sex <- as.factor(tab$sex)

g <- ggplot(tab,aes(y=TotalEarnings,x=sex))
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex),size=3) +
  geom_errorbar(aes(ymin = TotalEarnings - se,ymax = TotalEarnings + se,color=sex),
                linewidth=2,position=position_dodge(width=0.5)) +
  facet_wrap(~dataset) +
  ylab("Total Earnings") + xlab("Dataset") + 
  scale_x_discrete(labels = c("Female","Male")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))

# set up mixed_by models
setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
rm(decode_formula)
decode_formula <- NULL

# sex only
decode_formula[[1]] <- formula(~rt_lag_sc*sex + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# age by sex interaction
decode_formula[[4]] = formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[5]] = formula(~rt_lag_sc*sex + rt_lag_sc*race_white + rt_lag_sc*hispanic_yn + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# highest level of education
decode_formula[[6]] = formula(~rt_lag_sc*sex + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~rt_lag_sc*sex + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex with vmax
decode_formula[[8]] = formula(~rt_lag_sc*sex*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~rt_lag_sc*sex + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex with entropy
decode_formula[[10]] = formula(~rt_lag_sc*sex*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q3, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = "dataset",return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2)))
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Sex-clock-Study1-fmri-meg-pred-rt-',j,'.Rdata'))
}

# plot emtrends for sex
load(paste0(curr_date,'-Sex-clock-Study1-fmri-meg-pred-rt-4.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$sex <- as.factor(tab$sex)

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex))
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex),size=3) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,
                    ymax = rt_lag_sc.trend + SE,
                    color=sex),linewidth=2,position=position_dodge(width=0.5)) +
  facet_wrap(~dataset) +
  scale_x_discrete(labels = c("Female","Male")) + scale_y_reverse() +
  ylab("RT Swing (RT ~ previous RT)") + xlab("Dataset") + 
  theme(legend.position = "none",
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))

# look at differences by last_outcome
rm(decode_formula)
decode_formula <- NULL

# sex only
decode_formula[[1]] <- formula(~rt_lag_sc*sex*reward_lag_rec + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*age*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*race_white*reward_lag_rec + rt_lag_sc*hispanic_yn*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*edu*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q3, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = "dataset",return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc=c(-2,-1,0,1,2)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Sex-clock-Study1-fmri-meg-pred-rt-lastoutcome-',j,'.Rdata'))
}

load(paste0(curr_date,'-Sex-clock-Study1-fmri-meg-pred-rt-lastoutcome-5.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

# by sex and last outcome: won't work for models without outcome
tab <- as.data.frame(ddq$emtrends_list$RTxO) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset,reward_lag_rec) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$last_outcome <- ifelse(tab$reward_lag_rec==-0.5, "Omission","Reward")
tab$sex1 <- ifelse(tab$sex==0,"Female","Male")

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex1))
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex1),size=2) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,
                    ymax = rt_lag_sc.trend + SE,
                    color=sex1),linewidth=2,position=position_dodge(width=0.5)) +
  facet_grid(~dataset*last_outcome, labeller = as_labeller(c("fMRI" = "fMRI", "MEG" = "MEG",
                                                         "Omission" = "Omission", "Reward" = "Reward"),)) +
  ylab("RT Swing (RT ~ previous RT)") + xlab("Sex") + scale_y_reverse() +
  theme(legend.position = "none",
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16))


g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex1),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~last_outcome) + scale_y_reverse()

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=last_outcome),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~sex1) + 
  ylab("RT Swing (RT ~ previous RT)") + xlab("Prior Trial Outcome")

#getting out sensitivity details:
extract_sensitivity(file.path(rootdir, "fmri/mixed_by_output/"), "Study1",
                    file.path(rootdir, "fmri/mixed_by_output/sensitivity_study1.rds"))

#########################################################################################################
#########################################################################################################

# Now, we can do the same for age!
setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
rm(decode_formula)
decode_formula <- NULL

# age only
decode_formula[[1]] <- formula(~rt_lag_sc*age + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run))

# now adding sex
decode_formula[[3]] = formula(~rt_lag_sc*age + rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# age by sex interaction
decode_formula[[4]] = formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[5]] = formula(~rt_lag_sc*age + rt_lag_sc*race_white + rt_lag_sc*hispanic_yn + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# highest level of education
decode_formula[[6]] = formula(~rt_lag_sc*age + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~rt_lag_sc*age + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of age with vmax
decode_formula[[8]] = formula(~rt_lag_sc*age*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~rt_lag_sc*age + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of age with entropy
decode_formula[[10]] = formula(~rt_lag_sc*age*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q3, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = "dataset",return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:age), at = list(age=c(-1.5,1.5),rt_lag_sc=c(-2,-1,0,1,2)))
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:age), at=list(age=c(-1.5,1.5),rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-clock-Study1-fmri-meg-pred-rt-',j,'.Rdata'))
}

# plot emtrends for age
load(paste0(curr_date,'-Age-Sex-clock-Study1-fmri-meg-pred-rt-10.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(age,dataset) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=age),group=dataset)
g + geom_point(stat="identity",aes(color=dataset),size=3) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,
                    ymax = rt_lag_sc.trend + SE,
                    color=dataset),size=2) +
  geom_line(aes(color=dataset),linewidth=2) +
  scale_y_reverse() +
  ylab("RT Swing (RT ~ previous RT)") + xlab("Age (scaled)") +
  labs(color = "Dataset") +
  theme(
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12))

#########################################################################################################
#########################################################################################################
#########################################################################################################

# Softmax inverse temperature (beta)
#
# beta is one value per participant, so these are ordinary regressions on a
# participant-level frame -- no mixed models needed.
#
# NB two things about the earlier version of this block. It also regressed beta
# on v_max_wi and v_entropy_wi, which are trial-level quantities, using the
# trial-level frame: that repeats each participant's beta across hundreds of
# rows, inflating n from 71 to tens of thousands and rendering the standard
# errors meaningless. Those models are removed rather than corrected, because
# the question they addressed (does beta track within-trial value?) is not one a
# participant-level constant can answer. The age model is kept but computed on
# the participant-level frame, which the earlier version did not do.

bsub1 <- Q3 %>% filter(dataset == "fMRI") %>%
  distinct(id, sex, age, race_white, hispanic_yn, edu, beta)
cat("beta models, Study 1 fMRI: n =", nrow(bsub1), "participants\n")

# reported in the Supplement: no relationship with sex or age
summary(lm(beta ~ sex, data = bsub1))
summary(lm(beta ~ age, data = bsub1))

# covariate-adjusted, also reported
summary(lm(beta ~ sex + age, data = bsub1))
summary(lm(beta ~ sex * age, data = bsub1))
summary(lm(beta ~ sex + race_white + hispanic_yn, data = bsub1))
summary(lm(beta ~ sex + edu, data = bsub1))
