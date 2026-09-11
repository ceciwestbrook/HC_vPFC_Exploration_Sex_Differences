### 01b_build_mplus_data_study2.R — builds the Mplus .dat files for Study 2.
###
### This is the corrected version: protocol assignment (split_ksoc_bsoc / maxT)
### has been removed, so the vmPFC-hippocampus join no longer requires a trial
### count of exactly 240 or 300 and no longer silently drops participants with
### partial data. Yields 185 participants rather than 179.
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

######################################################################################################
######################################################################################################
#########################################   STUDY 2   ################################################
######################################################################################################
######################################################################################################

# set root directory - change for your needs
rootdir <- rootdir1
repo_directory <- file.path(rootdir,'BSOC')

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
vmPFC <- read_csv(file.path(repo_directory,'clock_aligned_bsocial_vmPFC.csv.gz')) %>%
  mutate(network = case_when(atlas_value %in% c(55, 56, 159, 160) ~ 'LIM',
                             atlas_value %in% c(65, 66, 67, 170, 171) ~ 'CTR',
                             atlas_value %in% c(84, 86, 88, 89, 161, 191, 192, 194) ~ 'DMN'),
         run = case_when(run == 'run1' ~ 1, run == 'run2' ~ 2),
         id = as.character(id)) %>%
  filter(evt_time > -5 & evt_time < 5) %>%
  select(id, run, trial, vmPFC_decon = decon_mean, evt_time, network)

# now load in HC data
load(file.path(repo_directory, 'BSOC_HC_clock_TRdiv2.Rdata'))   # loads `hc`

# Compress 12 bins to 2 by averaging the anterior 6 and posterior 6
hc <- hc %>%
  filter(evt_time > -5 & evt_time < 5) %>%
  group_by(id, run, trial, evt_time, HC_region) %>%
  summarize(HCwithin = mean(decon_mean, na.rm = TRUE), .groups = 'drop')

# Merge vmPFC and HC data
Q <- inner_join(vmPFC, hc, by = c("id", "run", "trial", "evt_time"))

#rm(hc,vmPFC)

# Get task behav data
behav <- read_csv(file.path(repo_directory, 'bsocial_clock_trial_df.csv')) %>%
  select(id, scanner_run, trial, iti_prev, iti_ideal, rt_csv, rt_lag) %>%
  rename(run = scanner_run)

Q$id <- as.numeric(Q$id)
Q <- inner_join(behav, Q, by = c("id", "run", "trial"))
n_distinct(Q$id)     # still expect 194

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;

Q1 <- Q %>% group_by(id,run,trial,network,HC_region, rt_csv, rt_lag) %>% 
  summarize(vmPFC_decon = mean(vmPFC_decon,na.rm=TRUE),HCwithin = mean(HCwithin,na.rm=TRUE)) %>% 
  ungroup()
Q1 <- Q1 %>% group_by(id,run,trial) %>%
  pivot_wider(values_from=c(vmPFC_decon,HCwithin),names_from = 'network') %>% 
  ungroup()

# add in age and sex variables
demo <- read_csv(file.path(repo_directory,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory,'2025-02-27-Partial-demo-pull-KSOC.csv'))
demo1$id <- demo1$registration_redcapid
demo <- demo %>% rename(sex=registration_birthsex,
                        gender=registration_gender,
                        group=registration_group) %>%
  select(id,group,age,sex)
demo1 <- demo1 %>% rename(sex=registration_birthsex,
                          gender=registration_gender,
                          group=registration_group) %>%
  select(id,group,age,sex)
demo2 <- rbind(demo,demo1)
Q1 <- inner_join(Q1,demo2,by=c('id'))
Q1$female <- ifelse(Q1$sex==1,1,0)
Q1 <- Q1 %>% select(!sex)
# recoding to make female the reference (0) group
Q1$sex <- ifelse(Q1$female==1,0,
                ifelse(Q1$female==0,1,NA))
Q1 <- Q1 %>% select(!female)
# Exclusions BEFORE scaling age. The original scaled first, which centred age on
# a sample that included the seven participants over 50 who are then removed.
# Harmless for these models, since age is not in USEVARIABLES, but wrong order.
Q1 <- apply_exclusions(Q1, "study2")   # participants over 50, and one failed run

Q1$age <- scale(Q1$age)

Q1_AH <- Q1 %>% filter(HC_region == "AH")
Q1_PH <- Q1 %>% filter(HC_region == "PH")

setwd(file.path(rootdir,"Mplus"))

save(Q1_AH,file=file.path(rootdir,'Mplus/Study2_HC_vmPFC_clock_AHforMplus.Rdata'))
save(Q1_PH,file=file.path(rootdir,'Mplus/Study2_HC_vmPFC_clock_PHforMplus.Rdata'))

# now to format for MPlus.
# for loading datasets that have already been prepared:
#load("Study2_HC_vmPFC_clock_AHforMplus.Rdata")
#load("Study2_HC_vmPFC_clock_PHforMplus.Rdata")

prepareMplusData(df = Q1_AH, filename = "Study2_HC_vmPFC_clock_AH_forMplus_taa.dat", overwrite = TRUE)
prepareMplusData(df = Q1_PH, filename = "Study2_HC_vmPFC_clock_PH_forMplus_taa.dat", overwrite = TRUE)

# now, after the Mplus models have been run to generate random slopes for vPFC-HC:

# we ran six models: 3 networks (DMN, CTR, LIM) x 2 HC regions (AH, PH).
fileheader <- "Study2_get_HC_vmPFC_randslopes_"
for(netname in c("CTR","DMN","LIM")){
  Q1_AH <- merge.data.frame(Q1_AH,getslopes("AH",netname,fileheader),by="id")
  Q1_PH <- merge.data.frame(Q1_PH,getslopes("PH",netname,fileheader),by="id")
}

# now saving these again
save(Q1_AH,file=file.path(rootdir,"Mplus/Study2_HC_vmPFC_clock_AHforMplus_withrs.Rdata"))
save(Q1_PH,file=file.path(rootdir,"Mplus/Study2_HC_vmPFC_clock_PHforMplus_withrs.Rdata"))

# create one omnibus file for the combined model
Q <- Q1_AH %>% select(id, run, trial, rt_csv, rt_lag, age, sex, AH_CTR_rs, AH_DMN_rs, AH_LIM_rs)
Q2 <- Q1_PH %>% select(id, run, trial, rt_csv, rt_lag, age, sex, PH_CTR_rs, PH_DMN_rs, PH_LIM_rs)
test <- merge.data.frame(Q,Q2,by=c("id","run","trial","rt_csv","rt_lag","age","sex"))

prepareMplusData(df = test, filename = "Study2_HC_vmPFC_clock_OMNIBUS_forMplus_taa.dat", overwrite = TRUE)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# extract new slopes including group and correlate with original slopes to look for signs of moderation
# by group

setwd(file.path(rootdir,"Mplus"))

# for loading datasets that have already been prepared:
load(file.path(rootdir,"Mplus/Study2_HC_vmPFC_clock_AHforMplus_withrs.Rdata"))
load(file.path(rootdir,"Mplus/Study2_HC_vmPFC_clock_PHforMplus_withrs.Rdata"))

Qc <- Q1[!(Q1$group %in% c("DNA","DEP")), ]
Qc$bpd <- factor(ifelse(Qc$group == "HC","HC","BPD"), levels = c("HC","BPD"))

# Qc already has bpd (HC/BPD) and already excludes DNA/DEP.
# Make everything numeric for Mplus (no strings allowed in .dat).
Qc2 <- Qc %>%
  mutate(
    id      = as.numeric(as.character(id)),
    sex_num = as.numeric(as.character(sex)),      # your coding: female=0, male=1
    bpd_num = ifelse(bpd == "BPD", 1, 0)          # HC=0, BPD=1
  )

# fixed column order; names are positional in Mplus so labels are arbitrary
mplus_vars <- c("id","run","trial","rt_csv","rt_lag",
                "vmPFC_decon_CTR","vmPFC_decon_DMN","vmPFC_decon_LIM",
                "HCwithin_CTR","HCwithin_DMN","HCwithin_LIM",
                "age","sex_num","bpd_num")

prepareMplusData(
  df = Qc2 %>% filter(HC_region == "AH") %>% select(all_of(mplus_vars)),
  filename = "Study2_vmPFC_clock_AH_groupadj.dat", overwrite = TRUE
)
prepareMplusData(
  df = Qc2 %>% filter(HC_region == "PH") %>% select(all_of(mplus_vars)),
  filename = "Study2_vmPFC_clock_PH_groupadj.dat", overwrite = TRUE
)

# create one omnibus file for the combined model
Q <- Q1_AH %>% select(id, run, trial, rt_csv, rt_lag, age, sex, AH_CTR_rs, AH_DMN_rs, AH_LIM_rs)
Q2 <- Q1_PH %>% select(id, run, trial, rt_csv, rt_lag, age, sex, PH_CTR_rs, PH_DMN_rs, PH_LIM_rs)
test <- merge.data.frame(Q,Q2,by=c("id","run","trial","rt_csv","rt_lag","age","sex"))

library(dplyr)
library(purrr)

# the six group-adjusted FSCORES output files
cells <- c("AH_CTR","AH_DMN","AH_LIM","PH_CTR","PH_DMN","PH_LIM")
fs_files <- setNames(
  paste0("Study2_get_HC_vmPFC_randslopes_adjustbygroup_", cells, ".dat"),
  cells
)

# readModels parses the .out + savedata and labels columns for you;
# point it at the .out files (same basename as your .inp)
out_files <- setNames(paste0("Study2_get_HC_vmPFC_randslopes_adjustbygroup_", cells, ".out"), cells)

read_rs <- function(outfile, cell) {
  m  <- readModels(outfile)
  sd <- m$savedata                      # data.frame with named columns
  # the factor-score mean column for the random slope is named "RS" (Mplus uppercases)
  rs_col <- grep("^RS$", names(sd), value = TRUE)          # exact "RS"
  if (length(rs_col) == 0) rs_col <- grep("^RS\\.Mean$|^RS_", names(sd), value = TRUE)[1]
  tibble(id = sd$ID, !!paste0(cell, "_new") := sd[[rs_col]]) %>%
    distinct(id, .keep_all = TRUE)      # one row per subject (FSCORES repeats per row)
}

new_slopes <- imap(out_files, read_rs) %>% reduce(full_join, by = "id")

# --- your ORIGINAL slopes (the ones in the Mplus MSEM file) ---
# replace `orig` with your object; needs id + the six *_rs columns
orig <- test %>% distinct(id, .keep_all = TRUE) %>%
  select(id, AH_CTR_rs, AH_DMN_rs, AH_LIM_rs, PH_CTR_rs, PH_DMN_rs, PH_LIM_rs)

new_slopes$id <- as.numeric(new_slopes$id)
orig$id       <- as.numeric(orig$id)
merged <- inner_join(orig, new_slopes, by = "id")

cors <- sapply(cells, function(cn)
  cor(merged[[paste0(cn, "_rs")]], merged[[paste0(cn, "_new")]], use = "complete.obs"))

data.frame(cell = cells, r = round(cors, 4), n = nrow(merged))

# for the covariate-adjusted group model:
# add person-level bpd to `test` (0 = HC, 1 = BPD), dropping DNA/DEP via the join
bpd_lookup <- Qc %>%
  distinct(id, bpd) %>%
  mutate(id = as.numeric(as.character(id)),
         bpd_num = ifelse(bpd == "BPD", 1, 0)) %>%
  select(id, bpd_num)

test2 <- test %>%
  mutate(id = as.numeric(as.character(id))) %>%
  inner_join(bpd_lookup, by = "id")

prepareMplusData(df = test2,
                 filename = "Study2_HC_vmPFC_clock_OMNIBUS_groupadj_forMplus_taa.dat",
                 overwrite = TRUE)
