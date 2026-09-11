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
### 03_study2_models.R — full Study 2 re-run, unattended
#
# Fixes the run_trial derivation, verifies trial alignment, then runs every
# Study 2 model that feeds the manuscript or supplement. Designed for overnight
# use: each model is wrapped so one failure cannot kill the run, completed
# outputs are skipped on restart, and progress goes to a timestamped log.
#
# Usage:  Rscript 03_study2_models.R           (or via the sbatch wrapper)
#         FORCE=1 Rscript 03_study2_models.R   to re-run models already on disk
#         VERIFY_ONLY=1 ...                  gates and diagnostics only
#         SAVE_DATA=1 ...                    also write Q_fmri (large)
#
# Requires, in rootdir1 or the repo directory (see source_first below):
#   mixed_by.R      model fitting
#   demographics.R  prep_demo_study2(), demo_summary(), study2_excluded_ids
#   plotting.R      subject_slopes(), for the Figure 2 input at stage 4
# Writes Mplus data to rootdir1/Mplus, where both studies' Mplus files live.
#
# UNTESTED against your data — the verification gate at stage 1 is there so a
# bad build stops the run instead of feeding 12 hours of compute.

## ---- library path -------------------------------------------------------
# Batch jobs do not always inherit the interactive session's library path.
# Prepend it explicitly so `library(fmri.pipeline)` resolves either way.
# Must match the running R version: this tree was built for R 4.6.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(lme4); library(lmerTest)
  library(emmeans); library(fmri.pipeline); library(parallel); library(doParallel)
})

## ---- paths and options ----------------------------------------------------
# rootdir1 set in config.R
# rootdir2 set in config.R          # note: fmri lowercase in S1, fMRI upper in S2
outdir <- file.path(rootdir2, 'fMRI/mixed_by_output')
logdir <- file.path(rootdir2, 'fMRI/rerun_logs')
dir.create(logdir, showWarnings = FALSE, recursive = TRUE)

# Sourced helpers live in rootdir1; a repo checkout is also searched in case
# the layout differs. Fails loudly rather than silently missing a function.
source_dirs <- c(rootdir1, file.path(rootdir1, 'HC_vPFC_repo/vPFC-HC-Clock'))
source_first <- function(fname) {
  hits <- file.path(source_dirs, fname)
  hit <- hits[file.exists(hits)]
  if (!length(hit)) stop('cannot find ', fname, ' in: ', paste(source_dirs, collapse = ', '))
  message('sourcing ', hit[1]); source(hit[1])
}
source_first('mixed_by.R')
source_first('demographics.R')

ncores   <- as.integer(Sys.getenv("NCORES", "20"))
FORCE    <- nzchar(Sys.getenv("FORCE"))
# VERIFY_ONLY=1 builds the data, runs both alignment gates, prints the Ns and
# demographics, saves the data objects, and stops before any models.
VERIFY_ONLY <- nzchar(Sys.getenv("VERIFY_ONLY"))
run_tag  <- strftime(Sys.time(), format = '%Y-%m-%d')
logfile  <- file.path(logdir, paste0('study2_rerun_', strftime(Sys.time(), '%Y%m%d_%H%M'), '.log'))

log_msg <- function(...) {
  line <- paste0('[', strftime(Sys.time(), '%H:%M:%S'), '] ', paste0(..., collapse = ''))
  cat(line, '\n', sep = ''); cat(line, '\n', sep = '', file = logfile, append = TRUE)
}
timed <- function(label, expr) {
  t0 <- Sys.time()
  out <- tryCatch(expr, error = function(e) { log_msg('FAILED  ', label, ' :: ', conditionMessage(e)); NULL })
  log_msg(if (is.null(out)) 'failed  ' else 'done    ', label,
          '  (', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min)')
  out
}

log_msg('=== Study 2 re-run starting | ncores=', ncores, ' | FORCE=', FORCE, ' ===')
log_msg('log: ', logfile)

## ===========================================================================
## STAGE 1 — build data with the run_trial fix, then verify alignment
## ===========================================================================

df <- read_csv(file.path(rootdir2, 'bsocial_clock_trial_df.csv'), show_col_types = FALSE)

split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300]
bsoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240]
log_msg('protocols: ', length(bsoc), ' bsocial (240 trials), ', length(ksoc), ' ksocial (300 trials)')

# --- verbatim from the original script ---
df_bsoc <- df %>% filter(id %in% bsoc) %>% mutate(block = case_when(trial <= 40 ~ 1,
                                                                    trial > 40 & trial <= 80 ~ 2,
                                                                    trial > 80 & trial <=120 ~ 3,
                                                                    trial > 120 & trial <=160 ~ 4,
                                                                    trial > 160 & trial <=200 ~ 5,
                                                                    trial > 200 & trial <=240 ~ 6))
df_bsoc <- df_bsoc %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial,
                                                     trial > 40 & trial <= 80 ~ trial-40,
                                                     trial > 80 & trial <=120 ~ trial-80,
                                                     trial > 120 & trial <=160 ~ trial-120,
                                                     trial > 160 & trial <=200 ~ trial-160,
                                                     trial > 200 & trial <=240 ~ trial-200))
df_bsoc <- df_bsoc %>% mutate(protocol = 'bsocial',
                              run_trial0_c = run_trial0-floor(run_trial0/40.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df_ksoc <- df %>% filter(id %in% ksoc) %>% mutate(block = case_when(trial <= 50 ~ 1,
                                                                    trial > 50 & trial <= 100 ~ 2,
                                                                    trial > 100 & trial <=150 ~ 3,
                                                                    trial > 150 & trial <=200 ~ 4,
                                                                    trial > 200 & trial <=250 ~ 5,
                                                                    trial > 250 & trial <=300 ~ 6))
df_ksoc <- df_ksoc %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial,
                                                     trial > 50 & trial <= 100 ~ trial-50,
                                                     trial > 100 & trial <=150 ~ trial-100,
                                                     trial > 150 & trial <=200 ~ trial-150,
                                                     trial > 200 & trial <=250 ~ trial-200,
                                                     trial > 250 & trial <=300 ~ trial-250))
df_ksoc <- df_ksoc %>% mutate(protocol = 'ksocial',
                              run_trial0_c = run_trial0-floor(run_trial0/50.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df <- rbind(df_bsoc,df_ksoc)

df <- df %>%
  group_by(id, scanner_run) %>%
  mutate(id = as.character(id),
         v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev),
         v_max_sc = scale(v_max),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_swing_sc = scale(rt_swing)) %>% ungroup()

behav <- df %>% select(id,scanner_run,trial,run_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag,rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,last_outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,total_earnings)
# --- end verbatim ---

# ---- THE FIX -------------------------------------------------------------
# Was: behav %>% rename(run = scanner_run) %>% select(!run_trial) %>%
#        mutate(run_trial = case_when(trial <= 150 ~ trial, trial > 150 ~ trial - 150))
# That hard-codes 150-trial scanner runs: correct for ksocial, wrong for bsocial's
# 120. The source `run_trial` still has to be dropped, because it is position
# within a 40/50-trial contingency block rather than within a scanner run — but it
# is now recomputed from each run's own minimum, which is protocol-agnostic and
# matches how the imaging side derives the same variable.
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>%
  group_by(id, run) %>% mutate(run_trial = trial - min(trial) + 1) %>% ungroup()
# --------------------------------------------------------------------------

demo2 <- prep_demo_study2(file.path(rootdir2, 'bsoc_clock_N171_dfx_demoonly.csv'),
                          file.path(rootdir2, '2025-02-27-Partial-demo-pull-KSOC.csv'))

## ---- imaging side ----
vmPFC <- read_csv(file.path(rootdir2, 'clock_aligned_bsocial_vmPFC.csv.gz'), show_col_types = FALSE)
vmPFC <- vmPFC %>% mutate(run1 = case_when(run == 'run1' ~ 1, run == 'run2' ~ 2)) %>%
  select(!run) %>% rename(run = run1) %>%
  group_by(id, run) %>% mutate(run_trial = trial - min(trial) + 1) %>% ungroup() %>%
  mutate(network = case_when(atlas_value %in% c(55,56,159,160) ~ 'LIM',
                             atlas_value %in% c(65,66,67,170,171) ~ 'CTR',
                             atlas_value %in% c(84,86,88,89,161,191,192,194) ~ 'DMN'),
         symmetry_group = case_when(atlas_value %in% c(65,66,170) ~ 1,
                                    atlas_value %in% c(86,161) ~ 2, atlas_value %in% c(56,159) ~ 3,
                                    atlas_value %in% c(84,191) ~ 4, atlas_value %in% c(88,192) ~ 5,
                                    atlas_value %in% c(67,171) ~ 6, atlas_value %in% c(89,194) ~ 7,
                                    atlas_value %in% c(55,160) ~ 8),
         id = as.character(id)) %>%
  filter(evt_time > -5, evt_time < 5) %>%
  # keep `trial` so it can serve as a redundant join key
  select(id, run, trial, run_trial, vmPFC_decon = decon_mean, atlas_value,
         evt_time, symmetry_group, network)

## ---- VERIFICATION GATE ----
log_msg('--- verifying trial alignment ---')
align <- inner_join(behav %>% distinct(id, run, run_trial, trial_behav = trial),
                    vmPFC %>% distinct(id, run, run_trial, trial_img = trial),
                    by = c('id','run','run_trial')) %>%
  mutate(offset = trial_behav - trial_img)
offs <- align %>% count(run, offset)
print(offs); capture.output(print(offs), file = logfile, append = TRUE)

if (any(align$offset != 0)) {
  log_msg('ABORT: non-zero trial offsets remain — alignment not fixed. Nothing was run.')
  quit(status = 1)
}
tr <- behav %>% distinct(id, run, run_trial) %>% count(id, name = 'n_trials') %>%
  left_join(split_ksoc_bsoc %>% mutate(id = as.character(id)), by = 'id') %>%
  group_by(protocol = ifelse(maxT == 300, 'ksocial', 'bsocial')) %>%
  summarize(n = n(), median_trials = median(n_trials), .groups = 'drop')
print(tr); capture.output(print(tr), file = logfile, append = TRUE)
log_msg('alignment OK (all offsets zero)')

## ---- hippocampus, merge, censor, exclusions ----
load(file.path(rootdir2, 'BSOC_HC_clock_TRdiv2.Rdata'))   # object: hc

## ---- SECOND VERIFICATION GATE: is hc's run_trial run position or block position? ----
# The behavioural source file's run_trial turned out to be position within a
# 40/50-trial contingency block, not within a scanner run. hc's run_trial is
# never recomputed, so if it is also block position, the vmPFC-hc join is
# misaligned too. Run position maxes at 120/150; block position at 40/50.
log_msg('--- verifying hippocampus run_trial ---')
hc_rng <- hc %>% group_by(run) %>%
  summarize(rt_min = min(run_trial, na.rm = TRUE), rt_max = max(run_trial, na.rm = TRUE),
            n_ids = n_distinct(id), .groups = 'drop')
print(hc_rng); capture.output(print(hc_rng), file = logfile, append = TRUE)

if (max(hc_rng$rt_max, na.rm = TRUE) <= 50) {
  log_msg('ABORT: hippocampus run_trial maxes at ', max(hc_rng$rt_max),
          ' — looks like block position, not run position. The vmPFC-hc join would be',
          ' misaligned. Nothing was run.')
  quit(status = 1)
}

# If hc carries `trial`, check the pairing directly rather than inferring from ranges.
if ('trial' %in% names(hc)) {
  hc_off <- inner_join(hc %>% distinct(id, run, run_trial, trial_hc = trial) %>%
                         mutate(id = as.character(id)),
                       vmPFC %>% distinct(id, run, run_trial, trial_img = trial),
                       by = c('id','run','run_trial')) %>%
    mutate(offset = trial_hc - trial_img) %>% count(run, offset)
  print(hc_off); capture.output(print(hc_off), file = logfile, append = TRUE)
  if (any(hc_off$offset != 0)) {
    log_msg('ABORT: vmPFC-hc trial offsets are non-zero. Nothing was run.')
    quit(status = 1)
  }
  log_msg('vmPFC-hc alignment OK')
} else {
  log_msg('hc has no `trial` column — range check only; treat as provisional')
}

hc <- hc %>% filter(evt_time > -5, evt_time < 5) %>%
  group_by(id, run, run_trial, evt_time, HC_region) %>%
  summarize(decon1 = mean(decon_mean, na.rm = TRUE), .groups = 'drop') %>%
  group_by(id, run) %>%
  mutate(HCwithin = scale(decon1), HCbetween = mean(decon1, na.rm = TRUE)) %>%
  ungroup() %>% mutate(id = as.character(id)) %>% select(!decon1)

Q_fmri <- inner_join(vmPFC, hc, by = c('id','run','run_trial','evt_time')) %>%
  inner_join(behav, by = c('id','run','trial','run_trial')) %>%     # `trial` as guard key
  inner_join(demo2, by = 'id') %>%
  apply_exclusions('study2')   # participant AND single-run exclusions

Q_fmri <- Q_fmri %>%
  mutate(across(c(vmPFC_decon, HCwithin, HCbetween),
                ~ replace(.x, evt_time > rt_csv + iti_ideal | evt_time < -iti_prev, NA)),
         age = as.vector(scale(age)))

Q_behav <- behav %>% inner_join(demo2, by = 'id') %>%
  apply_exclusions('study2') %>%
  filter(rt_csv < 4, rt_csv > 0.2) %>%
  mutate(age = as.vector(scale(age)),
         reward_lag_rec = case_when(last_outcome == 'Reward' ~ 0.5,
                                    last_outcome == 'Omission' ~ -0.5))
# Softmax beta: the ffx/fixedparams global_statistics file returns the GROUP
# estimate, constant across participants (2.602 for all 187), so regressions on
# it fit only floating-point noise. The mfx file carries per-subject estimates.
# NB: mfx values are empirical-Bayes, shrunk toward the group mean, so
# between-subject variance is attenuated — state this in STAR Methods.
modelfits <- read.csv(file.path(rootdir2,
                                'fMRIEmoClock_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv'))
modelfits$id <- sub('_1$', '', as.character(modelfits$id))
log_msg('beta fits: n = ', nrow(modelfits), ' | unique = ', length(unique(modelfits$beta)),
        ' | sd = ', signif(sd(modelfits$beta, na.rm = TRUE), 4))
if (length(unique(modelfits$beta)) < 5)
  log_msg('WARNING: beta looks constant — check that this is the mfx file')
n_before <- n_distinct(Q_behav$id)
# left_join, not inner: beta appears in no behavioural formula, so the fits file
# must not determine who is in the RT models. Participants without a fit get NA
# and drop only from the beta regressions.
Q_behav <- left_join(Q_behav, modelfits[, c('id','beta')], by = 'id')
n_nobeta <- Q_behav %>% distinct(id, beta) %>% filter(is.na(beta)) %>% nrow()
log_msg('behavioural N = ', n_before, ' (unchanged by the beta join); ',
        n_nobeta, ' participants without a beta estimate')

log_msg('N behavioral = ', n_distinct(Q_behav$id), ' | N fMRI = ', n_distinct(Q_fmri$id))

## ---- verification printout (before any save, so a crash cannot cost it) ----
log_msg('--- usable trials per participant, by protocol ---')
trial_chk <- Q_fmri %>% distinct(id, run, run_trial) %>% count(id, name = 'n_trials') %>%
  left_join(split_ksoc_bsoc %>% mutate(id = as.character(id)), by = 'id') %>%
  group_by(protocol = ifelse(maxT == 300, 'ksocial', 'bsocial')) %>%
  summarize(n = n(), median_trials = median(n_trials),
            min = min(n_trials), max = max(n_trials), .groups = 'drop')
print(trial_chk); capture.output(print(trial_chk), file = logfile, append = TRUE)
log_msg('expected after the fix: bsocial 240, ksocial 300')

gxs_f <- xtabs(~ sex1 + group, data = distinct(as.data.frame(Q_fmri), id, .keep_all = TRUE))
gxs_b <- xtabs(~ sex1 + group, data = distinct(as.data.frame(Q_behav), id, .keep_all = TRUE))
rw    <- distinct(as.data.frame(Q_behav), id, .keep_all = TRUE) %>% count(race_white)
print(gxs_f); print(gxs_b); print(rw)
capture.output(print(gxs_f), print(gxs_b), print(rw), file = logfile, append = TRUE)

## ---- save --------------------------------------------------------------
# Q_behav is small (trial-level for 187 participants) and is needed to build
# Figure 2's violins, so it is always written. Q_fmri is very large and its
# save needs a big contiguous allocation, so it stays behind SAVE_DATA=1.
timed('save Q_behav', {
  saveRDS(Q_behav, file.path(outdir, paste0(run_tag, '_Study2_Qbehav_fixed.rds')),
          compress = FALSE); TRUE })
if (nzchar(Sys.getenv('SAVE_DATA'))) {
  timed('save Q_fmri', {
    saveRDS(Q_fmri, file.path(outdir, paste0(run_tag, '_Study2_Qfmri_fixed.rds')),
            compress = FALSE); TRUE })
} else {
  log_msg('SAVE_DATA not set — Q_fmri not written (set SAVE_DATA=1 to keep it)')
}

if (VERIFY_ONLY) {
  log_msg('=== VERIFY_ONLY complete — gates passed, no models run ===')
  quit(status = 0, save = 'no')
}

## ===========================================================================
## STAGE 2 — helper for running and skipping
## ===========================================================================
splits <- c('evt_time','network','HC_region')

run_mb <- function(tag, i, data, formula, outcome, split_on = NULL,
                   emm = NULL, emt = NULL, keep_models = FALSE) {
  f <- file.path(outdir, paste0(run_tag, '_', tag, i, '.Rdata'))
  done <- list.files(outdir, pattern = paste0(tag, i, '\\.Rdata$'), full.names = TRUE)
  if (!FORCE && length(done) && any(file.info(done)$mtime > Sys.Date() - 1)) {
    log_msg('skip    ', tag, i, ' (already present today)'); return(invisible(NULL))
  }
  timed(paste0(tag, i), {
    # return_models = TRUE keeps every fitted lmer object: with 60 split cells and
    # 20 workers that exhausts memory and writes very large files. Downstream code
    # only uses the coefficients and the emmeans/emtrends grids, so it is off
    # unless a caller explicitly asks.
    ddf <- mixed_by(data, outcomes = outcome, rhs_model_formulae = formula,
                    split_on = split_on, padjust_by = 'term', padjust_method = 'fdr',
                    ncores = ncores, refit_on_nonconvergence = 3,
                    return_models = keep_models,
                    tidy_args = list(effects = c('fixed','ran_vals','ran_pars','ran_coefs'),
                                     conf.int = TRUE),
                    emmeans_spec = emm, emtrends_spec = emt)
    save(ddf, file = f)
    log_msg('  ', basename(f), ' = ', round(file.info(f)$size / 1e6, 1), ' MB')
    TRUE
  })
}

emm_sexHC <- list(Sex = list(outcome = 'vmPFC_decon', model_name = 'model1',
                             specs = formula(~sex:HCwithin), at = list(HCwithin = c(-1.5, 1.5))))
emm_ageHC <- list(Age = list(outcome = 'vmPFC_decon', model_name = 'model1',
                             specs = formula(~age:HCwithin), at = list(HCwithin = c(-1.5, 1.5))))
emm_rt <- list(RT = list(outcome = 'rt_csv_sc', model_name = 'model1',
                         specs = formula(~rt_lag_sc:sex), at = list(rt_lag_sc = c(-2,-1,0,1,2))))
emt_rt <- list(RT = list(outcome = 'rt_csv_sc', model_name = 'model1', var = 'rt_lag_sc',
                         specs = formula(~rt_lag_sc:sex), at = list(rt_lag_sc = c(-2,-1,0,1,2))))

## ===========================================================================
## STAGE 3 — fMRI models (sex series, age series, group sensitivity)
## ===========================================================================
log_msg('=== stage 3: fMRI models ===')

sexf <- list(
  formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*HCwithin + race_white*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)))
for (i in seq_along(sexf))
  run_mb('censored-vmPFC-HC-network-clock-Study2-', i, Q_fmri, sexf[[i]],
         'vmPFC_decon', splits, emm_sexHC)

agef <- list(
  formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + race_white*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)))
for (i in seq_along(agef))
  run_mb('censored-vmPFC-HC-network-clock-Study2-agemodels-', i, Q_fmri, agef[[i]],
         'vmPFC_decon', splits, emm_ageHC)

# group sensitivity (HC vs BPD; DEP and DNA dropped for this analysis only)
Qc <- Q_fmri %>% filter(!group %in% c('DNA','DEP')) %>%
  mutate(bpd = factor(ifelse(group == 'HC', 'HC', 'BPD'), levels = c('HC','BPD')))
grpf <- list(
  formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + bpd + (1 | id/run)),
  formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + bpd + (1 | id/run)),
  formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + bpd + (1 | id/run)),
  formula(~sex*bpd*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~age*HCwithin + sex*bpd*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)),
  formula(~sex*age*bpd*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)))
for (i in seq_along(grpf)) {
  emm <- list(Sex = list(outcome = 'vmPFC_decon', model_name = 'model1',
                         specs = if (i >= 4) formula(~sex:bpd:HCwithin) else formula(~sex:HCwithin),
                         at = list(HCwithin = c(-1.5, 1.5))))
  emt <- if (i >= 4) list(Coupling = list(outcome = 'vmPFC_decon', model_name = 'model1',
                                          var = 'HCwithin', specs = formula(~sex:bpd))) else NULL
  run_mb('censored-vmPFC-HC-network-clock-Study2-sensitivity_group', i, Qc, grpf[[i]],
         'vmPFC_decon', splits, emm, emt)
}

## ---- HC-only neural series -------------------------------------------------
# The Results report hippocampus-vPFC coupling within the healthy-control
# subsample. The existing files (2026-06-09, n = 52) predate the trial-alignment
# fix, so both series are re-run here. Uses the same formulas as the full-sample
# sex and age series.
Q_hcf <- Q_fmri %>% filter(group == 'HC')
log_msg('HC-only fMRI subsample: n = ', n_distinct(Q_hcf$id))
for (i in seq_along(sexf))
  run_mb('censored-vmPFC-HC-network-clock-Study2-HConly', i, Q_hcf, sexf[[i]],
         'vmPFC_decon', splits, emm_sexHC)
for (i in seq_along(agef))
  run_mb('censored-vmPFC-HC-network-clock-Study2-agemodels-HConly', i, Q_hcf, agef[[i]],
         'vmPFC_decon', splits, emm_ageHC)

## ===========================================================================
## STAGE 4 — behavioral models
## ===========================================================================
log_msg('=== stage 4: behavioral models ===')

bsex <- list(
  formula(~rt_lag_sc*sex + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*race_white + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)))
for (i in seq_along(bsex))
  run_mb('Age-Sex-clock-study2-pred-rt-', i, Q_behav, bsex[[i]], 'rt_csv_sc', NULL, emm_rt, emt_rt)

bage <- list(
  formula(~rt_lag_sc*age + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*race_white + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
  formula(~rt_lag_sc*age*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)))
emm_age <- list(RT = list(outcome = 'rt_csv_sc', model_name = 'model1',
                          specs = formula(~rt_lag_sc:age), at = list(age = c(-1.5,1.5), rt_lag_sc = c(-2,-1,0,1,2))))
emt_age <- list(RT = list(outcome = 'rt_csv_sc', model_name = 'model1', var = 'rt_lag_sc',
                          specs = formula(~rt_lag_sc:age), at = list(age = c(-1.5,1.5), rt_lag_sc = c(-2,-1,0,1,2))))
for (i in seq_along(bage))
  run_mb('Age-clock-Study2-pred-rt-', i, Q_behav, bage[[i]], 'rt_csv_sc', NULL, emm_age, emt_age)

# last-outcome series
blo <- list(
  formula(~rt_lag_sc*sex*reward_lag_rec + (1 | id/run)),
  formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)),
  formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*age*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)),
  formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*race_white*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)),
  formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*edu*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)))
emm_lo <- c(emm_rt, list(RTxO = list(outcome = 'rt_csv_sc', model_name = 'model1',
                                     specs = formula(~rt_lag_sc:reward_lag_rec:sex), at = list(rt_lag_sc = c(-2,-1,0,1,2)))))
emt_lo <- c(emt_rt, list(RTxO = list(outcome = 'rt_csv_sc', model_name = 'model1', var = 'rt_lag_sc',
                                     specs = formula(~rt_lag_sc:reward_lag_rec:sex), at = list(rt_lag_sc = c(-2,-1,0,1,2)))))
for (i in seq_along(blo))
  run_mb('Age-Sex-clock-Study2-pred-rt-lastoutcome-', i, Q_behav, blo[[i]], 'rt_csv_sc', NULL, emm_lo, emt_lo)

# HC-only. Try (1|id/run) first now that both runs are present; fall back to (1|id).
Q_hc <- Q_behav %>% filter(group == 'HC')
bhc <- list(formula(~rt_lag_sc*sex + (1 | id/run)),
            formula(~rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
            formula(~rt_lag_sc*sex + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
            formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)))
for (i in seq_along(bhc)) {
  ok <- run_mb('Age-Sex-clock-study2-pred-rt-HConly', i, Q_hc, bhc[[i]], 'rt_csv_sc', NULL, emm_rt, emt_rt)
  if (is.null(ok)) {
    log_msg('  retrying HConly', i, ' with (1 | id) only')
    f2 <- as.formula(gsub('1 \\| id/run', '1 | id', deparse1(bhc[[i]])))
    run_mb('Age-Sex-clock-study2-pred-rt-HConly-idonly', i, Q_hc, f2, 'rt_csv_sc', NULL, emm_rt, emt_rt)
  }
}

# behavioral group sensitivity
Qcb <- Q_behav %>% filter(!group %in% c('DNA','DEP')) %>%
  mutate(bpd = factor(ifelse(group == 'HC','HC','BPD'), levels = c('HC','BPD')))
bgrp <- list(formula(~rt_lag_sc*sex*bpd + (1 | id/run)),
             formula(~rt_lag_sc*sex*bpd + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
             formula(~rt_lag_sc*sex*bpd + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)),
             formula(~rt_lag_sc*sex*age*bpd + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)))
emm_g <- list(RT = list(outcome='rt_csv_sc', model_name='model1', specs=formula(~rt_lag_sc:sex:bpd),
                        at = list(rt_lag_sc = c(-2,-1,0,1,2))))
emt_g <- list(RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc',
                        specs=formula(~rt_lag_sc:sex:bpd)))
for (i in seq_along(bgrp))
  run_mb('Age-Sex-clock-study2-pred-rt-sensitivity_group-', i, Qcb, bgrp[[i]], 'rt_csv_sc', NULL, emm_g, emt_g)

## ---- Figure 2 input: per-participant RT-swing slopes -----------------------
# The descriptive layer of Figure 2 cannot be rebuilt from mixed_by output (the
# behavioural models are random-intercept only), so it is computed here from the
# same frame the models used and saved as a few-kilobyte file. Study 1's half is
# added in the session, where Q3 lives.
timed('Study 2 participant slopes', {
  ok <- tryCatch({ source_first('plotting.R'); TRUE },
                 error = function(e) { log_msg('  plotting.R unavailable: ',
                                               conditionMessage(e)); FALSE })
  if (!ok) stop('skipped: could not source plotting.R')
  subj2 <- subject_slopes(Q_behav) %>%
    left_join(distinct(Q_behav, id, sex), by = 'id') %>%
    mutate(panel = 'Experiment 2',
           Sex = factor(ifelse(sex == 1, 'Male', 'Female'),
                        levels = c('Female','Male'))) %>%
    select(panel, id, Sex, slope)
  log_msg('  n participants with a slope: ', nrow(subj2))
  saveRDS(subj2, file.path(outdir, paste0(run_tag, '_Study2_subject_slopes.rds')))
  TRUE
})

## ===========================================================================
## STAGE 5 — text-reported summaries (fast; all printed to the log)
## ===========================================================================
log_msg('=== stage 5: summaries ===')
sink(logfile, append = TRUE, split = TRUE)

cat('\n\n########## DEMOGRAPHICS (behavioral sample) ##########\n')
s2b <- demo_summary(distinct(as.data.frame(Q_behav), id, .keep_all = TRUE), 'Study 2 behavioral')
print_demo_summary(s2b)
saveRDS(s2b, file.path(outdir, paste0(run_tag, '_table1_study2_behav.rds')))

cat('\n\n########## DEMOGRAPHICS (fMRI sample) ##########\n')
s2f <- demo_summary(distinct(as.data.frame(Q_fmri), id, .keep_all = TRUE), 'Study 2 fMRI')
print_demo_summary(s2f)
saveRDS(s2f, file.path(outdir, paste0(run_tag, '_table1_study2_fmri.rds')))

cat('\n\n########## EARNINGS BY SEX ##########\n')
earn <- Q_behav %>% distinct(id, sex1, total_earnings)
print(t.test(total_earnings ~ sex1, data = earn))
print(earn %>% group_by(sex1) %>% summarize(mean = mean(total_earnings), sd = sd(total_earnings), n = n()))

cat('\n\n########## SOFTMAX BETA (subject-level mfx estimates) ##########\n')
bsub <- Q_behav %>% distinct(id, sex, sex1, age, race_white, edu, beta) %>% filter(!is.na(beta))
cat('n =', nrow(bsub), '| unique beta =', length(unique(bsub$beta)),
    '| sd =', signif(sd(bsub$beta, na.rm = TRUE), 4), '\n')
print(summary(bsub$beta))
for (nm in c('beta ~ sex','beta ~ sex + age','beta ~ sex*age',
             'beta ~ sex + race_white','beta ~ sex + edu',
             'beta ~ age','beta ~ age + sex','beta ~ age*sex')) {
  cat('\n--- ', nm, ' ---\n'); print(summary(lm(as.formula(nm), data = bsub)))
}
cat('\n--- focal sex effect, for the writeup ---\n')
mb <- lm(beta ~ sex, data = bsub)
print(cbind(summary(mb)$coefficients, confint(mb))['sex', , drop = FALSE])
cat('sex coded 0 = female, 1 = male; positive = higher beta (less stochastic choice) in males\n')
print(bsub %>% group_by(sex1) %>%
        summarize(mean_beta = mean(beta), sd = sd(beta), n = n(), .groups = 'drop'))

cat('\n\n########## GROUP x SEX CELL COUNTS ##########\n')
print(xtabs(~ sex1 + group, data = distinct(as.data.frame(Q_behav), id, .keep_all = TRUE)))
sink()

## ===========================================================================
## STAGE 6 — Mplus data prep for the MSEMs (Mplus itself run separately)
## ===========================================================================
log_msg('=== stage 6: Mplus data prep ===')
timed('mplus data prep', {
  library(MplusAutomation)
  # All Mplus files for both studies live in rootdir1/Mplus
  mplus_dir <- file.path(rootdir1, 'Mplus')
  if (!dir.exists(mplus_dir)) stop('no Mplus directory at ', mplus_dir)
  Q1 <- Q_fmri %>%
    group_by(id, run, trial, run_trial, network, HC_region, rt_csv, rt_lag, age, sex) %>%
    summarize(vmPFC_decon = mean(vmPFC_decon, na.rm = TRUE),
              HCwithin = mean(HCwithin, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(values_from = c(vmPFC_decon, HCwithin), names_from = 'network')
  for (rg in c('AH','PH')) {
    d <- Q1 %>% filter(HC_region == rg)
    prepareMplusData(df = as.data.frame(d),
                     filename = file.path(mplus_dir,
                                          paste0('bsocial_HC_vmPFC_clock_', rg, '_forMplus_', run_tag, '.dat')),
                     overwrite = TRUE)
  }
  TRUE
})
log_msg('NOTE: random-slope extraction and both MSEMs must be run in Mplus, then')
log_msg('      re-merged and re-run — see the notes accompanying this script.')

log_msg('=== run complete ===')