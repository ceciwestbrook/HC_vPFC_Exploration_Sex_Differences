### figure2_rt_swings.R — Figure 2: RT-swing slopes by sex, both studies.
###
###   source("figures/figure2_rt_swings.R")
###
### Three panels across the top (Study 1 fMRI, Study 1 MEG, Study 2), each a
### violin of participant-level slopes with the model-estimated group slope and
### its 95% CI overlaid; a difference row below showing the male − female
### interaction estimate per panel.
###
### This figure spans both studies, which is why it lives here rather than in
### either study's analysis script. The original code was embedded in the Study 2
### section of the combined analysis script; the plotting functions it used now
### live in R/plotting.R.
###
### Inputs: the participant-level frames and the model-2 outputs for both
### studies. If `Q3` and `Q_behav` are not in the session, the script says what
### to run rather than failing obscurely.

if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
source(file.path(repo_directory, "R/plotting.R"))

library(dplyr); library(tidyr); library(purrr); library(ggplot2); library(patchwork)

panel_lv <- c("Experiment 1\n(fMRI)", "Experiment 1\n(MEG)", "Experiment 2")

## ---- participant-level frames ---------------------------------------------
## Study 1: Q3 is built by analysis/02_study1_models.R. It is not saved there by
## default, so either keep the session or add a saveRDS() call (see README).
if (!exists("Q3")) {
  f <- try(newest(mbo1, "Study1_Q3\\.rds$"), silent = TRUE)
  if (inherits(f, "try-error"))
    stop("Q3 not found. Either run analysis/02_study1_models.R in this session, ",
         "or save it there with\n  saveRDS(Q3, file.path(mbo1, ",
         "paste0(curr_date, \"_Study1_Q3.rds\")))", call. = FALSE)
  Q3 <- readRDS(f)
}
## Study 2: the behavioural frame, saved by analysis/03_study2_models.R
if (!exists("Q")) {
  f <- try(newest(mbo2, "Study2_Qbehav_fixed\\.rds$"), silent = TRUE)
  if (inherits(f, "try-error"))
    stop("Study 2 behavioural frame not found. Run analysis/03_study2_models.R ",
         "with SAVE_DATA=1.", call. = FALSE)
  Q <- readRDS(f)
}

## ---- participant-level slopes ---------------------------------------------
subj <- bind_rows(
  map_dfr(c("fMRI", "MEG"), function(ds) {
    d <- Q3 %>% filter(dataset == ds)
    subject_slopes(d) %>%
      left_join(distinct(d, id, sex), by = "id") %>%
      mutate(panel = paste0("Experiment 1\n(", ds, ")"),
             Sex = ifelse(sex == 1, "Male", "Female"))
  }) %>% select(panel, Sex, slope),
  subject_slopes(Q) %>%
    left_join(distinct(Q, id, sex1), by = "id") %>%
    mutate(panel = "Experiment 2",
           Sex = ifelse(sex1 == "M", "Male", "Female")) %>%
    select(panel, Sex, slope)
) %>%
  mutate(Sex = factor(Sex, levels = c("Female", "Male")),
         panel = factor(panel, levels = panel_lv))

cat("participant-level slopes:\n")
print(subj %>% count(panel, Sex))

## ---- model-estimated group slopes and the sex difference -------------------
## Model 2 in each behavioural series: the reported specification.
d1 <- load_obj(newest(mbo1, "[-_]Sex-clock-Study1-fmri-meg-pred-rt-2\\.Rdata$"))
d2 <- load_obj(newest(mbo2, "[-_]Age-Sex-clock-study2-pred-rt-2\\.Rdata$"))

s1d <- as.data.frame(d1$emtrends_list$RT) %>% distinct(dataset, sex, .keep_all = TRUE)
s2d <- as.data.frame(d2$emtrends_list$RT) %>% distinct(sex1, .keep_all = TRUE)

est <- bind_rows(
  tibble(panel = paste0("Experiment 1\n(", s1d$dataset, ")"),
         Sex = ifelse(s1d$sex == 1, "Male", "Female"),
         slope = s1d$rt_lag_sc.trend, se = s1d$std.error),
  tibble(panel = "Experiment 2",
         Sex = ifelse(s2d$sex1 == "M", "Male", "Female"),
         slope = s2d$rt_lag_sc.trend, se = s2d$std.error)
) %>%
  mutate(Sex = factor(Sex, levels = levels(subj$Sex)),
         panel = factor(panel, levels = panel_lv),
         lo = slope - 1.96 * se, hi = slope + 1.96 * se)

## The difference row. sex_diff_table() handles the term-name difference between
## studies (`sex` in Study 1, `sex1M` in Study 2) and warns if a file looks stale.
dif <- sex_diff_table(d1, d2, panel_lv = panel_lv)
cat("\nsex difference per panel (expect three negative estimates):\n")
print(dif)

## ---- assemble --------------------------------------------------------------
p <- plot_subject_slopes(subj = subj, dif = dif)
save_fig(p, "Figure2_RTswing_sex_violin.pdf", figdir, width = 8, height = 5.5)

cat("\nReported values for the legend: Study 1 fMRI -0.13, MEG -0.19, Study 2 -0.06,\n",
    "all p < 0.001. Axes are reversed, so lower slopes mean larger RT swings.\n")
