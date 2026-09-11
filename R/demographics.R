### demographics.R — demographic recoding for both studies
# Sourced once at the top of the analysis script. Defines:
#   prep_demo_study1(path)          -> one row per participant
#   prep_demo_study2(path1, path2)  -> one row per participant
if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
source(file.path(repo_directory, "R/exclusions.R"))
# Both return: id, age, sex (0 = F, 1 = M), sex1 ("F"/"M"), race_cat, race_white,
#              n_races, edu, plus study-specific extras (Hispanic, group, gender).
#
# Race category order in both sources: American Indian or Alaska Native, Asian,
# Black or African American, Native Hawaiian or Other Pacific Islander, White.
# Endorsement is coded -1 in Study 1 and 1 in Study 2, so it is tested as nonzero.
# Study 2 additionally has registration_race___999 = declined to provide.

library(dplyr)

race_lv <- c("American Indian or Alaska Native","Asian","Black or African American",
             "Native Hawaiian or Other Pacific Islander","White",
             "More than one race","Not reported")

race_factor <- function(ind, declined = NULL) {
  end <- as.matrix(ind) != 0
  end[is.na(end)] <- FALSE
  n_end <- rowSums(end)
  lab <- ifelse(n_end > 1, "More than one race",
                ifelse(n_end == 1, race_lv[apply(end, 1, which.max)], "Not reported"))
  factor(lab, levels = race_lv)
}

collapse_white <- function(race_cat) {
  factor(ifelse(race_cat == "Not reported", NA,
                ifelse(race_cat == "White", "White", "Other or multiple")),
         levels = c("White","Other or multiple"))
}

# -8 / -9 are sentinel missing codes in the Study 1 registry
na_sentinel <- function(x, codes = c(-8, -9)) ifelse(x %in% codes, NA, x)

### Replace the two prep functions in demographics.R with these.
# Only change from before: race_1..race_5 (logical) are created and returned,
# so demo_summary() works identically for both studies.

prep_demo_study1 <- function(path) {
  d <- read.csv(path, header = TRUE)
  race_cols <- c("AmerIndianAlaskan","Asian","Black","HawaiianPacIsl","White")
  stopifnot(all(race_cols %in% names(d)))
  r <- d[, race_cols, drop = FALSE]
  
  d %>%
    mutate(id     = as.character(id),
           sex    = ifelse(female == 0, 1, ifelse(female == 1, 0, NA)),
           sex1   = ifelse(sex == 1, "M", "F"),
           race_1 = !is.na(r[[1]]) & r[[1]] != 0,
           race_2 = !is.na(r[[2]]) & r[[2]] != 0,
           race_3 = !is.na(r[[3]]) & r[[3]] != 0,
           race_4 = !is.na(r[[4]]) & r[[4]] != 0,
           race_5 = !is.na(r[[5]]) & r[[5]] != 0,
           race_cat   = race_factor(r),
           n_races    = rowSums(as.matrix(r) != 0, na.rm = TRUE),
           race_white = collapse_white(race_cat),
           hispanic_yn = factor(ifelse(is.na(Hispanic), NA,
                                       ifelse(Hispanic != 0, "Hispanic or Latino",
                                              "Not Hispanic or Latino")),
                                levels = c("Not Hispanic or Latino","Hispanic or Latino")),
           edu    = na_sentinel(edu)) %>%
    select(id, age, sex, sex1, race_1, race_2, race_3, race_4, race_5,
           race_cat, race_white, n_races, hispanic_yn, edu) %>%
    distinct(id, .keep_all = TRUE)
}

prep_demo_study2 <- function(path_bsoc, path_ksoc) {
  race_cols <- paste0("registration_race___", 1:5)
  
  tidy_one <- function(d) {
    stopifnot(all(race_cols %in% names(d)))
    r <- d[, race_cols, drop = FALSE]
    d %>%
      mutate(race_1 = !is.na(r[[1]]) & r[[1]] != 0,
             race_2 = !is.na(r[[2]]) & r[[2]] != 0,
             race_3 = !is.na(r[[3]]) & r[[3]] != 0,
             race_4 = !is.na(r[[4]]) & r[[4]] != 0,
             race_5 = !is.na(r[[5]]) & r[[5]] != 0,
             race_cat   = race_factor(r, declined = d$registration_race___999),
             n_races    = rowSums(as.matrix(r) != 0, na.rm = TRUE),
             race_white = collapse_white(race_cat),
             sex        = ifelse(registration_birthsex == 1, 0,
                                 ifelse(registration_birthsex == 2, 1, NA)),
             sex1       = ifelse(sex == 1, "M", "F"),
             edu        = registration_edu) %>%
      rename(gender = registration_gender, group = registration_group) %>%
      select(id, age, sex, sex1, race_1, race_2, race_3, race_4, race_5,
             race_cat, race_white, n_races, edu, gender, group)
  }
  
  d1 <- read_csv(path_bsoc, show_col_types = FALSE) %>% mutate(id = as.character(id))
  d2 <- read_csv(path_ksoc, show_col_types = FALSE) %>%
    mutate(id = as.character(registration_redcapid))
  
  bind_rows(tidy_one(d1), tidy_one(d2)) %>% distinct(id, .keep_all = TRUE)
}

# Participants excluded from Study 2 for age or acquisition problems
study2_excluded_ids <- excluded_ids("study2")

### Append to demographics.R
# Adds race_1..race_5 (logical, common names across studies) to the prep output,
# plus a study-agnostic summary function so each study can be summarized
# independently and the table assembled later from small saved files.

# --- add to BOTH prep functions, replacing the select() line ---
# Study 1:  r <- d[, c("AmerIndianAlaskan","Asian","Black","HawaiianPacIsl","White")]
# Study 2:  r <- d[, paste0("registration_race___", 1:5)]
# then in the mutate, before select():
#     race_1 = r[[1]] != 0, race_2 = r[[2]] != 0, race_3 = r[[3]] != 0,
#     race_4 = r[[4]] != 0, race_5 = r[[5]] != 0,
# and add race_1:race_5 to the select() list.
# (NA becomes FALSE via the coalesce in race_endorsed below.)

race_labels <- c("American Indian or Alaska Native","Asian","Black or African American",
                 "Native Hawaiian or Other Pacific Islander","White")

# subj: ONE ROW PER PARTICIPANT from the analytic Q, i.e.
#       distinct(as.data.frame(Q), id, .keep_all = TRUE)
# Returns a tidy list of summaries for a single study.
demo_summary <- function(subj, study) {
  rc <- paste0("race_", 1:5)
  stopifnot(all(rc %in% names(subj)))
  m <- as.matrix(subj[, rc]); m[is.na(m)] <- FALSE
  k <- rowSums(m); N <- nrow(subj)
  
  pct <- function(x) round(100 * x / N, 1)
  
  list(
    study = study,
    n     = N,
    continuous = tibble(
      variable = c("Age","Years of education"),
      mean = c(mean(subj$age, na.rm = TRUE), mean(subj$edu, na.rm = TRUE)),
      sd   = c(sd(subj$age,   na.rm = TRUE), sd(subj$edu,   na.rm = TRUE)),
      min  = c(min(subj$age,  na.rm = TRUE), min(subj$edu,  na.rm = TRUE)),
      max  = c(max(subj$age,  na.rm = TRUE), max(subj$edu,  na.rm = TRUE)),
      missing = c(sum(is.na(subj$age)), sum(is.na(subj$edu)))
    ) %>% mutate(across(mean:sd, ~ round(.x, 2))),
    sex = tibble(category = c("Female","Male"),
                 n   = c(sum(subj$sex == 0, na.rm = TRUE), sum(subj$sex == 1, na.rm = TRUE))) %>%
      mutate(pct = pct(n)),
    race = tibble(category = race_labels, n = colSums(m)) %>% mutate(pct = pct(n)),
    race_other = tibble(category = c("More than one category","Not reported"),
                        n = c(sum(k > 1), sum(k == 0))) %>% mutate(pct = pct(n)),
    ethnicity = if ("hispanic_yn" %in% names(subj)) {
      subj %>% count(hispanic_yn, .drop = FALSE) %>% mutate(pct = pct(n))
    } else NULL,
    group = if ("group" %in% names(subj)) {
      subj %>% count(group, .drop = FALSE) %>% mutate(pct = pct(n))
    } else NULL,
    group_by_sex = if ("group" %in% names(subj)) {
      as.data.frame(xtabs(~ sex1 + group, data = subj))
    } else NULL
  )
}

# print one study's summary in table order
print_demo_summary <- function(s) {
  cat("\n========== ", s$study, " (N = ", s$n, ") ==========\n", sep = "")
  print(s$continuous); print(s$sex); print(s$race); print(s$race_other)
  if (!is.null(s$ethnicity))    print(s$ethnicity)
  if (!is.null(s$group))      { print(s$group); print(s$group_by_sex) }
  invisible(s)
}