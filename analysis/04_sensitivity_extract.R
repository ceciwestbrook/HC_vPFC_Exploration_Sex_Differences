### 04_sensitivity_extract.R — race/education sensitivity extraction
# Sourced once. Each study calls extract_sensitivity() in its own section after
# its models finish, writing a small .rds. assemble_sensitivity() builds the
# report from whichever files exist, so the two studies never need to be run,
# or even finished, together.
#
# Model indices:  neural  base = 1, race = 4, edu = 5
#                 behav   base = 2, race = 5, edu = 6

library(dplyr); library(tidyr); library(purrr); library(tibble)

load_obj <- function(f){
  e <- new.env(); load(f, envir = e); o <- ls(e)
  n <- intersect(c("ddf","ddq"), o); e[[ if (length(n)) n[1] else o[1] ]]
}
# Pick a model file. Files are named <date>-<tag> by the original script and
# <date>_<tag> by the batch runner, so patterns allow either separator. With
# `date` given, that vintage is returned; otherwise the latest date prefix wins
# (by filename, not mtime, so a touched file cannot masquerade as current).
pick_file <- function(dir, rx, date = NULL){
  f <- list.files(dir, pattern = rx, full.names = TRUE)
  if (!length(f)) return(NA_character_)
  d <- sub("^(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", basename(f))
  if (!is.null(date)) { f <- f[d == date]; if (!length(f)) return(NA_character_); return(f[1]) }
  f[order(d, decreasing = TRUE)][1]
}
# fixed-effect rows whose term components match a set, order-independent
get_term <- function(ddx, comps){
  cd <- as_tibble(ddx$coef_df_reml) %>% filter(effect == "fixed")
  hit <- vapply(cd$term, function(t) setequal(strsplit(t, ":", fixed = TRUE)[[1]], comps), logical(1))
  cd[hit, ]
}
# whichever sex variable this file used
sex_term <- function(ddx){
  trm <- unique(unlist(strsplit(ddx$coef_df_reml$term, ":", fixed = TRUE)))
  hit <- intersect(c("sex","sex1M","sex1"), trm)
  if (!length(hit)) NA_character_ else hit[1]
}

# templates: one pattern per study x modality, either separator
patterns <- list(
  Study1 = list(neural = "censored-vmPFC-HC-network-clock-Study1-sexmodels%d\\.Rdata$",
                behav  = "[-_]Sex-clock-Study1-fmri-meg-pred-rt-%d\\.Rdata$"),
  Study2 = list(neural = "censored-vmPFC-HC-network-clock-Study2-%d\\.Rdata$",
                behav  = "[-_]Age-Sex-clock-study2-pred-rt-%d\\.Rdata$")
)
idx <- list(neural = c(base = 1, race = 4, edu = 5),
            behav  = c(base = 2, race = 5, edu = 6))

extract_one <- function(dir, study, modality, cov, i, date = NULL){
  f <- pick_file(dir, sprintf(patterns[[study]][[modality]], i), date)
  if (is.na(f)) { message("  [missing] ", study, "/", modality, "/", cov); return(NULL) }
  ddx <- load_obj(f)
  sx  <- sex_term(ddx)
  if (is.na(sx)) { warning(study,"/",modality,"/",cov,": no sex term in ", basename(f)); return(NULL) }
  
  trm <- unique(ddx$coef_df_reml$term)
  cov_ok <- switch(cov, base = TRUE,
                   race = any(grepl("race", trm)),
                   edu  = any(grepl("edu",  trm)))
  if (!cov_ok) warning(study,"/",modality,"/",cov,": covariate term absent from ", basename(f))
  
  comps <- if (modality == "neural") c(sx, "HCwithin") else c("rt_lag_sc", sx)
  g <- get_term(ddx, comps)
  if (!nrow(g)) { warning(study,"/",modality,"/",cov,": focal term ",
                          paste(comps, collapse = ":"), " absent"); return(NULL) }
  
  padjcol    <- intersect(c("padj_fdr_term","padj"), names(g))[1]
  split_cols <- intersect(c("network","HC_region","evt_time","dataset"), names(g))
  g %>% transmute(study, modality, cov, term,
                  across(all_of(split_cols)),
                  across(all_of(c("estimate","std.error","statistic",
                                  "p.value","conf.low","conf.high"))),
                  padj = if (!is.na(padjcol)) .data[[padjcol]] else NA_real_,
                  cov_ok, file = basename(f),
                  run_date = sub("^(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", basename(f)))
}

# Call this in each study's section. covs = c("base","race") to skip education.
# `date` pins a vintage; omit it for the latest.
extract_sensitivity <- function(dir, study, out_path = NULL,
                                modalities = c("neural","behav"),
                                covs = c("base","race","edu"), date = NULL){
  message("Extracting ", study, " from ", dir, if (!is.null(date)) paste0(" [", date, "]") else "")
  res <- map_dfr(modalities, function(md)
    map_dfr(covs, ~ extract_one(dir, study, md, .x, idx[[md]][[.x]], date)))
  if (!nrow(res)) { warning("nothing extracted for ", study); return(invisible(NULL)) }
  print(res %>% distinct(modality, cov, cov_ok, file))
  if (!is.null(out_path)) { saveRDS(res, out_path); message("Saved ", out_path) }
  invisible(res)
}

# ---------------------------------------------------------------------------
# Assembly: reads whatever exists, reports on that
# ---------------------------------------------------------------------------
assemble_sensitivity <- function(paths){
  have <- paths[file.exists(unlist(paths))]
  if (!length(have)) stop("no sensitivity files found")
  message("Using: ", paste(basename(unlist(have)), collapse = ", "))
  res <- map_dfr(have, readRDS)
  
  cat("\n===== integrity: covariate present where expected =====\n")
  print(res %>% distinct(study, modality, cov, cov_ok, file))
  
  cat("\n===== behavioral: rt_lag x sex across covariates =====\n")
  bh <- res %>% filter(modality == "behav")
  if (nrow(bh)) print(bh %>%
                        mutate(dataset = if ("dataset" %in% names(bh)) dataset else NA_character_) %>%
                        arrange(study, dataset, match(cov, c("base","race","edu"))) %>%
                        select(study, dataset, cov, estimate, std.error, p.value, conf.low, conf.high))
  
  nr <- res %>% filter(modality == "neural")
  if (nrow(nr)) {
    cat("\n===== neural: stability of sex x HCwithin =====\n")
    print(nr %>% select(study, cov, network, HC_region, evt_time, estimate) %>%
            pivot_wider(names_from = cov, values_from = estimate) %>%
            group_by(study) %>%
            summarize(across(any_of(c("race","edu")),
                             list(r = ~ cor(.x, base, use = "complete.obs"),
                                  max_abs_d = ~ max(abs(.x - base), na.rm = TRUE))),
                      n_cells = n(), .groups = "drop"))
    
    cat("\n--- FDR-significant cells per study x covariate ---\n")
    print(nr %>% group_by(study, cov) %>%
            summarize(n_sig = sum(padj < .05, na.rm = TRUE), n = n(), .groups = "drop") %>%
            arrange(study, match(cov, c("base","race","edu"))))
    
    cat("\n--- AH x LIM cells (mediation pathway) ---\n")
    print(nr %>% filter(network == "LIM", HC_region == "AH") %>%
            select(study, cov, evt_time, estimate, padj) %>%
            pivot_wider(names_from = cov, values_from = c(estimate, padj)) %>%
            arrange(study, evt_time))
  }
  invisible(res)
}