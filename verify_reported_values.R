### verify_reported_values.R — check the paper's headline numbers against the
### model output, without re-running anything.
###
###   source("verify_reported_values.R")
###
### Reads existing model output and Mplus .out files and compares a handful of
### values to what the manuscript reports. Every check prints the expected value,
### the value found, and whether they agree to the reported precision.
###
### The point is not to re-derive the analysis. It is to catch the failure mode
### this project actually hit: a script silently loading a stale file, or an
### edited number drifting out of step with the model it came from. Run it before
### submission, and again after any re-run.
###
### A check can fail for three reasons, in decreasing order of likelihood:
###   1. the file it read is not the one the manuscript was written from
###   2. the manuscript value is stale
###   3. the model genuinely changed
### The output names the file it read for each check, so 1 can be ruled out first.

if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
source(file.path(repo_directory, "R/plotting.R"))   # load_obj, newest, find_term

suppressPackageStartupMessages({
  library(dplyr); library(MplusAutomation)
})

## ---- reporting helpers ----------------------------------------------------
.results <- new.env(parent = emptyenv())
.results$rows <- list()

check <- function(label, expected, found, tol = NULL, source_file = "") {
  if (is.na(found)) {
    verdict <- "NOT FOUND"
  } else {
    if (is.null(tol)) {
      # match to the precision the manuscript reports
      dp <- nchar(sub("^[^.]*\\.?", "", format(expected, scientific = FALSE)))
      tol <- 10^(-dp) / 2
    }
    verdict <- if (abs(found - expected) <= tol) "ok" else "MISMATCH"
  }
  .results$rows[[length(.results$rows) + 1]] <-
    list(label = label, expected = expected, found = found, verdict = verdict,
         file = basename(source_file))
  cat(sprintf("  %-9s %-46s reported %-9s found %-9s %s\n",
              verdict, label,
              format(expected, scientific = FALSE),
              if (is.na(found)) "--" else format(round(found, 4), scientific = FALSE),
              if (nzchar(source_file)) paste0("[", basename(source_file), "]") else ""))
  invisible(verdict)
}

safe <- function(expr) tryCatch(expr, error = function(e) { message("    (", conditionMessage(e), ")"); NA_real_ })

## ===========================================================================
cat("\n=== Behavioural sex difference in RT swings (Results, Figure 2) ===\n")
## ===========================================================================
f <- safe(newest(mbo1, "[-_]Sex-clock-Study1-fmri-meg-pred-rt-2\\.Rdata$"))
if (!is.na(f)) {
  d <- load_obj(f); g <- find_term(d, c("rt_lag_sc", "sex"))
  check("Study 1 fMRI, rt_lag x sex", -0.13,
        safe(g$estimate[g$dataset == "fMRI"]), tol = 0.005, source_file = f)
  check("Study 1 MEG, rt_lag x sex", -0.19,
        safe(g$estimate[g$dataset == "MEG"]), tol = 0.005, source_file = f)
}
f <- safe(newest(mbo2, "[-_]Age-Sex-clock-study2-pred-rt-2\\.Rdata$"))
if (!is.na(f)) {
  d <- load_obj(f); g <- find_term(d, c("rt_lag_sc", "sex"))
  check("Study 2, rt_lag x sex", -0.06, safe(g$estimate[1]), tol = 0.005, source_file = f)
}

## ===========================================================================
cat("\n=== Neural sex difference in connectivity (Results) ===\n")
## ===========================================================================
f <- safe(newest(mbo2, "censored-vmPFC-HC-network-clock-Study2-1\\.Rdata$"))
if (!is.na(f)) {
  d <- load_obj(f); g <- as_tibble(find_term(d, c("sex", "HCwithin")))
  check("Study 2 sex x HCwithin: cells total", 102, nrow(g), tol = 0, source_file = f)
  check("Study 2 sex x HCwithin: FDR-significant", 102,
        sum(g$padj_fdr_term < .05, na.rm = TRUE), tol = 0, source_file = f)
}

## ===========================================================================
cat("\n=== Age x sex on connectivity: only AH-DMN replicates (Figure 5) ===\n")
## ===========================================================================
f1 <- safe(newest(mbo1, "censored-vmPFC-HC-network-clock-Study1-agemodels-?3\\.Rdata$"))
f2 <- safe(newest(mbo2, "censored-vmPFC-HC-network-clock-Study2-agemodels-3\\.Rdata$"))
if (!is.na(f1) && !is.na(f2)) {
  sig_signs <- function(f, net, reg) {
    g <- as_tibble(find_term(load_obj(f), c("age", "sex", "HCwithin"))) %>%
      filter(network == net, HC_region == reg, padj_fdr_term < .05)
    if (!nrow(g)) return("none")
    paste0(nrow(g), " cells, ",
           if (all(g$estimate > 0)) "all positive"
           else if (all(g$estimate < 0)) "all negative" else "mixed signs")
  }
  for (net in c("DMN", "CTR", "LIM")) for (reg in c("AH", "PH")) {
    cat(sprintf("  %-4s %s   Study 1: %-28s Study 2: %s\n", reg, net,
                sig_signs(f1, net, reg), sig_signs(f2, net, reg)))
  }
  cat("  Expected: AH-DMN significant and same-signed in both. AH-CTR significant\n",
      " in both but OPPOSITE signs, which is why the legend claims only one\n",
      " replicable effect.\n")
}

## ===========================================================================
cat("\n=== MSEM: the mediation (Results, Figure 6, Tables S2/S4) ===\n")
## ===========================================================================
msem <- function(file, section, param) {
  p <- file.path(mplusdir, file)
  if (!file.exists(p)) { message("    missing: ", file); return(c(NA_real_, NA_real_)) }
  v <- readModels(p)$parameters$unstandardized
  r <- v[v$paramHeader == section & v$param == param, ]
  if (!nrow(r)) return(c(NA_real_, NA_real_))
  c(r$est[1], r$pval[1])
}

x <- msem("Study2_omnibus_msem_simpres.out", "New.Additional.Parameters", "LIM_IND")
check("Study 2 AH-LIM indirect, estimate", -0.031, x[1], tol = 0.0005,
      source_file = "Study2_omnibus_msem_simpres.out")
check("Study 2 AH-LIM indirect, p", 0.006, x[2], tol = 0.0005,
      source_file = "Study2_omnibus_msem_simpres.out")

x <- msem("Study2_omnibus_msem_2f.out", "New.Additional.Parameters", "AH_IND")
check("Study 2 two-factor AH indirect, estimate", -0.023, x[1], tol = 0.0005,
      source_file = "Study2_omnibus_msem_2f.out")
check("Study 2 two-factor AH indirect, p", 0.031, x[2], tol = 0.0005,
      source_file = "Study2_omnibus_msem_2f.out")

x <- msem("Study2_omnibus_msem_simpres_groupadj.out", "New.Additional.Parameters", "LIM_IND")
check("Study 2 AH-LIM indirect, group-adjusted", -0.035, x[1], tol = 0.0005,
      source_file = "Study2_omnibus_msem_simpres_groupadj.out")

## ===========================================================================
cat("\n=== Sample sizes (Table 1, STAR Methods) ===\n")
## ===========================================================================
clusters <- function(file) {
  p <- file.path(mplusdir, file)
  if (!file.exists(p)) return(NA_real_)
  l <- grep("Number of clusters", readLines(p, warn = FALSE), value = TRUE)[1]
  if (is.na(l)) NA_real_ else as.numeric(gsub("\\D", "", l))
}
check("MSEM clusters, Study 2", 185, clusters("Study2_omnibus_msem_simpres.out"),
      tol = 0, source_file = "Study2_omnibus_msem_simpres.out")
check("MSEM clusters, Study 2 group-adjusted", 177,
      clusters("Study2_omnibus_msem_simpres_groupadj.out"), tol = 0,
      source_file = "Study2_omnibus_msem_simpres_groupadj.out")

f <- safe(newest(mbo2, "Study2_Qbehav_fixed\\.rds$"))
if (!is.na(f)) {
  check("Study 2 behavioural N", 187, n_distinct(readRDS(f)$id), tol = 0, source_file = f)
}

## ===========================================================================
cat("\n=== Summary ===\n")
## ===========================================================================
v <- vapply(.results$rows, function(r) r$verdict, character(1))
cat(sprintf("  %d checks: %d ok, %d mismatched, %d not found\n",
            length(v), sum(v == "ok"), sum(v == "MISMATCH"), sum(v == "NOT FOUND")))
if (any(v == "MISMATCH")) {
  cat("\n  MISMATCHES -- check, in this order, whether the file read is the one the\n",
      " manuscript was written from; whether the manuscript value is stale; or\n",
      " whether the model changed:\n")
  for (r in .results$rows[v == "MISMATCH"])
    cat(sprintf("    %-46s reported %s, found %s  [%s]\n", r$label,
                format(r$expected), format(round(r$found, 4)), r$file))
}
if (any(v == "NOT FOUND")) {
  cat("\n  NOT FOUND -- the output file or term is absent. Expected if you have not\n",
      " run that part of the pipeline; a problem if you have:\n")
  for (r in .results$rows[v == "NOT FOUND"]) cat("    ", r$label, "\n")
}
invisible(NULL)
