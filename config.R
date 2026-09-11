### config.R — every path the analysis needs, in one place.
###
### Edit the two roots below, then `source("config.R")` at the top of any
### script. Nothing else in the repository should contain an absolute path.
###
### Why two roots: Study 1 (MMClock) and Study 2 (BSOC) live in separate trees
### and the original scripts switched between them mid-file. They are named
### explicitly here so it is always clear which study a path belongs to.

## --- edit these two ---------------------------------------------------------
rootdir1 <- "/ix/cladouceur/DNPL_sex_differences_paper"        # Study 1 (MMClock)
rootdir2 <- "/ix/cladouceur/DNPL_sex_differences_paper/BSOC"   # Study 2 (BSOC)

## --- derived; no need to edit ----------------------------------------------
# The directory holding this config.R, i.e. the repository root. Resolved from
# the sourced path so scripts work whether run from the root or a subdirectory.
repo_directory <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE),
  error = function(e) normalizePath(".", mustWork = FALSE))

mbo1    <- file.path(rootdir1, "fmri/mixed_by_output")   # Study 1 model output
mbo2    <- file.path(rootdir2, "fMRI/mixed_by_output")   # Study 2 model output
mplusdir <- file.path(rootdir1, "Mplus")                 # .inp, .dat and .out
figdir  <- file.path(rootdir2, "figures")                # figure output
logdir  <- file.path(rootdir2, "fMRI/rerun_logs")        # batch logs

## Participant and run exclusions. Deliberately NOT in the repository: the list
## names individual participants in a clinical sample. It travels with the data.
## See exclusions_TEMPLATE.csv for the format, and R/exclusions.R for the reader.
exclusions_file <- file.path(rootdir1, "exclusions_for_repo/exclusions.csv")
id_remap_file   <- file.path(rootdir1, "exclusions_for_repo/id_remap.csv")

## Compute resources for mixed_by(). Override with the NCORES environment
## variable when submitting a batch job.
ncores <- as.integer(Sys.getenv("NCORES", unset = "8"))

## R library path, if you keep packages outside the system library.
## Comment out if not needed.
# .libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))

for (d in c(mbo1, mbo2, mplusdir)) {
  if (!dir.exists(d)) warning("config.R: directory not found: ", d, call. = FALSE)
}
