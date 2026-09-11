### exclusions.R — participant and run exclusions, read from a data-side file.
###
### The exclusion lists are NOT in this repository. They name individual
### participants, and this is a clinical sample; a public script stating which
### participants were dropped is more disclosure than the analysis requires.
### The lists live in a CSV that travels with the data instead.
###
### Expected file, at `exclusions_file` (see config.R):
###
###   study,id,run,reason
###   study2,100001,,age_over_50
###   study2,100003,1,preprocessing_failure
###
### `run` is blank for whole-participant exclusions and a run number when only
### one run is dropped. `reason` is a short non-clinical tag, used for reporting
### counts rather than for filtering.
###
### A template with fabricated ids is in `exclusions_TEMPLATE.csv`.

load_exclusions <- function(path = exclusions_file) {
  if (!file.exists(path)) {
    stop("exclusions file not found at:\n  ", path,
         "\n\nThis file is not distributed with the code because it names ",
         "individual participants.\nIt should accompany the data. See ",
         "exclusions_TEMPLATE.csv for the expected format.", call. = FALSE)
  }
  ex <- utils::read.csv(path, colClasses = "character", na.strings = c("", "NA"))
  need <- c("study", "id", "run", "reason")
  if (!all(need %in% names(ex)))
    stop("exclusions file must have columns: ", paste(need, collapse = ", "),
         call. = FALSE)
  ex$run <- suppressWarnings(as.integer(ex$run))
  ex
}

## Participant ids excluded outright, for one study.
excluded_ids <- function(study, path = exclusions_file) {
  ex <- load_exclusions(path)
  unique(ex$id[ex$study == study & is.na(ex$run)])
}

## Single runs excluded, as a data frame of id and run.
excluded_runs <- function(study, path = exclusions_file) {
  ex <- load_exclusions(path)
  ex[ex$study == study & !is.na(ex$run), c("id", "run")]
}

## Apply both to a trial-level frame. Reports what it removed, so the numbers
## in the paper can be checked against the log rather than taken on trust.
apply_exclusions <- function(dat, study, path = exclusions_file,
                             id_col = "id", run_col = "run", verbose = TRUE) {
  ids  <- excluded_ids(study, path)
  runs <- excluded_runs(study, path)
  n0 <- length(unique(as.character(dat[[id_col]])))

  keep <- !as.character(dat[[id_col]]) %in% ids
  if (nrow(runs) && run_col %in% names(dat)) {
    drop_run <- paste(as.character(dat[[id_col]]), dat[[run_col]]) %in%
                paste(runs$id, runs$run)
    keep <- keep & !drop_run
  }
  out <- dat[keep, , drop = FALSE]

  if (isTRUE(verbose)) {
    n1 <- length(unique(as.character(out[[id_col]])))
    message(sprintf("exclusions (%s): %d participants -> %d (%d removed), %d run(s) dropped",
                    study, n0, n1, n0 - n1, nrow(runs)))
    tab <- table(load_exclusions(path)$reason[load_exclusions(path)$study == study])
    for (r in names(tab)) message(sprintf("   %-24s %d", r, tab[[r]]))
  }
  out
}


### --- id remapping -----------------------------------------------------------
### A small number of participants were entered under incorrect ids. The mapping
### is data-side for the same reason as the exclusions.
###
### Expected file, at `id_remap_file` (see config.R):
###
###   from,to
###   900001,100024
###
### If the file is absent, ids pass through unchanged and a message says so —
### the analysis still runs, but the affected participants will not merge.

remap_ids <- function(ids, path = id_remap_file) {
  if (!exists("id_remap_file") || !file.exists(path)) {
    message("no id remap file at ", if (exists("id_remap_file")) path else "id_remap_file",
            "; ids passed through unchanged")
    return(ids)
  }
  map <- utils::read.csv(path, colClasses = "character")
  if (!all(c("from", "to") %in% names(map)))
    stop("id remap file must have columns: from, to", call. = FALSE)
  out <- as.character(ids)
  hit <- match(out, map$from)
  n <- sum(!is.na(hit))
  out[!is.na(hit)] <- map$to[hit[!is.na(hit)]]
  if (n) message("remapped ", n, " id(s)")
  out
}
