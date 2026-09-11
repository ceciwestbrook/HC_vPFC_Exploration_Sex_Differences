# Anterior hippocampal–prefrontal connectivity and sex differences in exploration

Analysis code for Westbrook et al., *Anterior hippocampal-prefrontal connectivity
mediates male bias for exploration in humans*.

## What this is

This repository documents the analysis rather than reproducing it end to end.
Study 1 data are public; Study 2 imaging data are access-controlled through the
NIMH Data Archive, so the pipeline cannot be run from scratch without applying
for access. See **Data availability** below.

## Setup

Edit the two paths at the top of `config.R`:

```r
rootdir1 <- "/path/to/study1"   # MMClock
rootdir2 <- "/path/to/study2"   # BSOC
```

Every other path in the repository derives from those two. No script should
contain an absolute path; if you find one, it is a bug.

Requires R (developed under 4.6.0) with `lme4`, `lmerTest`, `emmeans`,
`tidyverse`, `data.table`, `MplusAutomation`, `wesanderson`, `patchwork`, and
Mplus 8.7 for the structural equation models.

## Run order

One entry point per study. The multilevel SEMs alternate between R and Mplus, so
the sequence is not a single pipeline.

**Study 1 (MMClock)**

| Step | Script | What it does |
|---|---|---|
| 1 | `analysis/02_study1_models.R` | Multilevel models: behaviour, connectivity, age and sex effects, sensitivity analyses. |
| 2 | `analysis/01a_build_mplus_data_study1.R` | Writes the `.dat` files for the random-slope models. |
| 3 | *Mplus* | `mplus/Study1_get_HC_vmPFC_randslopes_{AH,PH}_{CTR,DMN,LIM}.inp` — six models. |
| 4 | `analysis/01a_build_mplus_data_study1.R` (second half) | Merges the extracted slopes, writes the omnibus `.dat`. |
| 5 | *Mplus* | `Study1_omnibus_msem_2f.inp` and `Study1_omnibus_msem_simpres.inp`. |

**Study 2 (BSOC)**

| Step | Script | What it does |
|---|---|---|
| 1 | `analysis/03_study2_models.R` | All Study 2 models, as a batch job. Submit via `slurm/run_study2.sbatch`. Runs alignment checks first and aborts if they fail; skips models whose output already exists. `VERIFY_ONLY=1` stops after the data build, `SAVE_DATA=1` also writes the model frames. |
| 2 | `analysis/01b_build_mplus_data_study2.R` | Writes the `.dat` files for the random-slope models. |
| 3 | *Mplus* | `mplus/Study2_get_HC_vmPFC_randslopes_{AH,PH}_{CTR,DMN,LIM}.inp` — six models. Upload the `.out` **and savedata** files. |
| 4 | `analysis/01b_build_mplus_data_study2.R` (second half) | Merges the slopes, writes the omnibus and group-adjusted `.dat` files. |
| 5 | *Mplus* | The four MSEMs: `Study2_omnibus_msem_2f`, `_simpres`, and the two `_groupadj` variants. `FBITERATIONS = 20000`. |

**Both studies**

| Step | Script | What it does |
|---|---|---|
| 6 | `analysis/04_sensitivity_extract.R` | Race and education sensitivity comparisons. |
| 7 | `figures/figure2_rt_swings.R` | Figure 2. Spans both studies, so it lives here rather than in either analysis script. |
| 8 | `figures/figures_medusa.R` | Figures 3, 4, 5 and S2, S4, S5. |
| 9 | `figures/figures_behavioral.R` | Figures S1, S6, S7. |

## Repository layout

```
config.R                     all paths, edit this
R/                           functions sourced by the analysis scripts
  mixed_by.R                 the mixed_by() model runner
  plotting.R                 figure functions (plot_term_ts, plot_emm, save_fig)
  demographics.R             participant preparation and summaries
  get_trial_data.R           task data loading
  get_random_slope.R         extracts random slopes from Mplus savedata
analysis/
  01a_build_mplus_data_study1.R   Mplus .dat files, Study 1
  01b_build_mplus_data_study2.R   Mplus .dat files, Study 2
  02_study1_models.R              Study 1 multilevel models
  03_study2_models.R              Study 2 multilevel models, batch
  04_sensitivity_extract.R        race and education sensitivity
figures/
  figure2_rt_swings.R             Figure 2, both studies
  figures_medusa.R                Figures 3, 4, 5, S2, S4, S5
  figures_behavioral.R            Figures S1, S6, S7
mplus/                       Mplus input files
slurm/                       batch submission scripts
```

## Checking the numbers

`verify_reported_values.R` compares the paper's headline values against the
model output without re-running anything:

```r
source("verify_reported_values.R")
```

It checks the behavioural sex differences, the neural cell counts, which age ×
sex effects replicate, the MSEM indirect effects, and the sample sizes. Each
check prints the reported value, the value found, and the file it came from.

It exists because of the failure mode this project actually hit: a script
silently loading a stale file, or an edited number drifting out of step with the
model it came from. Run it before submission and after any re-run. A mismatch
most often means the file read is not the one the manuscript was written from —
which is why the filename is printed alongside each result.

## Participant exclusions

The exclusion lists are not in this repository. They name individual
participants in a clinical sample, so they live in a CSV that travels with the
data:

```
exclusions.csv    participant and single-run exclusions
id_remap.csv      corrections for participants entered under the wrong id
```

Set their locations in `config.R`. `exclusions_TEMPLATE.csv` and
`id_remap_TEMPLATE.csv` show the format with fabricated ids. `R/exclusions.R`
reads them and provides `excluded_ids()`, `excluded_runs()`,
`apply_exclusions()` and `remap_ids()`.

`apply_exclusions()` reports what it removed, so the sample sizes in the paper
can be checked against the log rather than taken on trust. If the exclusions
file is missing the scripts stop with a clear message; if the remap file is
missing, ids pass through unchanged and the affected participants simply will
not merge.

## Two things worth knowing before you run anything

**Model indices are not stable across scripts.** The behavioural age series,
sex series and last-outcome series each number their formulae independently, and
the same index means different models in different series. The figure scripts
name the file they load, so check that rather than assuming.

## Data availability

Study 1 behavioural and imaging data: [Zenodo 10.5281/zenodo.3978642](https://doi.org/10.5281/zenodo.3978642)

SCEPTIC reinforcement learning model: [Zenodo 10.5281/zenodo.3978658](https://doi.org/10.5281/zenodo.3978658)

Study 2 behavioural data and clinical assessments: [NIMH Data Archive Collection 2739](https://doi.org/10.15154/r1g0-f618).
Study 2 imaging data are available from the lead contact on reasonable request,
subject to participant consent and IRB approval.

## Citation

This code is associated with a manuscript currently under review.
