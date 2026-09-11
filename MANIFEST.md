# File triage

A keep/archive decision for every script in `/ix/cladouceur/DNPL_sex_differences_paper`, with the
reasoning, so a coauthor can see why each call was made.

## Included, renamed for run order

| Repository path | Original | Why |
|---|---|---|
| `analysis/01a_build_mplus_data_study1.R` | `build_df_for_Mplus_redoStudy2.R`, Study 1 half | Split by study so each has one entry point. |
| `analysis/01b_build_mplus_data_study2.R` | same file, Study 2 half | **The corrected version** — protocol assignment removed. See the warning below. |
| `analysis/02_study1_models.R` | `code_for_paper.R`, lines 1–812 | Study 1 only. The original covered both studies in 2,288 lines with a mid-file root switch; split so each study has one entry point. |
| `analysis/03_study2_models.R` | `run_study2_all.R` | Study 2, all models, batch-submittable with alignment gates and skip logic. **Supersedes the Study 2 half of `code_for_paper.R`**, which is archived rather than published because it still contained the protocol-assignment code. Verified to cover every model output the old half wrote. |
| `analysis/04_sensitivity_extract.R` | `sensitivity_extract.R` | Race and education sensitivity. |
| `figures/figure2_rt_swings.R` | rescued from `code_for_paper.R` | **Figure 2.** The code was embedded in the Study 2 section of the combined script even though the figure spans both studies, so splitting by study would have dropped it. Rewritten to call the functions in `R/plotting.R`. |
| `figures/figures_medusa.R` | `Regenerate_medusa_panels.R` | Figures 3, 4, 5, S2, S4, S5. |
| `figures/figures_behavioral.R` | `makeplots_s1_s6_s7.R` | Figures S1, S6, S7. |
| `R/mixed_by.R` | `mixed_by.R` | Model runner, sourced by the analysis scripts. |
| `R/plotting.R` | `plotting.R` | Figure functions. |
| `R/demographics.R` | `demographics.R` | Participant preparation. |
| `R/get_trial_data.R` | `get_trial_data.R` | Task data loading. |
| `R/get_random_slope.R` | `get_random_slope.R` | Extracts slopes from Mplus savedata. |
| `R/exclusions.R` | new | Reads the data-side exclusion and id-remap files. Written for this repository. |
| `slurm/run_study2.sbatch` | same | Batch submission. |
| `verify_reported_values.R` | new | Checks the paper's headline numbers against the model output. Written for this repository. |

## Excluded — superseded

**The Study 2 half of `code_for_paper.R`** (lines 813–2288) is archived at
`archive_local/02_study2_half_SUPERSEDED.R`, outside the repository. It still
contains `split_ksoc_bsoc` and `maxT` matching, so running it would reproduce
the 179-participant sample rather than the published 185. Every model it wrote
is covered by `03_study2_models.R`.


**`code_for_paper_build_df_for_Mplus.R` — do not publish.** Despite the name,
this is *not* the version behind the paper. It still contains the
`split_ksoc_bsoc` protocol assignment, the `maxT` matching and the hard-coded
`440223` patch — the code that silently dropped six participants and produced
the 179-participant MSEM sample. The published results use 185. The file naming
is actively misleading and it is worth renaming locally to something like
`build_df_for_Mplus_PRE_FIX.R` so nobody reaches for it later.

| File | Why excluded |
|---|---|
| `code_for_paper_build_df_for_Mplus.R` | Pre-fix; produced the superseded 179-participant sample. |
| `build_df_for_Mplus.R` | Study 1 only, older `mmclock_*` output naming. |
| `build_df_for_Mplus_BSOC.R` | March 2025; old `BSOC_*` naming, protocol-assignment code. |
| `BSOC/build_df_forMplus_BSOC.R` | Same vintage, second copy in the BSOC subdirectory. |
| `BSOC/Build_mixed_by_vmPFC_HC_BSOC.R` | Superseded by `03_study2_models.R`. |
| `Build_mixed_by_vmPFC_HC.R`, `Build_mixed_by_vmPFC_HC_Papale.R` | Earlier model-building scripts. |
| `plot_mixed_by_vmPFC_HC.R` and its `_network_simplified` / `_simplified_networksymmetry` variants | Fully superseded by `R/plotting.R`, which provides a `plot_mixed_by_vmPFC_HC()` wrapper with the old signature. Every call in the analysis script passes `toprocess='network-by-HC'`, which the wrapper handles, so nothing is lost. |
| `plot_emmeans_subject_level_random_slopes.R`, `plot_emtrends_...`, `plot_subject_level_random_slopes.R` | Superseded by `plotting.R`. |
| `alex_vmPFC_HC_medusa.R`, `alex_vmPFC_HC_medusa_Ceciedits.R`, `Alex_plotting code.R` | Earlier collaborator code. |
| `regenerate_figures.R` | Superseded by `Regenerate_medusa_panels.R`. |
| `MSEMS_meta_regression.R`, `MSEM_with_lavaan.R` | Exploratory; not reported. |

## Changed, not just moved

**The beta (softmax inverse temperature) block in `02_study1_models.R` was
rewritten.** The original ran ten regressions of `beta` on trial-level
predictors (`v_max_wi`, `v_entropy_wi`) using the trial-level frame. Because
`beta` is one value per participant, this repeated each participant's beta
across hundreds of rows, inflating *n* from 71 to tens of thousands and making
the standard errors meaningless. Those models are removed rather than corrected:
a participant-level constant cannot answer whether beta tracks within-trial
value. The age model is kept but now computed on the participant-level frame,
which the original did not do. The values the Supplement reports (no relation
with sex or age) come from the participant-level models and are unaffected.

## Excluded — not code

`*.Rdata` (`bsocial_*`, `HC_clock_Aug2023`, `MMclock_clock_Aug2023`,
`mmclock_HC_vmPFC_clock`), `slurm-11148230.out`, the `tail` file (almost
certainly a mistyped redirect, safe to delete), and the `BSOC/`, `fmri/`, `meg/`,
`Mplus/`, `trial_data/` and `HC_vPFC_repo/` directories, which hold data and
model output.

## Excluded — throwaway diagnostics

`diag_msem_sample.R` and `diag_msem_sample.sbatch` were written to diagnose the
six missing MSEM participants. That question is resolved, so they document a
fixed bug rather than the analysis. Keep them locally if you want the record.

## Still to do

- [ ] Copy the Mplus `.inp` files for the paper into `mplus/`. **Not all of
      them**: the Mplus directory also holds `BSOC_*`, `bsocial_*`, `mmclock_*`,
      `_0to4_*` and `_evt_time0_*` variants from earlier time-window analyses.
      The paper uses `Study1_*` and `Study2_*` only.
- [ ] Choose a license. Without one, others cannot legally reuse the code.
- [x] ~~Decide about the hard-coded participant IDs.~~ **Done.** All real ids
      are out of the code. The exclusion list and the id remap now live in
      data-side CSVs (`exclusions.csv`, `id_remap.csv`), read by
      `R/exclusions.R`, gitignored, with templates in the repository showing the
      format using fabricated ids. `exclusions.csv` is written out for you and
      complete; `id_remap.csv` needs the two real pairs copied from the archived
      `get_trial_data.R`, which is deliberately not reproduced here.
- [ ] Fix the stale hard-coded dates in `04_sensitivity_extract.R` (F2), which
      will silently select the wrong model output once newer files exist.
- [x] ~~Redundant exclusion line (F5)~~ **Done.** `01b` applied the participant
      exclusions twice; now one `apply_exclusions()` call handles participants
      and runs together.
- [x] ~~Age scaled before exclusions~~ **Done.** `01b` scaled age and then
      removed the seven participants over 50, so age was centred on a sample
      that included them. Harmless here (age is not in `USEVARIABLES`) but the
      wrong order; exclusions now come first.
- [ ] `02_mixed_by_models.R` overwrites `basemodel_formula[[1]]` in the Study 1
      age block (F3), which makes the model numbering misleading.
- [x] ~~Save `Q3` (F4)~~ **Done.** `02_study1_models.R` now writes
      `<date>_Study1_Q3.rds` to the Study 1 output directory, so
      `figure2_rt_swings.R` and the Study 1 panels of Figure S6 can run without
      rebuilding it.
- [ ] Test from a clean checkout. Nothing here has been run — the restructuring
      was done without access to R or the data.
