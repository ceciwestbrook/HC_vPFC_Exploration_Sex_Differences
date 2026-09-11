# Mplus input files

Copy in only the files behind the published results:

```bash
cp /ix/cladouceur/DNPL_sex_differences_paper/Mplus/Study1_*.inp .
cp /ix/cladouceur/DNPL_sex_differences_paper/Mplus/Study2_*.inp .
```

**Not the whole directory.** It also contains `BSOC_*`, `bsocial_*`, `mmclock_*`,
`*_0to4_*` and `*_evt_time0_*` files from earlier time-window analyses that the
paper does not report.

Expected, per study: six random-slope models
(`get_HC_vmPFC_randslopes_{AH,PH}_{CTR,DMN,LIM}`), the two MSEMs
(`omnibus_msem_2f`, `omnibus_msem_simpres`), and for Study 2 the two
group-adjusted variants (`*_groupadj`).

The `.dat` files these read are built by `analysis/01_build_mplus_data.R` and are
gitignored, since they contain participant-level data.
