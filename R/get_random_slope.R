if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
get_random_slope <- function(ddf,effect_of_interest){
  qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term==effect_of_interest)
  
  qdf <- qdf %>% rename(id=level)
  qdf <- qdf %>% group_by(id,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()
  qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
  qdf$id <- as.character(qdf$id)
  #qdf <- qdf %>% select(!outcome)
  source(file.path(repo_directory, "R/get_trial_data.R"))
  df <- get_trial_data(repo_directory=rootdir1,dataset='mmclock_fmri')
  df <- df %>% select(rt_vmax_lag,iti_prev,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                      v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                      rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
    group_by(id, run) %>% 
    mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           rt_vmax_sc = scale(rt_vmax),
           v_max_wi_lag = lag(v_max_wi),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change))
  df$id <- as.character(df$id)
  Q2 <- inner_join(qdf,df,by='id')
  Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
  return(Q2)
}


