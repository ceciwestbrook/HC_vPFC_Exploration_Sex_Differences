### figures_S1_S6_S7.R — supplementary figures rebuilt from the re-run models
#
#   source("figures/figures_behavioral.R")
#
# The originals use ggplot's default two-colour hue palette and facet by
# condition with Sex on the x axis. That is deliberately NOT the FantasticFox1
# palette in plotting.R, so these helpers carry their own theme.
#
# Note: the last-outcome models carry TWO emtrends specs. `$RT` gives slopes by
# sex only (averaged over outcome); `$RTxO` splits by outcome as well. Use RTxO.

if (!exists("rootdir1")) {
  .cfg <- NULL
  for (.p in c("config.R", "../config.R", "../../config.R")) {
    if (file.exists(.p)) { .cfg <- .p; break }
  }
  if (is.null(.cfg)) stop("cannot find config.R; run from inside the repository")
  source(.cfg)
}
source(file.path(repo_directory, "R/plotting.R"))
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
library(dplyr); library(tidyr); library(ggplot2); library(patchwork)


hue2 <- c(Female = "#F8766D", Male = "#00BFC4")   # ggplot defaults, as in the originals

orig_theme <- function(base = 12) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = base - 1),
          axis.text = element_text(colour = "black"))
}

grab <- function(dir, tag) {
  rx <- paste0("^\\d{4}-\\d{2}-\\d{2}[-_].*", tag, "\\.Rdata$")
  f <- list.files(dir, pattern = rx, full.names = TRUE)
  if (!length(f)) stop("no file matching ", tag, " in ", dir)
  d <- sub("^(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", basename(f))
  chosen <- f[order(d, decreasing = TRUE)][1]
  message("  ", basename(chosen))
  load_obj(chosen)
}

# pull the outcome-split emtrends grid and label it for plotting
lo_grid <- function(ddf, study_label) {
  g <- ddf$emtrends_list$RTxO
  if (is.null(g)) stop("no RTxO emtrends spec in this model")
  g <- as.data.frame(g)
  # reward_lag_rec is +0.5 reward / -0.5 omission; sex is 0 female / 1 male
  g %>%
    mutate(Outcome = factor(ifelse(reward_lag_rec > 0, "Reward", "Omission"),
                            levels = c("Omission", "Reward")),
           Sex = factor(ifelse(sex == 1, "Male", "Female"),
                        levels = c("Female", "Male")),
           lo = rt_lag_sc.trend - 1.96 * std.error,
           hi = rt_lag_sc.trend + 1.96 * std.error,
           study = study_label,
           facet = if ("dataset" %in% names(.))
             paste0(dataset, "\n", Outcome) else as.character(Outcome)) %>%
    distinct(facet, Sex, .keep_all = TRUE)
}

# One panel: Sex on x, facetted by condition, reversed y.
# Point + capped error bar rather than a crossbar: a box invites the question
# "box of what?", whereas a point with whiskers reads as an estimate and its
# interval without explanation. The interval is named in the caption so the
# reader never has to guess.
#
# Colours come from sex_pal in plotting.R, matching Figure 2. Since both study
# halves are being regenerated together, they restyle consistently. Swap in
# `hue2` if you would rather keep the original salmon/teal.
slope_panel <- function(d, facet_levels = NULL, ylab = "RT swing (RT ~ previous RT)",
                        show_y = TRUE, pal = sex_pal) {
  if (!is.null(facet_levels)) d$facet <- factor(d$facet, levels = facet_levels)
  p <- ggplot(d, aes(x = Sex, y = rt_lag_sc.trend, colour = Sex)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.8) +
    geom_point(size = 3.2) +
    facet_wrap(~ facet, nrow = 1) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_y_reverse() +
    labs(x = NULL, y = if (show_y) ylab else NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 11),
          axis.text = element_text(colour = "black"))
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  p
}


## ==========================================================================
## FIGURE S1 — exploratory behaviour by diagnostic group
## ==========================================================================
## Matches the original layout, which differs from S7: Group on the x axis,
## Sex as colour, connecting lines, legend on the right. The lines belong here
## because the interaction is the point of the figure.
message("\n=== Figure S1 ===")
sg <- grab(mbo2, "clock-study2-pred-rt-sensitivity_group-2")
et <- as.data.frame(sg$emtrends_list$RT)
cat("  columns:", paste(names(et), collapse = ", "), "\n")

grpvar <- intersect(c("bpd", "group", "bpd_num"), names(et))[1]
if (is.na(grpvar)) stop("no group column in the emtrends grid; check the spec")

d_s1 <- et %>%
  mutate(Group = factor(ifelse(grepl("^BPD|^1$", as.character(.data[[grpvar]])),
                               "BPD", "HC"), levels = c("HC", "BPD")),
         Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
         lo = rt_lag_sc.trend - 1.96 * std.error,
         hi = rt_lag_sc.trend + 1.96 * std.error) %>%
  distinct(Group, Sex, .keep_all = TRUE)
print(d_s1 %>% select(Group, Sex, rt_lag_sc.trend, std.error, lo, hi))

p_s1 <- ggplot(d_s1, aes(x = Group, y = rt_lag_sc.trend, colour = Sex, group = Sex)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.8) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = sex_pal, name = "Sex") +
  scale_y_reverse() +
  labs(x = "Group", y = "RT swing (slope of RT on previous RT)")+
       #caption = "Points are estimated RT-swing slopes; whiskers are 95% confidence intervals. Lower slopes indicate larger RT swings.") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(colour = "black"),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey30"))
save_fig(p_s1, "FigureS1_behaviour_by_group.pdf", figdir, width = 5.4, height = 4.2)

cat("\n  Expected: HC 0.475 F / 0.325 M, BPD 0.420 F / 0.400 M.\n",
    " LEGEND CHECK: the female line moves (0.42 -> 0.475 in HC) while the male\n",
    " values hold, so the convergence is no longer symmetric and the claim that\n",
    " females are 'comparable across groups' (0.475 vs 0.420) no longer holds.\n",
    " The F test is now F(1, 42504) = 36.93, p = 1.2e-9, up from F(1, 21097) = 14.62.\n")

## ==========================================================================
## FIGURE S7 — RT swings by last outcome, both studies
## ==========================================================================
message("\n=== Figure S7 ===")
lo1 <- grab(mbo1, "clock-Study1-fmri-meg-pred-rt-lastoutcome-2")
lo2 <- grab(mbo2, "clock-Study2-pred-rt-lastoutcome-2")
d1 <- lo_grid(lo1, "Study 1")
d2 <- lo_grid(lo2, "Study 2")

cat("\nStudy 1 slopes:\n");  print(d1 %>% select(facet, Sex, rt_lag_sc.trend, std.error))
cat("\nStudy 2 slopes:\n");  print(d2 %>% select(facet, Sex, rt_lag_sc.trend, std.error))

# match the original facet order: fMRI Omission, fMRI Reward, MEG Omission, MEG Reward
lev1 <- c("fMRI\nOmission","fMRI\nReward","MEG\nOmission","MEG\nReward")
lev1 <- lev1[lev1 %in% d1$facet]
p1 <- slope_panel(d1, lev1) +
  labs(title = "a) Study 1") +
  theme(plot.title = element_text(face = "bold", size = 13))
p2 <- slope_panel(d2, c("Omission","Reward"), show_y = FALSE) +
  labs(title = "b) Study 2") +
  theme(plot.title = element_text(face = "bold", size = 13))

# share the y range across panels so the two halves are comparable
yr <- range(c(d1$lo, d1$hi, d2$lo, d2$hi))
p1 <- p1 + scale_y_reverse(limits = rev(yr))
p2 <- p2 + scale_y_reverse(limits = rev(yr))

p_s7 <- p1 + p2 +
  plot_layout(widths = c(length(lev1), 2)) +
  plot_annotation(
    caption = "Points are estimated RT-swing slopes; whiskers are 95% confidence intervals. Lower slopes indicate larger RT swings.",
    theme = theme(plot.caption = element_text(hjust = 0, size = 9, colour = "grey30")))
save_fig(p_s7, "FigureS7_last_outcome_both_studies.pdf", figdir,
         width = 9.5, height = 3.9)

cat("\n  Study 2's three-way is now 0.02, p = 0.25 \u2014 null. The legend's claim of an\n",
    " opposite-direction effect in Study 2 no longer holds. Study 1 is unchanged\n",
    " (fMRI -0.09, MEG -0.07, both p < 0.001).\n")

## ==========================================================================
## FIGURE S6 — total earnings by sex, both studies
## ==========================================================================
## Styled to match Figure S7: point plus capped whiskers, captioned so the
## interval needs no explanation. NB the original used +/- 1 SE; this uses 95%
## CIs for consistency with S7, so the bars are about twice as long.
## Study 1 needs Q3 in scope.
message("\n=== Figure S6 ===")

earn_tab <- function(d, label) {
  d %>% group_by(sex) %>%
    summarize(TotalEarnings = mean(total_earnings),
              se = sd(total_earnings) / sqrt(n()), n = n(), .groups = "drop") %>%
    mutate(sex = factor(ifelse(sex == 1, "Male", "Female"),
                        levels = c("Female", "Male")),
           dataset = label)
}

tabs <- list()
if (exists("Q3")) {
  for (ds in c("fMRI", "MEG")) {
    d <- Q3 %>% filter(dataset == ds) %>%
      distinct(id, sex, total_earnings) %>% drop_na()
    if (nrow(d)) {
      cat("\nStudy 1", ds, ":\n"); print(t.test(total_earnings ~ sex, data = d))
      tabs[[ds]] <- earn_tab(d, paste("Study 1", ds))
    }
  }
} else message("  Q3 not loaded - Study 1 panels skipped")

Q_behav <- readRDS(newest(mbo2, "Study2_Qbehav_fixed\\.rds$"))
d2e <- Q_behav %>% distinct(id, sex, total_earnings) %>% drop_na()
cat("\nStudy 2:\n"); print(t.test(total_earnings ~ sex, data = d2e))
tabs[["Study2"]] <- earn_tab(d2e, "Study 2")

tab <- bind_rows(tabs)
print(tab)

p_s6 <- ggplot(tab, aes(x = sex, y = TotalEarnings, colour = sex)) +
  geom_errorbar(aes(ymin = TotalEarnings - 1.96 * se, ymax = TotalEarnings + 1.96 * se),
                width = 0.18, linewidth = 0.8) +
  geom_point(size = 3.2) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = sex_pal, guide = "none") +
  labs(x = NULL, y = "Total earnings (points)") +
       #caption = "Points are mean total earnings; whiskers are 95% confidence intervals.") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 11),
        axis.text = element_text(colour = "black"),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey30"))
save_fig(p_s6, "FigureS6_earnings_both_studies.pdf", figdir,
         width = 2.9 * n_distinct(tab$dataset), height = 4)

cat("\n  NOTE: whiskers are now 95% CIs, matching Figure S7. The original used +/- 1 SE,\n",
    " so the bars will look wider than before - state the interval in the legend.\n",
    " Study 2 is now t(74.7) = 1.84, p = 0.069, females higher (9522 vs 9176).\n",
    " The Supplement's Study 1 MEG t was a transcription error: t(61) = 0.98, not\n",
    " 0.33; the p of 0.33 was correct.\n")

message("\nWritten to ", figdir)
message("Check facet labels and panel widths against the originals before replacing them.")