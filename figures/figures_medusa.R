### figures_medusa.R — all MEDuSA panels, both studies, one function
#
#   source("figures/figures_medusa.R")
#
# Why this exists: Study 2's panels were regenerated with plot_term_ts() while
# Study 1's came from the original plot_mixed_by_vmPFC_HC(), so they no longer
# match within a composite figure. Study 1's models were never affected by any
# of the bugs, so regenerating them here is purely restyling — the estimates do
# not change. Both halves now come from the same call with the same theme.
#
# flipy = FALSE throughout: these are connectivity estimates, where higher means
# more coupling. The default reverses the y axis, which is right for RT-swing
# slopes and wrong here.

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


## ---- panel geometry -------------------------------------------------------

## ---- why the two studies looked so different ------------------------------
# Beyond the TR count, three things made the panels mismatch, two of which are
# now fixed in plotting.R:
#  * Point sizes saturated. Every Study 2 cell is p < .001, so every point took
#    the maximum size and hid the lines and error bars, while Study 1's mixed
#    significance gave it small points. The size range is now narrow (1.6-3.4).
#  * Legend keys differed, because unused significance levels were dropped.
#    They are now retained, so both panels carry the same key.
#  * The y scales differ by up to an order of magnitude between studies. Each
#    panel now gets its OWN scale, so neither is flattened. Say in the legend
#    that the axes differ, since a reader comparing panel heights would
#    otherwise infer the effects are comparable in size when Study 1's age x sex
#    effect is roughly a tenth of Study 2's.

grab <- function(dir, tag) {
  rx <- paste0("^\\d{4}-\\d{2}-\\d{2}[-_]", tag, "\\.Rdata$")
  f  <- list.files(dir, pattern = rx, full.names = TRUE)
  if (!length(f)) { message("  [missing] ", tag); return(NULL) }
  d <- sub("^(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", basename(f))
  chosen <- f[order(d, decreasing = TRUE)][1]
  message("  ", basename(chosen))
  load_obj(chosen)
}

audit <- function(ddf, comps, label) {
  g <- find_term(ddf, comps)
  if (!nrow(g)) { message("  [", label, "] term not found"); return(invisible(NULL)) }
  g <- as_tibble(g)
  cat(sprintf("  %-34s %3d of %3d cells significant | %+.4f to %+.4f\n",
              label, sum(g$padj_fdr_term < .05, na.rm = TRUE), nrow(g),
              min(g$estimate), max(g$estimate)))
  invisible(g)
}

# One file per study, both at the same dimensions, matching how the published
# figures are built.
#
# The published figures give each study its own full-width row, assembled in
# Photoshop (Figure 3 panels c and d, for instance). So both studies are written
# at the SAME width, each with its own axes and legend, and the composite stacks
# them. Equal width is the point: it is what lets Study 2's 17 timepoints crowd
# together into overlapping points while Study 1's 9 stay separated, which is
# how the published figures distinguish the two samples at a glance.
PANEL_W <- 7.5      # same for both studies
PANEL_H <- 3.0

make_pair <- function(tag1, tag2, comps1, comps2, ylab, stem, share_y = FALSE) {
  message("\n=== ", stem, " ===")
  d1 <- grab(mbo1, tag1); d2 <- grab(mbo2, tag2)
  # share_y = FALSE by default: each study gets its own y range. The effects
  # differ by up to an order of magnitude between studies (Study 1's age x sex
  # effect spans about +/- 5e-04 against Study 2's -0.008 to 0), so a shared
  # scale flattens the smaller one into a line. Pass share_y = TRUE for a
  # specific figure if the two are close enough in magnitude to share an axis.
  yl <- if (share_y) shared_ylim(list(d1, d2), comps1) else NULL
  if (!is.null(yl)) cat(sprintf("  shared y range: %+.5f to %+.5f\n", yl[1], yl[2]))
  if (!is.null(d1)) audit(d1, comps1, "Study 1")
  if (!is.null(d2)) audit(d2, comps2, "Study 2")
  
  # single-study case: one panel, nothing to stack
  if (is.null(d1) || is.null(d2)) {
    d  <- if (is.null(d1)) d2 else d1
    cm <- if (is.null(d1)) comps2 else comps1
    save_fig(plot_term_ts(d, cm, flipy = FALSE, ylab = ylab, ylim = yl),
             paste0(stem, ".pdf"), figdir, width = PANEL_W, height = PANEL_H)
    return(invisible(list(s1 = d1, s2 = d2)))
  }
  
  # one file per study, same width and height, each with its own axes and
  # legend, ready to stack in the composite
  save_fig(plot_term_ts(d1, comps1, flipy = FALSE, ylab = ylab, ylim = yl),
           paste0(stem, "_Study1.pdf"), figdir, width = PANEL_W, height = PANEL_H)
  save_fig(plot_term_ts(d2, comps2, flipy = FALSE, ylab = ylab, ylim = yl),
           paste0(stem, "_Study2.pdf"), figdir, width = PANEL_W, height = PANEL_H)
  invisible(list(s1 = d1, s2 = d2))
}

## ==========================================================================
## FIGURE 3 — sex x HCwithin
## ==========================================================================
make_pair("censored-vmPFC-HC-network-clock-Study1-sexmodels1",
          "censored-vmPFC-HC-network-clock-Study2-1",
          c("sex", "HCwithin"), c("sex", "HCwithin"),
          "Male \u2212 Female [AU]", "Figure3_sex_HCwithin")

## ==========================================================================
## FIGURE 4 — age x HCwithin
## ==========================================================================
make_pair("censored-vmPFC-HC-network-clock-Study1-agemodels-?1",
          "censored-vmPFC-HC-network-clock-Study2-agemodels-1",
          c("age", "HCwithin"), c("age", "HCwithin"),
          "Age effect [AU]", "Figure4_age_HCwithin")

## ==========================================================================
## FIGURE 5 — age x sex x HCwithin
## ==========================================================================
f5 <- make_pair("censored-vmPFC-HC-network-clock-Study1-agemodels-?3",
                "censored-vmPFC-HC-network-clock-Study2-agemodels-3",
                c("age", "sex", "HCwithin"), c("age", "sex", "HCwithin"),
                "Age \u00d7 sex [AU]", "Figure5_agesex_HCwithin")
cat("\n  LEGEND CHECK for Figure 5: PH-CTR no longer replicates. Study 1 has 2\n",
    " significant cells, both positive (0 and +1 s); Study 2 has 7, all negative,\n",
    " at the window edges. AH-DMN is the only replicating age x sex effect.\n")

## ==========================================================================
## FIGURES S4 and S5 — Vmax and entropy moderation
## ==========================================================================
## Study 1 model indices for these are a guess: the sex series runs 1-10 in both
## studies, so sexmodels8 and sexmodels10 should be the Vmax and entropy
## moderation models. Check the term list printed below; if the terms are absent,
## the indices differ and need adjusting.
make_pair("censored-vmPFC-HC-network-clock-Study1-sexmodels8",
          "censored-vmPFC-HC-network-clock-Study2-8",
          c("sex", "v_max_wi", "HCwithin"), c("sex", "v_max_wi", "HCwithin"),
          "Sex \u00d7 Vmax [AU]", "FigureS4_vmax_moderation")
make_pair("censored-vmPFC-HC-network-clock-Study1-sexmodels10",
          "censored-vmPFC-HC-network-clock-Study2-10",
          c("sex", "v_entropy_wi", "HCwithin"), c("sex", "v_entropy_wi", "HCwithin"),
          "Sex \u00d7 entropy [AU]", "FigureS5_entropy_moderation")
cat("\n  LEGEND CHECK for S4 and S5: in Study 2, Vmax is significant in 6 of 102\n",
    " cells (1 in AH-LIM) and entropy in 29 of 102 (3 in AH-LIM), so the AH-LIM\n",
    " timing claims in the current legends are not supportable.\n")

## ==========================================================================
## FIGURE S2 — sex effect on coupling, by diagnostic group (Study 2 only)
## ==========================================================================
## Reproduces the original plot_coupling_sexeffect() exactly. Two things about
## it that differ from the other MEDuSA panels:
##  * It plots the Male - Female contrast computed WITHIN each group from the
##    Coupling emtrends grid, not the sex x group x HCwithin interaction term.
##    That is why the original shows large positive values in HC and near-zero
##    in BPD rather than a single interaction trace.
##  * Significance is |sexeff| > 1.96 * se_diff, uncorrected, with se_diff taken
##    as sqrt(se_M^2 + se_F^2). Two levels only, NS and p < .05.
## No Study 1 counterpart: Study 1 has no diagnostic groups.
message("\n=== FigureS2 sex effect on coupling by group (Study 2 only) ===")

plot_coupling_sexeffect <- function(ddf, save_as = NULL, flipy = FALSE) {
  if (is.null(ddf$emtrends_list) || is.null(ddf$emtrends_list$Coupling))
    stop("no Coupling emtrends in this model - it is only specified for group ",
         "models 4-6, so check the file")
  d <- as.data.frame(ddf$emtrends_list$Coupling)
  sexcol <- if ("sex1" %in% names(d)) "sex1" else "sex"
  d$Sex  <- ifelse(d[[sexcol]] %in% c("M", "1", 1), "Male", "Female")
  ycol   <- grep("\\.trend$", names(d), value = TRUE)[1]
  
  dw <- d %>%
    select(t = evt_time, network, HC_region, bpd, Sex,
           est = all_of(ycol), se = std.error) %>%
    pivot_wider(names_from = Sex, values_from = c(est, se)) %>%
    mutate(sexeff  = est_Male - est_Female,
           se_diff = sqrt(se_Male^2 + se_Female^2),
           sig     = factor(ifelse(abs(sexeff) > 1.96 * se_diff, "p < .05", "NS"),
                            levels = c("NS", "p < .05")),
           network2 = factor(network, levels = c("DMN", "CTR", "LIM")),
           Group    = factor(bpd, levels = c("HC", "BPD")))
  
  cat(sprintf("  cells: %d | significant: %d\n", nrow(dw), sum(dw$sig == "p < .05")))
  print(dw %>% group_by(Group, HC_region, network2) %>%
          summarize(n = n(), sig = sum(sig == "p < .05"),
                    mean_eff = mean(sexeff), .groups = "drop"), n = Inf)
  
  pd <- position_dodge(width = 0.33)
  p <- ggplot(dw, aes(x = t, y = sexeff, group = network2, color = network2)) +
    geom_point(aes(size = sig, alpha = sig), position = pd) +
    geom_hline(yintercept = 0, lty = "dashed", color = "#A9A9A9") +
    geom_vline(xintercept = 0, lty = "dashed", color = "#A9A9A9") +
    geom_line(position = pd, linewidth = 1) +
    geom_errorbar(aes(ymin = sexeff - se_diff, ymax = sexeff + se_diff),
                  position = pd, width = 0, color = "black") +
    facet_grid(HC_region ~ Group) +
    scale_color_manual(values = net_pal) +
    scale_size_manual(values = c("NS" = 2.5, "p < .05" = 5)) +
    scale_alpha_manual(values = c("NS" = .4, "p < .05" = 1)) +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
    xlab("Time relative to trial onset [s]") +
    ylab("Male > Female [AU]") +
    labs(color = NULL, size = NULL, alpha = NULL) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 13),
          strip.text = element_text(size = 13), panel.spacing = unit(1.2, "lines"))
  if (isTRUE(flipy)) p <- p + scale_y_reverse()
  if (!is.null(save_as)) ggsave(file.path(figdir, save_as), p, width = 9, height = 5)
  p
}

sg4 <- grab(mbo2, "censored-vmPFC-HC-network-clock-Study2-sensitivity_group4")
if (!is.null(sg4)) plot_coupling_sexeffect(sg4, "FigureS2_coupling_sexeffect_model4.pdf")

message("\nAll MEDuSA panels written to ", figdir)
message("Both study halves now come from plot_term_ts() with the same theme, so ",
        "they should drop into the composites without a visible style change.")
message("One file per study, both at ", PANEL_W, " x ", PANEL_H, " inches, to be ",
        "stacked in the composite.")
message("Equal width is deliberate: it is what makes Study 2's 17 timepoints ",
        "crowd into overlapping points while Study 1's 9 stay separated.")
message("Y axes are free per study, so neither panel is flattened. Note in the ",
        "legends that the axes differ between panels.")