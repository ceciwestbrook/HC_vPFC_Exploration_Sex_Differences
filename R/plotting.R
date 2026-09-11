### plotting.R — all figures for the vPFC-HC sex differences paper
#
# Source once alongside mixed_by.R and demographics.R. Four plotting functions,
# one shared aesthetic layer, one save helper.
#
#   plot_term_ts()        estimate +/- SE over time for a fixed effect, from
#                         coef_df_reml. Point size/alpha keyed to FDR p.
#                         -> Figures 3, 4, 5, and any neural effect to inspect.
#   plot_mixed_by_vmPFC_HC()  compatibility wrapper: same signature as the old
#                         function, loops over terms and writes one PDF each.
#   plot_emm()            emmeans or emtrends output, either over time (lines +
#                         ribbon) or by group (points + intervals).
#   subject_slopes()      per-participant RT-swing slopes from trial data.
#   plot_subject_slopes() violin + jitter + difference panel -> Figure 2.
#
# Conventions: sex is numeric (0 = female, 1 = male) in both studies. Network
# colours follow the manuscript: CTR yellow, DMN orange, LIM blue. y axes are
# reversed where lower values mean larger RT swings (flipy = TRUE).

library(dplyr); library(tidyr); library(purrr); library(ggplot2)
library(stringr); library(wesanderson); library(patchwork)

## ---------------------------------------------------------------- aesthetics
.pal3    <- wes_palette("FantasticFox1", 3, type = "discrete")
net_pal  <- c(CTR = .pal3[2], DMN = .pal3[1], LIM = .pal3[3])
sex_pal  <- setNames(wes_palette("FantasticFox1", 5)[c(3, 5)], c("Female", "Male"))
# Wide size range, matching the published figures: point size encodes
# significance, so in the larger Study 2 sample where every cell is p < .001 the
# points are all maximum size and overlap. That overlap is intended, not a
# defect — it is what makes the two studies visually distinguishable at a glance.
p_sizes  <- c("NS" = 2, "p < .05" = 4, "p < .01" = 5, "p < .001" = 6)
p_alphas <- c("NS" = .35, "p < .05" = .7, "p < .01" = .85, "p < .001" = 1)

fig_theme <- function(base = 13) {
  theme_bw(base_size = base) +
    theme(legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = base + 1),
          axis.text = element_text(size = base, colour = "black"),
          strip.text = element_text(size = base),
          panel.spacing = unit(1.2, "lines"))
}

save_fig <- function(p, filename, dir = ".", width = 9, height = 3.5) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir, filename)
  ggsave(f, p, width = width, height = height, device = cairo_pdf)
  message("wrote ", f); invisible(f)
}

## ------------------------------------------------------------ file helpers
# Load the single object out of a mixed_by .Rdata file, whatever it was named.
load_obj <- function(f) {
  if (!file.exists(f)) stop("no such file: ", f)
  e <- new.env(); load(f, envir = e); o <- ls(e)
  n <- intersect(c("ddf", "ddq"), o); e[[ if (length(n)) n[1] else o[1] ]]
}
# Most recently modified file matching a regex — avoids hard-coding run dates.
newest <- function(dir, rx) {
  f <- list.files(dir, pattern = rx, full.names = TRUE)
  if (!length(f)) stop("no file matching ", rx, " in ", dir)
  f[order(file.info(f)$mtime, decreasing = TRUE)][1]
}

# match a term by its components, so sex:HCwithin and HCwithin:sex both hit
find_term <- function(ddf, components) {
  cd <- as_tibble(ddf$coef_df_reml) %>% filter(effect == "fixed")
  hit <- vapply(cd$term, function(t)
    setequal(strsplit(t, ":", fixed = TRUE)[[1]], components), logical(1))
  cd[hit, ]
}

.p_level <- function(d) {
  pcol <- if ("padj_fdr_term" %in% names(d) && !all(is.na(d$padj_fdr_term)))
    d$padj_fdr_term else d$p.value
  factor(cut(pcol, breaks = c(-Inf, .001, .01, .05, Inf),
             labels = c("p < .001", "p < .01", "p < .05", "NS")),
         levels = c("NS", "p < .05", "p < .01", "p < .001"))
}

## ------------------------------------------------- 1. term over time (MEDUSA)
# ddf   : a mixed_by result
# term  : term name, or a character vector of its components (order-free)
# facet : split columns to facet on; defaults to whichever of HC_region /
#         dataset are present. Colour is by network when present.
# ylim: pass c(lo, hi) to fix the y range, so panels from the two studies can
# share a scale. Use shared_ylim() to compute a range covering both.
# NB flipy defaults to TRUE here; the neural calls pass flipy = FALSE.
# dodge_width: horizontal offset between networks at each timepoint. Defaults
# to 0, so points sit exactly on their timepoints. The original code used 0.33,
# but the network traces are well separated vertically in these figures, so the
# dodge solved a collision that does not occur while placing points slightly off
# their true x values. Set it above 0 if traces in a particular panel overlap.
plot_term_ts <- function(ddf, term, facet = NULL, ylab = NULL, flipy = TRUE,
                         xlab = "Time relative to trial onset [s]", title = NULL,
                         ylim = NULL, dodge_width = 0) {
  d <- if (length(term) > 1) find_term(ddf, term) else
    as_tibble(ddf$coef_df_reml) %>% filter(effect == "fixed", term == !!term)
  if (!nrow(d)) { warning("term not found: ", paste(term, collapse = ":")); return(invisible(NULL)) }
  
  term_label <- unique(d$term)[1]
  d <- d %>% mutate(t = evt_time, p_level = .p_level(.))
  if ("network" %in% names(d)) d$network <- factor(d$network, levels = c("DMN", "CTR", "LIM"))
  if (is.null(facet)) facet <- intersect(c("HC_region", "dataset"), names(d))
  if (is.null(ylab))  ylab  <- paste0(term_label, " [AU]")
  
  aes_col <- if ("network" %in% names(d)) aes(colour = network, group = network) else aes()
  # When dodging, the error bar layer needs the same `group` aesthetic as the
  # lines and points. Without it the layer falls into a single group,
  # position_dodge has nothing to dodge by, and the bars stay centred while the
  # points shift — which reads as points not centred on their intervals.
  aes_eb <- if ("network" %in% names(d))
    aes(ymin = estimate - std.error, ymax = estimate + std.error, group = network) else
      aes(ymin = estimate - std.error, ymax = estimate + std.error)
  dodge <- position_dodge(width = dodge_width)
  p <- ggplot(d, aes(x = t, y = estimate)) +
    geom_hline(yintercept = 0, lty = "dashed", colour = "#A9A9A9") +
    geom_vline(xintercept = 0, lty = "dashed", colour = "#A9A9A9") +
    geom_line(aes_col, position = dodge, linewidth = 1) +
    geom_point(modifyList(aes_col, aes(size = p_level, alpha = p_level)),
               position = dodge) +
    # error bars last, so they draw over the points as in the published figures
    geom_errorbar(aes_eb, position = dodge, width = 0, colour = "black") +
    scale_colour_manual(values = net_pal, drop = FALSE) +
    # drop = TRUE lists only the significance levels present, as in the
    # published figures
    scale_size_manual(values = p_sizes, drop = TRUE) +
    scale_alpha_manual(values = p_alphas, drop = TRUE) +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
    labs(x = xlab, y = ylab, title = title) + fig_theme()
  if (length(facet)) p <- p + facet_wrap(as.formula(paste("~", paste(facet, collapse = "+"))))
  if (isTRUE(flipy)) {
    p <- p + if (is.null(ylim)) scale_y_reverse() else scale_y_reverse(limits = rev(ylim))
  } else if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  p
}

# Range covering the same term across several models, including error bars, so
# panels from different studies can be put on one scale.
shared_ylim <- function(ddfs, components, pad = 0.04) {
  vals <- unlist(lapply(ddfs, function(d) {
    if (is.null(d)) return(NULL)
    g <- find_term(d, components)
    if (!nrow(g)) return(NULL)
    c(g$estimate - g$std.error, g$estimate + g$std.error)
  }))
  if (!length(vals)) return(NULL)
  r <- range(vals, na.rm = TRUE)
  r + c(-1, 1) * diff(r) * pad
}

## --------------------------------- 2. compatibility wrapper for existing calls
# Keeps the old call sites working: loops over every fixed effect in the model
# and writes one PDF per term, named as before.
plot_mixed_by_vmPFC_HC <- function(ddf, toalign = "clock", toprocess = "network-by-HC",
                                   totest = "", behavmodel = "", model_iter = 1,
                                   hc_LorR = "LR", CTRflag = "FALSE", flipy = "FALSE",
                                   save_dir = ".", dodge_width = 0) {
  if (!grepl("network-by-HC", toprocess))
    stop("this wrapper handles toprocess = 'network-by-HC'; use plot_term_ts() directly ",
         "for other layouts")
  xlab <- if (toalign == "feedback") "Time relative to feedback [s]" else
    "Time relative to trial onset [s]"
  terms <- unique(as_tibble(ddf$coef_df_reml) %>% filter(effect == "fixed") %>% pull(term))
  for (fe in terms) {
    p <- plot_term_ts(ddf, fe, flipy = as.logical(flipy), xlab = xlab,
                      dodge_width = dodge_width)
    if (is.null(p)) next
    fname <- paste0(behavmodel, "-", totest, "_", toalign, "_line_", toprocess, "_",
                    str_replace_all(fe, "[^[:alnum:]]", "_"), "-", hc_LorR, "-",
                    model_iter, ".pdf")
    save_fig(p, fname, save_dir)
  }
  invisible(terms)
}

## ------------------------------------------------- 3. emmeans / emtrends output
# ddf  : a mixed_by result
# spec : name in ddf$emtrends_list or ddf$emmeans_list
# x    : "evt_time" for a time course, or a grouping variable for points
# colour_by / facet : optional grouping
plot_emm <- function(ddf, spec, x = "evt_time", colour_by = NULL, facet = NULL,
                     which = c("emtrends", "emmeans"), ylab = NULL, flipy = TRUE,
                     ci = TRUE, ylim = NULL) {
  which <- match.arg(which)
  lst <- if (which == "emtrends") ddf$emtrends_list else ddf$emmeans_list
  if (is.null(lst[[spec]])) { warning("no ", which, " spec named ", spec); return(invisible(NULL)) }
  d <- as.data.frame(lst[[spec]])
  
  ycol <- grep("\\.trend$", names(d), value = TRUE)
  if (!length(ycol)) ycol <- intersect(c("emmean", "estimate", "prob"), names(d))
  ycol <- ycol[1]
  locol <- intersect(c("asymp.LCL", "lower.CL", "conf.low"), names(d))[1]
  hicol <- intersect(c("asymp.UCL", "upper.CL", "conf.high"), names(d))[1]
  if (is.null(ylab)) ylab <- ycol
  
  # label a numeric sex column if present
  if ("sex" %in% names(d) && is.numeric(d$sex))
    d$Sex <- factor(ifelse(d$sex == 1, "Male", "Female"), levels = c("Female", "Male"))
  if (is.null(colour_by) && "Sex" %in% names(d)) colour_by <- "Sex"
  
  aes_col <- if (!is.null(colour_by))
    aes(colour = .data[[colour_by]], group = .data[[colour_by]]) else aes()
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[ycol]])) +
    geom_hline(yintercept = 0, lty = "dashed", colour = "#A9A9A9")
  
  if (x == "evt_time") {                      # time course: line + ribbon
    p <- p + geom_vline(xintercept = 0, lty = "dashed", colour = "#A9A9A9") +
      geom_line(aes_col, linewidth = 1)
    if (ci && !is.na(locol)) p <- p +
        geom_ribbon(modifyList(aes_col, aes(ymin = .data[[locol]], ymax = .data[[hicol]],
                                            fill = .data[[colour_by]])),
                    alpha = 0.15, colour = NA)
  } else {                                    # by group: points + intervals
    lo <- if (!is.na(locol)) locol else NULL
    p <- p + geom_point(aes_col, size = 4, position = position_dodge(width = 0.4))
    if (ci) p <- p + geom_errorbar(
      modifyList(aes_col, if (!is.null(lo))
        aes(ymin = .data[[lo]], ymax = .data[[hicol]]) else
          aes(ymin = .data[[ycol]] - std.error, ymax = .data[[ycol]] + std.error)),
      width = 0.15, linewidth = 1, position = position_dodge(width = 0.4))
  }
  
  if (!is.null(colour_by) && identical(colour_by, "Sex")) {
    p <- p + scale_colour_manual(values = sex_pal) + scale_fill_manual(values = sex_pal)
  }
  p <- p + labs(x = if (x == "evt_time") "Time relative to trial onset [s]" else NULL,
                y = ylab) + fig_theme()
  if (length(facet)) p <- p + facet_wrap(as.formula(paste("~", paste(facet, collapse = "+"))))
  if (isTRUE(flipy)) {
    p <- p + if (is.null(ylim)) scale_y_reverse() else scale_y_reverse(limits = rev(ylim))
  } else if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  p
}

# Range covering the same term across several models, including error bars, so
# panels from different studies can be put on one scale.
shared_ylim <- function(ddfs, components, pad = 0.04) {
  vals <- unlist(lapply(ddfs, function(d) {
    if (is.null(d)) return(NULL)
    g <- find_term(d, components)
    if (!nrow(g)) return(NULL)
    c(g$estimate - g$std.error, g$estimate + g$std.error)
  }))
  if (!length(vals)) return(NULL)
  r <- range(vals, na.rm = TRUE)
  r + c(-1, 1) * diff(r) * pad
}

## ------------------------------------------------------- 4. Figure 2 (violins)
# Participant-level analogue of the reported model: slope of RT on prior RT at
# mean trial position. Descriptive layer for the violins.
subject_slopes <- function(dat, min_trials = 40) {
  dat %>% select(id, rt_csv_sc, rt_lag_sc, trial_neg_inv_sc) %>% drop_na() %>%
    group_by(id) %>% filter(n() >= min_trials) %>%
    group_modify(~ tibble(slope = coef(lm(rt_csv_sc ~ rt_lag_sc * trial_neg_inv_sc,
                                          data = .x))[["rt_lag_sc"]])) %>% ungroup()
}

# subj : panel, Sex, slope   (one row per participant; needs trial-level data)
# dif  : panel, estimate, conf.low, conf.high  (model-based sex difference)
# Either may be NULL: pass only `dif` to draw the difference panel from saved
# model files alone, or only `subj` for the violins without it. The violins
# cannot be reconstructed from mixed_by output, because the behavioural models
# are random-intercept only and so carry no per-participant slopes.
plot_subject_slopes <- function(subj = NULL, dif = NULL,
                                ylab = "RT-swing slope\n(RT on previous RT)",
                                diff_lab = "Difference\n(95% CI)") {
  stopifnot(!is.null(subj) || !is.null(dif))
  
  p_dif <- NULL
  if (!is.null(dif)) {
    p_dif <- ggplot(mutate(dif, x = 1), aes(x = x, y = estimate)) +
      geom_hline(yintercept = 0, lty = "dashed", colour = "grey60") +
      geom_linerange(aes(ymin = conf.low, ymax = conf.high), linewidth = 1.2) +
      geom_point(size = 3) + facet_wrap(~ panel, nrow = 1) +
      scale_x_continuous(breaks = 1, labels = "Male \u2212 Female", limits = c(0.4, 1.6)) +
      scale_y_reverse() + labs(x = NULL, y = diff_lab) +
      fig_theme(14) + theme(panel.grid.major.x = element_blank())
  }
  if (is.null(subj)) {                      # difference panel alone: keep strips
    return(p_dif + theme(strip.background = element_rect(fill = "grey92", colour = NA),
                         strip.text = element_text(size = 13, face = "bold")))
  }
  
  est <- subj %>% group_by(panel, Sex) %>%
    summarize(se = sd(slope) / sqrt(n()), slope = mean(slope), .groups = "drop") %>%
    mutate(lo = slope - 1.96 * se, hi = slope + 1.96 * se)
  
  p_raw <- ggplot(subj, aes(x = Sex, y = slope)) +
    geom_violin(aes(fill = Sex, colour = Sex), width = 0.85, alpha = 0.20,
                linewidth = 0.4, trim = FALSE, scale = "width") +
    geom_jitter(aes(colour = Sex), width = 0.075, height = 0, size = 1.1, alpha = 0.45) +
    geom_linerange(data = est, aes(ymin = lo, ymax = hi), linewidth = 0.9) +
    geom_point(data = est, size = 2.8) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = sex_pal, guide = "none") +
    scale_fill_manual(values = sex_pal, guide = "none") +
    scale_y_reverse() + labs(x = NULL, y = ylab) +
    fig_theme(14) + theme(panel.grid.major.x = element_blank(),
                          strip.background = element_rect(fill = "grey92", colour = NA),
                          strip.text = element_text(size = 13, face = "bold"))
  if (is.null(p_dif)) return(p_raw)
  
  p_raw / (p_dif + theme(strip.background = element_blank(),
                         strip.text = element_blank())) +
    plot_layout(heights = c(3, 1.3))
}

## ------------------------------- helper: sex difference table from saved models
# Reads the model-2 objects and returns the `dif` data frame for the panel above.
# Warns if either study is missing, since a short table silently drops a facet.
sex_diff_table <- function(d1, d2,
                           panel_lv = c("Experiment 1\n(fMRI)", "Experiment 1\n(MEG)",
                                        "Experiment 2")) {
  s1 <- find_term(d1, c("rt_lag_sc", "sex"))
  s2 <- find_term(d2, c("rt_lag_sc", "sex"))
  for (nm in c("d1", "d2")) {
    d <- get(nm); s <- if (nm == "d1") s1 else s2
    if (!nrow(s)) warning(nm, ": no rt_lag_sc x sex term. Terms present: ",
                          paste(unique(unlist(strsplit(d$coef_df_reml$term, ":"))),
                                collapse = ", "),
                          ". A `sex1M` term means this is a pre-fix file.")
  }
  out <- bind_rows(
    if (nrow(s1)) s1 %>% transmute(panel = paste0("Experiment 1\n(", dataset, ")"),
                                   estimate, conf.low, conf.high, p.value),
    if (nrow(s2)) s2 %>% transmute(panel = "Experiment 2",
                                   estimate, conf.low, conf.high, p.value)
  ) %>% mutate(panel = factor(panel, levels = panel_lv)) %>% arrange(panel)
  if (nrow(out) < length(panel_lv))
    warning("only ", nrow(out), " of ", length(panel_lv), " panels available: ",
            paste(setdiff(panel_lv, as.character(out$panel)), collapse = " | "))
  out
}