# 07_incomplete_data.R
# Incomplete contact information: hide 20% or 40% of edges when computing
# centrality and choosing whom to isolate, but simulate spread on the full network.
# The same missing-edge samples are shared by all strategies.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
n <- vcount(g); m <- ecount(g)
tiebreak <- readRDS(file.path(PATHS$proc, "tiebreak.rds"))
sel_full <- readRDS(file.path(PATHS$proc, "selected_sets.rds"))
lambdas  <- readRDS(file.path(PATHS$res, "lambdas.rds"))
set.seed(PARAMS$seed + 7)

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
runs <- list(); stab <- list()
t0 <- Sys.time()
for (miss in PARAMS$miss_levels) {
  for (smp in seq_len(PARAMS$n_miss_samples)) {
    hidden <- sample.int(m, round(miss * m))
    g_obs  <- delete_edges(g, hidden)          # all 75 people kept, retained weights unchanged
    cs_obs <- centrality_scores(g_obs)
    for (k in PARAMS$k_values) {
      for (s in STRATEGIES) {
        iso <- select_top(cs_obs[[s]], k, tiebreak)
        stab[[length(stab) + 1]] <- data.table(miss = miss, sample = smp, k = k, strategy = s,
                                               jaccard_vs_full = jac(iso, sel_full[[s]][[as.character(k)]]))
        for (lv in names(lambdas)) {
          inf <- sim_fixed(W, lambdas[[lv]], PARAMS$gamma, iso, PARAMS$n_runs_miss)
          runs[[length(runs) + 1]] <- data.table(miss = miss, sample = smp, level = lv,
                                                 strategy = s, k = k, mean_ar = mean(inf) / n,
                                                 p_major = mean(inf >= PARAMS$major_cutoff))
        }
      }
      for (lv in names(lambdas)) {   # random isolation does not use the network
        inf <- sim_random(W, lambdas[[lv]], PARAMS$gamma, k, PARAMS$n_runs_miss)
        runs[[length(runs) + 1]] <- data.table(miss = miss, sample = smp, level = lv,
                                               strategy = "random", k = k, mean_ar = mean(inf) / n,
                                               p_major = mean(inf >= PARAMS$major_cutoff))
      }
    }
  }
  message(sprintf("missing level %.0f%% done (%.1f min)", 100 * miss, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
runs <- rbindlist(runs); stab <- rbindlist(stab)
fwrite(runs, file.path(PATHS$res, "incomplete_sample_means.csv"))
fwrite(stab, file.path(PATHS$res, "incomplete_selection_stability.csv"))

# ---- Summaries ------------------------------------------------------------------------------
# Uncertainty: treat each missing-edge sample as one unit; the spread of sample means
# includes both epidemic randomness and variation across missing networks.
summ <- runs[, .(mean_ar = mean(mean_ar), se = sd(mean_ar) / sqrt(.N),
                 p_major = mean(p_major), n_samples = .N), by = .(miss, level, strategy, k)]
summ[, `:=`(lo = mean_ar - 1.96 * se, hi = mean_ar + 1.96 * se)]

# Add complete-data results (miss = 0) from script 06
full <- fread(file.path(PATHS$tab, "T10_removal_summary.csv"))[strategy != "none",
          .(miss = 0, level, strategy, k, mean_ar, se, p_major, n_samples = NA_integer_, lo, hi)]
summ <- rbind(full, summ)
summ[, level := factor(level, names(lambdas))]
setorder(summ, level, k, strategy, miss)
fwrite(summ, file.path(PATHS$tab, "T12_incomplete_summary.csv"))

# Changes and differences with intervals.
# (a) change from complete data: incomplete-data mean (between-sample SE) minus the
#     complete-data mean from script 06 (independent runs), SEs combined in quadrature.
# (b) difference from random isolation, paired by missing-edge sample (both are run within
#     the same sample), SE = sd of the per-sample differences / sqrt(30).
cmp_full <- summ[miss == 0, .(level, strategy, k, full_mean = mean_ar, full_se = se)]
chg <- merge(summ[miss > 0, .(miss, level, strategy, k, mean_ar, se)], cmp_full, by = c("level", "strategy", "k"))
chg[, `:=`(change_vs_complete = mean_ar - full_mean,
           change_lo = mean_ar - full_mean - 1.96 * sqrt(se^2 + full_se^2),
           change_hi = mean_ar - full_mean + 1.96 * sqrt(se^2 + full_se^2))]
pair <- merge(runs[strategy != "random"], runs[strategy == "random", .(miss, sample, level, k, r_ar = mean_ar)],
              by = c("miss", "sample", "level", "k"))
pair <- pair[, .(diff_vs_random = mean(mean_ar - r_ar), diff_se = sd(mean_ar - r_ar) / sqrt(.N)),
             by = .(miss, level, strategy, k)]
chg <- merge(chg, pair, by = c("miss", "level", "strategy", "k"), all.x = TRUE)
chg[, `:=`(diffr_lo = diff_vs_random - 1.96 * diff_se, diffr_hi = diff_vs_random + 1.96 * diff_se)]
chg[, level := factor(level, names(lambdas))]
setorder(chg, miss, level, k, strategy)
fwrite(chg, file.path(PATHS$tab, "T12b_incomplete_changes.csv"))

stab_summ <- stab[, .(mean_jaccard = mean(jaccard_vs_full)), by = .(miss, strategy, k)]
fwrite(stab_summ, file.path(PATHS$tab, "T13_selection_stability.csv"))
print(dcast(stab_summ[, .(miss, k, strategy, j = round(mean_jaccard, 2))], miss + k ~ strategy, value.var = "j"))

# ---- Figures ---------------------------------------------------------------------------------
lab <- sprintf("%s (λ = %s)", names(lambdas), lambdas); names(lab) <- names(lambdas)
pd <- summ[k %in% c(5, 10)]
pd[, strategy := factor(strategy, c("random", STRATEGIES))]
pd[, k_lab := factor(paste("k =", k), paste("k =", c(5, 10)))]
p1 <- ggplot(pd, aes(miss, mean_ar, colour = strategy)) +
  geom_line(linewidth = 0.7) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.25, position = position_dodge(width = 0.04)) +
  facet_grid(k_lab ~ level, labeller = labeller(level = lab), scales = "free_y") +
  scale_colour_manual(values = STRAT_COLS) +
  scale_x_continuous(breaks = c(0, PARAMS$miss_levels), labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Share of contacts hidden when choosing whom to isolate",
       y = "Mean attack rate (95% interval)", colour = "Selection rule",
       title = "Effectiveness when contact information is incomplete")
save_fig(p1, "F9_incomplete_data.png", width = 6.1, height = 3.5)

stab_summ[, strategy := factor(strategy, STRATEGIES)]
p2 <- ggplot(stab_summ, aes(factor(k), mean_jaccard, fill = strategy)) +
  geom_col(position = "dodge") +
  facet_wrap(~ paste0(100 * miss, "% of contacts hidden")) +
  scale_fill_manual(values = STRAT_COLS[STRATEGIES]) +
  labs(x = "Number isolated (k)", y = "Mean Jaccard overlap\nwith complete-data set",
       fill = "Selection rule", title = "How stable is the selected set when contacts are missing?")
save_fig(p2, "F10_selection_stability.png", width = 6, height = 3)
