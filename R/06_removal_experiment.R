# 06_removal_experiment.R
# Main comparison on the complete network: isolate the top k people by each
# centrality measure (or k random people) and simulate outbreaks.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
n <- vcount(g)
sel <- readRDS(file.path(PATHS$proc, "selected_sets.rds"))
lambdas <- readRDS(file.path(PATHS$res, "lambdas.rds"))
set.seed(PARAMS$seed + 6)

runs <- list()
for (lv in names(lambdas)) {
  lam <- lambdas[[lv]]
  runs[[length(runs) + 1]] <- data.table(level = lv, strategy = "none", k = 0,
    infected = sim_fixed(W, lam, PARAMS$gamma, integer(0), PARAMS$n_runs))
  for (k in PARAMS$k_values) {
    runs[[length(runs) + 1]] <- data.table(level = lv, strategy = "random", k = k,
      infected = sim_random(W, lam, PARAMS$gamma, k, PARAMS$n_runs))
    for (s in STRATEGIES) {
      runs[[length(runs) + 1]] <- data.table(level = lv, strategy = s, k = k,
        infected = sim_fixed(W, lam, PARAMS$gamma, sel[[s]][[as.character(k)]], PARAMS$n_runs))
    }
  }
}
runs <- rbindlist(runs)
runs[, ar := infected / n]
fwrite(runs, file.path(PATHS$res, "removal_runs.csv"))

# ---- Summaries ----------------------------------------------------------------------------
summ <- runs[, {
  ci <- mean_ci(ar)
  pm <- mean(infected >= PARAMS$major_cutoff)
  .(mean_ar = ci[["mean"]], se = ci[["se"]], lo = ci[["lo"]], hi = ci[["hi"]],
    p_major = pm, p_major_se = sqrt(pm * (1 - pm) / .N),
    p_major_10 = mean(infected >= PARAMS$major_sens[1]),
    p_major_20 = mean(infected >= PARAMS$major_sens[2]),
    n_runs = .N)
}, by = .(level, strategy, k)]

# Difference from random isolation at the same k (independent runs)
rand <- summ[strategy == "random", .(level, k, r_mean = mean_ar, r_se = se)]
diff <- merge(summ[!strategy %in% c("random", "none")], rand, by = c("level", "k"))
diff[, `:=`(diff_vs_random = mean_ar - r_mean,
            diff_lo = mean_ar - r_mean - 1.96 * sqrt(se^2 + r_se^2),
            diff_hi = mean_ar - r_mean + 1.96 * sqrt(se^2 + r_se^2))]
diff[, `:=`(r_mean = NULL, r_se = NULL)]

# Difference from degree-based isolation (the simple benchmark)
deg <- summ[strategy == "degree", .(level, k, d_mean = mean_ar, d_se = se)]
diff <- merge(diff, deg, by = c("level", "k"))
diff[, `:=`(diff_vs_degree = mean_ar - d_mean,
            diffd_lo = mean_ar - d_mean - 1.96 * sqrt(se^2 + d_se^2),
            diffd_hi = mean_ar - d_mean + 1.96 * sqrt(se^2 + d_se^2))]
diff[strategy == "degree", c("diff_vs_degree", "diffd_lo", "diffd_hi") := NA]
diff[, `:=`(d_mean = NULL, d_se = NULL)]

lv_order <- names(lambdas)
summ[, level := factor(level, lv_order)]
diff[, level := factor(level, lv_order)]
setorder(summ, level, k, strategy); setorder(diff, level, k, strategy)
fwrite(summ, file.path(PATHS$tab, "T10_removal_summary.csv"))
fwrite(diff, file.path(PATHS$tab, "T11_removal_differences.csv"))
print(dcast(summ[, .(level, k, strategy, ar = round(mean_ar, 3))], level + k ~ strategy, value.var = "ar"))

# ---- Figures -------------------------------------------------------------------------------
lab <- sprintf("%s (λ = %s)", names(lambdas), lambdas)
names(lab) <- names(lambdas)
base <- summ[strategy == "none"]
plot_d <- rbind(summ[strategy != "none"],
                rbindlist(lapply(c("random", STRATEGIES), function(s) copy(base)[, strategy := s])))
plot_d[, strategy := factor(strategy, c("random", STRATEGIES))]

p1 <- ggplot(plot_d, aes(k, mean_ar, colour = strategy)) +
  geom_line(linewidth = 0.7) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.25, position = position_dodge(width = 0.6)) +
  facet_wrap(~ level, labeller = labeller(level = lab), scales = "free_y") +
  scale_colour_manual(values = STRAT_COLS) +
  scale_x_continuous(breaks = c(0, PARAMS$k_values)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Mean attack rate (95% MC interval)",
       colour = "Selection rule", title = "Outbreak size by isolation strategy")
save_fig(p1, "F6_attack_rate_by_strategy.png", width = 6.9, height = 2.8)

p2 <- ggplot(diff, aes(factor(k), diff_vs_random, colour = strategy)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = diff_lo, ymax = diff_hi), position = position_dodge(width = 0.6), size = 0.3) +
  facet_wrap(~ level, labeller = labeller(level = lab)) +
  scale_colour_manual(values = STRAT_COLS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Attack rate minus random isolation",
       colour = "Selection rule", title = "Reduction relative to random isolation (below zero = better)")
save_fig(p2, "F7_difference_vs_random.png", width = 6.9, height = 3)

p3 <- ggplot(plot_d, aes(k, p_major, colour = strategy)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_wrap(~ level, labeller = labeller(level = lab)) +
  scale_colour_manual(values = STRAT_COLS) +
  scale_x_continuous(breaks = c(0, PARAMS$k_values)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "P(at least 15 infected)",
       colour = "Selection rule", title = "Probability of a major outbreak")
save_fig(p3, "F8_major_outbreak.png", width = 6.9, height = 3)
