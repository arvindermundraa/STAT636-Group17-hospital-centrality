# 05_pilot_lambda.R
# No-isolation pilot runs to choose three transmission intensities (lambda, per contact hour).
# The three values are fixed here, before any strategy is compared.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
n <- vcount(g)
set.seed(PARAMS$seed + 5)

grid <- exp(seq(log(0.01), log(10), length.out = 40))
pilot <- rbindlist(lapply(grid, function(lam) {
  inf <- sim_fixed(W, lam, PARAMS$gamma, integer(0), 1000)
  data.table(lambda = lam, mean_ar = mean(inf) / n, p_major = mean(inf >= PARAMS$major_cutoff))
}))

# Smooth monotone curve on log(lambda), then invert at the target attack rates
fit <- loess(mean_ar ~ log(lambda), data = pilot, span = 0.3)
fine <- exp(seq(log(0.01), log(10), length.out = 2000))
pred <- cummax(pmin(pmax(predict(fit, data.frame(lambda = fine)), 0), 1))
lam_sel <- sapply(PARAMS$target_ar, function(a) fine[which(pred >= a)[1]])
lam_sel <- signif(lam_sel, 2)

# Confirm with fresh runs
confirm <- rbindlist(lapply(names(lam_sel), function(nm) {
  inf <- sim_fixed(W, lam_sel[[nm]], PARAMS$gamma, integer(0), 4000)
  ci <- mean_ci(inf / n)
  data.table(level = nm, lambda = lam_sel[[nm]], target_ar = PARAMS$target_ar[[nm]],
             mean_ar = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]],
             p_major = mean(inf >= PARAMS$major_cutoff),
             mean_ar_if_major = mean(inf[inf >= PARAMS$major_cutoff]) / n)
}))
print(confirm)

saveRDS(setNames(lam_sel, names(PARAMS$target_ar)), file.path(PATHS$res, "lambdas.rds"))
fwrite(pilot, file.path(PATHS$tab, "T8_pilot_grid.csv"))
fwrite(confirm, file.path(PATHS$tab, "T9_chosen_lambdas.csv"))

p <- ggplot(pilot, aes(lambda, mean_ar)) +
  geom_point(colour = "grey40", size = 1.4) +
  geom_line(data = data.frame(lambda = fine, mean_ar = pred), colour = "black") +
  geom_point(data = pilot, aes(y = p_major), colour = "#d95f02", size = 1.2, shape = 17) +
  geom_vline(xintercept = lam_sel, linetype = 2, colour = "#2171b5") +
  annotate("text", x = lam_sel, y = 1.02, label = names(lam_sel), colour = "#2171b5", size = 3.5) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(x = "Transmission intensity λ (per contact hour, log scale)",
       y = "Mean attack rate (circles)\nP(major outbreak) (triangles)",
       title = "No-isolation pilot: choosing three transmission intensities")
save_fig(p, "F5_pilot_lambda.png", width = 5.6, height = 3.2)
