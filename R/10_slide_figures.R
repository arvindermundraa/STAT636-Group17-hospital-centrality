# 10_slide_figures.R
# Presentation versions of the main figures, redrawn from the saved tables (no new simulation).
# Each figure is sized for its box on a 1920 x 1080 slide, at 2x resolution, so that the
# smallest text (axis tick labels) is at least 24 px on the slide.
# Output: output/figures/slides/S*.png

source("R/00_setup.R")
dir.create(file.path(PATHS$fig, "slides"), showWarnings = FALSE)
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
lay <- readRDS(file.path(PATHS$proc, "layout.rds"))
lambdas <- readRDS(file.path(PATHS$res, "lambdas.rds"))
lab <- sprintf("%s (λ = %s)", names(lambdas), lambdas); names(lab) <- names(lambdas)
DPI <- 200
theme_slide <- theme_minimal(base_size = 24, base_family = FONT) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.margin = margin(10, 20, 10, 10))
# box_w, box_h: size of the image box on the slide in px; the PNG is 2x that
save_slide <- function(p, name, box_w, box_h) {
  ggsave(file.path(PATHS$fig, "slides", name), p, width = 2 * box_w / DPI, height = 2 * box_h / DPI,
         dpi = DPI, bg = "white", device = ragg::agg_png)
}
cols <- STRAT_COLS
strat_lab <- c(random = "Random", degree = "Degree", betweenness = "Betweenness",
               eigenvector = "Eigenvector", harmonic = "Harmonic")

# S1: network by role (box 600 x 557)
ragg::agg_png(file.path(PATHS$fig, "slides", "S1_network_by_role.png"), width = 1200, height = 1114, res = 150,
              background = "white")
par(mai = c(0.75, 0, 0, 0), family = FONT)  # legend row underneath
plot(g, layout = lay, vertex.color = ROLE_COLS[V(g)$role], vertex.frame.color = "white",
     vertex.size = 3 + 2.2 * sqrt(strength(g)), vertex.label = NA,
     edge.width = 0.3 + 4 * E(g)$weight / max(E(g)$weight),
     edge.color = adjustcolor("grey40", alpha.f = 0.25))
legend("bottom", legend = c("Admin", "Doctors", "Nurses", "Patients"), pt.bg = ROLE_COLS, ncol = 4,
       pch = 21, col = "white", pt.cex = 3, cex = 2.0, bty = "n", inset = c(0, -0.09), xpd = NA,
       x.intersp = 0.7, text.width = NA)  # 2.0 x 12 pt at 150 dpi, shown at half size = 25 px
dev.off()

# S2: pilot curve (box 820 x 500)
pilot <- fread(file.path(PATHS$tab, "T8_pilot_grid.csv"))
p <- ggplot(pilot, aes(lambda, mean_ar)) +
  geom_point(colour = "grey30", size = 3) +
  geom_vline(xintercept = lambdas, linetype = 2, colour = "#2171b5", linewidth = 1) +
  # labels placed where they cannot overlap each other or the curve
  annotate("text", x = lambdas * 0.9, y = 0.99, label = names(lambdas), angle = 90,
           colour = "#2171b5", size = 8, hjust = 1, vjust = 0) +
  scale_x_log10() + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "λ per contact hour (log scale)", y = "Attack rate, no isolation") + theme_slide
save_slide(p, "S2_pilot_lambda.png", 820, 500)

# S3: main comparison (box 1440 x 600)
summ <- fread(file.path(PATHS$tab, "T10_removal_summary.csv"))
base <- summ[strategy == "none"]
pd <- rbind(summ[strategy != "none"], rbindlist(lapply(names(strat_lab), function(s) copy(base)[, strategy := s])))
pd[, `:=`(strategy = factor(strategy, names(strat_lab)), level = factor(level, names(lambdas)))]
p <- ggplot(pd, aes(k, mean_ar, colour = strategy)) +
  geom_line(linewidth = 1.2) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 1, position = position_dodge(width = 0.6)) +
  facet_wrap(~ level, labeller = labeller(level = lab), scales = "free_y") +
  scale_colour_manual(values = cols, labels = strat_lab) +
  scale_x_continuous(breaks = c(0, PARAMS$k_values)) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Attack rate", colour = NULL) + theme_slide
save_slide(p, "S3_attack_rate.png", 1440, 600)

# S4: weighted measures minus degree, pooled (box 1240 x 520)
st <- fread(file.path(PATHS$tab, "T19_stability_vs_degree.csv"))[strategy %in% STRATEGIES[-1]]
st[, `:=`(strategy = factor(strategy, STRATEGIES[-1]), level = factor(level, names(lambdas)))]
p <- ggplot(st, aes(factor(k), pooled_diff, colour = strategy)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 1) +
  geom_pointrange(aes(ymin = pooled_lo, ymax = pooled_hi), size = 0.7, linewidth = 1.2,
                  position = position_dodge(width = 0.55)) +
  facet_wrap(~ level, labeller = labeller(level = lab)) +
  scale_colour_manual(values = cols, labels = strat_lab) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), round(100 * x), " pts")) +
  labs(x = "Number of people isolated (k)", y = "Minus degree", colour = NULL) + theme_slide
save_slide(p, "S4_vs_degree.png", 1240, 520)

# S5: incomplete data, k = 10 only (box 1040 x 600)
inc <- fread(file.path(PATHS$tab, "T12_incomplete_summary.csv"))[k == 10]
inc[, `:=`(strategy = factor(strategy, names(strat_lab)), level = factor(level, names(lambdas)))]
p <- ggplot(inc, aes(miss, mean_ar, colour = strategy)) +
  geom_line(linewidth = 1.2) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.5, linewidth = 1, position = position_dodge(width = 0.05)) +
  facet_wrap(~ level, labeller = labeller(level = lab)) +
  scale_colour_manual(values = cols, labels = strat_lab) +
  scale_x_continuous(breaks = c(0, PARAMS$miss_levels), labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Contacts hidden when ranking (k = 10)", y = "Attack rate", colour = NULL) + theme_slide
save_slide(p, "S5_incomplete.png", 1040, 600)

# S6: observed vs rewired (box 1056 x 560)
cmp <- fread(file.path(PATHS$tab, "T21_rewiring_comparison.csv"))
cmp[, `:=`(strategy = factor(strategy, names(strat_lab)),
           network = factor(ifelse(network == "Observed ward", "Observed ward", "Rewired (20 networks)"),
                            c("Observed ward", "Rewired (20 networks)")))]
p <- ggplot(cmp, aes(factor(k), mean_ar, colour = strategy)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 1, position = position_dodge(width = 0.6)) +
  facet_wrap(~ network) +
  scale_colour_manual(values = cols, labels = strat_lab) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Attack rate (λ = 0.75)", colour = NULL) + theme_slide
save_slide(p, "S6_rewiring.png", 1056, 560)

# S7: rank scatter, eigenvector vs harmonic and degree vs betweenness (box 560 x 620)
cs <- fread(file.path(PATHS$tab, "T4_centrality_scores.csv"))
ar <- function(v) rank(-v, ties.method = "average")  # tied scores share their average rank (as in F3)
rk <- rbind(data.table(pair = "Degree vs betweenness", x = ar(cs$degree), y = ar(cs$betweenness), role = cs$role),
            data.table(pair = "Eigenvector vs harmonic", x = ar(cs$eigenvector), y = ar(cs$harmonic), role = cs$role))
rho <- c("Degree vs betweenness" = cor(cs$degree, cs$betweenness, method = "spearman"),
         "Eigenvector vs harmonic" = cor(cs$eigenvector, cs$harmonic, method = "spearman"))
rk[, pair := factor(pair, names(rho), sprintf("%s\nρ = %.2f", names(rho), rho))]
p <- ggplot(rk, aes(x, y, colour = role)) + geom_point(size = 3) +
  facet_wrap(~ pair, ncol = 1) + scale_colour_manual(values = ROLE_COLS) +
  scale_x_continuous(breaks = c(1, 25, 50, 75)) + scale_y_continuous(breaks = c(1, 25, 50, 75)) +
  labs(x = "Rank", y = "Rank", colour = NULL) + theme_slide +
  theme(strip.text = element_text(size = 20))
save_slide(p, "S7_rank_pairs.png", 560, 620)
message("slide figures written to ", file.path(PATHS$fig, "slides"))
