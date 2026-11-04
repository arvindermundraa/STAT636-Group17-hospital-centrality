# 09_robustness.R
# Three checks on the main comparison:
#  A. Replication: repeat the complete-data comparison with 4 more independent seeds
#     and pool all 5 (15,000 runs per condition) to see which differences are stable.
#  B. Major-outbreak cutoff: P(at least 10 / 15 / 20 infected) by strategy.
#  C. Extension from the proposal: degree-preserving rewiring with randomly reassigned
#     weights, to ask whether the gaps between rules come from structure beyond degree.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
n <- vcount(g)
sel <- readRDS(file.path(PATHS$proc, "selected_sets.rds"))
tiebreak <- readRDS(file.path(PATHS$proc, "tiebreak.rds"))
lambdas <- readRDS(file.path(PATHS$res, "lambdas.rds"))
lab <- sprintf("%s (λ = %s)", names(lambdas), lambdas); names(lab) <- names(lambdas)
k_rob <- c(3, 5, 10)

# ---- A. Replication with independent seeds ------------------------------------------------
rep_runs <- list()
for (rep in 1:4) {
  set.seed(PARAMS$seed + 900 + rep)
  for (lv in names(lambdas)) for (k in k_rob) {
    rep_runs[[length(rep_runs) + 1]] <- data.table(rep = rep, level = lv, strategy = "random", k = k,
      infected = sim_random(W, lambdas[[lv]], PARAMS$gamma, k, PARAMS$n_runs))
    for (s in STRATEGIES) rep_runs[[length(rep_runs) + 1]] <- data.table(rep = rep, level = lv,
      strategy = s, k = k,
      infected = sim_fixed(W, lambdas[[lv]], PARAMS$gamma, sel[[s]][[as.character(k)]], PARAMS$n_runs))
  }
  message("replicate ", rep, " done")
}
orig <- fread(file.path(PATHS$res, "removal_runs.csv"))[k %in% k_rob & strategy != "none"]
orig[, rep := 0L]
all_runs <- rbind(orig[, .(rep, level, strategy, k, infected)], rbindlist(rep_runs))
all_runs[, ar := infected / n]

per_rep <- all_runs[, .(mean_ar = mean(ar)), by = .(rep, level, strategy, k)]
per_rep[, rank_in_rep := frank(mean_ar, ties.method = "average"), by = .(rep, level, k)]
pooled <- all_runs[, {ci <- mean_ci(ar); .(mean_ar = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]], n_runs = .N)},
                   by = .(level, strategy, k)]
# How often (out of 5 replicates) each rule beats degree
vs_deg <- merge(per_rep[strategy != "degree"], per_rep[strategy == "degree", .(rep, level, k, deg = mean_ar)],
                by = c("rep", "level", "k"))
stab <- vs_deg[, .(beats_degree_in = sum(mean_ar < deg), mean_diff = mean(mean_ar - deg)),
               by = .(level, k, strategy)]
# Pooled difference from degree with CI
pd <- all_runs[, .(m = mean(ar), v = var(ar) / .N), by = .(level, strategy, k)]
pd <- merge(pd[strategy != "degree"], pd[strategy == "degree", .(level, k, md = m, vd = v)], by = c("level", "k"))
pd[, `:=`(diff = m - md, lo = m - md - 1.96 * sqrt(v + vd), hi = m - md + 1.96 * sqrt(v + vd))]
stab <- merge(stab, pd[, .(level, k, strategy, pooled_diff = diff, pooled_lo = lo, pooled_hi = hi)],
              by = c("level", "k", "strategy"))
stab[, level := factor(level, names(lambdas))]
setorder(stab, level, k, strategy)
fwrite(pooled, file.path(PATHS$tab, "T18_pooled_15000_runs.csv"))
fwrite(stab, file.path(PATHS$tab, "T19_stability_vs_degree.csv"))
print(stab[, .(level, k, strategy, beats_degree_in, pooled_diff = round(pooled_diff, 3),
               lo = round(pooled_lo, 3), hi = round(pooled_hi, 3))])

stab_plot <- stab[strategy != "random"]
stab_plot[, strategy := factor(strategy, STRATEGIES[-1])]
p1 <- ggplot(stab_plot, aes(factor(k), pooled_diff, colour = strategy)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = pooled_lo, ymax = pooled_hi), position = position_dodge(width = 0.5), size = 0.3) +
  facet_wrap(~ level, labeller = labeller(level = lab)) +
  scale_colour_manual(values = STRAT_COLS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Attack rate minus degree-based isolation",
       colour = "Selection rule",
       title = "Weighted measures minus degree, 15,000 pooled runs (below zero = better)")
save_fig(p1, "F13_vs_degree_pooled.png", width = 6.9, height = 2.8)

# ---- B. Major-outbreak cutoff sensitivity ---------------------------------------------------
runs <- fread(file.path(PATHS$res, "removal_runs.csv"))
cut <- rbindlist(lapply(c(10, 15, 20), function(cc)
  runs[, .(p = mean(infected >= cc)), by = .(level, strategy, k)][, cutoff := cc]))
base <- cut[strategy == "none"]
cut <- rbind(cut[strategy != "none"],
             rbindlist(lapply(c("random", STRATEGIES), function(s) copy(base)[, strategy := s])))
cut[, `:=`(strategy = factor(strategy, c("random", STRATEGIES)), level = factor(level, names(lambdas)))]
fwrite(cut, file.path(PATHS$tab, "T20_major_outbreak_cutoffs.csv"))
p2 <- ggplot(cut[level == "medium"], aes(k, p, colour = strategy)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_wrap(~ paste("At least", cutoff, "infected")) +
  scale_colour_manual(values = STRAT_COLS) +
  scale_x_continuous(breaks = c(0, PARAMS$k_values)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Probability", colour = "Selection rule",
       title = sprintf("Major-outbreak probability under three cutoffs (medium, λ = %s)", lambdas[["medium"]]))
save_fig(p2, "F14_major_outbreak_cutoffs.png", width = 6.9, height = 3)

# ---- C. Degree-preserving rewiring (exploratory extension) -----------------------------------
# Rewire edges keeping every person's number of contacts, then reassign the observed
# edge weights at random. This keeps degrees and the weight distribution but largely removes
# role/community structure and the link between degree and contact time.
set.seed(PARAMS$seed + 950)
n_rewire <- 20; runs_rw <- 300; lam <- lambdas[["medium"]]
comm <- readRDS(file.path(PATHS$res, "communities.rds"))$cluster
edge_key <- function(gg) { el <- as_edgelist(gg); paste(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2])) }
rw <- list(); rw_diag <- list()
for (r in seq_len(n_rewire)) {
  gr <- rewire(g, keeping_degseq(niter = 20 * ecount(g)))
  # degree sequence preserved, still a simple graph with the same number of edges
  stopifnot(identical(degree(gr), degree(g)), is_simple(gr), ecount(gr) == ecount(g))
  # Weight-assignment rule: the observed multiset of edge weights is permuted uniformly
  # at random over the rewired edges (old weights attached by rewire() are discarded).
  E(gr)$weight <- sample(E(g)$weight)
  Wr <- as_adjacency_matrix(gr, attr = "weight", sparse = FALSE)
  # How much structure the rewiring removes (no random numbers used here)
  rw_diag[[r]] <- data.table(net = r,
    share_edges_kept = mean(edge_key(gr) %in% edge_key(g)),
    modularity_obs_communities = modularity(gr, comm, weights = E(gr)$weight),
    cor_degree_strength = cor(degree(gr), strength(gr)))
  csr <- centrality_scores(gr)
  for (k in k_rob) {
    for (s in STRATEGIES) {
      inf <- sim_fixed(Wr, lam, PARAMS$gamma, select_top(csr[[s]], k, tiebreak), runs_rw)
      rw[[length(rw) + 1]] <- data.table(net = r, strategy = s, k = k, mean_ar = mean(inf) / n)
    }
    inf <- sim_random(Wr, lam, PARAMS$gamma, k, runs_rw)
    rw[[length(rw) + 1]] <- data.table(net = r, strategy = "random", k = k, mean_ar = mean(inf) / n)
  }
}
rw <- rbindlist(rw)
rw_diag <- rbind(data.table(net = 0, share_edges_kept = 1,
                            modularity_obs_communities = modularity(g, comm, weights = E(g)$weight),
                            cor_degree_strength = cor(degree(g), strength(g))),
                 rbindlist(rw_diag))
fwrite(rw_diag, file.path(PATHS$tab, "T21b_rewiring_diagnostics.csv"))   # net 0 = observed ward
print(rw_diag[, lapply(.SD, mean), by = .(observed = net == 0)])
rw_summ <- rw[, .(mean_ar = mean(mean_ar), se = sd(mean_ar) / sqrt(.N)), by = .(strategy, k)]
rw_summ[, `:=`(lo = mean_ar - 1.96 * se, hi = mean_ar + 1.96 * se, network = "Rewired (degree-preserving)")]
obs <- fread(file.path(PATHS$tab, "T10_removal_summary.csv"))[level == "medium" & k %in% k_rob & strategy != "none",
        .(strategy, k, mean_ar, se, lo, hi, network = "Observed ward")]
cmp <- rbind(obs, rw_summ)
fwrite(cmp, file.path(PATHS$tab, "T21_rewiring_comparison.csv"))
print(dcast(cmp[, .(network, k, strategy, ar = round(mean_ar, 3))], network + k ~ strategy, value.var = "ar"))

cmp[, strategy := factor(strategy, c("random", STRATEGIES))]
p3 <- ggplot(cmp, aes(factor(k), mean_ar, colour = strategy)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.3) +
  facet_wrap(~ network) +
  scale_colour_manual(values = STRAT_COLS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of people isolated (k)", y = "Mean attack rate (95% interval)", colour = "Selection rule",
       title = sprintf("Observed ward vs rewired networks (medium, λ = %s)", lam))
save_fig(p3, "F15_rewiring.png", width = 6.1, height = 2.7)
