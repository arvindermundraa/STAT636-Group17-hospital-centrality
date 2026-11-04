# 03_centrality.R
# Compute the four centrality measures, compare their rankings, and fix the
# isolation sets used in the complete-data experiment.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))

cs <- centrality_scores(g)

# One fixed random tie-breaking vector, reused everywhere (reproducible)
set.seed(PARAMS$seed)
tiebreak <- runif(vcount(g))
saveRDS(tiebreak, file.path(PATHS$proc, "tiebreak.rds"))

for (s in STRATEGIES) cs[[paste0("rank_", s)]] <- order(order(-cs[[s]], tiebreak))
fwrite(cs, file.path(PATHS$tab, "T4_centrality_scores.csv"))

# ---- Agreement between measures -------------------------------------------------------
spearman <- cor(cs[, STRATEGIES], method = "spearman")
fwrite(data.table(measure = rownames(spearman), round(spearman, 3)),
       file.path(PATHS$tab, "T5_rank_correlation.csv"))
print(round(spearman, 2))

# Overlap (Jaccard) of the selected sets for each k
sel <- lapply(setNames(STRATEGIES, STRATEGIES), function(s)
  lapply(setNames(PARAMS$k_values, PARAMS$k_values), function(k) select_top(cs[[s]], k, tiebreak)))
saveRDS(sel, file.path(PATHS$proc, "selected_sets.rds"))

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
overlap <- rbindlist(lapply(PARAMS$k_values, function(k) {
  pr <- combn(STRATEGIES, 2)
  data.table(k = k, pair = paste(pr[1, ], pr[2, ], sep = " vs "),
             jaccard = apply(pr, 2, function(p) jac(sel[[p[1]]][[as.character(k)]], sel[[p[2]]][[as.character(k)]])))
}))
fwrite(overlap, file.path(PATHS$tab, "T6_selection_overlap.csv"))

# Who gets selected (top 10 by each measure)
top10 <- rbindlist(lapply(STRATEGIES, function(s) {
  idx <- sel[[s]][["10"]]
  data.table(strategy = s, rank = 1:10, person = cs$name[idx], role = cs$role[idx])
}))
fwrite(top10, file.path(PATHS$tab, "T7_top10_selected.csv"))
print(dcast(top10, rank ~ strategy, value.var = "role"))

# ---- Figures ---------------------------------------------------------------------------
# Scatterplot matrix of ranks (1 = most central)
# Ranks with tied scores sharing their average rank, so the Spearman correlation printed in each
# panel is exactly the one in T5 (the tie-broken selection ranks would give slightly different values,
# because many people have betweenness 0).
rk <- as.data.frame(sapply(STRATEGIES, function(s) rank(-cs[[s]], ties.method = "average")))
# Printed at 4.6 in (Figure A3). pairs() shrinks text by 0.66 for a 4 x 4 grid, so pointsize 15
# gives tick labels of about 9.9 pt.
open_png(file.path(PATHS$fig, "F3_centrality_rank_matrix.png"), width = 4.6, height = 4.6, pointsize = 15)
pairs(rk, pch = 21, bg = ROLE_COLS[cs$role], col = "white", cex = 1.1,
      upper.panel = function(x, y, ...) {
        usr <- par("usr"); on.exit(par(usr = usr)); par(usr = c(0, 1, 0, 1))
        text(0.5, 0.5, sprintf("ρ = %.2f", cor(x, y, method = "spearman")), cex = 1.1)
      },
      main = "Centrality ranks (1 = most central)", cex.main = 1, font.main = 1)
dev.off()

p_ov <- ggplot(overlap, aes(factor(k), pair, fill = jaccard)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#2171b5", limits = c(0, 1)) +
  labs(x = "Number isolated (k)", y = NULL, fill = "Jaccard\noverlap",
       title = "How often two measures pick the same people")
save_fig(p_ov, "F4_selection_overlap.png", width = 6.5, height = 3.4)
