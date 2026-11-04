# 02_describe_network.R
# Descriptive statistics of the ward network: degree, strength, density, components,
# clustering, and mixing by role.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
role <- V(g)$role

# ---- Whole-network summary --------------------------------------------------------
summary_tab <- data.frame(
  measure = c("People", "Edges", "Density", "Connected components",
              "Mean degree", "Median degree", "Mean strength (contact h/day)",
              "Median strength (contact h/day)", "Global clustering (transitivity)",
              "Mean local clustering", "Diameter (unweighted)", "Mean shortest path (unweighted)",
              "Role assortativity (nominal, unweighted)"),
  value = c(vcount(g), ecount(g), edge_density(g), components(g)$no,
            mean(degree(g)), median(degree(g)), mean(strength(g)), median(strength(g)),
            transitivity(g, type = "global"), transitivity(g, type = "average"),
            diameter(g, weights = NA), mean_distance(g, weights = NA),
            assortativity_nominal(g, as.integer(factor(role)), directed = FALSE))
)
summary_tab$value <- signif(summary_tab$value, 4)
fwrite(summary_tab, file.path(PATHS$tab, "T1_network_summary.csv"))
print(summary_tab)

# ---- By role ------------------------------------------------------------------------
node_tab <- data.table(name = V(g)$name, role = role, degree = degree(g), strength = strength(g))
role_tab <- node_tab[, .(n = .N, mean_degree = mean(degree), mean_strength = mean(strength),
                         max_strength = max(strength)), by = role][order(role)]
fwrite(role_tab, file.path(PATHS$tab, "T2_role_summary.csv"))
print(role_tab)

# Mixing by role.
# P[r1, r2] = sum of w_ij over i in r1, j in r2. Row r1 of P sums to the total strength
# (person-level daily contact hours) of role r1, so each row of `share` is the share of
# that role's members' contact time spent with each role ("patients spend x% of their
# contact time with nurses"). A within-role contact counts for both people, as it should.
roles_f <- factor(role, levels = names(ROLE_COLS))
P <- sapply(levels(roles_f), function(r2) sapply(levels(roles_f), function(r1)
  sum(W[roles_f == r1, roles_f == r2])))
stopifnot(isTRUE(all.equal(unname(rowSums(P)), unname(as.numeric(tapply(strength(g), roles_f, sum))))))
share <- sweep(P, 1, rowSums(P), "/")
# Pair-level hours: each contacting pair counted once (within-role sums halved)
M <- P; diag(M) <- diag(M) / 2
fwrite(data.table(from = rownames(M), round(M, 2)), file.path(PATHS$tab, "T3_role_mixing_hours.csv"))
fwrite(data.table(from = rownames(share), round(share, 3)), file.path(PATHS$tab, "T3b_role_mixing_share.csv"))
print(round(share, 2))

# ---- Figures ----------------------------------------------------------------------------
set.seed(PARAMS$seed)
lay <- layout_with_fr(g, weights = E(g)$weight)
# Printed at 2.85 in wide in the report (Figure 1a); no title, the caption describes it
open_png(file.path(PATHS$fig, "F1_network_by_role.png"), width = 2.85, height = 2.65)
par(mai = c(0.28, 0, 0, 0))  # room for a one-line legend underneath
plot(g, layout = lay,
     vertex.color = ROLE_COLS[role], vertex.frame.color = "white",
     vertex.size = 3 + 2.2 * sqrt(strength(g)), vertex.label = NA,
     edge.width = 0.15 + 1.6 * E(g)$weight / max(E(g)$weight),
     edge.color = adjustcolor("grey40", alpha.f = 0.25))
legend("bottom", legend = c("Admin", "Doctors", "Nurses/aides", "Patients"), ncol = 4,
       pt.bg = ROLE_COLS, pch = 21, col = "white", pt.cex = 1.3, cex = 0.8, bty = "n",
       inset = c(0, -0.1), xpd = NA, x.intersp = 0.6, text.width = NA)
dev.off()

p_deg <- ggplot(node_tab, aes(degree, fill = role)) +
  geom_histogram(binwidth = 4, colour = "white") +
  scale_fill_manual(values = ROLE_COLS) +
  labs(x = "Degree (distinct contacts)", y = "People", fill = "Role")
p_str <- ggplot(node_tab, aes(role, strength, colour = role)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.1) +
  scale_colour_manual(values = ROLE_COLS, guide = "none") +
  labs(x = NULL, y = "Strength (contact hours per day)")
save_fig(p_deg + p_str + plot_layout(guides = "collect"), "F2_degree_strength.png", width = 6.3, height = 2.7)

saveRDS(lay, file.path(PATHS$proc, "layout.rds"))
