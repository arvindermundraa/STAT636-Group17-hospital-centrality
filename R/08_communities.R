# 08_communities.R
# Normalized Laplacian spectral clustering, eigengap choice of K, sensitivity to
# nearby K, ARI against professional roles, and cross-community contact of the
# people each strategy selects.

source("R/00_setup.R")
g <- readRDS(file.path(PATHS$proc, "hospital_graph.rds"))
W <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
role <- V(g)$role
sel <- readRDS(file.path(PATHS$proc, "selected_sets.rds"))
lay <- readRDS(file.path(PATHS$proc, "layout.rds"))

eig <- norm_laplacian_eigen(W)
ev <- eig$values[1:12]
gaps <- diff(ev)
K_star <- which.max(gaps[2:9]) + 1     # largest gap after the trivial first eigenvalue, K in 2..9
cat("Smallest eigenvalues:", round(ev, 3), "\nEigengap choice K =", K_star, "\n")

set.seed(PARAMS$seed + 8)
K_set <- sort(unique(pmax(2, c(K_star - 1, K_star, K_star + 1, 4))))   # 4 = number of roles
clus <- lapply(setNames(K_set, K_set), function(K) spectral_cluster(eig, K))
ari <- data.table(K = K_set, ARI_vs_roles = sapply(clus, function(cl) mclust::adjustedRandIndex(cl, role)),
                  modularity = sapply(clus, function(cl) modularity(g, cl, weights = E(g)$weight)))
print(ari)
fwrite(ari, file.path(PATHS$tab, "T14_spectral_ari.csv"))

cl <- clus[[as.character(K_star)]]
V(g)$community <- cl
xt <- as.data.frame.matrix(table(community = cl, role = role))
fwrite(data.table(community = rownames(xt), xt), file.path(PATHS$tab, "T15_community_by_role.csv"))
print(xt)

# Cross-community contact: share of each person's daily contact hours spent with
# people outside their own community
cross_share <- sapply(seq_len(vcount(g)), function(i) sum(W[i, cl != cl[i]]) / sum(W[i, ]))
node_comm <- data.table(name = V(g)$name, role = role, community = cl,
                        strength = strength(g), cross_share = cross_share)
fwrite(node_comm, file.path(PATHS$tab, "T16_node_communities.csv"))

sel_cross <- rbindlist(lapply(STRATEGIES, function(s) rbindlist(lapply(PARAMS$k_values, function(k) {
  idx <- sel[[s]][[as.character(k)]]
  data.table(strategy = s, k = k,
             mean_cross_share = mean(cross_share[idx]),
             mean_cross_hours = mean(cross_share[idx] * strength(g)[idx]),
             communities_covered = length(unique(cl[idx])))
}))))
sel_cross <- rbind(sel_cross, data.table(strategy = "everyone", k = NA,
                                         mean_cross_share = mean(cross_share),
                                         mean_cross_hours = mean(cross_share * strength(g)),
                                         communities_covered = K_star))
fwrite(sel_cross, file.path(PATHS$tab, "T17_selected_cross_community.csv"))
print(sel_cross)

# ---- Figures -----------------------------------------------------------------------------
p1 <- ggplot(data.frame(i = seq_along(ev), value = ev), aes(i, value)) +
  geom_point(size = 2) + geom_line() +
  geom_vline(xintercept = K_star + 0.5, linetype = 2, colour = "#d95f02") +
  annotate("text", x = K_star + 0.6, y = max(ev) * 0.95, hjust = 0, colour = "#d95f02",
           label = paste("largest gap: K =", K_star)) +
  scale_x_continuous(breaks = seq_along(ev)) +
  labs(x = "Index", y = "Eigenvalue", title = "Normalized Laplacian: smallest eigenvalues")
save_fig(p1, "F11_eigengap.png", width = 4.8, height = 2.9)

comm_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#1f78b4")
shapes <- c(ADM = "square", MED = "circle", NUR = "circle", PAT = "csquare")
# Printed at 2.85 in wide next to F1 (Figure 1b), same layout, so roles can be read off F1.
# (Role letters inside the nodes would be far below 9 pt at this size, so they are not drawn.)
open_png(file.path(PATHS$fig, "F12_network_communities.png"), width = 2.85, height = 2.65)
par(mai = c(0.28, 0, 0, 0))  # room for a one-line legend underneath
plot(g, layout = lay, vertex.color = comm_cols[cl], vertex.frame.color = "white",
     vertex.size = 3 + 2.2 * sqrt(strength(g)), vertex.label = NA,
     edge.width = 0.15 + 1.6 * E(g)$weight / max(E(g)$weight),
     edge.color = adjustcolor("grey40", alpha.f = 0.2))
legend("bottom", legend = paste("Community", seq_len(K_star)), ncol = K_star,
       pt.bg = comm_cols[seq_len(K_star)], pch = 21, col = "white", pt.cex = 1.3, cex = 0.8, bty = "n",
       inset = c(0, -0.1), xpd = NA, x.intersp = 0.6, text.width = NA)
dev.off()

saveRDS(list(K = K_star, cluster = cl, eigenvalues = ev), file.path(PATHS$res, "communities.rds"))
