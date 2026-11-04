# functions.R
# Shared functions: centrality scores, selection, and the weighted SIR simulator.

# ---- Centrality ---------------------------------------------------------------
# Degree counts distinct contacts. Eigenvector uses contact weights.
# Betweenness and harmonic use 1/weight as path length, so heavy contact = short distance.
centrality_scores <- function(g) {
  w <- E(g)$weight
  data.frame(
    name        = V(g)$name,
    role        = V(g)$role,
    degree      = degree(g),
    betweenness = betweenness(g, weights = 1 / w, normalized = TRUE),
    eigenvector = eigen_centrality(g, weights = w)$vector,
    harmonic    = harmonic_centrality(g, weights = 1 / w, normalized = TRUE),
    stringsAsFactors = FALSE
  )
}

# Top-k node indices for one score. Ties are broken by a fixed random vector.
select_top <- function(score, k, tiebreak) {
  if (k == 0) return(integer(0))
  order(-score, tiebreak)[seq_len(k)]
}

# ---- Weighted SIR -------------------------------------------------------------
# W: n x n matrix of daily contact hours. One step = one day.
# A susceptible j escapes each infected i with prob exp(-lambda * w_ij), so
# P(j infected today) = 1 - exp(-lambda * sum_i w_ij).
# People infected at the start of a day transmit, then recover with prob gamma.
# Newly infected people start transmitting the next day.
# isolated: node indices removed before the epidemic (never infected, never transmit).
sir_run <- function(W, lambda, gamma, active, seed_node, return_state = FALSE) {
  n <- nrow(W)
  stopifnot(active[seed_node])        # the first case must be a non-isolated person
  S <- active                         # isolated people are never susceptible
  I <- logical(n)
  S[seed_node] <- FALSE
  I[seed_node] <- TRUE
  ever <- 1L
  days <- 0L
  while (any(I)) {
    days <- days + 1L
    # I is the set infected at the START of the day; only they transmit today.
    pressure <- colSums(W[I, , drop = FALSE])
    new_inf  <- S & (runif(n) < 1 - exp(-lambda * pressure))
    # Recovery is drawn for start-of-day infecteds only, so people infected today
    # cannot recover today and start transmitting tomorrow.
    recover  <- I & (runif(n) < gamma)
    I[recover] <- FALSE
    S[new_inf] <- FALSE
    I[new_inf] <- TRUE
    ever <- ever + sum(new_inf)
  }
  out <- c(infected = ever, days = days)
  if (return_state) return(list(summary = out, ever_infected = active & !S))
  out
}

# Run n_runs epidemics with a fixed isolation set; the first case is drawn
# uniformly from the people who are not isolated.
# sample.int() on the candidate count (not sample(cand)) avoids R's sample(x) pitfall
# when only one candidate is left.
sim_fixed <- function(W, lambda, gamma, isolated, n_runs, return_seeds = FALSE) {
  n <- nrow(W)
  active <- rep(TRUE, n)
  active[isolated] <- FALSE
  cand <- which(active)
  seeds <- cand[sample.int(length(cand), n_runs, replace = TRUE)]
  inf <- vapply(seeds, function(s) sir_run(W, lambda, gamma, active, s)[["infected"]], numeric(1))
  if (return_seeds) return(list(infected = inf, seeds = seeds))
  inf
}

# Random isolation: a new random set of k people for every run.
sim_random <- function(W, lambda, gamma, k, n_runs, return_sets = FALSE) {
  n <- nrow(W)
  sets <- vector("list", n_runs); seeds <- integer(n_runs)
  inf <- vapply(seq_len(n_runs), function(r) {
    iso <- sample.int(n, k)
    active <- rep(TRUE, n)
    active[iso] <- FALSE
    cand <- which(active)
    s <- cand[sample.int(length(cand), 1)]
    if (return_sets) { sets[[r]] <<- iso; seeds[r] <<- s }
    sir_run(W, lambda, gamma, active, s)[["infected"]]
  }, numeric(1))
  if (return_sets) return(list(infected = inf, sets = sets, seeds = seeds))
  inf
}

# Mean with a normal-approximation 95% Monte Carlo interval.
mean_ci <- function(x) {
  m  <- mean(x)
  se <- sd(x) / sqrt(length(x))
  c(mean = m, se = se, lo = m - 1.96 * se, hi = m + 1.96 * se)
}

# ---- Spectral clustering --------------------------------------------------------
# Normalized Laplacian L = I - D^-1/2 W D^-1/2 (Ng, Jordan and Weiss version):
# take the K smallest eigenvectors, normalize rows, then k-means.
norm_laplacian_eigen <- function(W) {
  d <- rowSums(W)
  Dm <- diag(1 / sqrt(d))
  L <- diag(nrow(W)) - Dm %*% W %*% Dm
  e <- eigen(L, symmetric = TRUE)
  ord <- order(e$values)
  list(values = e$values[ord], vectors = e$vectors[, ord])
}

spectral_cluster <- function(eig, K, nstart = 200) {
  U <- eig$vectors[, seq_len(K), drop = FALSE]
  U <- U / sqrt(rowSums(U^2))
  kmeans(U, centers = K, nstart = nstart, iter.max = 100)$cluster
}
