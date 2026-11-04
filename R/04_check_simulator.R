# 04_check_simulator.R
# Checks of the SIR simulator against exact answers before using it for comparisons.
# Every check compares simulation output with a value derived independently of the
# simulator (closed-form or brute-force enumeration), so none of them is a tautology.

source("R/00_setup.R")
set.seed(PARAMS$seed + 4)
results <- list(); details <- list()

# 1. No transmission: only the first case is ever infected
W3 <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3)
x <- replicate(200, sir_run(W3, lambda = 0, gamma = 0.2, active = rep(TRUE, 3), seed_node = 1)[["infected"]])
results$no_transmission <- all(x == 1)
details$no_transmission <- "lambda = 0 gives exactly one case"

# 2. One-day infection probability on a single edge matches 1 - exp(-lambda * w)
W2 <- matrix(c(0, 2, 2, 0), 2); lam <- 0.3
# with gamma = 1 the first case transmits for exactly one day
emp <- mean(replicate(1e5, sir_run(W2, lam, gamma = 1, c(TRUE, TRUE), 1)[["infected"]] == 2))
results$edge_probability <- abs(emp - (1 - exp(-lam * 2))) < 0.01
details$edge_probability <- sprintf("empirical %.4f vs exact %.4f", emp, 1 - exp(-lam * 2))

# 3. Infectious period: mean number of days a case stays infected is 1/gamma = 5
dur <- replicate(5000, sir_run(matrix(0, 1, 1), lambda = 1, gamma = 0.2, active = TRUE, seed_node = 1)[["days"]])
results$mean_infectious_days <- abs(mean(dur) - 5) < 0.2
details$mean_infectious_days <- sprintf("mean infectious days %.2f (expected 5)", mean(dur))

# 4. Isolated people are never infected (checked on the infection state of the isolated
#    person itself, not only on the total count)
W4 <- matrix(5, 4, 4); diag(W4) <- 0
act <- c(TRUE, TRUE, FALSE, TRUE)
st <- replicate(300, sir_run(W4, lambda = 2, gamma = 0.2, active = act, seed_node = 1,
                             return_state = TRUE)$ever_infected)
results$isolation_respected <- !any(st[3, ]) && all(st[c(1, 2, 4), ])   # others all infected (p ~ 1)
details$isolation_respected <- sprintf("isolated node infected in %d of 300 runs; other three infected in %d of 300",
                                       sum(st[3, ]), sum(colSums(st[c(1, 2, 4), ]) == 3))

# 5. Two-person chain: probability the second person is ever infected
#    = sum_t (1-g)^(t-1) * p * (1-p)^(t-1) = p / (1 - (1-g)(1-p))
p <- 1 - exp(-lam * 2); gam <- 0.2
theory <- p / (1 - (1 - gam) * (1 - p))
sim <- mean(replicate(20000, sir_run(W2, lam, gam, c(TRUE, TRUE), 1)[["infected"]] == 2))
results$two_person_final_size <- abs(sim - theory) < 0.015
details$two_person_final_size <- sprintf("P(second infected) simulated %.4f vs exact %.4f", sim, theory)

# 6. Final-size distribution on a 4-person weighted graph against brute-force enumeration
#    of the exact Markov chain (all subsets of new infections and recoveries each day).
exact_final_size <- function(W, lambda, gamma, active, seed_node) {
  n <- nrow(W); memo <- new.env()
  key <- function(S, I) paste(c(as.integer(S), as.integer(I)), collapse = "")
  subsets <- function(v) {  # all logical sub-vectors of the TRUE positions of v
    idx <- which(v); out <- list()
    for (m in 0:(2^length(idx) - 1)) {
      s <- logical(length(v)); s[idx[bitwAnd(m, 2^(seq_along(idx) - 1)) > 0]] <- TRUE
      out[[length(out) + 1]] <- s
    }
    out
  }
  dist <- function(S, I) {   # distribution of the number never infected among active people
    k <- key(S, I)
    if (!is.null(memo[[k]])) return(memo[[k]])
    if (!any(I)) { d <- numeric(n + 1); d[sum(S) + 1] <- 1; memo[[k]] <- d; return(d) }
    q <- 1 - exp(-lambda * colSums(W[I, , drop = FALSE]))
    acc <- numeric(n + 1); p_self <- 0
    for (A in subsets(S)) for (B in subsets(I)) {
      pr <- prod(q[A]) * prod(1 - q[S & !A]) * gamma^sum(B) * (1 - gamma)^(sum(I) - sum(B))
      if (pr == 0) next
      S2 <- S & !A; I2 <- (I & !B) | A
      if (identical(S2, S) && identical(I2, I)) { p_self <- p_self + pr; next }
      acc <- acc + pr * dist(S2, I2)
    }
    d <- acc / (1 - p_self); memo[[k]] <- d; d
  }
  S0 <- active; S0[seed_node] <- FALSE; I0 <- logical(n); I0[seed_node] <- TRUE
  d <- dist(S0, I0)                        # P(number never infected = s), s = 0..n
  infected <- sum(active) - (0:n)
  ok <- infected >= 1
  setNames(d[ok], infected[ok])[order(infected[ok])]
}
W6 <- matrix(0, 4, 4)
W6[1, 2] <- 1.5; W6[2, 3] <- 0.4; W6[3, 4] <- 2.0; W6[1, 3] <- 0.7
W6 <- W6 + t(W6)
ex <- exact_final_size(W6, lambda = 0.8, gamma = 0.3, active = rep(TRUE, 4), seed_node = 1)
sm <- replicate(40000, sir_run(W6, 0.8, 0.3, rep(TRUE, 4), 1)[["infected"]])
emp6 <- sapply(as.integer(names(ex)), function(v) mean(sm == v))
se_mean <- sd(sm) / sqrt(length(sm))
# the enumerator itself reproduces the closed-form two-person answer of check 5
ex2 <- exact_final_size(W2, lam, gam, c(TRUE, TRUE), 1)
results$exact_final_size_4node <- abs(ex2[["2"]] - theory) < 1e-10 &&
  max(abs(emp6 - ex)) < 0.012 &&
  abs(mean(sm) - sum(as.integer(names(ex)) * ex)) < 4 * se_mean
details$exact_final_size_4node <- sprintf(
  "mean final size simulated %.4f vs exact %.4f; max |P sim - P exact| over sizes 1-4 = %.4f; enumerator reproduces check 5 (%.4f)",
  mean(sm), sum(as.integer(names(ex)) * ex), max(abs(emp6 - ex)), ex2[["2"]])

# 7. Only one person left un-isolated: sample() edge case. sim_fixed must always seed that
#    person (sample(x) with length(x) == 1 would draw from 1:x instead), and nobody else.
iso_all_but_3 <- c(1, 2, 4)
y <- sim_fixed(W4, lambda = 2, gamma = 0.2, isolated = iso_all_but_3, n_runs = 200)
results$single_candidate_seed <- all(y == 1)
details$single_candidate_seed <- "with 3 of 4 isolated, every run has exactly one case (the remaining person)"

# 8. sim_fixed draws the first case uniformly over the non-isolated people (chi-square
#    test), and sim_random draws a fresh isolation set every run with the first case
#    always outside that run's set.
W8 <- matrix(0.1, 10, 10); diag(W8) <- 0
fx <- sim_fixed(W8, lambda = 0.1, gamma = 0.5, isolated = c(2, 5), n_runs = 6000, return_seeds = TRUE)
cand <- setdiff(1:10, c(2, 5))
chi <- chisq.test(table(factor(fx$seeds, levels = cand)))
rd <- sim_random(W8, lambda = 0.1, gamma = 0.5, k = 3, n_runs = 500, return_sets = TRUE)
n_sets <- length(unique(vapply(rd$sets, function(s) paste(sort(s), collapse = "-"), "")))
seed_outside <- all(mapply(function(s, iso) !(s %in% iso), rd$seeds, rd$sets))
results$uniform_seed_and_random_redraw <- all(fx$seeds %in% cand) && chi$p.value > 0.001 &&
  n_sets > 50 && seed_outside
details$uniform_seed_and_random_redraw <- sprintf(
  "fixed sets: seeds never isolated, chi-square p = %.3f; random sets: %d distinct sets in 500 runs, seed outside set in all runs",
  chi$p.value, n_sets)

out <- data.frame(check = names(results), passed = unlist(results), detail = unlist(details[names(results)]))
print(out)
fwrite(out, file.path(PATHS$tab, "T0_simulator_checks.csv"))
if (!all(out$passed)) stop("Simulator check failed")
