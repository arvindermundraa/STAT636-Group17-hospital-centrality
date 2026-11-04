# run_all.R
# Runs the whole analysis in order (about 5 minutes on a laptop).
# Open STAT636_Group17.Rproj in RStudio first so paths start at the project root.

scripts <- c(
  "R/01_get_data.R",           # download SocioPatterns data, build weighted network
  "R/02_describe_network.R",   # descriptive statistics, role mixing, network figure
  "R/03_centrality.R",         # four centralities, agreement, isolation sets
  "R/04_check_simulator.R",    # sanity checks on the SIR simulator
  "R/05_pilot_lambda.R",       # choose three transmission intensities
  "R/06_removal_experiment.R", # main comparison on complete data
  "R/07_incomplete_data.R",    # 20% / 40% hidden contacts
  "R/08_communities.R",        # spectral communities, ARI, cross-community contact
  "R/09_robustness.R",         # 4 extra replicate seeds, outbreak cutoffs, rewiring extension
  "R/10_slide_figures.R"       # presentation versions of the main figures (from saved tables)
)
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  source(s, echo = FALSE, local = new.env())
  message(sprintf("done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
# Record the R version and package versions used for this run
writeLines(capture.output(sessionInfo()), "output/results/session_info.txt")
