# 00_setup.R
# Packages, folders, and the fixed settings used by every script.
# Run scripts from the project root (open STAT636_Group17.Rproj in RStudio).

needed <- c("igraph", "data.table", "ggplot2", "patchwork", "mclust", "ragg", "systemfonts")
missing_pkgs <- needed[!needed %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

suppressPackageStartupMessages({
  library(igraph)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(mclust)
})

PATHS <- list(
  raw  = "data/raw",
  proc = "data/processed",
  fig  = "output/figures",
  tab  = "output/tables",
  res  = "output/results"
)
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Settings fixed in the proposal
PARAMS <- list(
  seed           = 636017,
  gamma          = 0.2,             # daily recovery probability (mean infectious period 5 days)
  k_values       = c(1, 3, 5, 10),  # number of people isolated
  n_runs         = 3000,            # runs per complete-data condition (proposal started at 500;
                                    # raised to 3000 because 500 could not separate the centralities)
  major_cutoff   = 15,              # major outbreak: at least 15 people ever infected
  major_sens     = c(10, 20),       # sensitivity cutoffs
  miss_levels    = c(0.2, 0.4),     # share of edges hidden when ranking
  n_miss_samples = 30,              # missing-edge samples per level
  n_runs_miss    = 100,             # epidemic runs per sample and condition
  target_ar      = c(low = 0.25, medium = 0.50, high = 0.75), # pilot targets for lambda
  allow_rfid_fallback = FALSE       # stop if sociopatterns.org cannot be reached (see 01_get_data.R)
)

STRATEGIES <- c("degree", "betweenness", "eigenvector", "harmonic")
ROLE_COLS  <- c(ADM = "#8c6bb1", MED = "#e6550d", NUR = "#3182bd", PAT = "#31a354")
STRAT_COLS <- c(random = "grey45", degree = "#1b9e77", betweenness = "#d95f02",
                eigenvector = "#7570b3", harmonic = "#e7298a")

# ---- Fonts ------------------------------------------------------------------------------
# All figure text is Times New Roman, the font of the report and slides. ragg draws with the
# system fonts that systemfonts finds (on macOS: /System/Library/Fonts/Supplemental/Times New Roman.ttf).
# If the font is missing (e.g. a Linux machine without msttcorefonts), we fall back to the default
# serif font and say so; the tables do not depend on the font.
FONT <- "Times New Roman"
if (!FONT %in% systemfonts::system_fonts()$family) {
  message("NOTE: font '", FONT, "' was not found on this system, so figures will use the default ",
          "serif font instead. Install Times New Roman to reproduce the figures exactly ",
          "(tables and results are unaffected).")
  FONT <- "serif"
}

# Report figures are drawn at the width they are printed in the report (1:1), so a base size
# of 12 pt gives axis, legend and strip text of 9.6 pt on paper (at least 9 pt).
theme_set(theme_minimal(base_size = 12, base_family = FONT) +
            theme(panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))
update_geom_defaults("text", list(family = FONT))
update_geom_defaults("label", list(family = FONT))

# width/height in inches = printed size in the report
save_fig <- function(plot, name, width = 6.9, height = 3) {
  ggsave(file.path(PATHS$fig, name), plot, width = width, height = height, dpi = 300, bg = "white",
         device = ragg::agg_png)
}
# Base-graphics and igraph figures: open a ragg PNG device and set the font family
open_png <- function(file, width, height, res = 300, pointsize = 12) {
  ragg::agg_png(file, width = width, height = height, units = "in", res = res, pointsize = pointsize,
                background = "white")
  par(family = FONT)
}

source("R/functions.R")
