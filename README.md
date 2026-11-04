# Centrality-Based Epidemic Control in a Hospital Ward Contact Network

STAT 636 (Multivariate Analysis and Statistical Learning for Network Data), Texas A&M University,
Fall 2026, Prof. Jesús Arroyo. **Group 17:** Arvinder Singh Mundra, Hasini Kavuru, Satvika Reddy Akati.

**Question.** When a hospital ward can isolate only a few people, which centrality-based rule for
choosing them most reduces the expected outbreak size, and does it still work when contact information
is incomplete? We compare degree, weighted betweenness, weighted eigenvector and harmonic centrality
with random isolation, using a weighted SIR model simulated on the SocioPatterns hospital ward contact
network (75 people, 4 days of sensor data).

## Key results

All numbers come from the CSV tables in `output/tables/`, written by the scripts below.

- With 3 or more people isolated, every centrality rule beats random isolation at all three
  transmission levels (every 95% interval for the difference excludes zero). With 10 isolated, the
  attack rate falls from 17.5% (random) to 4.8–7.0% at low transmission, from 38.6% to 19.4–24.4% at
  medium and from 61.2% to 49.7–52.8% at high: 37–73% fewer infections at low and medium transmission.
  Isolating one person changes the attack rate by at most 3.1 points.
- The best rule depends on capacity. With 3–5 slots, at low and medium transmission, the contact-time
  measures beat degree by 0.8–3.0 points (15,000 pooled runs); degree spends a pick on an administrator
  with 61 brief contacts. With 10 slots, degree has the lowest attack rate, and harmonic is 2.3–4.6
  points worse.
- Hiding 20% or 40% of contacts changes who is selected (30–92% overlap with the complete-data choice)
  but changes attack rates by only −2.8 to +3.1 points; the rules stay ahead of random in all 72
  conditions with k ≥ 3.
- Spectral clustering finds 3 mixed-role communities (adjusted Rand index with roles 0.17–0.22). In
  degree-preserving rewired networks, outbreaks are much larger and degree is best or tied.

## Documents

The final report, slides and speaker script are in `docs/` of the GitHub repository (the same files
as the course submission):

- [Report (PDF)](docs/17%20-%20Centrality-Based%20Epidemic%20Control%20in%20a%20Hospital%20Ward%20Contact%20Network%20-%20Report.pdf)
- [Slides (PDF)](docs/17%20-%20Centrality-Based%20Epidemic%20Control%20in%20a%20Hospital%20Ward%20Contact%20Network%20-%20Slides.pdf)
- [Speaker script (PDF)](docs/Speaker_Script.pdf)

## How to run

1. Unzip the folder and double-click `STAT636_Group17.Rproj` to open it in RStudio. Scripts use paths
   relative to the project root.
2. Open `run_all.R` and click **Source**, or from a terminal in the project root run
   `Rscript run_all.R`. The full analysis takes about 5 minutes on a laptop (R 4.5.2).
3. Missing packages (igraph, data.table, ggplot2, patchwork, mclust, ragg, systemfonts) are installed
   on the first run. The versions used for the report are listed in `package_versions.txt`, and each run
   records its own in `output/results/session_info.txt`.
4. Figures use Times New Roman, drawn with `ragg`. On macOS and Windows the system copy is found
   automatically. If the font is missing (for example on Linux without `ttf-mscorefonts`), the scripts
   print a note and use the default serif font; the tables are unaffected.

Every random step is seeded from `PARAMS$seed = 636017` (`R/00_setup.R`). A rerun from a clean copy
reproduces every table exactly.

## Data source

`R/01_get_data.R` downloads the contact list from the SocioPatterns dataset page
https://sociopatterns.org/datasets/hospital-ward-dynamic-contact-network/
(file https://sociopatterns.org/assets/data/hospital_lyon_contacts.dat.gz; Vanhems et al., 2013,
*PLOS ONE* 8(9): e73970). It then checks the file against the published counts (32,424 contact
records, 75 people: 8 ADM, 11 MED, 27 NUR, 29 PAT) and against the MD5 checksum of the file used for
the report. `output/results/data_provenance.txt` records the URL, checksum, retrieval time,
observation period (4.022 days) and edge count (1,139).

If sociopatterns.org cannot be reached, the script stops. The CRAN copy `igraphdata::rfid` has
different node IDs, which changes the tie-breaking and all random draws, so it is used only if you
set `PARAMS$allow_rfid_fallback = TRUE`.

The data are distributed by the SocioPatterns collaboration (www.sociopatterns.org) under a Creative
Commons Attribution-NonCommercial-ShareAlike 3.0 licence.

## Scripts and outputs

| Script | What it does | Outputs |
|---|---|---|
| `R/00_setup.R` | Packages, folders, all fixed settings (`PARAMS`), Times New Roman plot theme and PNG devices | |
| `R/functions.R` | Centrality scores, top-k selection, weighted SIR simulator, spectral clustering | |
| `R/01_get_data.R` | Downloads and checks the data, builds the weighted network (weight = contact hours per day) | `data/raw/`, `data/processed/hospital_graph.rds`, `nodes.csv`, `edges.csv`, `output/results/data_provenance.txt` |
| `R/02_describe_network.R` | Density, components, degree and strength, clustering, role mixing | T1 network summary, T2 by role, T3 pair-level role-mixing hours, T3b share of each role's contact time by partner role; F1 network, F2 degree/strength |
| `R/03_centrality.R` | Degree, weighted betweenness, weighted eigenvector, harmonic; rank correlations; top-k sets; tie-break vector | T4 scores and ranks, T5 Spearman, T6 Jaccard overlap, T7 top 10; F3 rank matrix, F4 overlap |
| `R/04_check_simulator.R` | Eight checks of the SIR simulator against exact answers (stops if any fails) | T0 |
| `R/05_pilot_lambda.R` | No-isolation pilot; fixes low/medium/high λ (target attack rates 25/50/75%) before any comparison | T8 pilot grid, T9 chosen λ values, `output/results/lambdas.rds`; F5 |
| `R/06_removal_experiment.R` | Main comparison: k = 0, 1, 3, 5, 10, 3,000 runs per condition | T10 summary, T11 differences from random and degree; `output/results/removal_runs.csv`; F6, F7, F8 |
| `R/07_incomplete_data.R` | 20% or 40% of contacts hidden when ranking; 30 shared samples x 100 runs | T12 summary, T12b change from complete data and paired difference from random, T13 selection stability; F9, F10 |
| `R/08_communities.R` | Normalized-Laplacian spectral clustering, eigengap, ARI vs roles, cross-community contact | T14 ARI, T15 community x role, T16 node communities, T17 cross-community share of selected people; F11, F12 |
| `R/09_robustness.R` | 4 extra independent seeds (pooled 15,000 runs), major-outbreak cutoffs 10/15/20, degree-preserving rewiring extension | T18 pooled means, T19 pooled differences from degree, T20 cutoffs, T21 rewiring comparison, T21b rewiring diagnostics; F13, F14, F15 |
| `R/10_slide_figures.R` | Presentation versions of the main figures, redrawn from the saved tables with slide-size fonts (no new simulation) | `output/figures/slides/S1`–`S7` |

Tables are written to `output/tables/` (T*.csv) and figures to `output/figures/` (F*.png, 300 dpi,
drawn at the size they are printed in the report so that the smallest text is at least 9 pt).
Report figures: Figure 1 = F1 + F12, Figure 2 = F6, Figure 3 = F13, Figure 4 = F9, Figure 5 = F15,
Table 1 = T10. Appendix: A1 = F2, A2 = F5, A3 = F3, A4 = F7, A5 = F8, A6 = F14, A7 = F10, A8 = F11,
Table A1 = T0.

## Model, as specified in the proposal

- **Network:** 75 people. An edge means at least one recorded contact. The weight is total contact
  hours (20 s per record) divided by the 4.022 days covered by the records.
- **Centrality:** degree is unweighted. Eigenvector uses the weights as strengths. Betweenness and
  harmonic use 1/weight as path length. Rankings are computed once per observed network, and ties
  are broken by one fixed random vector (`data/processed/tiebreak.rds`).
- **SIR:** one step is one day. P(infection of j today) = 1 − exp(−λ Σ_i w_ij) over people i
  infected at the start of the day, which is the same as independent per-edge probabilities
  1 − exp(−λ w_ij). An infected person transmits and then recovers with probability 0.2.
  A newly infected person transmits from the next day. Isolated people are never infected and never
  transmit. The first case is drawn uniformly from the people not isolated.
- **Outcomes:** attack rate = everyone ever infected (including the first case) / 75. A major
  outbreak is 15 or more infected, with 10 and 20 as sensitivity cutoffs. Intervals are 95% Monte
  Carlo intervals.

## Deviations from the proposal

- The complete-data comparison uses **3,000 runs per condition instead of the planned 500**. The
  proposal allowed more replication when needed. With 500 runs the intervals were too wide to separate
  the centrality measures from one another.
- Added: four extra independent replicates of the comparison for k = 3, 5, 10 (`09_robustness.R`), to
  check which differences between rules are stable.
- Implementation choice not fixed by the proposal: the three λ values target mean attack rates of
  25%, 50% and 75% in the no-isolation pilot.
- Done as the optional extension: degree-preserving rewiring, with the observed weights permuted at
  random over the rewired edges (medium transmission only).
- Not done (extensions in the proposal): adaptive rankings and other missing-data mechanisms.

## Team and work division

| Member | Owned part | Scripts | Report | Slides |
|---|---|---|---|---|
| Arvinder Singh Mundra | Framing, data, network construction, integration | 00–02, 10, `run_all.R` | §1, §2, §5 | 1–3, close of 10 |
| Hasini Kavuru | Centrality measures, SIR simulator, choice of λ, main experiment | 03–06, `functions.R` | §3, §4.1–4.2 | 4–6, open questions on 10 |
| Satvika Reddy Akati | Isolation capacity, incomplete data, communities, robustness, rewiring | 07–09 | §4.3–4.5, Appendix | 7–9 |

## Licence

Code: MIT licence (see `LICENSE` in the repository). Data: **not included** in this repository or the
code archive. `R/01_get_data.R` downloads them from SocioPatterns. They are distributed by the
SocioPatterns collaboration under CC BY-NC-SA 3.0. If you use them, cite Vanhems et al. (2013) and
acknowledge the SocioPatterns collaboration (www.sociopatterns.org).
