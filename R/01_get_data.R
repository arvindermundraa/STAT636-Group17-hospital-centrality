# 01_get_data.R
# Download the SocioPatterns hospital ward contacts and build the weighted network.
# Output: data/processed/hospital_graph.rds, nodes.csv, edges.csv, data_provenance.txt

source("R/00_setup.R")

# ---- 1. Get the raw contact list ------------------------------------------------
# Primary source: SocioPatterns (Vanhems et al., 2013).
# File columns: t (end of 20-second interval, seconds), i, j, role of i, role of j.
dataset_pages <- c(
  "https://sociopatterns.org/datasets/hospital-ward-dynamic-contact-network/"
)
known_urls <- c(
  "https://sociopatterns.org/assets/data/hospital_lyon_contacts.dat.gz",   # current link (Nov 2026)
  "http://www.sociopatterns.org/wp-content/uploads/2013/09/detailed_list_of_contacts_Hospital.dat_.gz"
)
raw_file <- file.path(PATHS$raw, "hospital_lyon_contacts.dat.gz")
url_file <- file.path(PATHS$raw, "hospital_lyon_contacts.source_url.txt")  # records where raw_file came from
expected_md5 <- "c8c5534a4db1f24db5a86fbd3d1403ff"   # file downloaded from sociopatterns.org on 2026-11-03

# Read the dataset page and pick up the link to the contact list, in case it moved
links_from_pages <- function() {
  found <- character(0)
  for (pg in dataset_pages) {
    html <- tryCatch(suppressWarnings(paste(readLines(pg, warn = FALSE), collapse = " ")), error = function(e) "")
    hrefs <- regmatches(html, gregexpr('href="[^"]+"', html))[[1]]
    hrefs <- gsub('^href="|"$', "", hrefs)
    hrefs <- ifelse(startsWith(hrefs, "/"), paste0("https://sociopatterns.org", hrefs), hrefs)
    found <- c(found, hrefs[grepl("hospital", hrefs, ignore.case = TRUE) & grepl("\\.gz$", hrefs)])
  }
  unique(found)
}

used_url <- NA_character_
download_source <- function() {
  if (file.exists(raw_file) && file.size(raw_file) > 1000) {
    # Reuse a file downloaded earlier, but only with a record of its source URL
    used_url <<- if (file.exists(url_file)) readLines(url_file, n = 1) else "unknown (no source record)"
    return(TRUE)
  }
  src_urls <- unique(c(links_from_pages(), known_urls))
  for (u in src_urls) {
    ok <- tryCatch({
      download.file(u, raw_file, mode = "wb", quiet = TRUE)
      file.size(raw_file) > 1000
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) { used_url <<- u; writeLines(u, url_file); return(TRUE) }
  }
  if (file.exists(raw_file)) file.remove(raw_file)
  FALSE
}

if (download_source()) {
  ct <- read.table(gzfile(raw_file), header = FALSE,
                   col.names = c("t", "i", "j", "Si", "Sj"), stringsAsFactors = FALSE)
  ct <- as.data.table(ct)
  file_md5 <- unname(tools::md5sum(raw_file))
  if (file_md5 != expected_md5)
    warning("Raw file checksum differs from the file used for the report (", file_md5, ").")
  source_used <- paste0("SocioPatterns download: ", used_url, " (md5 ", file_md5, ")")
} else if (!isTRUE(PARAMS$allow_rfid_fallback)) {
  # The report's numbers depend on the SocioPatterns node IDs (node order drives the
  # tie-breaking vector and all random draws), so a silent fallback is not acceptable.
  stop("Could not download the contact list from sociopatterns.org. Check the dataset page ",
       dataset_pages[1], " for the current link, or set PARAMS$allow_rfid_fallback = TRUE ",
       "to use igraphdata::rfid (same contacts, different node IDs, so results will differ).")
} else {
  # Optional fallback: the same data as packaged on CRAN in igraphdata::rfid
  message("sociopatterns.org could not be reached; using igraphdata::rfid (same dataset, different node IDs).")
  rfid <- NULL
  if (requireNamespace("igraphdata", quietly = TRUE)) {
    data("rfid", package = "igraphdata", envir = environment())
  } else if (file.exists(file.path(PATHS$raw, "rfid.rda"))) {
    load(file.path(PATHS$raw, "rfid.rda"))
  } else {
    stop("No data source available. Check your internet connection.")
  }
  rfid <- upgrade_graph(rfid)
  el <- as_edgelist(rfid, names = FALSE)
  st <- V(rfid)$Status
  ct <- data.table(t = E(rfid)$Time, i = el[, 1], j = el[, 2],
                   Si = st[el[, 1]], Sj = st[el[, 2]])
  source_used <- "igraphdata::rfid (CRAN copy of the Vanhems et al. 2013 SocioPatterns data)"
}

# ---- 2. Check against the published description ----------------------------------
roles <- unique(rbind(data.table(id = ct$i, role = ct$Si), data.table(id = ct$j, role = ct$Sj)))
stopifnot(!any(duplicated(roles$id)))   # every person has one role
checks <- c(
  contact_records = nrow(ct),
  people          = nrow(roles),
  ADM = sum(roles$role == "ADM"), MED = sum(roles$role == "MED"),
  NUR = sum(roles$role == "NUR"), PAT = sum(roles$role == "PAT")
)
expected <- c(contact_records = 32424, people = 75, ADM = 8, MED = 11, NUR = 27, PAT = 29)
print(rbind(found = checks, expected = expected))
if (!all(checks == expected)) warning("Data do not match the published counts. Check the download.")

# ---- 3. Aggregate into a weighted, undirected network ----------------------------
# Weight = total contact hours / days observed = average daily contact hours.
# Each record covers the 20-second interval [t - 20, t] (SocioPatterns file description),
# so the observed window runs from min(t) - 20 to max(t): its length is max - min + 20 s.
# Each record therefore contributes 20/3600 contact hours.
obs_days <- (max(ct$t) - min(ct$t) + 20) / 86400
stopifnot(!any(ct$i == ct$j), all(ct$t %% 20 == 0))
ct[, `:=`(a = pmin(i, j), b = pmax(i, j))]
stopifnot(!anyDuplicated(ct[, .(t, a, b)]))   # no pair recorded twice in one interval
edges <- ct[, .(n_records = .N), by = .(a, b)]
edges[, hours := n_records * 20 / 3600]
edges[, weight := hours / obs_days]

setorder(roles, id)
roles[, name := as.character(id)]
g <- graph_from_data_frame(
  d = data.frame(from = as.character(edges$a), to = as.character(edges$b),
                 n_records = edges$n_records, hours = edges$hours, weight = edges$weight),
  directed = FALSE,
  vertices = data.frame(name = roles$name, role = roles$role)
)
stopifnot(is_simple(g), vcount(g) == 75)

saveRDS(g, file.path(PATHS$proc, "hospital_graph.rds"))
fwrite(roles[, .(id, role)], file.path(PATHS$proc, "nodes.csv"))
fwrite(edges[, .(a, b, n_records, hours, weight)], file.path(PATHS$proc, "edges.csv"))

writeLines(c(
  paste("Source:", source_used),
  paste("Retrieved:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste("Contact records:", nrow(ct)),
  paste("People:", nrow(roles)),
  sprintf("Observation period: %.3f days", obs_days),
  paste("Network edges (pairs with any contact):", ecount(g)),
  "Edge weight: total contact hours / observation days (average daily contact hours)"
), file.path(PATHS$res, "data_provenance.txt"))

cat(sprintf("\nBuilt network: %d people, %d edges, %.2f days observed.\n", vcount(g), ecount(g), obs_days))
