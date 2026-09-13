#!/usr/bin/env Rscript
# AR11_CollecTRI_fetch_regulons.R -- fetch CollecTRI (literature-curated
# regulatory interactions) and DoRothEA confidence A+B directly from the
# OmniPath REST API, for cross-validating AR06's chromatin-accessibility-
# inferred regulatory links against independent literature evidence
# (AR11b onward).
#
# decoupleR::get_collectri()/get_dorothea() are broken in the ks_1_2_r conda
# env: the main OmniPath server is reachable, but an unrelated species-lookup
# endpoint (omabrowser.org) intermittently 502s, which pushes OmnipathR onto
# its "static table" fallback -- and that fallback hits a real bug
# (unnest_evidences: "argument is of length zero") unrelated to this setup.
# The plain REST /interactions endpoint works fine, so we query it directly
# and reconstruct the signed edge list ourselves.
#
# Sign convention: weight = +1 (activating) if is_stimulation & !is_inhibition,
# -1 (repressive) if is_inhibition & !is_stimulation, NA (ambiguous/dual-mode,
# dropped) otherwise. Uses resource-filtered is_stimulation/is_inhibition, not
# consensus_*, so the sign reflects CollecTRI's/DoRothEA's own curation, not a
# cross-resource consensus that could include other databases' calls.
#
# Requires network access to omnipathdb.org.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/AR11_CollecTRI_fetch_regulons.R
# Self-checkpointing: skips entirely if both output RDS files exist.

suppressPackageStartupMessages({ library(jsonlite); library(dplyr); library(yaml) })

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[AR11] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

RDS_DIR <- "results/RDS"
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
COLLECTRI_PATH <- file.path(RDS_DIR, "AR11_collectri_signed.rds")
DOROTHEA_PATH  <- file.path(RDS_DIR, "AR11_dorothea_ab_signed.rds")

if (file.exists(COLLECTRI_PATH) && file.exists(DOROTHEA_PATH)) {
  cat("[AR11] Already complete (", COLLECTRI_PATH, " and ", DOROTHEA_PATH, " exist). Skipping. Delete to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

fetch_omnipath <- function(query) {
  url <- paste0("https://omnipathdb.org/interactions?", query,
                "&genesymbols=1&organisms=9606&format=json")
  df <- fromJSON(url)
  as_tibble(df)
}

add_sign <- function(df) {
  df %>%
    mutate(weight = case_when(
             is_stimulation & !is_inhibition ~ 1,
             is_inhibition & !is_stimulation ~ -1,
             TRUE ~ NA_real_)) %>%
    filter(!is.na(weight)) %>%
    transmute(source = source_genesymbol, target = target_genesymbol, weight, sources, curation_effort)
}

cat("[AR11] Fetching CollecTRI (literature-curated regulatory interactions) from OmniPath REST API ...\n")
collectri_raw <- fetch_omnipath("resources=CollecTRI&fields=sources,curation_effort")
collectri <- add_sign(collectri_raw) %>% distinct(source, target, .keep_all = TRUE)
cat(sprintf("[AR11] CollecTRI: %d raw rows -> %d signed, deduplicated edges (%d TFs, %d targets)\n",
            nrow(collectri_raw), nrow(collectri), length(unique(collectri$source)), length(unique(collectri$target))))
saveRDS(collectri, COLLECTRI_PATH)

cat("[AR11] Fetching DoRothEA confidence A+B (fallback for TFs uncovered by CollecTRI) ...\n")
dorothea_raw <- fetch_omnipath("datasets=dorothea&dorothea_levels=A,B&fields=sources,curation_effort")
dorothea <- add_sign(dorothea_raw) %>% distinct(source, target, .keep_all = TRUE)
cat(sprintf("[AR11] DoRothEA A+B: %d raw rows -> %d signed, deduplicated edges (%d TFs, %d targets)\n",
            nrow(dorothea_raw), nrow(dorothea), length(unique(dorothea$source)), length(unique(dorothea$target))))
saveRDS(dorothea, DOROTHEA_PATH)

cat("[DONE] AR11_CollecTRI_fetch_regulons complete -- saved", COLLECTRI_PATH, "and", DOROTHEA_PATH, "\n")
