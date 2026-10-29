# Rebuild only the satellite handoff; preserve all diver release fingerprints.
# Rscript code/R/export_landsat_correction.R /path/to/analysis-repo
local({
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(jsonlite); library(here)})
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  analysis <- normalizePath(args[1], mustWork = TRUE)
  source(here::here("code/R/01_utils.R"), local = TRUE)
  source(here::here("code/R/landsat_reference_utils.R"), local = TRUE)
  raw_path <- here::here("data/LANDSAT/MPA_Runs_new.csv")
  corrections_path <- here::here("data/landsat_empty_reference_corrections.csv")
  raw <- read.csv(raw_path, check.names = FALSE)
  clean <- remove_empty_landsat_references(raw, read.csv(corrections_path))
  years <- grep("^[0-9]{4}$", names(clean), value = TRUE)
  annual <- clean %>% pivot_longer(all_of(years), names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year)) %>% group_by(MPA, status, year) %>%
    summarise(value = if (all(!is.finite(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = status, values_from = value) %>%
    filter(is.finite(mpa), is.finite(reference)) %>%
    pivot_longer(c(mpa, reference), names_to = "status", values_to = "value") %>%
    mutate(CA_MPA_Name_Short = MPA, y = "Macrocystis pyrifera")
  annual <- calculate_proportions(annual, value_col = "value", correction_method = "adaptive")
  revised <- annual %>% select(CA_MPA_Name_Short, year, status, PropCorr) %>%
    pivot_wider(names_from = status, values_from = PropCorr) %>%
    mutate(Diff = mpa/reference, lnDiff = log(Diff))
  dest <- file.path(analysis, "data/harmonized")
  old_path <- file.path(dest, "harmonized_landsat_rr.csv")
  index_path <- file.path(dest, "harmonized_landsat_annual_index.csv")
  marker_path <- file.path(dest, "landsat_reference_correction.json")
  previous <- if (file.exists(marker_path)) fromJSON(marker_path, simplifyVector = FALSE) else NULL
  if (!is.null(previous) && (!isTRUE(previous$applied) ||
      previous$source_md5 != unname(tools::md5sum(raw_path)) ||
      previous$registry_md5 != unname(tools::md5sum(corrections_path)))) {
    stop("Source or correction registry changed; review before rebuilding.")
  }
  old <- read.csv(old_path)
  key <- function(x) paste(x$CA_MPA_Name_Short, x$year, sep = "|")
  stopifnot(!anyDuplicated(key(old)), !anyDuplicated(key(revised)),
            all(key(revised) %in% key(old)))
  idx <- match(key(revised), key(old))
  updated <- old[idx, , drop = FALSE]
  changed_columns <- c("mpa", "reference", "Diff", "lnDiff")
  updated[changed_columns] <- revised[changed_columns]
  affected <- read.csv(corrections_path)$MPA
  unchanged <- !updated$CA_MPA_Name_Short %in% affected
  stopifnot(max(abs(as.matrix(updated[unchanged, changed_columns]) -
                   as.matrix(old[idx[unchanged], changed_columns]))) < 1e-10)
  index <- annual %>% transmute(CA_MPA_Name_Short, year, status, index = Prop,
                                log_offset = PropCorr - Prop)
  check <- index %>% mutate(value = index + log_offset) %>%
    select(CA_MPA_Name_Short, year, status, value) %>% pivot_wider(names_from=status, values_from=value)
  stopifnot(max(abs(check$mpa - revised$mpa)) < 1e-10,
            max(abs(check$reference - revised$reference)) < 1e-10)
  before <- list(rr = unname(tools::md5sum(old_path)), index = unname(tools::md5sum(index_path)))
  diver_names <- c("harmonized_raw_responses.csv", "harmonized_response_ratios.csv", "harmonized_site_metadata.csv")
  diver_before <- tools::md5sum(file.path(dest, diver_names))
  write.csv(updated, old_path, row.names = FALSE)
  write.csv(index, index_path, row.names = FALSE)
  marker <- list(release = "landsat_empty_reference_slots_2026_10_16", applied = TRUE,
    removed_slots = read.csv(corrections_path), source_md5 = unname(tools::md5sum(raw_path)),
    registry_md5 = unname(tools::md5sum(corrections_path)),
    source_rows = nrow(raw), retained_source_rows = nrow(clean), affected_mpas = affected,
    old_rr_rows = nrow(old), new_rr_rows = nrow(updated), old_md5 = before,
    files = setNames(lapply(c(old_path, index_path), function(p) list(md5 = unname(tools::md5sum(p)),
          rows = nrow(read.csv(p)))), basename(c(old_path, index_path))),
    diver_inputs_unchanged = identical(diver_before, tools::md5sum(file.path(dest, diver_names))))
  stopifnot(marker$diver_inputs_unchanged)
  if (!is.null(previous)) {
    marker$old_rr_rows <- previous$old_rr_rows
    marker$old_md5 <- previous$old_md5
  }
  write_json(marker, marker_path, auto_unbox = TRUE, pretty = TRUE)
  csvs <- list.files(dest, pattern = "\\.csv$", full.names = TRUE)
  writeLines(paste(tools::md5sum(csvs), basename(csvs)), file.path(dest, "checksums.md5"))
  # Mirror the corrected satellite handoff, not the independently revised diver inputs.
  data_out <- here::here("output/harmonized")
  dir.create(data_out, recursive = TRUE, showWarnings = FALSE)
  file.copy(c(old_path, index_path, marker_path), data_out, overwrite = TRUE)
  cat("Removed six empty reference slots; Landsat pairs:", nrow(old), "->", nrow(updated),
      "; diver inputs unchanged.\n")
})
