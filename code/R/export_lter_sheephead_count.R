# Focused correction of the 2024 LTER sheephead biomass export; no model fits.
# Rscript code/R/export_lter_sheephead_count.R BASELINE_DIR TARGET_DIR ... [--apply]
# Default is read-only. Source this file to use the functions without exporting.

lter_count_reconstruct <- function(data_repo) {
  fish_file <- file.path(data_repo, "data/LTER/Annual_fish_comb_20240307.csv")
  sites_file <- file.path(data_repo, "data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv")
  if (unname(tools::md5sum(fish_file)) != "ddf0cb97132a35aede43a8057709ffee") {
    stop("Raw fish file does not match the verified knb-lter-sbc.17.39 release.")
  }
  fish <- utils::read.csv(fish_file)
  sites <- utils::read.csv(sites_file)
  sheephead <- fish[fish$SCIENTIFIC_NAME == "Semicossyphus pulcher", , drop = FALSE]
  stopifnot(nrow(fish) == 93811L, nrow(sheephead) == 1499L)
  if (anyDuplicated(sites$site_id[sites$site_id %in% sheephead$SITE])) {
    stop("Sheephead site lookup would multiply source rows.")
  }
  helpers <- new.env(parent = globalenv())
  sys.source(file.path(data_repo, "code/R/00c_analysis_constants.R"), envir = helpers)
  sys.source(file.path(data_repo, "code/R/01_utils.R"), envir = helpers)
  z <- merge(sheephead, sites, by.x = "SITE", by.y = "site_id")
  z <- z[z$CA_MPA_Name_Short %in% c("Campus Point SMCA", "Naples SMCA") &
           !(z$SITE %in% helpers$EXCLUDED_LTER_SITE_CODES) & z$COUNT != -99999, , drop = FALSE]
  if (any(!is.finite(z$COUNT) | z$COUNT < 0 | z$COUNT != floor(z$COUNT)) ||
      any(z$COUNT > 0 & (!is.finite(z$SIZE) | z$SIZE <= 0)) ||
      any(z$AREA != helpers$LTER_FISH_SURVEY_AREA_M2)) {
    stop("Unexpected count, positive-count size, or survey area; review before correction.")
  }
  # Preserve the existing zero-fill rule and all hierarchical averaging choices.
  z$old_biomass <- suppressWarnings(helpers$bio_sheephead(z$SIZE))
  z$old_biomass[is.na(z$old_biomass)] <- 0
  z$new_biomass <- z$COUNT * z$old_biomass
  transects <- dplyr::summarise(dplyr::group_by(z, YEAR, SITE, TRANSECT, status,
    CA_MPA_Name_Short, ChannelIsland, MPA_Start),
    COUNT = sum(COUNT), old_biomass = sum(old_biomass),
    new_biomass = sum(new_biomass), .groups = "drop")
  site_years <- dplyr::summarise(dplyr::group_by(transects, YEAR, SITE, status,
    CA_MPA_Name_Short, ChannelIsland, MPA_Start),
    COUNT = mean(COUNT), old_biomass = mean(old_biomass),
    new_biomass = mean(new_biomass), .groups = "drop")
  annual <- dplyr::summarise(dplyr::group_by(site_years, status, CA_MPA_Name_Short,
    MPA_Start, YEAR), COUNT = mean(COUNT), old_biomass = mean(old_biomass),
    new_biomass = mean(new_biomass), .groups = "drop")
  names(annual)[names(annual) == "YEAR"] <- "year"
  annual$y <- "SPUL"
  build <- function(biomass_col) {
    # Scale across both statuses, before dropping unpaired years, as in module 06.
    props <- helpers$calculate_proportions(annual, biomass_col)
    rr <- tidyr::pivot_wider(props[c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")],
      names_from = "status", values_from = "PropCorr")
    rr <- as.data.frame(helpers$calculate_log_response_ratio(rr))
    values <- annual[c("CA_MPA_Name_Short", "year", "status", biomass_col)]
    names(values)[names(values) == biomass_col] <- "value"
    values$value <- values$value / helpers$LTER_FISH_SURVEY_AREA_M2
    paired <- tidyr::pivot_wider(values, names_from = "status", values_from = "value")
    paired <- paired[stats::complete.cases(paired), , drop = FALSE]
    raw <- as.data.frame(tidyr::pivot_longer(paired, c("mpa", "reference"),
      names_to = "status", values_to = "value"))
    list(rr = rr, raw = raw)
  }
  list(old = build("old_biomass"), corrected = build("new_biomass"),
    retained_rows = nrow(z), repeated_individual_rows = sum(z$COUNT > 1),
    annual = annual, inputs = c(fish_file, sites_file,
      file.path(data_repo, "code/R/00c_analysis_constants.R"),
      file.path(data_repo, "code/R/01_utils.R")))
}

lter_count_prepare <- function(file, reconstruction, kind) {
  stopifnot(kind %in% c("raw", "rr"))
  bytes <- readBin(file, "raw", n = file.info(file)$size)
  # These canonical CSV exports use LF and one physical line per record.
  if (any(bytes == as.raw(13)) || tail(bytes, 1) != as.raw(10)) {
    stop("Expected LF-terminated canonical CSV: ", file)
  }
  lines <- readLines(file, warn = FALSE)
  data <- utils::read.csv(file, check.names = FALSE)
  if (length(lines) != nrow(data) + 1L) stop("Multiline CSV records are not supported: ", file)
  taxon_col <- if (kind == "raw") "taxon_name" else "y"
  required <- c("source", "resp", taxon_col, "CA_MPA_Name_Short", "year")
  stopifnot(all(required %in% names(data)))
  target <- which(data$source == "LTER" & data$resp == "Bio" &
                    data[[taxon_col]] == "Semicossyphus pulcher")
  keys <- c("CA_MPA_Name_Short", "year", if (kind == "raw") "status")
  fields <- if (kind == "raw") "value" else c("mpa", "reference", "Diff", "lnDiff")
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  old <- reconstruction$old[[kind]]
  corrected <- reconstruction$corrected[[kind]]
  expected_keys <- key(old)
  target_keys <- key(data[target, , drop = FALSE])
  if (anyDuplicated(expected_keys) || anyDuplicated(target_keys) ||
      !setequal(expected_keys, target_keys) || !identical(key(corrected), expected_keys)) {
    stop("Target key mismatch or duplicate; refusing to add/drop/reclassify rows: ", file)
  }
  m <- match(target_keys, expected_keys)
  before <- as.matrix(data[target, fields, drop = FALSE])
  expected <- as.matrix(old[m, fields, drop = FALSE])
  after <- as.matrix(corrected[m, fields, drop = FALSE])
  tolerance <- 1e-10
  if (any(!is.finite(c(before, expected, after))) || max(abs(before - expected)) > tolerance) {
    stop("Export does not reproduce the uncorrected canonical history (or already corrected): ", file)
  }
  result <- data
  changed_cells <- abs(after - before) > tolerance
  replacement <- before
  replacement[changed_cells] <- after[changed_cells]
  result[target, fields] <- replacement
  changed_rows <- target[rowSums(changed_cells) > 0]
  # Use the CSV writer for changed records, retaining every other original line.
  for (i in changed_rows) {
    con <- textConnection("record", "w", local = TRUE)
    utils::write.table(result[i, , drop = FALSE], con, sep = ",", row.names = FALSE,
      col.names = FALSE, quote = TRUE, na = "NA", qmethod = "double")
    close(con)
    if (length(record) != 1L) stop("Unexpected multiline replacement record.")
    lines[i + 1L] <- record
  }
  check <- utils::read.csv(text = paste(lines, collapse = "\n"), check.names = FALSE)
  if (!identical(data[-target, , drop = FALSE], check[-target, , drop = FALSE]) ||
      !identical(data[target, setdiff(names(data), fields), drop = FALSE],
                 check[target, setdiff(names(data), fields), drop = FALSE]) ||
      max(abs(as.matrix(check[target, fields, drop = FALSE]) - after)) > tolerance) {
    stop("Non-target values/columns changed or correction failed round-trip: ", file)
  }
  original <- readLines(file, warn = FALSE)
  unchanged_lines <- setdiff(seq_along(lines), changed_rows + 1L)
  stopifnot(identical(lines[unchanged_lines], original[unchanged_lines]))
  list(file = normalizePath(file), before_md5 = unname(tools::md5sum(file)),
    lines = lines, target_rows = length(target), changed_rows = changed_rows,
    changed_cells = sum(changed_cells), non_target_rows = nrow(data) - length(target),
    max_baseline_difference = max(abs(before - expected)), kind = kind)
}

lter_count_export <- function(data_repo, baseline_dir, target_dirs, apply = FALSE) {
  baseline_dir <- normalizePath(baseline_dir, mustWork = TRUE)
  target_dirs <- unique(normalizePath(target_dirs, mustWork = TRUE))
  if (baseline_dir %in% target_dirs) stop("Baseline cannot be an export target.")
  names <- c(raw = "harmonized_raw_responses.csv", rr = "harmonized_response_ratios.csv")
  baseline_files <- file.path(baseline_dir, names)
  baseline_md5 <- unname(tools::md5sum(baseline_files))
  if (anyNA(baseline_md5)) stop("Missing baseline input files.")
  reconstructed <- lter_count_reconstruct(data_repo)
  prepared <- list()
  # Validate EVERY destination before writing ANY destination.
  for (directory in target_dirs) for (kind in names(names)) {
    file <- file.path(directory, names[[kind]])
    if (!identical(unname(tools::md5sum(file)), baseline_md5[match(kind, names(names))])) {
      stop("Destination differs from preserved baseline: ", file)
    }
    prepared[[length(prepared) + 1L]] <- lter_count_prepare(file, reconstructed, kind)
  }
  summary <- do.call(rbind, lapply(prepared, function(x) data.frame(file = x$file,
    target_rows = x$target_rows, changed_rows = length(x$changed_rows),
    changed_cells = x$changed_cells, non_target_rows_preserved = x$non_target_rows,
    max_baseline_difference = x$max_baseline_difference)))
  print(summary, row.names = FALSE)
  cat("Retained raw sheephead rows:", reconstructed$retained_rows,
      "; COUNT > 1 rows:", reconstructed$repeated_individual_rows, "\n")
  if (!apply) {
    cat("DRY RUN: no files written; add --apply only after preserving the baseline.\n")
    return(invisible(list(summary = summary, prepared = prepared, reconstruction = reconstructed)))
  }
  if (!identical(unname(tools::md5sum(baseline_files)), baseline_md5) ||
      any(vapply(prepared, function(x) unname(tools::md5sum(x$file)) != x$before_md5, logical(1)))) {
    stop("Inputs changed during validation; no exports written.")
  }
  staged <- vapply(prepared, function(x) {
    file <- tempfile("lter-count-", tmpdir = dirname(x$file))
    writeLines(x$lines, file, useBytes = TRUE)
    stopifnot(identical(readLines(file, warn = FALSE), x$lines))
    file
  }, character(1))
  on.exit(unlink(staged), add = TRUE)
  for (i in seq_along(prepared)) {
    if (!file.rename(staged[i], prepared[[i]]$file)) {
      stop("Export rename failed: ", prepared[[i]]$file, "; restore from baseline if needed.")
    }
  }
  cat("Applied only LTER sheephead Bio corrections. All non-target lines preserved.\n")
  invisible(list(summary = summary, prepared = prepared, reconstruction = reconstructed))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  apply <- "--apply" %in% args
  args <- args[args != "--apply"]
  if (length(args) < 2L) stop("Usage: Rscript export_lter_sheephead_count.R BASELINE_DIR TARGET_DIR ... [--apply]")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  data_repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
  lter_count_export(data_repo, args[1], args[-1], apply = apply)
}
