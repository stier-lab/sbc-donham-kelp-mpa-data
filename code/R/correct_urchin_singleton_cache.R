# Rscript code/R/correct_urchin_singleton_cache.R BACKUP_DIRECTORY [--apply]
# With no --apply, reconstruct and validate the repair without writing files.

urchin_repair_definitions <- function(root) {
  wanted <- c("EXCLUDED_KFM_SITES", "EXCLUDED_REFERENCE_SITES",
              "bio_redurch", "bio_purpurch")
  env <- new.env(parent = baseenv())
  expressions <- parse(file.path(root, "code/R/01_utils.R"))
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir = env,
                      inherits = FALSE)))
  env
}

urchin_repair_key <- function(x) {
  cols <- c("CA_MPA_Name_Short", "site_status", "year", "y")
  stopifnot(!any(vapply(x[cols], function(v)
    any(grepl("\r", as.character(v)), na.rm = TRUE), logical(1))))
  key <- do.call(paste, c(x[cols], sep = "\r"))
  key[!complete.cases(x[cols])] <- NA_character_
  key
}

urchin_repair_pool_summary <- function(sizes) {
  stopifnot(!anyNA(sizes$count), all(sizes$count >= 0),
            all(sizes$count == floor(sizes$count)))
  sizes$y <- sizes$classcode
  sizes$.key <- urchin_repair_key(sizes)
  groups <- split(seq_len(nrow(sizes)), sizes$.key)
  result <- lapply(names(groups), function(key) {
    rows <- groups[[key]]
    a <- as.numeric(na.omit(rep(sizes$size[rows], sizes$count[rows])))
    data.frame(key = key, pool_n = length(a),
               singleton_size = if (length(a) == 1L) a else NA_real_)
  })
  do.call(rbind, result)
}

urchin_repair_pools <- function(root, definitions) {
  s <- read.csv(file.path(root, "data/ALL_sizefreq_2024.csv"))
  s <- subset(s,
    (campus %in% c("VRG", "UCSB") & classcode %in% c("STRPUR", "MESFRA")) |
    (campus == "KFM" & classcode %in%
       c("Strongylocentrotus purpuratus", "Strongylocentrotus franciscanus")) |
    (campus == "LTER" & classcode %in% c("SFL", "SFS", "SPS", "SPL")))
  s <- s[!s$site %in% c(definitions$EXCLUDED_KFM_SITES,
                        definitions$EXCLUDED_REFERENCE_SITES), ]
  s$classcode[s$classcode %in%
    c("STRPUR", "Strongylocentrotus purpuratus", "SPS", "SPL")] <- "STRPURAD"
  s$classcode[s$classcode %in%
    c("MESFRA", "Strongylocentrotus franciscanus", "SFS", "SFL")] <- "MESFRAAD"
  sites <- read.csv(file.path(root, "data/PISCO/master_site_table_Emilyedit.csv"))
  sites <- sites[!duplicated(sites[c("site", "CA_MPA_Name_Short")]), ]
  s <- merge(s, sites, by = c("site", "site_new", "campus"), all.x = TRUE)
  geo <- read.csv(file.path(root, "data/MBON",
    "SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv"))
  geo <- unique(geo[c("site_id", "CA_MPA_Name_Short", "status")])
  gi <- match(s$site, geo$site_id)
  fill <- !is.na(gi) & is.na(s$CA_MPA_Name_Short)
  s$CA_MPA_Name_Short[fill] <- geo$CA_MPA_Name_Short[gi[fill]]
  fill <- !is.na(gi) & is.na(s$site_status)
  s$site_status[fill] <- geo$status[gi[fill]]
  list(pisco = urchin_repair_pool_summary(subset(s, size >= 25)),
       lter = urchin_repair_pool_summary(s[complete.cases(s$size), ]))
}

urchin_repair_plan <- function(before, pools, program, definitions) {
  matched <- match(urchin_repair_key(before), pools$key)
  rows <- which(before$count > 0 & pools$pool_n[matched] == 1L)
  expected <- c(pisco = 50L, lter = 26L)
  stopifnot(program %in% names(expected), length(rows) == expected[[program]])
  sizes <- pools$singleton_size[matched[rows]]
  stopifnot(all(is.finite(sizes)), all(sizes >= 1),
            all(is.finite(before$biomass[rows])),
            all(before$y[rows] %in% c("MESFRAAD", "STRPURAD")))
  mass <- ifelse(before$y[rows] == "MESFRAAD",
                 definitions$bio_redurch(sizes), definitions$bio_purpurch(sizes))
  after <- before
  after$biomass[rows] <- before$count[rows] * mass
  log <- cbind(data.frame(program = program, cache_row = rows),
               before[rows, setdiff(names(before), "biomass"), drop = FALSE],
               data.frame(pool_n = 1L, singleton_size_mm = sizes,
                          before_biomass = before$biomass[rows],
                          after_biomass = after$biomass[rows]))
  urchin_repair_verify(before, after, rows, after$biomass[rows])
  list(before = before, after = after, rows = rows, changes = log)
}

urchin_repair_verify <- function(before, after, rows, expected) {
  stopifnot(identical(after$biomass[rows], expected))
  # Restoring only the authorized cells must recover the complete original object.
  restored <- after
  restored$biomass[rows] <- before$biomass[rows]
  stopifnot(identical(restored, before),
            identical(is.na(before), is.na(after)))
  invisible(TRUE)
}

run_urchin_singleton_repair <- function(archive, apply = FALSE,
                                       root = here::here()) {
  definitions <- urchin_repair_definitions(root)
  pools <- urchin_repair_pools(root, definitions)
  paths <- file.path(root, "data/cache",
                     paste0(c("pisco", "lter"), "_urchin_bootstrap.rds"))
  names(paths) <- c("pisco", "lter")
  plans <- lapply(names(paths), function(program)
    urchin_repair_plan(readRDS(paths[[program]]), pools[[program]],
                       program, definitions))
  names(plans) <- names(paths)
  log_columns <- unique(unlist(lapply(plans, function(plan) names(plan$changes))))
  changes <- do.call(rbind, lapply(plans, function(plan) {
    log <- plan$changes
    for (column in setdiff(log_columns, names(log))) log[[column]] <- NA
    log[log_columns]
  }))
  rownames(changes) <- NULL
  stopifnot(nrow(changes) == 76L, all(is.finite(changes$after_biomass)))
  print(data.frame(program = names(plans),
    singleton_rows = vapply(plans, function(p) length(p$rows), integer(1)),
    changed_cells = vapply(plans, function(p)
      sum(p$changes$before_biomass != p$changes$after_biomass), integer(1))))
  if (!apply) return(invisible(plans))
  if (!all(changes$before_biomass != changes$after_biomass)) {
    stop("Expected 76 unrepaired cells; refusing a partial or repeated repair")
  }
  backups <- file.path(archive, paste0("before_", basename(paths)))
  csv <- file.path(archive, "urchin_singleton_cache_corrections.csv")
  if (file.exists(csv)) stop("Change CSV already exists; refusing to overwrite: ", csv)
  for (i in seq_along(backups)) {
    if (file.exists(backups[i]) &&
        !identical(readRDS(backups[i]), plans[[i]]$before)) {
      stop("Existing backup differs from current source cache: ", backups[i])
    }
  }
  dir.create(archive, recursive = TRUE, showWarnings = FALSE)
  staged <- vapply(dirname(paths), function(directory)
    tempfile("urchin_singleton_", tmpdir = directory, fileext = ".rds"), character(1))
  on.exit(unlink(staged[file.exists(staged)]), add = TRUE)
  for (i in seq_along(paths)) {
    saveRDS(plans[[i]]$after, staged[i])
    stopifnot(identical(readRDS(staged[i]), plans[[i]]$after))
    if (!file.exists(backups[i])) {
      if (!file.copy(paths[i], backups[i], overwrite = FALSE)) {
        stop("Could not create backup: ", backups[i])
      }
    }
    stopifnot(identical(readRDS(backups[i]), plans[[i]]$before))
  }
  write.csv(changes, csv, row.names = FALSE)
  logged <- read.csv(csv)
  stopifnot(nrow(logged) == 76L,
            isTRUE(all.equal(logged$after_biomass, changes$after_biomass)))
  # Check for concurrent cache changes before replacing either source cache.
  for (i in seq_along(paths)) {
    stopifnot(identical(readRDS(paths[i]), plans[[i]]$before))
  }
  for (i in seq_along(paths)) {
    if (!file.rename(staged[i], paths[i])) stop("Could not install cache: ", paths[i])
    urchin_repair_verify(readRDS(backups[i]), readRDS(paths[i]),
                         plans[[i]]$rows, plans[[i]]$changes$after_biomass)
  }
  cat("PASS: 76 deterministic biomass repairs; all other cells, keys, attributes,",
      "and missingness identical.\nChanges:", csv, "\n")
  invisible(plans)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!(length(args) %in% 1:2) ||
      (length(args) == 2L && args[2] != "--apply")) {
    stop("Usage: Rscript correct_urchin_singleton_cache.R BACKUP_DIRECTORY [--apply]")
  }
  run_urchin_singleton_repair(args[1], apply = length(args) == 2L)
}
