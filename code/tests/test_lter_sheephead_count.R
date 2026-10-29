# Run from any directory: Rscript /path/to/code/tests/test_lter_sheephead_count.R BASELINE_DIR
# Only temporary copies are written. No production export or model fits.
local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply the preserved harmonized baseline directory.")
  baseline <- normalizePath(args[1], mustWork = TRUE)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo, "code/R/export_lter_sheephead_count.R"), envir = env)
  expect_error <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  }
  files <- c("harmonized_raw_responses.csv", "harmonized_response_ratios.csv")
  before_md5 <- tools::md5sum(file.path(baseline, files))
  temporary <- tempfile("lter-count-test-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  targets <- file.path(temporary, c("upstream", "analysis"))
  for (target in targets) {
    dir.create(target)
    stopifnot(all(file.copy(file.path(baseline, files), target)))
  }
  dry <- env$lter_count_export(repo, baseline, targets)
  stopifnot(identical(unname(tools::md5sum(file.path(targets[1], files))), unname(before_md5)),
    identical(dry$summary$target_rows, c(90L, 45L, 90L, 45L)),
    identical(dry$summary$changed_rows, c(23L, 45L, 23L, 45L)),
    identical(dry$summary$non_target_rows_preserved, c(6190L, 3296L, 6190L, 3296L)),
    dry$reconstruction$retained_rows == 912L,
    dry$reconstruction$repeated_individual_rows == 106L)
  reconstructed <- dry$reconstruction
  rr <- merge(reconstructed$old$rr, reconstructed$corrected$rr,
    by = c("CA_MPA_Name_Short", "year"), suffixes = c("_old", "_new"))
  delta <- rr$lnDiff_new - rr$lnDiff_old
  stopifnot(sum(abs(delta) > 1e-10) == 19L,
    sum(abs(delta[rr$CA_MPA_Name_Short == "Campus Point SMCA"]) > 1e-10) == 1L,
    abs(delta[rr$CA_MPA_Name_Short == "Campus Point SMCA" & rr$year == 2020] +
        0.2385135) < 1e-7)
  old <- reconstructed$old$raw
  corrected <- reconstructed$corrected$raw
  campus <- old$CA_MPA_Name_Short == "Campus Point SMCA" & old$status == "mpa"
  stopifnot(identical(old$value[campus], corrected$value[campus]))
  i <- old$CA_MPA_Name_Short == "Campus Point SMCA" & old$status == "reference" & old$year == 2020
  stopifnot(abs(old$value[i] - 2.95515610) < 1e-8,
    abs(corrected$value[i] - 3.75887078) < 1e-8)
  # Exercise the actual module-06 assignment, not just the reconstruction formula.
  assignment <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], quote(lter.fish.sub.site$biomass)),
    as.list(parse(file.path(repo, "code/R/06_lter_processing.R"))))
  stopifnot(length(assignment) == 1L)
  sys.source(file.path(repo, "code/R/00c_analysis_constants.R"), envir = env)
  sys.source(file.path(repo, "code/R/01_utils.R"), envir = env)
  env$lter.fish.sub.site <- data.frame(SIZE = c(25, 25, 25), COUNT = c(1, 2, 0))
  eval(assignment[[1]], env)
  stopifnot(identical(env$lter.fish.sub.site$biomass,
    c(env$bio_sheephead(25), 2 * env$bio_sheephead(25), 0)))
  for (kind in c("raw", "rr")) {
    index <- if (kind == "raw") 1L else 2L
    bad <- reconstructed
    bad$old[[kind]] <- bad$old[[kind]][-1, , drop = FALSE]
    expect_error(env$lter_count_prepare(file.path(targets[1], files[index]), bad, kind), "key mismatch")
  }
  expect_error(env$lter_count_export(repo, baseline, baseline), "Baseline cannot")
  applied <- env$lter_count_export(repo, baseline, targets, apply = TRUE)
  stopifnot(identical(unname(tools::md5sum(file.path(targets[1], files))),
    unname(tools::md5sum(file.path(targets[2], files)))),
    identical(tools::md5sum(file.path(baseline, files)), before_md5))
  for (x in applied$prepared) {
    original <- readLines(file.path(baseline, basename(x$file)))
    output <- readLines(x$file)
    unchanged <- setdiff(seq_along(original), x$changed_rows + 1L)
    stopifnot(identical(original[unchanged], output[unchanged]))
  }
  after_md5 <- tools::md5sum(file.path(targets[1], files))
  expect_error(env$lter_count_export(repo, baseline, targets, apply = TRUE), "differs from preserved baseline")
  expect_error(env$lter_count_prepare(file.path(targets[1], files[1]), reconstructed, "raw"), "already corrected")
  stopifnot(identical(tools::md5sum(file.path(targets[1], files)), after_md5))
  cat("PASS: release fingerprint, original-history concordance, module-06 count assignment,\n",
    "19 annual contrast changes, unchanged Campus MPA values, no production writes during tests,\n",
    "two-target parity, all non-target lines/columns, immutable baseline, missing-key failure,\n",
    "baseline-write refusal and repeat/double-correction refusal.\n", sep = "")
})
