# Run before repair, and with BACKUP_DIRECTORY after repair to check installed caches.
root <- here::here()
source(file.path(root, "code/R/correct_urchin_singleton_cache.R"))
definitions <- urchin_repair_definitions(root)
pools <- urchin_repair_pools(root, definitions)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
if (!length(args)) {
  dry_run <- run_urchin_singleton_repair(tempdir(), apply = FALSE, root = root)
  stopifnot(length(dry_run) == 2L)
}
for (program in c("pisco", "lter")) {
  name <- paste0(program, "_urchin_bootstrap.rds")
  installed <- readRDS(file.path(root, "data/cache", name))
  before <- if (length(args)) readRDS(file.path(args[1], paste0("before_", name))) else installed
  plan <- urchin_repair_plan(before, pools[[program]], program, definitions)
  if (length(args)) stopifnot(identical(installed, plan$after))
  again <- urchin_repair_plan(plan$after, pools[[program]], program, definitions)
  stopifnot(identical(again$after, plan$after))
  row <- plan$rows[1]
  bad <- plan$after
  bad$count[row] <- bad$count[row] + 1
  stopifnot(inherits(try(urchin_repair_verify(before, bad, plan$rows,
    plan$changes$after_biomass), silent = TRUE), "try-error"))
}

# One unique size with repeated individuals is not a singleton; NA rows can be.
x <- data.frame(CA_MPA_Name_Short = "test", site_status = "mpa", year = 2000,
                classcode = "MESFRAAD", size = c(70, NA), count = c(7, 3))
stopifnot(urchin_repair_pool_summary(x)$pool_n == 7L)
x$count[1] <- 1
stopifnot(urchin_repair_pool_summary(x)$pool_n == 1L,
          urchin_repair_pool_summary(x)$singleton_size == 70)

# Cross-check deterministic totals against the parent's fixed sampling helper.
expressions <- parse(file.path(root, "code/R/01_utils.R"))
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("bootstrap_biomass"))) eval(expr)
}
for (size in c(24.9, 40, 65, 70, 90, 100)) {
  for (fun in list(definitions$bio_redurch, definitions$bio_purpurch)) {
    result <- bootstrap_biomass(7, 1, data.frame(size = size, count = 1),
                                 fun, n_resamples = 10)
    stopifnot(isTRUE(all.equal(result$biomass, 7 * fun(size))), result$se == 0)
  }
}
cat("PASS: exact 50/26 scope, deterministic fixed-helper agreement, full-object",
    "preservation, unchanged missingness, idempotent calculations, and pool controls\n")
