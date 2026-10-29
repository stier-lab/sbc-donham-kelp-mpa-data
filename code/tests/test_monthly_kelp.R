source("code/R/01_utils.R")

expect_error <- function(expression) {
  failed <- tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  stopifnot(failed)
}
stopifnot(identical(names(MACRO_MONTHLY_SLOPES), as.character(1:12)))
stopifnot(all.equal(bio_macro(rep(1, 12), 1:12),
                    unname(MACRO_MONTHLY_SLOPES) * 1000, tolerance = 1e-12))
stopifnot(all.equal(bio_macro(c(0, 2), 7), c(0, 204.68), tolerance = 1e-12))
expect_error(bio_macro(1, 0))
expect_error(bio_macro(1, 13))
expect_error(bio_macro(1, 1.5))
expect_error(bio_macro(1, NA))
expect_error(bio_macro(c(1, 2), c(1, 2, 3)))
expect_error(bio_macro(-1, 7))
stopifnot(is.na(bio_macro(NA_real_, 7)))
expect_error(read_macro_monthly_slopes(tempfile()))

set.seed(102)
single <- bootstrap_macro_monthly(2, 7, 1, 1000)
next_draw <- runif(5)
set.seed(102)
for (i in 1:1000) invisible(sample(2, 1, replace = TRUE))
stopifnot(identical(next_draw, runif(5)), single$count == 2,
          abs(single$biomass - 204.68) < 1e-10)

stipes <- c(2, 4, 2)
months <- c(5, 7, 10)
set.seed(205)
paired <- bootstrap_macro_monthly(stipes, months, 3, 1000)
set.seed(205)
expected <- replicate(1000, {
  indices <- sample.int(3, 3, replace = TRUE)
  c(sum(stipes[indices]), sum(bio_macro(stipes[indices], months[indices])))
})
stopifnot(all.equal(unlist(paired), setNames(rowMeans(expected), c("count", "biomass"))))
set.seed(205)
legacy_count <- mean(replicate(1000, sum(sample(stipes, 3, replace = TRUE))))
stopifnot(paired$count == legacy_count)
expect_error(bootstrap_macro_monthly(2, 7, 0, 10))
cat("PASS: monthly coefficients, dates, zero/missing values, paired draws and single-record RNG guards\n")

# Historical source remains reproducible when requested explicitly.
legacy <- read_macro_monthly_slopes(file.path(dirname(here::here()), "Donham-Stier-CA-MPA-2026", "data/external/sbc11215/table5.csv"))
stopifnot(abs(bio_macro(2,7,legacy)-209.9)<1e-10)
