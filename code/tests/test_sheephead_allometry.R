env <- new.env(parent = globalenv())
sys.source("code/R/00c_analysis_constants.R", envir = env)
expressions <- as.list(parse("code/R/01_utils.R"))
for (name in c("bio_sheephead", "apply_sheephead_allometry")) {
  x <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
              identical(x[[2]], as.name(name)), expressions)
  stopifnot(length(x) == 1L)
  eval(x[[1]], env)
}
path <- "data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv"
fingerprint <- tools::md5sum(path)
original <- read.csv(path)
adopted <- env$apply_sheephead_allometry(original)
cells <- adopted
attr(cells, "sheephead_override") <- NULL
i <- which(original$pisco_classcode == "SPUL")
stopifnot(identical(cells[-i, ], original[-i, ]),
          adopted$WL_a[i] == 0.0144, adopted$WL_b[i] == 3.04,
          adopted$WL_W_units[i] == "g", adopted$WL_L_units[i] == "cm",
          adopted$WL_input_length[i] == "TL", identical(tools::md5sum(path), fingerprint))
stopifnot(isTRUE(all.equal(env$bio_sheephead(c(10, 30, 50)),
                          0.0144 * c(10, 30, 50)^3.04)))
# Exercise the actual module-04 count-weighted TL assignment.
x <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("for")) &&
            grepl("Fish.merged\\$biomass", paste(deparse(x), collapse = "")),
            as.list(parse("code/R/04_pisco_processing.R")))
stopifnot(length(x) == 1L)
env$Fish.merged <- data.frame(WL_input_length = "TL", WL_a = adopted$WL_a[i],
  WL_b = adopted$WL_b[i], count = c(1, 2, 0), fish_tl = 30, biomass = NA_real_)
eval(x[[1]], env)
stopifnot(isTRUE(all.equal(env$Fish.merged$biomass, c(1, 2, 0) * 0.0144 * 30^3.04)))
for (bad in list(original[-i, ], rbind(original, original[i, ]),
                 original[, setdiff(names(original), "WL_a")])) {
  error <- tryCatch(env$apply_sheephead_allometry(bad), error = identity)
  stopifnot(inherits(error, "error"))
}
cat("PASS: verified units/formula, actual count-weighted assignment, other species and raw file preserved; missing/duplicate SPUL rejected.\n")
