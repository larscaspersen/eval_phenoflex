# Run from package root against the installed candidate package selected by R_LIBS.
source("tests/regression/cases.R")
destination <- "tests/testthat/fixtures/weather-regression"
cat("Testing evalpheno", as.character(packageVersion("evalpheno")), "from", find.package("evalpheno"), "\n")
provenance <- read.csv(file.path(destination, "baseline_sources.csv"))
stopifnot(identical(unname(tools::md5sum(file.path(destination, "inputs.rds"))),
                    provenance$md5[basename(provenance$file) == "inputs.rds"]))
expected <- read.csv(file.path(destination, "baseline.csv"), colClasses = c(season = "character"),
                     na.strings = "NA")
# Empty labels must remain strings, not missing values.
actual <- run_regression_cases(readRDS(file.path(destination, "inputs.rds")))
comparison <- compare_regression(expected, actual)
dir.create("development/regression-results", recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, "development/regression-results/comparison.csv", row.names = FALSE)
if (!all(comparison$passed)) {
  print(comparison[!comparison$passed, ])
  stop("Regression differences found; baseline was not modified.")
}
cat("PASS:", nrow(comparison), "values match the saved baseline.\n")
