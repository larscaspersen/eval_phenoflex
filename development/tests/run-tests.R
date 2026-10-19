# Run from the package root: Rscript development/tests/run-tests.R
if (!file.exists("DESCRIPTION") || !dir.exists("development/tests/testthat"))
  stop("Run this script from the evalpheno package root.")

pkgload::load_all(".", quiet = TRUE)
testthat::test_dir("tests/testthat", package = "evalpheno", reporter = "summary")
testthat::test_dir("development/tests/testthat", package = "evalpheno", reporter = "summary")
