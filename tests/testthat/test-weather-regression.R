source(test_path("..", "regression", "cases.R"), local = TRUE)

test_that("station wrappers and evaluations match the frozen baseline", {
  path <- test_path("fixtures", "weather-regression")
  expected <- read.csv(file.path(path, "baseline.csv"), colClasses = c(season = "character"), na.strings = "NA")
  provenance <- read.csv(file.path(path, "baseline_sources.csv"))
  expect_identical(unname(tools::md5sum(file.path(path, "inputs.rds"))),
    provenance$md5[basename(provenance$file) == "inputs.rds"])
  actual <- run_regression_cases(readRDS(file.path(path, "inputs.rds")))
  comparison <- compare_regression(expected, actual)
  failures <- comparison[!comparison$passed, ]
  expect_equal(nrow(failures), 0L, info = paste(capture.output(head(failures, 20)), collapse = "\n"))
  expect_true(all(is.finite(actual$value[actual$case_id == "ordinary"])))
  expect_true(all(is.na(actual$value[actual$case_id == "no_bloom"])))
  expect_equal(actual$value[actual$case_id == "ordinary"], actual$value[actual$case_id == "characteristic"])
  expect_true(all(c("calendar_dec31", "calendar_jan01", "calendar_feb29",
                    "combined", "three_stages", "sequential", "parallel", "partial_overlap") %in% actual$case_id))
})

test_that("regression comparison detects changed values and missing predictions", {
  x <- data.frame(case_id = "example", status = "finite", value = 10)
  changed <- x; changed$value <- 11
  expect_false(compare_regression(x, changed)$passed)
  changed$value <- NA_real_
  expect_false(compare_regression(x, changed)$passed)
})
