legacy_combined_inputs <- function() {
  list(
    x = c(rep(40, 2), rep(300, 2), rep(0.5, 2), 25, 286.1, 47.7, 28, 4, 4, 1.6),
    modelfn = function(weather, par) weather,
    bloomJDays = list(c(100, 110), 120),
    SeasonList = list(list(102, 107), list(124)),
    ncult = 2
  )
}

test_that("legacy combined fitting accepts paired cultivar lists of different lengths", {
  inputs <- legacy_combined_inputs()
  result <- suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning")
  expect_named(result, c("F", "g"))
  expect_equal(result$F, 29)
  inputs$return_pred <- TRUE
  expect_equal(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"), c(102, 107, 124))
  inputs <- legacy_combined_inputs()
  inputs$bloomJDays <- lapply(inputs$bloomJDays, as.integer)
  expect_equal(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning")$F, 29)
})

test_that("legacy combined fitting rejects invalid cultivar counts before evaluation", {
  for (count in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "2", matrix(2))) {
    inputs <- legacy_combined_inputs()
    inputs$ncult <- count
    expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"), "ncult.*positive integer")
  }
})

test_that("legacy combined fitting requires outer lists matching ncult", {
  for (observations in list(c(100, 120), data.frame(a = 100, b = 120),
                           list(100), list(100, 110, 120))) {
    inputs <- legacy_combined_inputs()
    inputs$bloomJDays <- observations
    expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"), "bloomJDays.*length equal to ncult")
  }
  for (seasons in list(c(100, 120), data.frame(a = 100, b = 120),
                      list(list(100)), list(list(100), list(110), list(120)))) {
    inputs <- legacy_combined_inputs()
    inputs$SeasonList <- seasons
    expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"), "SeasonList.*length equal to ncult")
  }
})

test_that("legacy combined fitting checks observation vectors and nested season lists", {
  for (observations in list(c("100", "110"), list(100, 110), matrix(c(100, 110)))) {
    inputs <- legacy_combined_inputs()
    inputs$bloomJDays[[1]] <- observations
    expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"),
                 "bloomJDays\\[\\[1\\]\\].*numeric vector")
  }
  for (seasons in list(c(102, 107), data.frame(Temp = c(102, 107)))) {
    inputs <- legacy_combined_inputs()
    inputs$SeasonList[[1]] <- seasons
    expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"),
                 "SeasonList\\[\\[1\\]\\].*list of seasons")
  }
})

test_that("legacy combined fitting rejects unmatched seasons for each cultivar", {
  for (cultivar in 1:2) {
    for (return_pred in c(FALSE, TRUE)) {
      inputs <- legacy_combined_inputs()
      inputs$bloomJDays[[cultivar]] <- rep(100, 3)
      inputs$return_pred <- return_pred
      inputs$x <- numeric() # Validation must precede parameter conversion.
      inputs$modelfn <- function(...) stop("must not predict")
      expect_error(suppressWarnings(do.call(eval_phenoflex_combined, inputs), classes = "deprecatedWarning"),
                   paste0("Cultivar ", cultivar, ".*must have the same length"))
    }
  }
})
