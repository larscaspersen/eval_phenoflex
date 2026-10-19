rss_example <- function() {
  model <- pheno_model("sequential", chill_dynamic("kinetic"))
  parameters <- model$parameters
  parameters[c("yc", "zc")] <- c(0.1, 5)
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  warmer <- weather
  warmer$Temp <- 10
  list(model = model, parameters = parameters, seasons = list(weather, warmer))
}

test_that("RSS sums residuals across seasons using the supplied candidate parameters", {
  x <- rss_example()
  predicted <- vapply(x$seasons, function(weather) {
    predict_phenology(x$model, weather, x$parameters)
  }, numeric(1))
  expect_true(all(is.finite(predicted)))
  expect_equal(phenology_rss(x$parameters, x$model, x$seasons, predicted), 0)
  expect_equal(phenology_rss(x$parameters, x$model, x$seasons, predicted + c(2, -3)), 13)
  # Stored parameters yield no bloom; the explicit candidate must be used.
  expect_true(is.na(predict_phenology(x$model, x$seasons[[1]])))
  unchanged <- x$model$parameters
  phenology_rss(x$parameters, x$model, x$seasons, predicted)
  expect_identical(x$model$parameters, unchanged)
})

test_that("RSS excludes candidates that do not predict bloom", {
  x <- rss_example()
  x$parameters["zc"] <- 1e12
  expect_identical(phenology_rss(x$parameters, x$model, x$seasons, c(55, 56)), Inf)
})

test_that("RSS rejects unmatched or missing observations and invalid weather", {
  x <- rss_example()
  for (observed in list(55, numeric(), c(55, NA), c(55, Inf), c("55", "56"),
                        matrix(c(55, 56), ncol = 1))) {
    expect_error(phenology_rss(x$parameters, x$model, x$seasons, observed),
                 "one finite numeric phenology date per season")
  }
  expect_error(phenology_rss(x$parameters, x$model, list(), numeric()), "non-empty list")
  expect_error(phenology_rss(x$parameters, x$model, x$seasons[[1]], 55), "non-empty list")
  expect_error(phenology_rss(x$parameters, x$model, list(data.frame(Temp = 8)), 55))
  expect_error(phenology_rss(x$parameters, population_pheno_model(), x$seasons, c(55, 56)),
               "single pheno_model")
})
