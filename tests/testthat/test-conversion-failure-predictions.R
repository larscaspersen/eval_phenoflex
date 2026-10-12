test_that("single-cultivar evaluators preserve their output types on conversion failure", {
  testthat::local_mocked_bindings(
    .fit_chill_parameters = function(...) list(termcd = 3L), .package = "evalpheno"
  )
  parameters <- c(40, 190, 0.5, 25, 286.1, 47.7, 28, 4, 4, 1.6)
  seasons <- list(NULL, NULL, NULL)
  for (evaluator in list(eval_phenoflex_single, eval_phenoflex_single_twofixed)) {
    expect_identical(evaluator(parameters, function(...) stop("must not predict"),
                               c(100, 110, 120), seasons, return_pred = TRUE),
                     rep(NA_real_, 3))
    expect_identical(evaluator(parameters, function(...) stop("must not predict"),
                               c(100, 110, 120), seasons),
                     list(F = 1e6, g = rep(1e6, 5)))
  }
})

test_that("combined evaluators preserve prediction lengths and shapes on conversion failure", {
  testthat::local_mocked_bindings(
    .fit_chill_parameters = function(...) list(termcd = 3L), .package = "evalpheno"
  )
  seasons <- list(list(NULL, NULL), list(NULL))
  observed <- list(c(100, 110), 120)
  reduced <- c(rep(40, 2), rep(190, 2), rep(0.5, 2), 25, 286.1, 47.7, 28, 4, 4, 1.6)
  full <- c(rep(40, 2), rep(190, 2), rep(0.5, 2), 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  model <- function(...) stop("must not predict")
  expect_identical(eval_phenoflex_combined(reduced, model, observed, seasons,
                                          ncult = 2, return_pred = TRUE),
                   rep(NA_real_, 3))
  expect_identical(eval_phenoflex_combined(reduced, model, observed, seasons, ncult = 2),
                   list(F = 1e6, g = rep(1e6, 5)))
  observed_frames <- lapply(observed, function(dates) data.frame(pheno = dates))
  expect_identical(evaluation_function_meigo_nonlinear_combined(
    full, model, observed_frames, seasons, n_cult = 2, return_bloom_days = TRUE
  ), list(rep(NA_real_, 2), NA_real_))
  expect_identical(evaluation_function_meigo_nonlinear_combined(
    full, model, observed_frames, seasons, n_cult = 2
  ), list(F = 1e6, g = rep(1e6, 5)))
})
