test_that("weighted summaries handle failure votes and use weighted descriptive spread", {
  mat <- rbind(c(10, 20, 30), c(10, NA, 30), c(NA, NA, 30), c(NA, NA, NA))
  colnames(mat) <- c("a", "b", "c")
  w <- c(c = .5, a = .2, b = .3)
  out <- summarise_ensemble_predictions(mat, w)
  expect_equal(out$predicted[1], 23)
  expect_equal(out$ensemble_sd[1], sqrt(.2 * 13^2 + .3 * 3^2 + .5 * 7^2))
  expect_equal(out$predicted[2], (10 * .2 + 30 * .5) / .7)
  expect_true(all(is.na(out$predicted[3:4])))
  expect_equal(out$n_fail, 0:3)
  expect_equal(out$failure_fraction, (0:3) / 3)
  expect_equal(out$ensemble_sd[3], 0)
  expect_true(is.na(out$ensemble_sd[4]))
  # Count and weight votes differ when the two failing members have low weight.
  mat <- matrix(c(NA, NA, 30), 1, dimnames = list(NULL, c("a", "b", "c")))
  expect_equal(summarise_ensemble_predictions(mat, c(.1, .1, .8), failure_vote = "weight")$predicted, 30)
  # Surviving predictions with no positive weight must never produce day zero.
  expect_true(is.na(summarise_ensemble_predictions(matrix(c(NA, 10), 1), c(1, 0),
                                                   failure_threshold = 1)$predicted))
  expect_true(is.na(summarise_ensemble_predictions(matrix(c(NA, 10), 1))$predicted))
})

test_that("caps redistribute repeatedly and document infeasible survivor caps", {
  # One-pass redistribution would leave the second member at .45 above the cap.
  mat <- matrix(c(0, 100, 200), 1)
  out <- summarise_ensemble_predictions(mat, c(.8, .15, .05), max_weight = .4)
  expect_equal(out$predicted, 80) # capped weights .4, .4, .2
  expect_false(out$cap_relaxed)
  out <- summarise_ensemble_predictions(matrix(c(10, NA, NA), 1), c(.8, .15, .05),
                                        failure_threshold = 1, max_weight = .4)
  expect_equal(out$predicted, 10)
  expect_true(out$cap_relaxed)
  expect_equal(summarise_ensemble_predictions(matrix(c(10, 20), 1), max_weight = .5)$predicted, 15)
})

test_that("matrix and long predictions agree including singleton dimensions", {
  mat <- matrix(c(10, 20, 30, 40, NA, 60), 3, dimnames = list(c("z", "y", "x"), c("b", "a")))
  long <- data.frame(case_id = rep(rownames(mat), 2), member_id = rep(colnames(mat), each = 3),
                     predicted = as.vector(mat))
  expect_equal(summarise_ensemble_predictions(long), summarise_ensemble_predictions(mat))
  expect_equal(summarise_ensemble_predictions(matrix(c(10, 20, NA), 3, 1))$predicted,
               c(10, 20, NA))
  expect_equal(summarise_ensemble_predictions(matrix(10, 1, 1))$ensemble_sd, 0)
  expect_error(summarise_ensemble_predictions(long[-1, ]), "exactly one")
  expect_error(summarise_ensemble_predictions(rbind(long, long[1, ])), "exactly one")
  expect_error(summarise_ensemble_predictions(matrix(Inf, 1)), "finite or NA")
  for (w in list(c(0, 0), c(-1, 2), c(NA, 1), c(Inf, 1), 1))
    expect_error(summarise_ensemble_predictions(matrix(1:2, 1), w), "weights must")
  expect_error(summarise_ensemble_predictions(mat, c(c = 1, d = 2)), "names must match")
  expect_error(summarise_ensemble_predictions(mat, failure_threshold = 0), "failure_threshold")
  expect_error(summarise_ensemble_predictions(mat, max_weight = 0), "max_weight")
})

test_that("ensembles predict new data with ordinary date shapes and optional details", {
  x <- cv_example()
  a <- x$model
  b <- a
  b$parameters["zc"] <- 10
  models <- list(a = a, b = b)
  ensemble <- pheno_ensemble(models, weights = c(b = .75, a = .25))
  predictions <- predict_phenology(models, x$seasons)
  expected <- predictions$a * .25 + predictions$b * .75
  expect_equal(predict_phenology(ensemble, x$seasons), expected)
  expect_equal(predict_phenology(ensemble, x$seasons[[1]]), unname(expected[1]))
  details <- predict_phenology(ensemble, x$seasons, basic_output = FALSE)
  expect_equal(details$summary$predicted, unname(expected))
  expect_equal(nrow(details$individual_predictions), 16)
  expect_equal(predict_phenology_ensemble(ensemble, x$seasons), details$summary)
  expect_equal(ensemble$weighting, "explicit")
  expect_equal(pheno_ensemble(models)$weighting, "equal")
  expect_equal(models, list(a = a, b = b))
  expect_error(predict_phenology(ensemble, x$seasons, a$parameters), "overrides")
  expect_error(pheno_ensemble(models, weighting = "inverse_mse"), "requires")
  expect_error(pheno_ensemble(list(a, population_pheno_model())), "Ensemble entries")
  expect_error(pheno_ensemble(list(a, stage_pheno_models(a, c(5, 10)))), "same number")
})

test_that("CV weighting uses held-out MSE and matches member IDs", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2)
  cv <- fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric())
  cv$calibration$scores$mse <- c(0, 3)
  ensemble <- pheno_ensemble(cv, weighting = "inverse_mse", epsilon = 1)
  expect_equal(unname(ensemble$weights), c(.8, .2))
  expect_equal(ensemble$weighting, "inverse_mse")
  cv$calibration$scores$mse <- c(Inf, 3)
  expect_equal(unname(pheno_ensemble(cv, weighting = "inverse_mse")$weights), c(0, 1))
})
