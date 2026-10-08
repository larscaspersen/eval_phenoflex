test_that("single fitted objects predict directly and separate settings from diagnostics", {
  x <- cv_example()
  fit <- fit_phenology(x$model, phenology_rss, numeric(), numeric(),
                       seasons = x$seasons, observed = x$observed)
  expect_identical(predict_phenology(fit, x$seasons), predict_phenology(fit$model, x$seasons))
  expect_identical(predict_phenology(fit, x$seasons[[1]], basic_output = FALSE),
                   predict_phenology(fit$model, x$seasons[[1]], basic_output = FALSE))
  override <- fit$par
  override["zc"] <- 10
  expect_identical(predict_phenology(fit, x$seasons, override),
                   predict_phenology(fit$model, x$seasons, override))
  expect_identical(fit$predictions$observed, x$observed)
  expect_identical(fit$predictions$predicted, x$observed)
  expect_named(fit, c("par", "value", "model", "predictions", "calibration_settings", "diagnostics"))
  expect_equal(fit$calibration_settings$parameters, x$model$parameters)
  expect_equal(fit$calibration_settings$lower, x$model$parameters)
  expect_equal(fit$diagnostics$evaluations, 1)
  expect_equal(fit$diagnostics$termination, "fixed")
  expect_null(fit$diagnostics$optimizer_result)
  expect_equal(model_parameters(fit), fit$par)
  expect_equal(as_pheno_models(fit), list(fit$model))
})

test_that("fitted models store aligned calibration dates without evaluating unobserved weather", {
  x <- cv_example()
  observed <- as.list(x$observed)
  observed[2] <- list(NULL)
  x$seasons[[2]] <- NULL
  # Keep an explicit invalid slot; NULL observation ensures it is skipped.
  seasons <- append(x$seasons, list(NULL), after = 1)
  fit <- fit_phenology(x$model, phenology_rss, numeric(), numeric(),
                       seasons = seasons, observed = observed)
  expect_true(is.na(fit$predictions$predicted[2]))
  expect_equal(unname(fit$predictions$predicted[-2]), unname(unlist(observed[-2])))
  expect_identical(fit$predictions$observed, observed)
  expect_equal(fit$diagnostics$evaluations, 1)
})

test_that("custom evaluators retain observations and can supply a prediction callback", {
  x <- cv_example()
  calls <- 0L
  callback <- function(model) {
    calls <<- calls + 1L
    predict_phenology(model, x$seasons)
  }
  loss <- function(p, m, observed) sum(p)
  fit <- fit_phenology(x$model, loss, numeric(), numeric(),
                       observed = x$observed, prediction = callback)
  expect_equal(calls, 1)
  expect_identical(fit$predictions$predicted, x$observed)
  expect_identical(fit$predictions$observed, x$observed)
  fit <- fit_phenology(x$model, loss, numeric(), numeric(), observed = x$observed)
  expect_null(fit$predictions$predicted)
  expect_identical(fit$predictions$observed, x$observed)
  fit <- fit_phenology(x$model, phenology_rss, numeric(), numeric(),
                       seasons = x$seasons, observed = x$observed, prediction = FALSE)
  expect_null(fit$predictions$predicted)
  expect_error(fit_phenology(x$model, loss, prediction = TRUE), "prediction must")
})

test_that("combined, stage and population fit wrappers retain prediction shapes", {
  x <- cv_example()
  stages <- stage_pheno_models(x$model, c(early = 5, late = 10))
  observations <- predict_phenology(stages, x$seasons)
  fit <- fit_phenology(stages, phenology_rss_stages, numeric(), numeric(),
                       seasons = x$seasons, observed = observations)
  expect_equal(predict_phenology(fit, x$seasons), observations)
  expect_equal(fit$predictions$predicted, observations)
  expect_identical(fit$predictions$observed, observations)
  combined <- combined_pheno_model(x$model, 2)
  nested <- list(a = x$seasons[1:3], b = x$seasons[3:1])
  observations <- predict_phenology(combined, nested)
  fit <- fit_phenology(combined, phenology_rss_combined, numeric(), numeric(),
                       seasons = nested, observed = observations)
  expect_equal(predict_phenology(fit, nested), observations)
  expect_equal(fit$predictions$predicted, observations)
  expect_identical(fit$predictions$observed, observations)
  population <- population_pheno_model("sequential", chill_dynamic("kinetic"), n = 2)
  population$model <- x$model
  fit <- fit_phenology(population, function(p, m) sum(p), numeric(), numeric())
  expect_equal(predict_phenology(fit, x$seasons[[1]]), predict_phenology(population, x$seasons[[1]]))
})

test_that("CV results predict directly using the stored ensemble and retain paired dates", {
  x <- cv_example()
  data <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_indices = 8)
  cv <- fit_phenology_cv(x$model, data, phenology_rss, numeric(), numeric(), refit = TRUE)
  expect_s3_class(cv, "phenology_cv")
  expect_s3_class(cv$model, "pheno_ensemble")
  expect_equal(predict_phenology(cv, x$seasons), predict_phenology(cv$model, x$seasons))
  expect_equal(predict_phenology(cv, x$seasons, basic_output = FALSE),
               predict_phenology(cv$model, x$seasons, basic_output = FALSE))
  expect_equal(validate_phenology(cv), validate_phenology(cv$model, data))
  expect_equal(names(cv$calibration$validation$weather), "8")
  expect_null(cv$calibration$data)
  expect_null(cv$calibration$diagnostics)
  expect_true(all(c("predicted", "observed") %in% names(cv$predictions)))
  expect_false(any(cv$predictions$case_id %in% .cv_validation_indices(data)))
  expect_true(all(cv$predictions$set == "assessment"))
  expect_named(cv, c("model", "predictions", "calibration"))
  for (id in names(cv$model$models)) {
    expect_equal(predict_phenology(cv$model$models[[id]], x$seasons), x$observed)
  }
  expect_equal(predict_phenology(cv$calibration$refit, x$seasons), x$observed)
  weighted <- fit_phenology_cv(x$model, data, phenology_rss, numeric(), numeric(),
                               weighting = "inverse_mse", failure_threshold = 1, max_weight = .5)
  expect_equal(weighted$model$weighting, "inverse_mse")
  expect_equal(weighted$model$failure_threshold, 1)
  expect_equal(weighted$model$max_weight, .5)
  expect_error(predict_phenology(cv, x$seasons, x$model$parameters), "overrides")
})

test_that("CV wrappers preserve combined/stage output shapes for direct prediction", {
  x <- cv_example()
  stages <- stage_pheno_models(x$model, c(early = 5, late = 10))
  observed <- predict_phenology(stages, x$seasons)
  data <- prepare_phenology_cv(x$seasons, observed, layout = "stages", v = 2,
                              validation_indices = 8)
  cv <- fit_phenology_cv(stages, data, phenology_rss_stages, numeric(), numeric())
  expect_length(cv$calibration$validation$weather, 1)
  expect_equal(nrow(cv$calibration$validation$cases), 2)
  expect_equal(predict_phenology(cv, x$seasons), observed)
  expect_equal(validate_phenology(cv)$metrics$mse, 0)
  combined <- combined_pheno_model(x$model, 2)
  nested <- list(a = x$seasons[1:4], b = x$seasons[4:1])
  observed <- predict_phenology(combined, nested)
  data <- prepare_phenology_cv(nested, observed, layout = "combined", v = 2,
                              validation_indices = list(a = 1, b = 2))
  cv <- fit_phenology_cv(combined, data, phenology_rss_combined, numeric(), numeric())
  expect_equal(predict_phenology(cv, nested), observed)
  expect_equal(validate_phenology(cv)$metrics$mse, 0)
})

test_that("training pairs and diagnostics are optional and summaries use assessment only", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2, repeats = 4, validation_indices = 8)
  lean <- fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric(), refit = TRUE)
  detailed <- fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric(), refit = TRUE,
                               keep_diagnostics = TRUE, keep_training_predictions = TRUE)
  expect_equal(lean$model, detailed$model)
  assessment <- detailed$predictions[detailed$predictions$set == "assessment", ]
  rownames(assessment) <- NULL
  expect_equal(lean$predictions, assessment)
  for (id in names(detailed$model$models)) {
    training <- detailed$predictions[detailed$predictions$set == "training" &
                                       detailed$predictions$split_id == id, ]
    expect_equal(training$case_id, .cv_splits(p)[[id]]$train)
    expect_equal(training$predicted, training$observed)
    expect_false(any(training$case_id %in% .cv_validation_indices(p)))
  }
  expect_equal(nrow(detailed$calibration$diagnostics$restart_metrics), 9)
  expect_named(detailed$calibration$diagnostics$selected, c(names(detailed$model$models), "full"))
  expect_null(lean$calibration$diagnostics)
  expect_equal(summary(detailed), summary(lean))
  # Training errors cannot change the pooled held-out score.
  detailed$predictions$predicted[detailed$predictions$set == "training"] <- NA_real_
  expect_equal(summary(detailed)$metrics, summary(lean)$metrics)
  expect_equal(summary(lean)$metrics$n_observed, 28)
  expect_equal(summary(lean)$metrics$n_failed, 0)
  expect_equal(summary(lean)$metrics$rmse, 0)
  expect_equal(nrow(summary(lean)$scores), 8)
  expect_lt(length(capture.output(print(lean))), 12)
  expect_lt(length(capture.output(print(p))), 12)
  expect_output(print(lean), "2 more rows")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(lean, path)
  restored <- readRDS(path)
  expect_equal(predict_phenology(restored, x$seasons), x$observed)
  expect_equal(validate_phenology(restored), validate_phenology(lean$model, p))
  for (option in c("keep_diagnostics", "keep_training_predictions"))
    expect_error(do.call(fit_phenology_cv, c(list(x$model, p, phenology_rss, numeric(), numeric()),
                                           setNames(list(NA), option))), "single logical")
})
