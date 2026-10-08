test_that("preparation exposes balanced entry assignments and independent validation", {
  x <- cv_example()
  set.seed(901)
  state <- .Random.seed
  prepared <- prepare_phenology_cv(x$seasons, x$observed, v = 3, repeats = 2,
                                   validation_indices = 7:8, seed = 17)
  expect_identical(.Random.seed, state)
  expect_s3_class(prepared, "phenology_cv_data")
  expect_named(prepared, c("data", "assignments", "settings"))
  expect_named(prepared$assignments, c("case_id", "member", "season", "year", "set", "repeat_id", "fold"))
  expect_false("observed" %in% names(prepared$data$cases))
  expect_length(.cv_splits(prepared), 6)
  expect_equal(.cv_cases(prepared)$year[.cv_validation_indices(prepared)], 2007:2008)
  expect_equal(as.vector(table(prepared$assignments$repeat_id, prepared$assignments$fold)), rep(2L, 6))
  for (s in .cv_splits(prepared)) {
    expect_length(s$train, 4)
    expect_length(s$assessment, 2)
    expect_length(intersect(s$train, s$assessment), 0)
    expect_length(intersect(c(s$train, s$assessment), .cv_validation_indices(prepared)), 0)
  }
  for (r in 1:2) {
    assessed <- unlist(lapply(.cv_splits(prepared), function(s)
      if (s$repeat_id == r) s$assessment else integer()))
    expect_equal(sort(unname(assessed)), 1:6)
  }
  again <- prepare_phenology_cv(x$seasons, x$observed, v = 3, repeats = 2,
                                validation_indices = 7:8, seed = 17)
  expect_identical(again, prepared)
  expect_output(print(prepared), "independent validation entries")
})

test_that("preparation groups all observations of a year and preserves NULL dates", {
  x <- cv_example()
  nested <- list(a = x$seasons, b = x$seasons[8:1])
  obs <- list(a = as.list(x$observed), b = rev(x$observed))
  obs$a[2] <- list(NULL)
  p <- prepare_phenology_cv(nested, obs, layout = "combined", v = 3,
                            validation_indices = list(a = 8, b = 1),
                            groups = list(x$years, rev(x$years)), seed = 9)
  expect_true(is.na(.cv_cases(p)$observed[2]))
  for (s in .cv_splits(p)) {
    train <- .cv_cases(p)$year[s$train]
    assessment <- .cv_cases(p)$year[s$assessment]
    expect_false(any(train %in% assessment))
    expect_true(all(table(.cv_cases(p)$year, .cv_cases(p)$case_id %in% s$assessment) %in% c(0, 2)))
    subset <- evalpheno:::.cv_subset(p, s$train)
    if (2002 %in% train) expect_true(any(vapply(subset$observed$a, is.null, logical(1))))
  }
})

test_that("validation proportion, seed restoration and bad preparation inputs", {
  x <- cv_example()
  a <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_prop = 0.25, seed = 7)
  b <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_prop = 0.25, seed = 7)
  expect_identical(a, b)
  expect_length(.cv_validation_indices(a), 2)
  expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 1), "at least two")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 8, validation_indices = 8), "at least v")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, validation_indices = 999), "within range")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, validation_prop = 1), "strictly")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, validation_prop = .2,
                                   validation_indices = 8), "not both")
  expect_error(prepare_phenology_cv(x$seasons, x$observed[-1]), "one finite numeric")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, years = 1:2), "per season")
  expect_error(prepare_phenology_cv(x$seasons, rep(list(NULL), 8), v = 2), "observed date")
  no_year <- lapply(x$seasons, function(s) {s$Year <- NULL; s})
  expect_error(prepare_phenology_cv(no_year, x$observed), "Supply years")
  expect_s3_class(prepare_phenology_cv(no_year, x$observed, years = x$years), "phenology_cv_data")
})

test_that("fitting consumes inspected splits and does not use reserved years", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 3, repeats = 2, validation_indices = 7:8)
  # Make independent observations deliberately unrelated; they must not affect CV.
  p$data$observed[7:8] <- 200
  calls <- list()
  loss <- function(parameters, model, seasons, observed, offset) {
    calls[[length(calls) + 1L]] <<- vapply(seasons, function(s) max(s$Year), numeric(1))
    expect_false(any(calls[[length(calls)]] %in% 2007:2008))
    expect_equal(length(seasons), length(observed))
    phenology_rss(parameters, model, seasons, observed) + offset
  }
  state <- .Random.seed
  cv <- fit_phenology_cv(x$model, p, loss, lower = numeric(), upper = numeric(),
                         restarts = 2, refit = TRUE, seed = 17, offset = 3, keep_diagnostics = TRUE)
  expect_identical(.Random.seed, state)
  expect_identical(cv$calibration$settings$preparation, p$settings)
  expect_length(cv$model$models, 6)
  expect_length(calls, 14)
  expect_equal(cv$calibration$scores$training_loss, rep(3, 6))
  expect_equal(cv$calibration$scores$mse, rep(0, 6))
  expect_equal(nrow(cv$predictions), 12)
  expect_equal(cv$calibration$refit_loss, 3)
  expect_equal(nrow(cv$calibration$diagnostics$restart_metrics), 14)
  expect_false(anyDuplicated(cv$calibration$diagnostics$restart_metrics$seed) > 0)
  expect_output(print(cv), "fold models")
  validation <- validate_phenology(pheno_ensemble(cv), p)
  expect_equal(validation$predictions$year, 2007:2008)
  expect_gt(validation$metrics$mse, 0)
  expect_equal(validation$metrics$n_observed, 2)
  expect_equal(validate_phenology(cv$calibration$refit, p), validate_phenology(x$model, p))
  again <- fit_phenology_cv(x$model, p, loss, lower = numeric(), upper = numeric(),
                            restarts = 2, seed = 17, offset = 3)
  expect_identical(again$predictions, cv$predictions)
})

test_that("restarts are selected using training loss and held-out failures score Inf", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2)
  cv <- fit_phenology_cv(x$model, p, function(p, m, seasons, observed) runif(1),
                         lower = numeric(), upper = numeric(), restarts = 3, seed = 5, keep_diagnostics = TRUE)
  for (id in names(cv$model$models)) expect_equal(cv$calibration$scores$training_loss[cv$calibration$scores$split_id == id],
    min(cv$calibration$diagnostics$restart_metrics$training_loss[cv$calibration$diagnostics$restart_metrics$split_id == id]))
  failure <- x$model
  failure$parameters["zc"] <- 1e12
  cv <- fit_phenology_cv(failure, p, function(p, m, seasons, observed) 0,
                         lower = numeric(), upper = numeric())
  expect_true(all(is.infinite(cv$calibration$scores$mse)))
  expect_true(all(cv$calibration$scores$n_failed == 4))
  expect_error(pheno_ensemble(cv, weighting = "inverse_mse"), "positive total")
  expect_error(validate_phenology(x$model, p), "No independent")
})

test_that("stage and combined fits preserve within-fit sharing and separate outputs", {
  x <- cv_example()
  stages <- stage_pheno_models(x$model, c(early = 5, late = 10))
  obs <- predict_phenology(stages, x$seasons)
  p <- prepare_phenology_cv(x$seasons, obs, layout = "stages", v = 2, validation_indices = 8)
  cv <- fit_phenology_cv(stages, p, phenology_rss_stages, numeric(), numeric())
  expect_equal(cv$calibration$scores$mse, c(0, 0))
  expect_equal(nrow(cv$predictions), 14)
  ensemble <- pheno_ensemble(cv)
  expect_equal(predict_phenology(ensemble, x$seasons), obs)
  expect_equal(validate_phenology(ensemble, p)$metrics$mse, 0)
  details <- predict_phenology_ensemble(ensemble, x$seasons, return_members = TRUE)
  expect_equal(unique(details$summary$output), c("early", "late"))
  expect_equal(nrow(details$individual_predictions), 32)
  # Unnamed observations still match named stages by position.
  p <- prepare_phenology_cv(x$seasons, unname(obs), layout = "stages", v = 2)
  expect_s3_class(fit_phenology_cv(stages, p, phenology_rss_stages, numeric(), numeric()), "phenology_cv")
  nested <- list(a = x$seasons, b = x$seasons[8:1])
  combined <- combined_pheno_model(x$model, 2)
  obs <- predict_phenology(combined, nested)
  p <- prepare_phenology_cv(nested, obs, layout = "combined", v = 2,
                            validation_indices = list(a = 8, b = 1))
  cv <- fit_phenology_cv(combined, p, phenology_rss_combined, numeric(), numeric())
  expect_equal(cv$calibration$scores$mse, c(0, 0))
  expect_equal(predict_phenology(pheno_ensemble(cv), nested), obs)
  expect_equal(validate_phenology(pheno_ensemble(cv), p)$metrics$mse, 0)
})

test_that("modified splits cannot leak assessment or independent years", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_indices = 8)
  bad <- p
  bad$assignments <- rbind(bad$assignments, transform(bad$assignments[bad$assignments$set == "validation", ],
                                                   set = "cv", repeat_id = 1L, fold = 1L))
  expect_error(fit_phenology_cv(x$model, bad, phenology_rss, numeric(), numeric()), "Invalid CV split")
  bad <- p
  bad$assignments <- rbind(bad$assignments, bad$assignments[1, ])
  expect_error(fit_phenology_cv(x$model, bad, phenology_rss, numeric(), numeric()), "Invalid CV split")
  expect_error(fit_phenology_cv(x$model, p, lower = numeric(), upper = numeric()), "evaluation")
  expect_error(fit_phenology_cv(population_pheno_model(), p, phenology_rss), "layout")
  expect_error(fit_phenology_cv(x$model, p, phenology_rss, seasons = x$seasons), "come from data")
  expect_error(fit_phenology_cv(x$model, p, function(...) stop("bad loss"), numeric(), numeric()), "restart 1: bad loss")
})

test_that("prepared CV works with both actual optimizers and partial bounds", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_indices = 8)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 20, trace = FALSE) else list(max.call = 50)
    cv <- fit_phenology_cv(x$model, p, phenology_rss,
                           lower = c(yc = .09, zc = 1), upper = c(yc = .11, zc = 20),
                           optimizer = optimizer, seed = 17, max_iterations = 2, control = control,
                           keep_diagnostics = TRUE)
    expect_equal(cv$calibration$scores$mse, c(0, 0))
    for (id in names(cv$model$models)) {
      expect_gt(cv$calibration$diagnostics$selected[[id]]$diagnostics$evaluations, 1)
      expect_false(is.null(cv$calibration$diagnostics$selected[[id]]$diagnostics$optimizer_result))
      par <- model_parameters(cv$model$models[[id]])
      fixed <- setdiff(names(par), c("yc", "zc"))
      expect_identical(par[fixed], x$model$parameters[fixed])
    }
  }
})

test_that("edited assignments remain authoritative and grouping prevents leakage", {
  x <- cv_example()
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2, repeats = 2,
                            groups = rep(1:4, each = 2))
  # A coherent edit changes the assessment rows consumed by fitting.
  p$assignments$fold[p$assignments$repeat_id == 1] <-
    3L - p$assignments$fold[p$assignments$repeat_id == 1]
  cv <- fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric())
  expect_equal(cv$predictions$case_id[cv$predictions$split_id == "repeat1_fold1"],
               p$assignments$case_id[p$assignments$repeat_id == 1 & p$assignments$fold == 1])
  p$assignments$fold[1] <- 3L - p$assignments$fold[1]
  expect_error(fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric()), "Invalid CV split")
})

test_that("same-year locations can be withheld by index without withholding the year", {
  x <- cv_example()
  x$seasons <- lapply(x$seasons, function(s) {s$Year <- 2020; s})
  names(x$seasons) <- paste0("location", seq_along(x$seasons))
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 3, repeats = 2,
                            validation_indices = 2, seed = 17)
  expect_identical(.cv_validation_indices(p), 2L)
  expect_equal(p$settings$validation_indices, 2)
  expect_equal(unique(.cv_cases(p)$year), 2020)
  expect_equal(p$assignments$season[p$assignments$set == "validation"], 2)
  for (s in .cv_splits(p)) {
    expect_false(2 %in% c(s$train, s$assessment))
    expect_equal(sort(c(s$train, s$assessment)), c(1, 3:8))
  }
  calls <- list()
  loss <- function(parameters, model, seasons, observed) {
    calls[[length(calls) + 1L]] <<- names(seasons)
    expect_false("location2" %in% names(seasons))
    phenology_rss(parameters, model, seasons, observed)
  }
  cv <- fit_phenology_cv(x$model, p, loss, numeric(), numeric(), refit = TRUE)
  expect_length(calls, 7)
  expect_equal(cv$calibration$scores$mse, rep(0, 6))
  validation <- validate_phenology(pheno_ensemble(cv), p)
  expect_equal(validation$predictions$season_name, "location2")
  expect_equal(validation$metrics$mse, 0)
  a <- prepare_phenology_cv(x$seasons, x$observed, v = 2, validation_prop = .25, seed = 7)
  expect_length(.cv_validation_indices(a), 2)
  expect_equal(sum(a$assignments$set == "cv"), 6)
})

test_that("combined validation uses local indices and proportions count observed entries", {
  x <- cv_example()
  nested <- list(a = x$seasons, b = x$seasons)
  obs <- list(a = x$observed, b = x$observed)
  p <- prepare_phenology_cv(nested, obs, layout = "combined", v = 3,
                            validation_indices = list(a = 2, b = c(1, 3)))
  expect_equal(.cv_validation_indices(p), c(2, 9, 11))
  expect_equal(.cv_cases(p)$member[.cv_validation_indices(p)], c(1, 2, 2))
  expect_equal(.cv_cases(p)$season[.cv_validation_indices(p)], c(2, 1, 3))
  model <- combined_pheno_model(x$model, 2)
  cv <- fit_phenology_cv(model, p, phenology_rss_combined, numeric(), numeric(), refit = TRUE)
  expect_equal(validate_phenology(cv$calibration$refit, p)$metrics$n_observed, 3)
  obs$a <- as.list(obs$a)
  obs$b <- as.list(obs$b)
  obs$a[7:8] <- rep(list(NULL), 2)
  obs$b[7:8] <- rep(list(NULL), 2)
  # Twelve observed entries, regardless of only eight distinct years.
  p <- prepare_phenology_cv(nested, obs, layout = "combined", v = 2,
                            validation_prop = .25, seed = 8)
  expect_length(.cv_validation_indices(p), 3)
  expect_false(any(is.na(.cv_cases(p)$observed[.cv_validation_indices(p)])))
  expect_error(prepare_phenology_cv(nested, obs, layout = "combined", v = 2,
                                   validation_indices = 2), "list in member order")
  expect_error(prepare_phenology_cv(nested, obs, layout = "combined", v = 2,
                                   validation_indices = list(a = 2, b = 9)), "within range")
  expect_error(prepare_phenology_cv(nested, obs, layout = "combined", v = 2,
                                   validation_indices = list(b = 2, a = 3)), "member order")
})

test_that("stage validation expands common-season indices and proportions to all stages", {
  x <- cv_example()
  stages <- stage_pheno_models(x$model, c(early = 5, late = 10))
  obs <- predict_phenology(stages, x$seasons)
  p <- prepare_phenology_cv(x$seasons, obs, layout = "stages", v = 2,
                            validation_indices = c(2, 5))
  expect_equal(.cv_validation_indices(p), c(2, 5, 10, 13))
  expect_equal(sort(unique(.cv_cases(p)$season[.cv_validation_indices(p)])), c(2, 5))
  p <- prepare_phenology_cv(x$seasons, obs, layout = "stages", v = 2,
                            validation_prop = .25, seed = 4)
  expect_length(.cv_validation_indices(p), 4)
  expect_length(unique(.cv_cases(p)$season[.cv_validation_indices(p)]), 2)
  expect_true(all(table(.cv_cases(p)$entry_id[.cv_validation_indices(p)]) == 2))
  for (s in .cv_splits(p))
    expect_false(any(.cv_cases(p)$entry_id[s$train] %in% .cv_cases(p)$entry_id[s$assessment]))
  bad <- p
  bad$assignments <- bad$assignments[!(bad$assignments$set == "validation" &
                                       bad$assignments$case_id == .cv_validation_indices(p)[1]), ]
  expect_error(fit_phenology_cv(stages, bad, phenology_rss_stages, numeric(), numeric()),
               "All stages")
})

test_that("optional groups control CV folds and bad indices fail clearly", {
  x <- cv_example()
  groups <- rep(c("site1-year1", "site2-year1", "site1-year2", "site2-year2"), each = 2)
  p <- prepare_phenology_cv(x$seasons, x$observed, v = 2, groups = groups, validation_indices = 1)
  expect_identical(.cv_validation_indices(p), 1L)
  # An explicit index reserves exactly that entry, even within a supplied group.
  expect_true(2 %in% .cv_splits(p)[[1]]$train || 2 %in% .cv_splits(p)[[1]]$assessment)
  for (s in .cv_splits(p)) expect_false(any(groups[s$train] %in% groups[s$assessment]))
  expect_s3_class(fit_phenology_cv(x$model, p, phenology_rss, numeric(), numeric()), "phenology_cv")
  for (indices in list(0, -1, 9, NA_real_, Inf, 1.5, c(2, 2), TRUE, matrix(2)))
    expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 2,
                                     validation_indices = indices), "validation_indices")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 2, groups = 1:2), "per season")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 2, groups = rep(NA, 8)), "per season")
  expect_error(prepare_phenology_cv(x$seasons, x$observed, v = 2, groups = rep("one", 8)), "at least v")
})
