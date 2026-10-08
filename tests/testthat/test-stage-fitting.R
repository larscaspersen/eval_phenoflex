stage_example <- function() {
  base <- pheno_model("sequential", chill_dynamic("kinetic"))
  base$parameters["yc"] <- 0.1
  models <- stage_pheno_models(base, c(early = 5, middle = 10, late = 15))
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  warmer <- weather
  warmer$Temp <- 10
  list(base = base, models = models, seasons = list(cold = weather, warm = warmer))
}

test_that("stage collections contain ordinary models with shared parameters", {
  x <- stage_example()
  expect_s3_class(x$models, "pheno_model_list")
  expect_named(x$models, c("early", "middle", "late"))
  expect_true(all(vapply(x$models, inherits, logical(1), what = "pheno_model")))
  expect_identical(vapply(x$models, function(m) m$parameters[["zc"]], numeric(1)),
                    c(early = 5, middle = 10, late = 15))
  p <- model_parameters(x$models)
  expect_identical(p[c("yc", "zc1", "zc2", "zc3")], c(yc = 0.1, zc1 = 5, zc2 = 10, zc3 = 15))
  expect_true(validate_parameters(x$models, rev(p)))
  expect_identical(model_parameters(x$base), x$base$parameters)
  for (i in 1:3) {
    single <- x$base
    single$parameters["zc"] <- c(5, 10, 15)[i]
    expect_identical(x$models[[i]], single)
  }
  expected <- lapply(x$models, function(m) {
    vapply(x$seasons, function(w) predict_phenology(m, w), numeric(1))
  })
  expect_identical(predict_phenology(x$models, x$seasons), expected)
  expect_true(all(vapply(1:2, function(year) {
    all(diff(vapply(expected, function(stage) stage[year], numeric(1))) >= 0)
  }, logical(1))))
  for (structure in c("parallel", "phenoflex")) {
    expect_true(validate_parameters(stage_pheno_models(pheno_model(structure), c(50, 100)),
                                    model_parameters(stage_pheno_models(pheno_model(structure), c(50, 100)))))
  }
})

test_that("NULL observations skip only their own stage and season", {
  x <- stage_example()
  p <- model_parameters(x$models)
  predicted <- predict_phenology(x$models, x$seasons)
  observed <- list(list(predicted[[1]][1] + 2, NULL),
                   list(NULL, predicted[[2]][2] - 3),
                   predicted[[3]] + c(4, -5))
  expect_equal(phenology_rss_stages(p, x$models, x$seasons, observed), 54)
  # Later stage never occurs, but none of its seasons is observed.
  p["zc3"] <- 1e12
  observed[[3]] <- rep(list(NULL), 2)
  expect_equal(phenology_rss_stages(p, x$models, x$seasons, observed), 13)
  observed[[3]][1] <- list(60)
  expect_identical(phenology_rss_stages(p, x$models, x$seasons, observed), Inf)
  observed <- rep(list(rep(list(NULL), 2)), 3)
  expect_error(phenology_rss_stages(p, x$models, x$seasons, observed), "At least one")
})

test_that("unobserved stage-year pairs are not predicted", {
  x <- stage_example()
  p <- model_parameters(x$models)
  calls <- 0L
  testthat::local_mocked_bindings(predict_phenology = function(model, weather, parameters, ...) {
    calls <<- calls + 1L
    expect_equal(weather$Temp[1], 8)
    parameters[["zc"]]
  })
  observed <- list(list(6, NULL), list(NULL, NULL), list(18, NULL))
  expect_equal(phenology_rss_stages(p, x$models, x$seasons, observed), 10)
  expect_identical(calls, 2L)
})

test_that("ordinary model lists can be calibrated without a combined wrapper", {
  x <- stage_example()
  simple <- as_pheno_models(x$models)
  models <- pheno_model_list(simple, shared = setdiff(names(simple[[1]]$parameters), "zc"), ordered = "zc")
  expect_identical(model_parameters(models), model_parameters(x$models))
  expect_false(inherits(simple, "pheno_model_list"))
  expect_identical(as_pheno_models(models), simple)
  observed <- predict_phenology(models, x$seasons)
  fit <- fit_phenology(simple, phenology_rss_combined,
                       lower = numeric(), upper = numeric(),
                       seasons = rep(list(x$seasons), 3), observed = observed)
  expect_equal(fit$value, 0)
  expect_s3_class(fit$model, "pheno_model_list")
})

test_that("combined specifications convert to full single models using fitted parameters", {
  base <- pheno_model("phenoflex", chill_dynamic("kinetic"))
  combined <- combined_pheno_model(base, 3, cultivar_specific = c("yc", "zc"))
  p <- combined$parameters
  p[c("yc1", "yc2", "yc3")] <- c(10, 20, 30)
  p["Tb"] <- 2
  simple <- as_pheno_models(combined, rev(p))
  expect_true(all(vapply(simple, inherits, logical(1), what = "pheno_model")))
  for (i in 1:3) expect_identical(simple[[i]]$parameters, cultivar_parameters(combined, i, p))
  expect_identical(combined$model$parameters[["Tb"]], 4)
  fit <- fit_phenology(combined, function(p, m) sum(p), numeric(), numeric(), parameters = p)
  expect_identical(as_pheno_models(fit), simple)
  expect_identical(model_parameters(fit), fit$par)
  fit$model$parameters["Tb"] <- 3
  expect_identical(as_pheno_models(fit), simple)
  expect_identical(as_pheno_models(base), list(base))
})

test_that("stage order and observation shapes are checked", {
  x <- stage_example()
  for (heat in list(numeric(), c(10, 5), c(5, NA), c(5, Inf), c(0, 5), matrix(c(5, 10))))
    expect_error(stage_pheno_models(x$base, heat), "positive, finite and non-decreasing")
  expect_error(stage_pheno_models(pheno_model("partial_overlap"), c(50, 100)), "zc parameter")
  p <- model_parameters(x$models)
  p[c("zc1", "zc2")] <- c(20, 10)
  expect_error(validate_parameters(x$models, p), "non-decreasing")
  expect_error(as_pheno_models(x$models, p), "non-decreasing")
  models <- as_pheno_models(x$models)
  models[[2]]$parameters["yc"] <- 0.2
  expect_error(pheno_model_list(models, shared = "yc"), "equal initial values")
  expect_error(pheno_model_list(models, shared = "unknown"), "model parameter names")
  expect_error(pheno_model_list(models, shared = "zc", ordered = "zc"), "unshared")
  expect_error(pheno_model_list(list(x$base, pheno_model("parallel"))), "identical structure")
  observed <- predict_phenology(x$models, x$seasons)
  expect_error(phenology_rss_stages(model_parameters(x$models), x$models, x$seasons, observed[-1]), "per cultivar")
  observed[[2]] <- list(55)
  expect_error(phenology_rss_stages(model_parameters(x$models), x$models, x$seasons, observed), "Cultivar 2")
})

test_that("both optimizers fit ordered stage thresholds with sparse observations", {
  x <- stage_example()
  target <- model_parameters(x$models)
  target[c("zc1", "zc2", "zc3")] <- c(6, 12, 18)
  observed <- predict_phenology(x$models, x$seasons, target)
  observed[[1]] <- list(observed[[1]][1], NULL)
  observed[[2]] <- list(NULL, observed[[2]][2])
  baseline <- phenology_rss_stages(model_parameters(x$models), x$models, x$seasons, observed)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 30, trace = FALSE) else list(max.call = 300)
    fit <- fit_phenology(x$models, phenology_rss_stages,
                         lower = c(zc = c(1, 8, 14)), upper = c(zc = c(8, 14, 22)),
                         seasons = x$seasons, observed = observed,
                         optimizer = optimizer, seed = 17, max_iterations = 20, control = control)
    expect_lt(fit$value, baseline)
    expect_equal(fit$value, phenology_rss_stages(fit$par, x$models, x$seasons, observed))
    expect_true(all(diff(fit$par[c("zc1", "zc2", "zc3")]) >= 0))
    expect_identical(model_parameters(fit$model), fit$par)
    fixed <- setdiff(names(target), c("zc1", "zc2", "zc3"))
    expect_identical(fit$par[fixed], target[fixed])
    for (i in 1:3) expect_identical(fit$model[[i]]$parameters[["zc"]], fit$par[[paste0("zc", i)]])
    expect_equal(predict_phenology(fit$model, x$seasons), predict_phenology(x$models, x$seasons, fit$par))
    expect_identical(as_pheno_models(fit), as_pheno_models(fit$model))
  }
})
