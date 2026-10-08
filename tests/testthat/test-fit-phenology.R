fit_example <- function() {
  model <- pheno_model(structure = "sequential", chill = chill_dynamic("kinetic"))
  lower <- upper <- model$parameters
  lower[c("yc", "zc")] <- c(10, 100)
  upper[c("yc", "zc")] <- c(80, 400)
  loss <- function(parameters, model, target) {
    sum(((parameters[c("yc", "zc")] - target) / c(1, 10))^2)
  }
  list(model = model, lower = lower, upper = upper, loss = loss,
       target = c(yc = 45, zc = 250))
}

test_that("both optimizers minimize the same named loss and retain fixed parameters", {
  x <- fit_example()
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 30, trace = FALSE) else
      list(max.call = 3000, simple.function = TRUE)
    fit <- fit_phenology(x$model, x$loss, x$lower[c("zc", "yc")], x$upper[c("yc", "zc")],
                         optimizer = optimizer, seed = 17, control = control,
                         max_iterations = 80, target = x$target)
    expect_s3_class(fit, "phenology_fit")
    expect_identical(names(fit$par), parameter_schema(x$model)$name)
    expect_lt(fit$value, 0.1)
    expect_equal(fit$value, x$loss(fit$par, x$model, x$target))
    fixed <- x$lower == x$upper
    expect_identical(fit$par[fixed], x$model$parameters[fixed])
    expect_identical(fit$model$parameters, fit$par)
    expect_identical(x$model$parameters, default_parameters(x$model))
    expect_gt(fit$diagnostics$evaluations, 1L)
    expect_identical(fit$diagnostics$termination, "optimizer")
    expect_false(is.null(fit$diagnostics$optimizer_result))
    expect_equal(fit$calibration_settings$control[[if (optimizer == "DEoptim") "itermax" else "maxit"]], 80)
  }
})

test_that("seeded optimization is reproducible and restores the caller's RNG", {
  x <- fit_example()
  set.seed(901)
  before <- .Random.seed
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 20, itermax = 3, trace = FALSE) else
      list(maxit = 3, max.call = 100, seed = -999)
    run <- function() fit_phenology(x$model, x$loss, x$lower, x$upper,
                                    optimizer = optimizer, seed = 0, control = control,
                                    target = x$target)
    a <- run()
    expect_identical(.Random.seed, before)
    b <- run()
    expect_identical(a$par, b$par)
    expect_identical(a$value, b$value)
    expect_identical(a$diagnostics$evaluations, b$diagnostics$evaluations)
    expect_identical(.Random.seed, before)
    if (optimizer == "GenSA") expect_identical(a$calibration_settings$control$seed, -1L)
  }
})

test_that("GenSA enforces its native time limit and returns native diagnostics", {
  skip_if_not_installed("GenSA")
  x <- fit_example()
  calls <- 0L
  slow_loss <- function(parameters, model) {
    calls <<- calls + 1L
    Sys.sleep(0.02)
    (parameters[["yc"]] - 45)^2
  }
  fit <- fit_phenology(x$model, slow_loss, x$lower, x$upper,
                       optimizer = "GenSA", max_iterations = 1000,
                       max_seconds = 0.08, control = list(max.call = 1000))
  expect_identical(fit$diagnostics$termination, "optimizer")
  expect_gte(fit$diagnostics$elapsed, 0.08)
  expect_lt(calls, 1000)
  expect_equal(fit$value, (fit$par[["yc"]] - 45)^2)
  expect_equal(fit$diagnostics$evaluations, calls)
  expect_false(is.null(fit$diagnostics$optimizer_result))
})

test_that("GenSA receives the full time budget even after a slow baseline", {
  skip_if_not_installed("GenSA")
  x <- fit_example()
  received_time <- NULL
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    received_time <<- control$max.time
    list(par = par, value = fn(par))
  }, .package = "GenSA")
  fit <- fit_phenology(x$model, function(p, m) {Sys.sleep(0.02); sum(p)},
                       x$lower, x$upper, optimizer = "GenSA", max_seconds = 0.01,
                       control = list(max.time = 99))
  expect_identical(received_time, 0.01)
  expect_false(is.null(fit$diagnostics$optimizer_result))
  expect_identical(fit$diagnostics$evaluations, 2L)
  fit_phenology(x$model, function(p, m) sum(p), x$lower, x$upper,
                optimizer = "GenSA", control = list(max.time = 0.05))
  expect_identical(received_time, 0.05)
})

test_that("DEoptim rejects time limits it does not support", {
  x <- fit_example()
  expect_error(fit_phenology(x$model, x$loss, x$lower, x$upper, max_seconds = 1),
               "DEoptim has no native time limit")
  expect_error(fit_phenology(x$model, x$loss, x$lower, x$upper,
                             control = list(max.time = 1)),
               "DEoptim has no native time limit")
})

test_that("fixed models and population specifications return fitted models", {
  model <- population_pheno_model(seed = 17)
  p <- model$model$parameters
  fit <- fit_phenology(model, function(parameters, model) sum(parameters), p, p)
  expect_identical(fit$diagnostics$termination, "fixed")
  expect_identical(fit$diagnostics$evaluations, 1L)
  expect_null(fit$diagnostics$optimizer_result)
  expect_identical(fit$par, p)
  expect_identical(fit$model, model)
  skip_if_not_installed("DEoptim")
  lower <- c(yc = 10)
  upper <- c(yc = 80)
  fitted <- fit_phenology(model, function(parameters, model) (parameters[["yc"]] - 45)^2,
                          lower, upper, max_iterations = 3,
                          control = list(NP = 10, trace = FALSE))
  expect_identical(fitted$model$model$parameters, fitted$par)
  expect_identical(fitted$model$n, model$n)
  expect_identical(fitted$model$sd, model$sd)
  expect_identical(fitted$model$seed, model$seed)
})

test_that("omitted parameters use initial values throughout both optimizer searches", {
  x <- fit_example()
  initial <- x$model$parameters
  initial["Tf"] <- 2
  constants <- initial[!names(initial) %in% c("yc", "zc")]
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    constants_preserved <- TRUE
    complete_names <- TRUE
    evaluation <- function(parameters, model, target) {
      constants_preserved <<- constants_preserved &&
        identical(parameters[names(constants)], constants)
      complete_names <<- complete_names && identical(names(parameters), names(initial))
      sum(((parameters[c("yc", "zc")] - target) / c(1, 10))^2)
    }
    control <- if (optimizer == "DEoptim") list(NP = 20, trace = FALSE) else
      list(max.call = 100)
    fit <- fit_phenology(x$model, evaluation, lower = c(yc = 10, zc = 100),
                         upper = c(zc = 400, yc = 80), parameters = rev(initial),
                         optimizer = optimizer, control = control,
                         max_iterations = 3, target = x$target)
    expect_true(constants_preserved)
    expect_true(complete_names)
    expect_identical(fit$par[names(constants)], constants)
    expect_identical(fit$calibration_settings$lower[names(constants)], constants)
    expect_identical(fit$calibration_settings$upper[names(constants)], constants)
    expect_identical(fit$model$parameters, fit$par)
    expect_identical(x$model$parameters[["Tf"]], 4)
    native_parameters <- if (optimizer == "DEoptim") fit$diagnostics$optimizer_result$optim$bestmem else
      fit$diagnostics$optimizer_result$par
    expect_length(native_parameters, 2L)
  }
})

test_that("empty bounds and equal subset bounds keep every parameter constant", {
  x <- fit_example()
  initial <- x$model$parameters
  initial["Tf"] <- 2
  for (bounds in list(numeric(), initial["yc"])) {
    fit <- fit_phenology(x$model, function(parameters, model) sum(parameters),
                         lower = bounds, upper = bounds, parameters = initial)
    expect_identical(fit$par, initial)
    expect_identical(fit$calibration_settings$lower, initial)
    expect_identical(fit$calibration_settings$upper, initial)
    expect_identical(fit$diagnostics$evaluations, 1L)
    expect_identical(fit$diagnostics$termination, "fixed")
    expect_null(fit$diagnostics$optimizer_result)
  }
})

test_that("seed handling preserves absent RNG state and permits unseeded evaluation", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) before <- .Random.seed
  on.exit({
    if (had_seed) assign(".Random.seed", before, .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  model <- pheno_model()
  p <- model$parameters
  fit_phenology(model, function(p, m) runif(1), p, p, seed = 17)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  set.seed(901)
  state <- .Random.seed
  fit_phenology(model, function(p, m) runif(1), p, p, seed = NULL)
  expect_false(identical(.Random.seed, state))
})

test_that("bad bounds, starting points, limits and controls fail clearly", {
  x <- fit_example()
  run <- function(lower = x$lower, upper = x$upper, ...)
    fit_phenology(x$model, x$loss, lower, upper, target = x$target, ...)
  expect_error(run(lower = unname(x$lower)), "named numeric vector")
  expect_error(run(lower = c(unknown = 1), upper = c(unknown = 2)), "model parameter names")
  expect_error(run(lower = c(yc = 10), upper = c(zc = 400)), "same subset")
  expect_error(run(lower = c(yc = 10), upper = numeric()), "same subset")
  bad <- x$lower
  names(bad)[1] <- names(bad)[2]
  expect_error(run(lower = bad), "named numeric vector")
  bad <- x$lower
  bad[1] <- Inf
  expect_error(run(lower = bad), "finite named")
  expect_error(run(lower = x$upper, upper = x$lower), "lower must not exceed")
  bad <- x$upper
  bad["yc"] <- 20
  expect_error(run(upper = bad), "Initial parameters")
  for (limit in list(0, -1, NA_real_, Inf, c(1, 2))) {
    expect_error(run(max_iterations = limit), "max_iterations")
    expect_error(run(max_seconds = limit), "max_seconds")
  }
  expect_error(run(max_iterations = 1.5), "integer")
  expect_error(run(seed = -1), "seed")
  expect_error(run(control = list(3)), "named list")
  expect_error(run(control = list(trace = FALSE, trace = TRUE)), "duplicate")
  expect_error(run(control = list(cluster = "invalid")), "parallel cluster object")
  expect_error(run(control = list(itermax = 0)), "control\\$itermax")
  expect_error(run(optimizer = "GenSA", control = list(maxit = 0)), "control\\$maxit")
})

test_that("invalid domain candidates are excluded but evaluator errors propagate", {
  skip_if_not_installed("DEoptim")
  x <- fit_example()
  # The box includes yc <= 0, outside the model's domain.
  x$lower["yc"] <- -20
  fit <- fit_phenology(x$model, function(parameters, model) {
    expect_gt(parameters[["yc"]], 0)
    (parameters[["yc"]] - 45)^2
  }, x$lower, x$upper, seed = 17, max_iterations = 4,
  control = list(NP = 20, trace = FALSE))
  expect_true(is.finite(fit$value))
  for (bad in list(NA_real_, NaN, -Inf, c(1, 2), list(f = 1), "1")) {
    expect_error(fit_phenology(x$model, function(p, m) bad, x$lower, x$upper),
                 "one numeric loss")
  }
  set.seed(901)
  before <- .Random.seed
  expect_error(fit_phenology(x$model, function(p, m) stop("broken evaluator"),
                             x$lower, x$upper), "broken evaluator")
  expect_identical(.Random.seed, before)
  p <- x$model$parameters
  expect_error(fit_phenology(x$model, function(p, m) Inf, p, p), "No candidate")
})

test_that("both optimizers receive RSS from the explicitly supplied example evaluator", {
  model <- pheno_model("sequential", chill_dynamic("kinetic"))
  model$parameters[c("yc", "zc")] <- c(0.1, 5)
  lower <- upper <- model$parameters
  lower["yc"] <- 0.09
  upper["yc"] <- 0.11
  lower["zc"] <- 1
  upper["zc"] <- 20
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  target <- model$parameters
  target["zc"] <- 15
  observed <- predict_phenology(model, weather, target)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 20, trace = FALSE) else
      list(max.call = 200)
    fit <- fit_phenology(model, evaluation = phenology_rss,
                         lower = lower, upper = upper, seasons = list(weather),
                         observed = observed, optimizer = optimizer,
                         seed = 17, max_iterations = 20, control = control)
    expect_lt(fit$value, phenology_rss(model$parameters, model, list(weather), observed))
    expect_equal(fit$value, phenology_rss(fit$par, model, list(weather), observed))
    expect_equal(predict_phenology(fit$model, weather), observed, tolerance = 1 / 24)
  }
})

test_that("the fitter requires an explicitly supplied evaluation function", {
  model <- pheno_model()
  p <- model$parameters
  expect_error(fit_phenology(model, lower = p, upper = p),
               "evaluation must be supplied as a function")
})
