test_that("DEoptim uses its returned optimum and counts without reevaluating it", {
  skip_if_not_installed("DEoptim")
  model <- pheno_model()
  calls <- 0L
  testthat::local_mocked_bindings(DEoptim = function(fn, lower, upper, control) {
    list(optim = list(bestmem = c(yc = 45), bestval = 0, nfeval = 20L))
  }, .package = "DEoptim")
  fit <- fit_phenology(model, function(p, m) {
    calls <<- calls + 1L
    (p[["yc"]] - 45)^2
  }, c(yc = 10), c(yc = 80))
  expect_identical(fit$par[["yc"]], 45)
  expect_identical(fit$value, 0)
  expect_identical(fit$diagnostics$evaluations, 21L)
  expect_identical(calls, 1L)
  expect_identical(fit$model$parameters, fit$par)
})

test_that("PSOCK DEoptim matches serial fitting and preserves caller clusters and RNG", {
  skip_if_not_installed("DEoptim")
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl))
  initialize_worker <- function() .libPaths(.Library)
  environment(initialize_worker) <- baseenv()
  parallel::clusterCall(cl, initialize_worker)
  model <- pheno_model("sequential", chill_dynamic("kinetic"))
  initial <- model$parameters
  initial["Tf"] <- 2
  target <- c(yc = 45, zc = 250)
  loss <- function(parameters, model, observed, scale) {
    sum(((parameters[c("yc", "zc")] - observed) / scale)^2)
  }
  run <- function(control) fit_phenology(model, loss, lower = c(zc = 100, yc = 10),
    upper = c(yc = 80, zc = 400), parameters = initial, observed = target,
    scale = c(1, 10), prediction = FALSE, seed = 17, max_iterations = 5,
    control = c(list(NP = 20, trace = FALSE), control))
  set.seed(901)
  before <- .Random.seed
  serial <- run(list())
  fit <- run(list(cluster = cl, parallelType = "parallel"))
  expect_identical(.Random.seed, before)
  expect_identical(fit$par, serial$par)
  expect_identical(fit$value, serial$value)
  expect_lt(fit$value, loss(initial, model, target, c(1, 10)))
  expect_equal(fit$diagnostics$evaluations, serial$diagnostics$evaluations)
  expect_equal(fit$diagnostics$evaluations, 1 + fit$diagnostics$optimizer_result$optim$nfeval)
  expect_identical(fit$model$parameters, fit$par)
  expect_identical(fit$predictions$observed, target)
  fixed <- !names(initial) %in% c("yc", "zc")
  expect_identical(fit$par[fixed], initial[fixed])
  expect_identical(run(list(cluster = cl))$par, fit$par)
  expect_length(parallel::clusterCall(cl, Sys.getpid), 2)
})

test_that("parallel DEoptim retains a better baseline and handles infinite losses", {
  skip_if_not_installed("DEoptim")
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl))
  model <- pheno_model()
  control <- list(cluster = cl, NP = 10, trace = FALSE)
  run <- function(loss) fit_phenology(model, loss, c(yc = 10), c(yc = 80),
                                     max_iterations = 2, control = control)
  for (other_loss in c(1, Inf)) {
    fit <- run(function(p, m) if (p[["yc"]] == m$parameters[["yc"]]) 0 else other_loss)
    expect_identical(fit$par, model$parameters)
    expect_identical(fit$value, 0)
    expect_gt(fit$diagnostics$evaluations, 1)
  }
  fit <- run(function(p, m) if (p[["yc"]] == m$parameters[["yc"]]) Inf else
    (p[["yc"]] - 45)^2)
  expect_true(is.finite(fit$value))
  expect_equal(fit$value, (fit$par[["yc"]] - 45)^2)
  expect_error(run(function(p, m) Inf), "No candidate")
})

test_that("PSOCK workers execute the native phenology evaluator and predictions", {
  skip_if_not_installed("DEoptim")
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl))
  model <- pheno_model("sequential", chill_dynamic("kinetic"),
                       parameters = c(yc = 0.1, zc = 5))
  weather <- data.frame(Temp = rep(8, 480), Year = 2008, JDay = rep(50:69, each = 24))
  target <- model$parameters
  target["zc"] <- 15
  observed <- predict_phenology(model, weather, target)
  fit <- fit_phenology(model, phenology_rss, lower = c(zc = 1), upper = c(zc = 20),
    seasons = list(weather), observed = observed, seed = 17, max_iterations = 5,
    control = list(cluster = cl, NP = 10, trace = FALSE))
  expect_lt(fit$value, phenology_rss(model$parameters, model, list(weather), observed))
  expect_equal(fit$value, phenology_rss(fit$par, model, list(weather), observed))
  expect_equal(fit$predictions$predicted, predict_phenology(fit, list(weather)))
})

test_that("parallel workers exclude invalid domains and propagate evaluator errors", {
  skip_if_not_installed("DEoptim")
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl))
  model <- pheno_model()
  control <- list(cluster = cl, NP = 10, trace = FALSE,
                  initialpop = matrix(c(-10, seq(10, 80, length.out = 9)), ncol = 1))
  run <- function(loss) fit_phenology(model, loss, c(yc = -20), c(yc = 80),
                                     max_iterations = 1, control = control)
  fit <- run(function(p, m) {
    if (p[["yc"]] <= 0) stop("invalid candidate was evaluated")
    (p[["yc"]] - 45)^2
  })
  expect_true(is.finite(fit$value))
  expect_equal(fit$diagnostics$evaluations, 1 + fit$diagnostics$optimizer_result$optim$nfeval)
  for (bad in list(NA_real_, -Inf, c(1, 2), matrix(1))) {
    expect_error(run(function(p, m) {
      if (p[["yc"]] == m$parameters[["yc"]]) return(1)
      bad
    }), "one numeric loss")
  }
  expect_error(run(function(p, m) {
    if (p[["yc"]] == m$parameters[["yc"]]) return(1)
    stop("broken worker evaluator")
  }), "broken worker evaluator")
  expect_length(parallel::clusterCall(cl, Sys.getpid), 2)
})

test_that("automatic parallel clusters are stopped on success and worker errors", {
  skip_if_not_installed("DEoptim")
  skip_if_not_installed("parallelly")
  cl <- NULL
  create <- parallel::makePSOCKcluster
  testthat::local_mocked_bindings(makePSOCKcluster = function(...) {
    cl <<- create(2)
    cl
  }, .package = "parallel")
  model <- pheno_model()
  run <- function(loss, type) fit_phenology(model, loss, c(yc = 10), c(yc = 80),
    max_iterations = 1, control = list(parallelType = type, NP = 10, trace = FALSE))
  for (type in list("parallel", "auto", 1)) {
    fit <- run(function(p, m) (p[["yc"]] - 45)^2, type)
    expect_equal(fit$value, (fit$par[["yc"]] - 45)^2)
    expect_error(parallel::clusterCall(cl, Sys.getpid))
  }
  expect_error(run(function(p, m) {
    if (p[["yc"]] == m$parameters[["yc"]]) return(1)
    stop("broken worker evaluator")
  }, "parallel"), "broken worker evaluator")
  expect_error(parallel::clusterCall(cl, Sys.getpid))
  expect_error(fit_phenology(model, function(p, m) 1, c(yc = 10), c(yc = 80),
    max_iterations = 1, control = list(parallelType = "parallel", trace = FALSE,
                                     packages = "evalpheno_missing_worker_package")),
    "evalpheno_missing_worker_package")
  expect_error(parallel::clusterCall(cl, Sys.getpid))
})

test_that("PSOCK controls load custom packages and export global evaluator dependencies", {
  skip_if_not_installed("DEoptim")
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl))
  # These names intentionally live outside the evaluator closure.
  loss <- eval(quote(function(p, m) .evalpheno_parallel_helper(p)), .GlobalEnv)
  helper <- eval(quote(function(p) {
    purrr::map_dbl(as.list(p["yc"]), function(z) (z - .evalpheno_parallel_target)^2)[[1]]
  }), .GlobalEnv)
  exports <- list(.evalpheno_parallel_helper = helper, .evalpheno_parallel_target = 45)
  existed <- vapply(names(exports), exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
  previous <- mget(names(exports)[existed], envir = .GlobalEnv)
  on.exit({
    rm(list = names(exports), envir = .GlobalEnv)
    list2env(previous, envir = .GlobalEnv)
  }, add = TRUE)
  list2env(exports, envir = .GlobalEnv)
  model <- pheno_model()
  fit <- fit_phenology(model, loss, c(yc = 10), c(yc = 80), max_iterations = 3,
    control = list(cluster = cl, NP = 10, trace = FALSE, packages = "purrr",
                   parVar = c(".evalpheno_parallel_helper", ".evalpheno_parallel_target")))
  expect_equal(fit$value, (fit$par[["yc"]] - 45)^2)
  expect_lt(fit$value, loss(model$parameters, model))
  expect_true(all(unlist(parallel::clusterCall(cl, function() "package:purrr" %in% search()))))
})

test_that("foreach DEoptim uses a registered PSOCK backend", {
  skip_if_not_installed("DEoptim")
  skip_if_not_installed("doParallel")
  cl <- parallel::makePSOCKcluster(2)
  on.exit({foreach::registerDoSEQ(); parallel::stopCluster(cl)})
  initialize_worker <- function(paths) .libPaths(paths)
  environment(initialize_worker) <- baseenv()
  parallel::clusterCall(cl, initialize_worker, .libPaths())
  doParallel::registerDoParallel(cl)
  model <- pheno_model()
  run <- function(type) fit_phenology(model, function(p, m, target) (p[["yc"]] - target)^2,
    c(yc = 10), c(yc = 80), target = 45, seed = 17, max_iterations = 3,
    control = list(parallelType = type, NP = 10, trace = FALSE))
  serial <- run("none")
  fit <- run("foreach")
  expect_identical(fit$par, serial$par)
  expect_identical(fit$value, serial$value)
  expect_equal(fit$diagnostics$evaluations, serial$diagnostics$evaluations)
})
