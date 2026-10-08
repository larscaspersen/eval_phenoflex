test_that("observations are explicit and forwarded with other evaluation arguments", {
  expect_true("observed" %in% names(formals(fit_phenology)))
  model <- pheno_model()
  dates <- c(year_a = 38, year_b = 43)
  calls <- 0L
  loss <- function(parameters, model, seasons, observed, offset) {
    calls <<- calls + 1L
    expect_identical(observed, dates)
    expect_identical(seasons, c("year_a", "year_b"))
    sum((parameters[["yc"]] - observed)^2) + offset
  }
  fit <- fit_phenology(model, loss, numeric(), numeric(),
                       seasons = c("year_a", "year_b"), observed = dates, offset = 5)
  expect_identical(fit$value, 18)
  expect_identical(calls, 1L)
})

test_that("nested observations and explicit NULL are passed without alteration", {
  model <- pheno_model()
  dates <- list(early = list(90, NULL), late = list(NULL, 100))
  fit <- fit_phenology(model, function(p, m, observed) {
    expect_identical(observed, dates)
    7
  }, numeric(), numeric(), observed = dates)
  expect_identical(fit$value, 7)
  fit <- fit_phenology(model, function(p, m, observed) {
    expect_false(missing(observed))
    expect_null(observed)
    3
  }, numeric(), numeric(), observed = NULL)
  expect_identical(fit$value, 3)
  # Evaluators that capture observations still work with the existing interface.
  fit <- fit_phenology(model, function(p, m) sum((p[["yc"]] - dates$early[[1]])^2),
                       numeric(), numeric())
  expect_identical(fit$value, 2500)
})

test_that("both optimizer objectives receive observations on every candidate", {
  model <- pheno_model()
  dates <- c(30, 40)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    run <- function(fn) {
      fn(c(yc = 30))
      value <- fn(c(yc = 35))
      list(optim = list(bestmem = c(yc = 35), bestval = value, nfeval = 2L))
    }
    testthat::local_mocked_bindings(DEoptim = function(fn, lower, upper, control) run(fn),
                                    .package = "DEoptim")
    testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) run(fn),
                                    .package = "GenSA")
    calls <- 0L
    fit <- fit_phenology(model, function(p, m, observed) {
      calls <<- calls + 1L
      expect_identical(observed, dates)
      sum((p[["yc"]] - observed)^2)
    }, c(yc = 10), c(yc = 80), optimizer = optimizer, observed = dates, max_iterations = 1)
    expect_identical(fit$value, 50)
    expect_identical(fit$par[["yc"]], 35)
    expect_identical(calls, 3L)
  }
})
