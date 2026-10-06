test_that("stable residuals match independent high-precision reference equations", {
  references <- jsonlite::read_json(test_path("fixtures", "chill_residual_reference.json"))
  for (r in references) {
    actual <- solve_nle(unlist(r$x), unlist(r$params))
    expect_true(all(is.finite(actual)), info = r$name)
    expect_equal(actual, unlist(r$residuals), tolerance = 1e-10, info = r$name)
  }
})

test_that("ordinary residuals and converted parameters retain legacy results", {
  legacy <- new.env(parent = baseenv())
  sys.source(test_path("fixtures", "solve_nle_legacy.R"), legacy)
  sys.source(test_path("fixtures", "convert_parameters_legacy.R"), legacy)
  for (p in list(c(279, 286.1, 47.7, 28), c(279, 287, 28, 26), c(277, 285, 40, 30))) {
    expect_equal(solve_nle(c(500, 15000), p), legacy$solve_nle(c(500, 15000), p),
                 tolerance = 1e-10)
    par <- c(40, 190, .5, 25, p, 4, 36, 4, 1.6)
    old <- suppressWarnings(legacy$convert_parameters(par))
    expect_false(is.list(old))
    new <- convert_parameters(par)
    expect_false(is.list(new))
    # Relative checks per coefficient (A1 is many orders larger than E0).
    expect_equal(new[5:8] / old[5:8], rep(1, 4), tolerance = 1e-6)
    expect_equal(new[-(5:8)], par[-(5:8)])
    expect_lt(max(abs(solve_nle(new[5:6], p))), 1e-8)
  }
})

test_that("transformed solver enforces representable ordered energies", {
  p <- c(279, 286.1, 47.7, 28)
  expect_equal(solve_nle_transformed(log(c(500, 14500)), p),
               solve_nle(c(500, 15000), p), tolerance = 1e-12)
  for (z in list(c(1000, 1), c(-1000, 1), c(1, -1000), c(700, -700)))
    expect_identical(solve_nle_transformed(z, p), c(Inf, Inf))
  fit <- .fit_chill_parameters(p)
  expect_equal(fit$termcd, 1L)
  expect_gt(fit$x[1], 0)
  expect_gt(fit$x[2], fit$x[1])
  expect_equal(fit$x, c(exp(fit$z[1]), exp(fit$z[1]) + exp(fit$z[2])))
  expect_lt(max(abs(fit$fvec)), 1e-8)
})

test_that("domain errors are explicit and calibration failure policies are preserved", {
  p <- c(279, 286.1, 47.7, 28)
  for (x in list(c(0, 1), c(2, 1), c(1, 1), c(NA, 2), c(1, Inf), 1))
    expect_error(solve_nle(x, p), "energies")
  for (bad in list(c(279, 279, 20, 20), c(279, 297, 20, 20),
                   c(279, 300, 20, 20), c(279, 286, 0, 20), c(279, 286, 20, NA))) {
    expect_error(solve_nle(c(500, 15000), bad), "params")
    par <- c(40, 190, .5, 25, bad, 4, 36, 4, 1.6)
    expect_identical(convert_parameters(par), list(F = 1e6, g = rep(1e6, 5)))
    expect_true(all(is.na(convert_parameters(par, "NA")[5:8])))
  }
  expect_error(solve_nle_transformed(c(NA, 1), p), "z must")
  expect_error(convert_parameters(1:4), "length 12")
})

test_that("calibration adapters receive physical parameters and preserve penalties", {
  par <- c(40, 190, .5, 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  expected <- convert_parameters(par)
  callback <- function(x, par) {
    expect_equal(par, expected, tolerance = 1e-10)
    100
  }
  result <- evaluation_function_meigo_nonlinear(par, callback, 102, list(NULL))
  expect_equal(result$F, 4)
  expect_equal(eval_all_daoptim(par, callback, 102, list(NULL), intermed_chill = TRUE), 4)
  reduced <- par[-c(5, 10)]
  expect_equal(eval_phenoflex_single(reduced, callback, 102, list(NULL))$F, 4)
  expect_equal(eval_phenoflex_single_twofixed(reduced, callback, 102, list(NULL))$F, 4)
  expect_equal(evaluation_function_gensa(par, callback, 102, list(NULL)), 4)
  expect_equal(eval_phenoflex_combined(reduced, callback, list(102),
                                      list(list(NULL)), ncult = 1)$F, 4)
  expect_equal(evaluation_function_meigo_nonlinear_combined(par, callback,
    list(data.frame(pheno = 102)), list(list(NULL)), n_cult = 1)$F, 4)
  expect_equal(evaluation_function_meigo_nonlinear_fixed(par[c(1:4, 6, 7, 9)], callback,
    list(data.frame(pheno = 102)), list(list(NULL)), n_cult = 1,
    theta_star = 279, pie_c = 28, Tc = 36, Tb = 4, slope = 1.6)$F, 4)
  stage_x <- c(40, 100, 190, 300, .5, 25, 286.1, 47.7, 28, 4, 4, 1.6)
  expect_error(eval_phenoflex_three_stages(stage_x, function(x, par) {
    expect_equal(par, expected, tolerance = 1e-10)
    stop("conversion reached")
  }, data.frame(), list(NULL)), "conversion reached")
  par[6] <- 300
  never <- function(...) stop("Model should not run for invalid conversion")
  expect_equal(evaluation_function_gensa(par, never, 102, list(NULL)), 1e10)
  expect_equal(eval_all_daoptim(par, never, 102, list(NULL), intermed_chill = TRUE), Inf)
  expect_equal(evaluation_function_meigo_nonlinear(par, never, 102, list(NULL))$F, 1e6)
})

test_that("step stagnation is not mistaken for a root", {
  testthat::local_mocked_bindings(nleqslv = function(...) {
    list(x = log(c(500, 14500)), fvec = c(0, 0), termcd = 2L)
  }, .package = "nleqslv")
  fit <- .fit_chill_parameters(c(279, 286.1, 47.7, 28))
  expect_equal(fit$termcd, 3L)
  expect_gt(max(abs(fit$fvec)), 1e-8)
})
