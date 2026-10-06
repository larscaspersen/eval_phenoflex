test_that("date conversion handles leap years and repeated days in different years", {
  expect_equal(return_JDay(12, rep(60, 24), rep(2024, 24)), 60)
  expect_equal(return_JDay(24, rep(60, 24), rep(2024, 24)), 60.5)
  expect_equal(return_JDay(1, c(366, 1), c(2024, 2025)), 0)
  expect_equal(return_JDay(1, c(365, 1), c(2023, 2024)), 0)
  expect_equal(return_JDay(1, c(60, 60), c(2023, 2024)), -305)
  expect_true(is.na(return_JDay(0, 1, 2024)))
  expect_true(is.na(return_JDay(NA_real_, 1, 2024)))
  expect_error(return_JDay(3, 1:2, c(2024, 2024)), "index")
})

test_that("optimizer adapters preserve scalar RSS and reject silent recycling", {
  model <- function(x, par) x
  expect_equal(eval_all_daoptim(1, model, c(10, 12), list(11, 14)), 5)
  expect_equal(eval_fixed_daoptim(1, model, c(10, 12), list(11, 14)), 5)
  expect_equal(eval_all_daoptim(1, model, c(10, 12), list(NA, 14), na_penalty = 20), 104)
  expect_error(eval_all_daoptim(1, model, 10, list(11, 14)), "one finite")
  expect_error(eval_all_daoptim(1, model, 10, list(c(11, 14))), "one numeric")
  expect_equal(eval_all_daoptim(rep(1, 12), function(...) stop("must not run"),
                               10, list(11), check_constrain = TRUE), Inf)
})

test_that("GenSA rejects excessive E0 Q10", {
  par <- c(40, 190, 0.5, 25, log(4) * 297 * 279 / 10,
           log(2) * 297 * 279 / 10, 6319.5, 5.939917e13, 4, 36, 4, 1.6)
  expect_equal(evaluation_function_gensa(par, function(...) stop("must not run"),
    100, list(NULL), intermed_par = FALSE), 1e10)
})

test_that("the detailed PhenoFlex wrapper preserves its list interface", {
  weather <- data.frame(Temp = rep(10, 48), JDay = rep(1:2, each = 24), Year = 2024)
  par <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e13,
           4, 36, 4, 1.6)
  out <- custom_PhenoFlex_GDHwrapper_v2(weather, par)
  expect_named(out, c("JDay", "chill_heat"))
  expect_equal(dim(out$chill_heat), c(48L, 4L))
  expect_equal(out$JDay, custom_PhenoFlex_GDHwrapper(weather, par))
})
