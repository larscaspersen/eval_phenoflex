test_that("both PhenoFlex wrappers reject either energy outside the Q10 limits", {
  testthat::local_mocked_bindings(
    PhenoFlex = function(...) stop("must not predict with invalid Q10"),
    .package = "chillR"
  )
  parameters <- c(40, 190, 0.5, 25, rep(log(2) * 297 * 279 / 10, 2),
                  6319.5, 5.939917e13, 4, 36, 4, 1.6)
  for (wrapper in list(custom_PhenoFlex_GDHwrapper, custom_PhenoFlex_GDHwrapper_v2)) {
    for (energy in 5:6) {
      for (q10 in c(1.4, 3.6)) {
        invalid <- parameters
        invalid[energy] <- log(q10) * 297 * 279 / 10
        expect_identical(wrapper(NULL, invalid, constraints = TRUE), NA)
      }
    }
  }
})

test_that("both PhenoFlex wrappers accept the Q10 boundaries and preserve the opt-out", {
  calls <- 0L
  testthat::local_mocked_bindings(
    PhenoFlex = function(temp, ...) {
      calls <<- calls + 1L
      list(bloomindex = 12L, y = rep(0, length(temp)),
           z = c(rep(0, 11), rep(200, length(temp) - 11)))
    },
    .package = "chillR"
  )
  weather <- data.frame(Temp = rep(10, 48), JDay = rep(1:2, each = 24), Year = 2024)
  parameters <- c(40, 190, 0.5, 25, rep(log(2) * 297 * 279 / 10, 2),
                  6319.5, 5.939917e13, 4, 36, 4, 1.6)
  for (wrapper in list(custom_PhenoFlex_GDHwrapper, custom_PhenoFlex_GDHwrapper_v2)) {
    for (energy in 5:6) {
      for (q10 in c(1.5, 2, 3.5)) {
        valid <- parameters
        valid[energy] <- log(q10) * 297 * 279 / 10
        before <- calls
        result <- wrapper(weather, valid, constraints = TRUE)
        expect_identical(calls, before + 1L)
        if (is.list(result)) {
          expect_equal(result$JDay, 1)
        } else {
          expect_equal(result, 1)
        }
      }
    }
    invalid <- parameters
    invalid[5:6] <- rep(log(4) * 297 * 279 / 10, 2)
    before <- calls
    expect_no_error(wrapper(weather, invalid, constraints = FALSE))
    expect_identical(calls, before + 1L)
  }
})
