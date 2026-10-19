test_that("sequential heat starts in the interval that reaches the chill threshold", {
  # At 8 C (281 K in the legacy kernel), choose k1=1 and xs=2.
  # With Tf=8 the converted fraction is 1/2, so the first interval
  # produces 1-exp(-1) chill units, exceeding yc=0.5.
  A1 <- exp(2 / 281)
  args <- list(temp = rep(8, 4), times = 0:3, yc = 0.5, zc = 1,
               E0 = 1, E1 = 2, A0 = 2 * A1 * exp(-1 / 281), A1 = A1,
               Tf = 8, Tb = 4, Tu = 25, Tc = 36, basic_output = FALSE)
  expected_chill <- 1 - exp(-1)
  expected_heat <- (25 - 4) / 2 * (1 - cos(pi * (8 - 4) / (25 - 4)))
  out <- do.call(seq_model, args)
  expect_equal(out$y[2], expected_chill, tolerance = 1e-12)
  expect_equal(out$z[1], 0)
  expect_equal(out$z[2], expected_heat, tolerance = 1e-12)
  expect_equal(out$bloomindex, 2L)

  # The same crossing must be detected even at the final available interval.
  args$temp <- rep(8, 2); args$times <- 0:1
  expect_equal(do.call(seq_model, args)$bloomindex, 2L)

  # Equality is sufficient; no additional chill-producing interval is needed.
  args$yc <- out$y[2]
  expect_equal(do.call(seq_model, args)$bloomindex, 2L)

  # Retain zero heat before chill completion.
  args$yc <- 1e9
  no_chill <- do.call(seq_model, args)
  expect_equal(no_chill$z, c(0, 0))
  expect_equal(no_chill$bloomindex, 0L)

  # Reaching chill does not create heat below the heat base temperature.
  args$yc <- 0.5; args$Tb <- 10
  no_heat <- do.call(seq_model, args)
  expect_equal(no_heat$y[2], expected_chill, tolerance = 1e-12)
  expect_equal(no_heat$z, c(0, 0))
  expect_equal(no_heat$bloomindex, 0L)
})
