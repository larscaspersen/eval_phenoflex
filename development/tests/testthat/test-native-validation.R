test_that("native structure kernels reject unsafe weather and non-hourly input", {
  for (fn in list(seq_model, parallel_model, po_model)) {
    expect_error(fn(numeric(0), numeric(0)), "at least two")
    expect_error(fn(8, 0), "at least two")
    expect_error(fn(c(8, 8, 8), 0:1), "same length")
    expect_error(fn(c(8, 8), 0:2), "same length")
    for (bad in c(NA_real_, NaN, Inf, -Inf)) {
      expect_error(fn(c(8, bad), 0:1), "temp.*finite")
      expect_error(fn(c(8, 8), c(0, bad)), "times.*finite")
    }
    for (times in list(c(0, 0), c(1, 0), c(0, 2), c(0, 0.5), c(0, 1, 3)))
      expect_error(fn(rep(8, length(times)), times), "consecutive hourly")
    expect_error(fn(c(8, -273), 0:1), "temp.*greater")
    expect_error(fn(c(281, 0), 0:1, deg_celsius = FALSE), "temp.*greater")

    # Fractional origins are valid: only spacing is constrained.
    expected <- fn(rep(8, 4), 0:3, basic_output = FALSE)
    expect_equal(fn(rep(8, 4), (0:3) + 100.25, basic_output = FALSE), expected)
    expect_equal(fn(rep(8, 4), (0:3) + 0.1, basic_output = FALSE), expected,
                 tolerance = 1e-12)
    # Match the legacy offset of 273 and convert all temperature thresholds.
    kelvin <- fn(rep(281, 4), 0:3, Tf = 277, Tb = 277, Tu = 298, Tc = 309,
                 deg_celsius = FALSE, basic_output = FALSE)
    expect_equal(kelvin, expected)
  }
})

test_that("native structure kernels validate shared and structure parameters", {
  for (fn in list(seq_model, parallel_model, po_model)) {
    run <- function(...) fn(rep(8, 4), 0:3, ...)
    for (name in c("yc", "A0", "A1", "E0", "E1", "slope", "Delta")) {
      for (bad in c(0, -1, Inf, NA_real_))
        expect_error(do.call(run, setNames(list(bad), name)), paste0(name, ".*positive"))
    }
    expect_error(run(E0 = 2, E1 = 1), "E0 < E1")
    for (name in c("Tf", "Tb", "Tu", "Tc")) {
      expect_error(do.call(run, setNames(list(NA_real_), name)), paste0(name, ".*finite"))
      expect_error(do.call(run, setNames(list(-273), name)), paste0(name, ".*greater"))
      expect_error(do.call(run, c(setNames(list(0), name), list(deg_celsius = FALSE))),
                   paste0(name, ".*greater"))
    }
    expect_error(run(Tb = 25), "Tb < Tu < Tc")
    expect_error(run(Tc = 25), "Tb < Tu < Tc")
  }
  for (fn in list(seq_model, parallel_model))
    expect_error(fn(rep(8, 4), 0:3, zc = 0), "zc.*positive")
  for (bad in c(-0.1, 1.1, Inf, NA_real_))
    expect_error(parallel_model(rep(8, 4), 0:3, kmin = bad), "kmin")
  expect_type(parallel_model(rep(8, 4), 0:3, kmin = 0), "list")
  expect_type(parallel_model(rep(8, 4), 0:3, kmin = 1), "list")
  expect_error(po_model(rep(8, 4), 0:3, b1 = 0), "b1.*positive")
  for (name in c("b2", "b3", "ol")) {
    for (bad in c(-1, Inf, NA_real_))
      expect_error(do.call(po_model, c(list(temp = rep(8, 4), times = 0:3),
                                      setNames(list(bad), name))), paste0(name, ".*non-negative"))
  }
  expect_type(po_model(rep(8, 4), 0:3, b2 = 0, b3 = 0, ol = 0), "list")
})
