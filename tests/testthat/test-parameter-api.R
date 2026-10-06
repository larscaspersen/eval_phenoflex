test_that("parameter names retain the legacy vocabulary and data interface", {
  expect_identical(phenoflex_parnames_old, phenoflex_parnames_kinetic)
  expect_identical(phenoflex_parnames_new, phenoflex_parnames_characteristic)
  expect_identical(phenoflex_parnames_kinetic[5:8], c("E0", "E1", "A0", "A1"))
  expect_identical(phenoflex_parnames_characteristic[5:8], c("theta_star", "theta_c", "tau", "pie_c"))
  expect_identical(phenoflex_parnames_old[-(5:8)], phenoflex_parnames_new[-(5:8)])
  for (name in c("phenoflex_parnames_old", "phenoflex_parnames_new",
                 "phenoflex_parnames_kinetic", "phenoflex_parnames_characteristic")) {
    env <- new.env()
    utils::data(list = name, package = "evalpheno", envir = env)
    expect_identical(env[[name]], get(name))
  }
})

test_that("both conversion directions round trip and preserve other parameters", {
  legacy <- new.env(parent = baseenv())
  sys.source(test_path("fixtures", "kinetic_to_characteristic_legacy.R"), legacy)
  for (p in list(c(279, 286.1, 47.7, 28), c(279, 287, 28, 26), c(277, 285, 40, 30))) {
    characteristic <- c(40, 190, .5, 25, p, 4, 36, 4, 1.6)
    kinetic <- characteristic_to_kinetic(characteristic)
    restored <- kinetic_to_characteristic(kinetic)
    expect_equal(restored, characteristic, tolerance = 1e-6)
    expect_equal(characteristic_to_kinetic(restored) / kinetic, rep(1, 12), tolerance = 1e-6)
    expect_identical(restored[-(5:8)], characteristic[-(5:8)])
    expect_identical(kinetic[-(5:8)], characteristic[-(5:8)])
    expect_equal(restored, legacy$convert_parameters_old_to_new(kinetic), tolerance = 1e-12)
    expect_identical(convert_parameters(characteristic), kinetic)
    expect_identical(convert_parameters_old_to_new(kinetic), restored)
    names(characteristic) <- phenoflex_parnames_characteristic
    named_kinetic <- characteristic_to_kinetic(characteristic)
    expect_identical(names(named_kinetic), phenoflex_parnames_kinetic)
    expect_equal(kinetic_to_characteristic(named_kinetic), characteristic, tolerance = 1e-6)
  }
})

test_that("inverse conversion rejects invalid kinetic inputs", {
  expect_error(kinetic_to_characteristic(1:4), "twelve")
  for (chill in list(c(1, 1, 1, 1), c(0, 2, 1, 1), c(1, 2, -1, 1), c(1, Inf, 1, 1)))
    expect_error(kinetic_to_characteristic(c(rep(1, 4), chill, rep(1, 4))), "finite")
})
