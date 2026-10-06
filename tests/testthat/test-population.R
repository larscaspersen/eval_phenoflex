population_par <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e13,
                    4, 36, 4, 1.6)
population_weather <- data.frame(Temp = c(rep(6, 1200), rep(20, 1200)),
                                  JDay = rep(1:100, each = 24))

test_that("population sampling handles zero variation and preserves RNG", {
  expect_equal(sample_bud_population(3)$yc_pop, rep(40, 3))
  expect_equal(sample_bud_population(1)$zc_pop, 190)
  set.seed(17)
  before <- .Random.seed
  a <- sample_bud_population(100, yc_sd = 2, zc_sd = 3)
  expect_identical(.Random.seed, before)
  expect_identical(a, sample_bud_population(100, yc_sd = 2, zc_sd = 3))
  expect_false(identical(a, sample_bud_population(100, yc_sd = 2, zc_sd = 3, seed = 8)))
  expect_error(sample_bud_population(n = 0), "n")
  expect_error(sample_bud_population(yc_sd = -1), "yc_sd")
  expect_error(sample_bud_population(dist_chill = "typo"), "arg")
})

test_that("skewed sampling honours seeds independently for both requirements", {
  skip_if_not_installed("sn")
  a <- sample_bud_population(200, yc_sd = 1, zc_sd = 1,
    dist_chill = "normal_skewed", dist_heat = "normal_skewed", skew = c(2, 2))
  expect_false(isTRUE(all.equal(a$yc_pop - 40, a$zc_pop - 190)))
  expect_identical(a, sample_bud_population(200, yc_sd = 1, zc_sd = 1,
    dist_chill = "normal_skewed", dist_heat = "normal_skewed", skew = c(2, 2)))
})

test_that("single-bud model agrees with chillR PhenoFlex", {
  p <- population_par
  actual <- phenoflex_population(population_weather, p, n = 1, stop_at_zc = FALSE)
  reference <- chillR::PhenoFlex(temp = population_weather$Temp,
    times = seq_len(nrow(population_weather)), yc = p[1], zc = p[2], s1 = p[3],
    Tu = p[4], E0 = p[5], E1 = p[6], A0 = p[7], A1 = p[8], Tf = p[9],
    Tc = p[10], Tb = p[11], slope = p[12], stopatzc = FALSE, basic_output = FALSE)
  expect_equal(actual$y, reference$y)
  expect_equal(as.vector(actual$z), reference$z)
  expect_equal(actual$bloomindex, reference$bloomindex)
})

test_that("constant populations replicate a single bud and compatibility helper", {
  single <- phenoflex_population(population_weather, population_par, n = 1)
  pop <- phenoflex_population(population_weather, population_par, n = 3)
  expect_equal(pop$bloomindex, rep(single$bloomindex, 3))
  expect_equal(pop$z[, 2], as.vector(single$z))
  expect_identical(pop, helper_run_pop_model(population_par, 0, 0, NULL,
                                            population_weather, n = 3))
  basic <- phenoflex_population(population_weather, population_par, n = 3,
                                basic_output = TRUE)
  expect_equal(basic$bloomindex, pop$bloomindex)
  expect_null(basic[["z"]])
})

test_that("forcing output and supplied populations retain their dimensions", {
  pop <- list(yc_pop = c(30, 40), zc_pop = c(150, 190))
  out <- phenoflex_population(population_weather, population_par, n = 2,
                              population = pop, jday_cut = c(20, 80))
  expect_equal(dim(out$exp), c(2L, 2L))
  expect_equal(out$yc_pop, pop$yc_pop)
  full <- phenoflex_population(population_weather, population_par, n = 2,
                               population = pop, jday_cut = c(20, 80), stop_at_zc = FALSE)
  expect_identical(out, full)
  expect_error(phenoflex_population(population_weather, population_par,
                                    jday_cut = 200), "cutting day")
  expect_error(phenoflex_population(population_weather, population_par,
                                    population = pop), "population")
  expect_error(phenoflex_population(population_weather[1, ], population_par), "two valid")
})
