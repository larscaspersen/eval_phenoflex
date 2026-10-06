modular_population_weather <- data.frame(Temp = rep(8, 480), Year = 2008,
  JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))

test_that("population specifications reuse the single-model parameter schema", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    model <- population_pheno_model(structure, n = 3)
    expect_s3_class(model, "population_pheno_model")
    expect_true(validate_model_spec(model))
    expect_identical(parameter_schema(model), parameter_schema(model$model))
    expect_identical(default_parameters(model), default_parameters(model$model))
    population <- sample_population_parameters(model)
    expect_equal(nrow(population), 3)
    expect_identical(colnames(population), names(model$sd))
    expect_equal(population[1, ], default_parameters(model)[colnames(population)])
    expect_equal(population[1, ], population[3, ])
  }
  expect_error(population_pheno_model(n = 0), "n")
  expect_error(population_pheno_model(n = 1.5), "n")
  expect_error(population_pheno_model(sd = c(E0 = 1)), "structure parameters")
  expect_error(population_pheno_model(sd = c(yc = -1)), "non-negative")
  expect_error(population_pheno_model(sd = c(yc = 1, yc = 2)), "structure parameters")
  expect_error(population_pheno_model(skew = 1), "named vector")
  expect_error(population_pheno_model(seed = -1), "seed")
})

test_that("sampling is reproducible, independent and preserves caller RNG state", {
  model <- population_pheno_model(n = 100, sd = c(yc = 1, zc = 2), seed = 17)
  set.seed(9)
  before <- .Random.seed
  a <- sample_population_parameters(model)
  expect_identical(.Random.seed, before)
  expect_identical(a, sample_population_parameters(model))
  expect_false(identical(a[, "yc"] - 40, (a[, "zc"] - default_parameters(model)[["zc"]]) / 2))
  model$seed <- 18
  expect_false(identical(a, sample_population_parameters(model)))
  model$seed <- NULL
  sample_population_parameters(model)
  expect_false(identical(.Random.seed, before))
  # Zero-dispersion single buds are deterministic, as are larger populations.
  single <- population_pheno_model(n = 1)
  expect_equal(unname(sample_population_parameters(single)[1, "yc"]), 40)
  # Deliberately invalid draws are rejected, never truncated or resampled.
  bad <- population_pheno_model(n = 100, sd = c(yc = 1e6), seed = 17)
  expect_error(sample_population_parameters(bad), "Invalid population bud")
})

test_that("skew sampling shares the existing sampler and keeps actual moments", {
  skip_if_not_installed("sn")
  model <- population_pheno_model(n = 10000, sd = c(yc = 1, zc = 2),
    distribution = "normal_skewed", skew = c(yc = 2, zc = -2), seed = 17)
  a <- sample_population_parameters(model)
  expect_identical(a, sample_population_parameters(model))
  p <- default_parameters(model)
  expect_equal(mean(a[, "yc"]), p[["yc"]], tolerance = .05)
  expect_equal(stats::sd(a[, "yc"]), 1, tolerance = .05)
  expect_equal(mean(a[, "zc"]), p[["zc"]], tolerance = .05)
  expect_equal(stats::sd(a[, "zc"]), 2, tolerance = .05)
})

test_that("constant populations match single buds for all structures and scalings", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    for (scaling in c("scaled", "unscaled")) {
      model <- population_pheno_model(structure, chill_dynamic("kinetic"), heat_gdh(scaling), n = 2)
      p <- default_parameters(model)
      p["yc"] <- .1
      if (structure == "partial_overlap") p[c("b1", "b2")] <- c(20, 30) else p["zc"] <- 30
      if (scaling == "unscaled") {
        names_heat <- intersect(c("zc", "b1", "b2"), names(p))
        p[names_heat] <- p[names_heat] / (p[["Tu"]] - p[["Tb"]])
      }
      for (stop in c(TRUE, FALSE)) {
        reference <- predict_phenology(model$model, modular_population_weather, p,
                                       stopatzc = stop, basic_output = FALSE)
        actual <- predict_population_phenology(model, modular_population_weather, p,
                                               stopatzc = stop, basic_output = FALSE)
        expect_gt(reference$bloomindex, 0)
        expect_equal(actual$bloomindex, rep(reference$bloomindex, 2))
        expect_equal(actual$z[, 1], reference$z)
        expect_equal(actual$z[, 2], reference$z)
        expect_equal(actual$bloom_jday, rep(return_JDay(reference$bloomindex,
          modular_population_weather$JDay, modular_population_weather$Year), 2))
        expect_identical(actual, predict_phenology(model, modular_population_weather, p,
                                                  stopatzc = stop, basic_output = FALSE))
      }
    }
  }
})

test_that("explicit populations preserve bud identity and allow correlated traits", {
  model <- population_pheno_model("phenoflex", chill_dynamic("kinetic"), n = 3)
  p <- default_parameters(model)
  population <- rbind(c(yc = .1, zc = 1, s1 = .5),
                      c(yc = .2, zc = 2, s1 = .5),
                      c(yc = 1e12, zc = 1e12, s1 = .5))
  set.seed(9); before <- .Random.seed
  out <- predict_population_phenology(model, modular_population_weather, p,
                                      population = population, basic_output = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(out$population, population)
  expect_equal(dim(out$z), c(480L, 3L))
  expect_equal(dim(out$chill), c(480L, 3L))
  expect_equal(out$bloomindex[3], 0)
  expect_true(is.na(out$bloom_jday[3]))
  for (i in 1:3) {
    bud <- p; bud[colnames(population)] <- population[i, ]
    reference <- predict_phenology(model$model, modular_population_weather, bud, basic_output = FALSE)
    expect_equal(out$bloomindex[i], reference$bloomindex)
    expect_equal(out$z[, i], reference$z)
  }
  reversed <- population[, rev(colnames(population)), drop = FALSE]
  expect_identical(predict_population_phenology(model, modular_population_weather, rev(p),
    population = reversed, basic_output = FALSE), out)
  expect_error(predict_population_phenology(model, modular_population_weather,
    population = population[-1, ]), "n rows")
  invalid <- population; invalid[1, "yc"] <- 0
  expect_error(predict_population_phenology(model, modular_population_weather,
    population = invalid), "Invalid population bud 1")
})

test_that("shared chill and potential heat are calculated once per population", {
  chill_fn <- calculate_chill_dynamic
  heat_fn <- calculate_heat_gdh_unscaled
  counts <- c(chill = 0, heat = 0)
  local_mocked_bindings(
    calculate_chill_dynamic = function(...) { counts["chill"] <<- counts["chill"] + 1; chill_fn(...) },
    calculate_heat_gdh_unscaled = function(...) { counts["heat"] <<- counts["heat"] + 1; heat_fn(...) },
    .package = "evalpheno")
  model <- population_pheno_model(n = 4)
  predict_population_phenology(model, modular_population_weather)
  expect_equal(unname(counts), c(1, 1))
})

test_that("forcing retains field heat and freezes chill with correct elapsed hours", {
  for (structure in c("sequential", "parallel", "phenoflex", "partial_overlap")) {
    model <- population_pheno_model(structure, chill_dynamic("kinetic"), heat_gdh("scaled"), n = 2)
    p <- default_parameters(model)
    p["yc"] <- .1
    if (structure == "partial_overlap") {
      p[c("b1", "b2", "b3", "ol")] <- c(5000, 3000, .1, 100)
    } else p["zc"] <- 5000
    cuts <- c(240, 480)
    out <- predict_population_phenology(model, modular_population_weather, p,
      stopatzc = TRUE, basic_output = FALSE, cut_indices = cuts,
      forcing_temperature = p[["Tu"]], max_hours_forcing = 500)
    full <- predict_population_phenology(model, modular_population_weather, p,
      stopatzc = FALSE, basic_output = FALSE, cut_indices = cuts,
      forcing_temperature = p[["Tu"]], max_hours_forcing = 500)
    expect_identical(out$forcing, full$forcing)
    expect_equal(dim(out$forcing$hours_to_bloom), c(2L, 2L))
    y <- full$chill[, "y"]
    heat_increment <- p[["Tu"]] - p[["Tb"]] # scaled GDH at its optimum
    for (j in seq_along(cuts)) {
      cut <- cuts[j]
      expect_gt(y[cut], p[["yc"]])
      effectiveness <- switch(structure,
        sequential = 1, partial_overlap = 1,
        parallel = 1,
        phenoflex = stats::plogis(p[["s1"]] * p[["yc"]] * (y[cut] - p[["yc"]]) / y[cut]))
      increment <- heat_increment * effectiveness
      requirement <- if (structure == "partial_overlap") {
        first <- which(y >= p[["yc"]])[1]
        p[["b1"]] + p[["b2"]] * exp(-p[["b3"]] * (y[cut] - y[first]))
      } else p[["zc"]]
      start <- full$z[cut, 1]
      expected <- max(0, ceiling((requirement - start) / increment))
      expect_lte(expected, 500)
      expect_equal(out$forcing$hours_to_bloom[j, ], rep(expected, 2))
      expect_equal(out$forcing$z_at_cut[j, ], rep(start, 2))
      expect_equal(out$forcing$z[[j]][, 1], start + (0:500) * increment, tolerance = 1e-10)
    }
    # A cut before any chill cannot complete a sequential/overlap/PhenoFlex bud.
    early <- predict_population_phenology(model, modular_population_weather, p,
      cut_indices = 1, forcing_temperature = p[["Tu"]], max_hours_forcing = 2)
    expect_true(all(is.na(early$forcing$hours_to_bloom)))
    p[intersect(c("zc", "b1", "b2"), names(p))] <- .01
    bloomed <- predict_population_phenology(model, modular_population_weather, p,
      cut_indices = 480, max_hours_forcing = 1)
    expect_equal(as.vector(bloomed$forcing$hours_to_bloom), c(0, 0))
  }
})

test_that("forcing validates cuts, duration and output flags", {
  model <- population_pheno_model(n = 1)
  for (cuts in list(0, 481, c(1, 1), NA_real_, 1.5))
    expect_error(predict_population_phenology(model, modular_population_weather,
      cut_indices = cuts), "cut_indices")
  expect_error(predict_population_phenology(model, modular_population_weather,
    cut_indices = 1, forcing_temperature = -273), "forcing_temperature")
  expect_error(predict_population_phenology(model, modular_population_weather,
    cut_indices = 1, max_hours_forcing = 0), "max_hours_forcing")
  expect_error(predict_population_phenology(model, modular_population_weather,
    stopatzc = NA), "Output flags")
})
