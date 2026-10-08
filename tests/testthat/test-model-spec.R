test_that("sequential specifications expose exact, validated parameter names", {
  model <- pheno_model("sequential", heat = heat_gdh("scaled"))
  expect_s3_class(model, "pheno_model")
  expect_true(validate_model_spec(model))
  p <- default_parameters(model)
  expect_identical(names(p), c("yc", "zc", "theta_star", "theta_c", "tau", "pie_c",
                               "Tf", "slope", "Tb", "Tu", "Tc"))
  expect_identical(parameter_schema(model)$component,
                   c(rep("structure", 2), rep("chill", 6), rep("heat", 3)))
  expect_identical(names(p), parameter_schema(model)$name)
  expect_true(validate_parameters(model, p))
  expect_true(validate_parameters(model, rev(p)))
  
  #parallel
  parallel <- pheno_model("parallel")
  parallel_p <- default_parameters(parallel)
  expect_length(parallel_p, 12)
  expect_equal(parallel_p[["kmin"]], 0.1)
  expect_true(validate_parameters(parallel, parallel_p))
  ##
  
  
  expect_error(pheno_model(chill = list(name = "utah")), "chill_dynamic")
  expect_error(pheno_model(heat = list(name = "gdd")), "heat_gdh")
  expect_error(validate_parameters(model, unname(p)), "schema names")
  expect_error(validate_parameters(model, p[-1]), "schema names")
  expect_error(validate_parameters(model, c(p, extra = 1)), "schema names")
  duplicate <- p; names(duplicate)[1] <- names(duplicate)[2]
  expect_error(validate_parameters(model, duplicate), "duplicates")
  for (name in c("yc", "zc", "slope")) {
    bad <- p; bad[name] <- 0
    expect_error(validate_parameters(model, bad), "positive")
  }
  bad <- p; bad["Tu"] <- bad["Tb"]
  expect_error(validate_parameters(model, bad), "Tb < Tu < Tc")
  bad <- p; bad["theta_c"] <- 300
  expect_error(validate_parameters(model, bad), "theta_c")
  bad <- p; bad[1] <- Inf
  expect_error(validate_parameters(model, bad), "finite")
  kinetic <- pheno_model(chill = chill_dynamic("kinetic"))
  p <- default_parameters(kinetic)
  expect_true(validate_parameters(kinetic, p))
  p["E1"] <- p["E0"]
  expect_error(validate_parameters(kinetic, p), "E0 < E1")
})

test_that("named sequential predictions reproduce every frozen sequential case", {
  path <- test_path("fixtures", "weather-regression")
  fixtures <- readRDS(file.path(path, "inputs.rds"))
  baseline <- read.csv(file.path(path, "baseline.csv"))
  for (representation in c("characteristic", "kinetic")) {
    model <- pheno_model("sequential", chill = chill_dynamic(representation),
                         heat = heat_gdh("scaled"))
    p <- default_parameters(model)
    p[c("yc", "zc")] <- c(20, 100)
    if (representation == "kinetic")
      p[c("E0", "E1", "A0", "A1")] <- fixtures$parameters$kinetic[5:8]
    for (station in names(fixtures$seasons)) for (year in names(fixtures$seasons[[station]])) {
      weather <- fixtures$seasons[[station]][[year]]
      expected <- baseline$value[baseline$case_id == "sequential" &
                                   baseline$station == station & baseline$season == year]
      expect_length(expected, 1)
      result <- predict_phenology(model, weather, p)
      expect_equal(result,
                   expected, tolerance = 1e-10)
      expect_identical(predict_phenology(model, weather, rev(p)), predict_phenology(model, weather, p))
      no_bloom <- p; no_bloom[c("yc", "zc")] <- 1e12
      expect_true(is.na(predict_phenology(model, weather, no_bloom)))
    }
  }
})

test_that("sequential adapter rejects invalid hourly weather", {
  model <- pheno_model("sequential", heat = heat_gdh("scaled"))
  weather <- data.frame(Temp = rep(10, 48), Year = 2008,
                        JDay = rep(59:60, each = 24), Hour = rep(0:23, 2))
  expect_type(predict_phenology(model, weather), "double")
  expect_error(predict_phenology(model, weather[-1, ]), "complete hourly")
  bad <- weather; bad$Temp[1] <- NA
  expect_error(predict_phenology(model, bad), "finite numeric")
  bad <- weather; bad$JDay[25:48] <- 61
  expect_error(predict_phenology(model, bad), "date gaps")
  bad <- weather; bad$Hour[1] <- 1
  expect_error(predict_phenology(model, bad), "0 to 23")
  bad <- weather; bad$JDay <- 367
  expect_error(predict_phenology(model, bad), "invalid calendar")
})



test_that("model parameters default, validate and drive predictions", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    for (representation in c("characteristic", "kinetic")) {
      for (scaling in c("scaled", "unscaled")) {
        model <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling))
        expect_identical(model$parameters, default_parameters(model))
        p <- model$parameters
        p["yc"] <- 45
        custom <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling),
                              parameters = rev(p))
        expect_identical(custom$parameters, rev(p))
        expect_identical(default_parameters(custom), default_parameters(model))
        partial <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling),
                                parameters = p[-1])
        expect_identical(partial$parameters, default_parameters(model))
        partial <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling),
                                parameters = c(Tf = 5, yc = 50))
        expected <- default_parameters(model)
        expected[c("Tf", "yc")] <- c(5, 50)
        expect_identical(partial$parameters, expected)
      }
    }
  }
  model <- pheno_model("sequential", chill_dynamic("kinetic"))
  p <- model$parameters
  p[c("yc", "zc")] <- c(1, 1)
  model <- pheno_model("sequential", chill_dynamic("kinetic"), parameters = p)
  weather <- data.frame(Temp = rep(10, 24 * 60), Year = 2008,
                        JDay = rep(1:60, each = 24))
  expect_true(is.finite(predict_phenology(model, weather)))
  expect_identical(predict_phenology(model, weather),
                   predict_phenology(model, weather, parameters = p))
  override <- p
  override[c("yc", "zc")] <- 1e12
  expect_true(is.na(predict_phenology(model, weather, parameters = override)))
  expect_identical(model$parameters, p)
  bad <- p; bad["yc"] <- -1
  expect_error(pheno_model("sequential", chill_dynamic("kinetic"), parameters = bad), "positive")
})

test_that("partial constructor parameters fill defaults and validate the completed model", {
  model <- pheno_model("parallel", parameters = c(yc = 50, zc = 180))
  expected <- default_parameters(model)
  expected[c("yc", "zc")] <- c(50, 180)
  expect_identical(model$parameters, expected)
  expect_identical(pheno_model(parameters = numeric())$parameters,
                   default_parameters(pheno_model()))
  expect_error(pheno_model(parameters = c(50, 180)), "named numeric")
  expect_error(pheno_model(parameters = c(yc = 50, unknown = 1)), "schema names")
  expect_error(pheno_model(parameters = c(yc = 50, yc = 60)), "duplicates")
  expect_error(pheno_model(parameters = setNames(50, NA_character_)), "schema names")
  expect_error(pheno_model(parameters = setNames(50, "")), "schema names")
  expect_error(pheno_model(parameters = c(yc = "50")), "named numeric")
  expect_error(pheno_model(parameters = matrix(50, dimnames = list("yc", NULL))),
                "named numeric")
  for (value in c(NA_real_, NaN, Inf, -Inf))
    expect_error(pheno_model(parameters = c(yc = value)), "finite")
  expect_error(pheno_model(parameters = c(yc = 0)), "positive")
  expect_error(pheno_model(parameters = c(Tb = 30)), "Tb < Tu < Tc")
  expect_error(pheno_model(chill = chill_dynamic("kinetic"), parameters = c(E0 = 20000)),
                "E0 < E1")
  expect_error(pheno_model("parallel", parameters = c(kmin = 2)), "between 0 and 1")
  # Explicit prediction overrides still require a complete vector.
  expect_error(validate_parameters(model, c(yc = 50, zc = 180)), "schema names")
})
