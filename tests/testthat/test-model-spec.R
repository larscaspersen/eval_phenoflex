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
      expect_equal(return_JDay(result$bloomindex, weather$JDay, weather$Year),
                   expected, tolerance = 1e-10)
      expect_identical(predict_phenology(model, weather, rev(p)), predict_phenology(model, weather, p))
      no_bloom <- p; no_bloom[c("yc", "zc")] <- 1e12
      expect_equal(predict_phenology(model, weather, no_bloom)$bloomindex, 0)
    }
  }
})

test_that("sequential adapter rejects invalid hourly weather", {
  model <- pheno_model("sequential", heat = heat_gdh("scaled"))
  weather <- data.frame(Temp = rep(10, 48), Year = 2008,
                        JDay = rep(59:60, each = 24), Hour = rep(0:23, 2))
  expect_named(predict_phenology(model, weather), "bloomindex")
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

