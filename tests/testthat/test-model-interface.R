test_that("each structure exposes its own schema under both heat scalings", {
  structure_names <- list(
    sequential = c("yc", "zc"),
    parallel = c("yc", "zc", "kmin"),
    partial_overlap = c("yc", "b1", "b2", "b3", "ol"),
    phenoflex = c("yc", "zc", "s1")
  )
  for (structure in names(structure_names)) {
    for (representation in c("kinetic", "characteristic")) {
      for (scaling in c("scaled", "unscaled")) {
        model <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling))
        schema <- parameter_schema(model)
        chill_names <- if (representation == "kinetic")
          c("E0", "E1", "A0", "A1", "Tf", "slope") else
          c("theta_star", "theta_c", "tau", "pie_c", "Tf", "slope")
        expect_identical(schema$name,
                         c(structure_names[[structure]], chill_names, "Tb", "Tu", "Tc"))
        expect_identical(schema$component,
                         c(rep("structure", length(structure_names[[structure]])),
                           rep("chill", 6), rep("heat", 3)))
        expect_identical(schema$representation[schema$component == "chill"],
                         c(rep(representation, 4), "native", "native"))
        p <- default_parameters(model)
        expect_identical(names(p), schema$name)
        expect_equal(unname(p), schema$default)
        expect_true(validate_parameters(model, p))
        expect_true(validate_parameters(model, rev(p)))
      }
    }
    scaled <- default_parameters(pheno_model(structure, heat = heat_gdh("scaled")))
    unscaled <- default_parameters(pheno_model(structure, heat = heat_gdh("unscaled")))
    requirements <- intersect(c("zc", "b1", "b2"), names(scaled))
    factor <- scaled[["Tu"]] - scaled[["Tb"]]
    expect_equal(unscaled[requirements], scaled[requirements] / factor)
    other <- setdiff(names(scaled), requirements)
    expect_identical(unscaled[other], scaled[other])
  }
})

test_that("heat specifications reject malformed and unsupported selections", {
  for (scaling in c("scaled", "unscaled"))
    expect_true(validate_model_spec(pheno_model(heat = heat_gdh(scaling))))
  expect_error(heat_gdh("unknown"), "arg")
  expect_error(pheno_model(heat = list(name = "gdh")), "exactly name and scaling")
  expect_error(pheno_model(heat = list(name = "gdh", scaling = "scaled", extra = 1)),
               "exactly name and scaling")
  expect_error(pheno_model(heat = list(name = "gdd", scaling = "scaled")), "GDH heat model")
  for (bad in list(NULL, NA_character_, c("scaled", "unscaled"), "unknown"))
    expect_error(pheno_model(heat = list(name = "gdh", scaling = bad)), "GDH scaling")
  expect_error(pheno_model("unknown"), "Supported structures")
})

test_that("structure domains are checked without imposing calibration bounds", {
  positive <- list(sequential = c("yc", "zc"), parallel = c("yc", "zc"),
                   partial_overlap = c("yc", "b1"), phenoflex = c("yc", "zc", "s1"))
  for (structure in names(positive)) {
    model <- pheno_model(structure)
    p <- default_parameters(model)
    for (name in positive[[structure]]) for (value in c(0, -1)) {
      bad <- p; bad[name] <- value
      expect_error(validate_parameters(model, bad), "positive")
    }
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      bad <- p; bad[1] <- value
      expect_error(validate_parameters(model, bad), "finite")
    }
    expect_error(validate_parameters(model, matrix(p, nrow = 1)), "named numeric vector")
    bad <- p; names(bad)[1] <- NA_character_
    expect_error(validate_parameters(model, bad), "schema names")
    bad <- p; bad["Tb"] <- -273
    expect_error(validate_parameters(model, bad), "Tb must exceed")
  }
  model <- pheno_model("parallel")
  p <- default_parameters(model)
  for (value in c(0, 1)) {
    p["kmin"] <- value
    expect_true(validate_parameters(model, p))
  }
  for (value in c(-0.01, 1.01)) {
    p["kmin"] <- value
    expect_error(validate_parameters(model, p), "kmin")
  }
  model <- pheno_model("partial_overlap")
  p <- default_parameters(model)
  # b2 is an additive heat requirement, so it need not exceed b1.
  p[c("b2", "b3", "ol")] <- 0
  expect_true(validate_parameters(model, p))
  p["b2"] <- p[["b1"]] / 2
  p["ol"] <- 2
  expect_true(validate_parameters(model, p))
  for (name in c("b2", "b3", "ol")) {
    bad <- p; bad[name] <- -1
    expect_error(validate_parameters(model, bad), "non-negative")
  }
  p[c("b1", "ol")] <- c(.Machine$double.xmax, 2)
  expect_error(validate_parameters(model, p), "b1 \\* ol must be finite")
})

test_that("named predictions dispatch all modules, representations and output flags", {
  fixtures <- readRDS(test_path("fixtures", "weather-regression", "inputs.rds"))
  weather <- fixtures$seasons[[1]][[1]]
  for (representation in c("kinetic", "characteristic")) {
    for (scaling in c("scaled", "unscaled")) {
      for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
        model <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling))
        p <- default_parameters(model)
        p["yc"] <- 20
        if (structure == "partial_overlap") {
          p[c("b1", "b2", "b3", "ol")] <- c(150, 100, .1, .5)
        } else p["zc"] <- 100
        if (scaling == "unscaled") {
          requirements <- intersect(c("zc", "b1", "b2"), names(p))
          p[requirements] <- p[requirements] / (p[["Tu"]] - p[["Tb"]])
        }
        kinetic <- if (representation == "kinetic") p[c("E0", "E1", "A0", "A1")] else {
          converted <- characteristic_to_kinetic(
            c(rep(0, 4), unname(p[c("theta_star", "theta_c", "tau", "pie_c")]), rep(0, 4)),
            failure_return = "NA")
          setNames(converted[5:8], c("E0", "E1", "A0", "A1"))
        }
        times <- seq_along(weather$Temp)
        chill <- calculate_chill_dynamic(weather$Temp, times,
          E0 = kinetic[["E0"]], E1 = kinetic[["E1"]],
          A0 = kinetic[["A0"]], A1 = kinetic[["A1"]],
          Tf = p[["Tf"]], slope = p[["slope"]])
        heat_fn <- if (scaling == "scaled") calculate_heat_gdh else calculate_heat_gdh_unscaled
        heat <- heat_fn(weather$Temp, times, Tb = p[["Tb"]], Tu = p[["Tu"]], Tc = p[["Tc"]])
        structure_fn <- switch(structure,
          sequential = apply_sequential_structure, parallel = apply_parallel_structure,
          partial_overlap = apply_partial_overlap_structure, phenoflex = apply_phenoflex_structure)
        structure_args <- as.list(p[parameter_schema(model)$component == "structure"])
        for (stop in c(TRUE, FALSE)) {
          expected <- do.call(structure_fn, c(list(chill = chill, heat = heat),
            structure_args, list(stopatzc = stop, basic_output = FALSE)))
          actual <- predict_phenology(model, weather, p, stopatzc = stop, basic_output = FALSE)
          expect_equal(actual, expected, tolerance = 1e-10)
          expect_named(actual, c("bloomindex", "chill", "z"))
          expect_gt(actual$bloomindex, 0)
          expect_equal(dim(actual$chill), c(nrow(weather), 3L))
          expect_length(actual$z, nrow(weather))
          basic <- predict_phenology(model, weather, p, stopatzc = stop)
          expect_named(basic, "bloomindex")
          expect_equal(basic$bloomindex, actual$bloomindex)
          expect_identical(predict_phenology(model, weather, rev(p), stopatzc = stop), basic)
        }
        no_bloom <- p
        no_bloom[intersect(c("zc", "b1", "b2"), names(p))] <- 1e12
        result <- predict_phenology(model, weather, no_bloom)
        expect_equal(result$bloomindex, 0)
        expect_true(is.na(return_JDay(result$bloomindex, weather$JDay, weather$Year)))
      }
    }
  }
})

test_that("switching GDH scaling preserves bloom and rescales full heat trajectories", {
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    scaled <- pheno_model(structure, chill_dynamic("kinetic"), heat_gdh("scaled"))
    unscaled <- pheno_model(structure, chill_dynamic("kinetic"), heat_gdh("unscaled"))
    p <- default_parameters(scaled)
    # Non-default temperature thresholds ensure dispatch uses supplied parameters.
    p[c("Tb", "Tu", "Tc", "yc")] <- c(2, 24, 38, .1)
    if (structure == "partial_overlap") p[c("b1", "b2")] <- c(20, 30) else p["zc"] <- 30
    q <- p
    requirements <- intersect(c("zc", "b1", "b2"), names(p))
    factor <- p[["Tu"]] - p[["Tb"]]
    q[requirements] <- p[requirements] / factor
    for (stop in c(TRUE, FALSE)) {
      a <- predict_phenology(scaled, weather, p, stopatzc = stop, basic_output = FALSE)
      b <- predict_phenology(unscaled, weather, q, stopatzc = stop, basic_output = FALSE)
      expect_gt(a$bloomindex, 0)
      expect_equal(a$bloomindex, b$bloomindex)
      expect_equal(a$chill, b$chill)
      expect_equal(a$z, factor * b$z, tolerance = 1e-10)
    }
  }
})

test_that("named predictions retain legacy station predictions in their native scaling", {
  fixtures <- readRDS(test_path("fixtures", "weather-regression", "inputs.rds"))
  baseline <- read.csv(test_path("fixtures", "weather-regression", "baseline.csv"))
  for (structure in c("parallel", "partial_overlap", "phenoflex")) {
    scaling <- if (structure == "phenoflex") "unscaled" else "scaled"
    model <- pheno_model(structure, chill_dynamic("kinetic"), heat_gdh(scaling))
    p <- default_parameters(model)
    p[c("E0", "E1", "A0", "A1")] <- fixtures$parameters$kinetic[5:8]
    p["yc"] <- 20
    if (structure != "partial_overlap") p["zc"] <- 100 else
      p[c("b1", "b2", "b3", "ol")] <- c(150, 100, .1, .5)
    for (station in names(fixtures$seasons)) for (year in names(fixtures$seasons[[station]])) {
      weather <- fixtures$seasons[[station]][[year]]
      case <- if (structure == "phenoflex") "ordinary" else structure
      expected <- baseline$value[baseline$case_id == case &
                                  baseline$station == station & baseline$season == year]
      expect_length(expected, 1)
      actual <- predict_phenology(model, weather, p)$bloomindex
      expect_equal(return_JDay(actual, weather$JDay, weather$Year), expected, tolerance = 1e-10)
    }
  }
})
