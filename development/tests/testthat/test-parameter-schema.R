test_that("model and schema printing use decimal notation without changing values", {
  old_options <- options(scipen = -9)
  on.exit(options(old_options), add = TRUE)
  model <- pheno_model(parameters = c(s1 = 0.00000001))
  parameters <- model$parameters
  schema <- parameter_schema(model)
  defaults <- schema$default

  model_output <- capture.output(result <- withVisible(print(model)))
  expect_false(any(grepl("[0-9][eE][+-][0-9]", model_output)))
  expect_false(any(grepl("[.][0-9]{3}", model_output)))
  expect_true(any(grepl("2567000000000000000", model_output, fixed = TRUE)))
  expect_identical(result$value, model)
  expect_false(result$visible)
  expect_identical(model$parameters, parameters)
  expect_equal(getOption("scipen"), -9)

  schema_output <- capture.output(result <- withVisible(print(schema)))
  expect_false(any(grepl("[0-9][eE][+-][0-9]", schema_output)))
  expect_true(any(grepl("2567000000000000000", schema_output, fixed = TRUE)))
  expect_identical(result$value, schema)
  expect_false(result$visible)
  expect_identical(schema$default, defaults)
  expect_type(schema$default, "double")
  expect_equal(getOption("scipen"), -9)

  expect_false(any(grepl("[.][0-9]{3}", schema_output)))
})

test_that("parameter displays round to two places and omit trailing zeros", {
  model <- pheno_model(parameters = c(yc = 1.57, s1 = 1.5701, Tb = -0.001))
  parameters <- model$parameters
  output <- capture.output(print(model))
  expect_true(any(grepl("1.57", output, fixed = TRUE)))
  expect_false(any(grepl("1.570", output, fixed = TRUE)))
  expect_false(any(grepl("[.][0-9]{3}", output)))
  expect_false(any(grepl("(^|[[:space:]])-0([[:space:]]|$)", output)))
  expect_identical(model$parameters, parameters)

  schema <- parameter_schema(model)
  schema$default[1:5] <- c(1.57, 1.5701, 1.5, 40, -0.001)
  original <- schema$default
  display <- format(schema)
  expect_identical(as.character(trimws(display$default[1:5])),
                   c("1.57", "1.57", "1.5", "40", "0"))
  expect_identical(schema$default, original)
  expect_type(schema$default, "double")
})

test_that("decimal schema printing is inherited by indexed and population schemas", {
  model <- pheno_model()
  models <- list(combined_pheno_model(model, 2),
    pheno_model_list(list(model, model)),
    stage_pheno_models(model, c(200, 300)),
    population_pheno_model(n = 2))
  for (model in models) {
    schema <- parameter_schema(model)
    expect_s3_class(schema, "pheno_parameter_schema")
    expect_s3_class(schema, "data.frame")
    expect_type(schema$default, "double")
    expect_false(any(grepl("[0-9][eE][+-][0-9]", capture.output(print(schema)))))
  }
})

test_that("chill requirement units follow the selected submodel", {
  chills <- list(dynamic = chill_dynamic(), hours = chill_hours(),
    utah = chill_utah(), positive = chill_positive_utah())
  expected_units <- c(dynamic = "Chill Portions", hours = "Chilling Hours",
    utah = "Utah Chill Units", positive = "Utah Chill Units")
  for (name in names(chills)) {
    for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
      schema <- parameter_schema(pheno_model(structure, chill = chills[[name]]))
      units <- setNames(schema$unit, schema$name)
      expect_identical(unname(units["yc"]), unname(expected_units[name]))
      expect_true(all(!is.na(schema$description) & nzchar(schema$description)))
      expect_true(all(!is.na(schema$unit) & nzchar(schema$unit)))
      for (parameter in intersect(c("s1", "b3"), schema$name))
        expect_identical(unname(units[parameter]), paste0("1 / (", expected_units[name], ")"))
    }
  }
})

test_that("heat requirement units follow GDH scaling, including heat-only models", {
  for (scaling in c("scaled", "unscaled")) {
    expected <- if (scaling == "scaled") "degree C h" else "h"
    for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
      for (chill in list(chill_dynamic(), chill_hours(), chill_none())) {
        schema <- parameter_schema(pheno_model(structure, chill, heat_gdh(scaling)))
        requirements <- schema$name %in% c("zc", "b1", "b2")
        expect_true(any(requirements))
        expect_true(all(schema$unit[requirements] == expected))
      }
    }
  }
})

test_that("descriptions distinguish structural roles and characteristic units", {
  descriptions <- lapply(c("sequential", "parallel", "phenoflex"), function(structure) {
    schema <- parameter_schema(pheno_model(structure))
    schema$description[schema$name == "yc"]
  })
  expect_match(descriptions[[1]], "Chill requirement", fixed = TRUE)

  model <- pheno_model()
  schema <- parameter_schema(model, "characteristic")
  units <- setNames(schema$unit, schema$name)
  expect_identical(unname(units[c("theta_star", "theta_c", "tau", "pie_c")]),
    c("K", "K", "h", "h"))
  expect_false("representation" %in% names(schema))
  expect_identical(names(model$parameters)[4:7], c("E0", "E1", "A0", "A1"))
})

test_that("indexed schemas preserve descriptions and units", {
  for (scaling in c("scaled", "unscaled")) {
    base <- pheno_model(chill = chill_hours(), heat = heat_gdh(scaling))
    reference <- parameter_schema(base)
    models <- list(
      combined_pheno_model(base, 2),
      pheno_model_list(list(base, base)),
      pheno_model_list(list(base, base), shared = character()),
      stage_pheno_models(base, c(200, 300)))
    for (model in models) {
      schema <- parameter_schema(model)
      index <- match(schema$base_name, reference$name)
      expect_identical(schema$description, reference$description[index])
      expect_identical(schema$unit, reference$unit[index])
      expect_false("representation" %in% names(schema))
    }
    population <- population_pheno_model(chill = chill_hours(), heat = heat_gdh(scaling), n = 2)
    expect_identical(parameter_schema(population), reference)
  }
})

test_that("custom parameter units remain unknown even for built-in parameter names", {
  fun <- function(weather, parameters) matrix(0, nrow(weather), 1)
  for (defaults in list(numeric(), c(E0 = 2, slope = 1, width = 5))) {
    schema <- parameter_schema(pheno_model(chill = chill_custom(fun, defaults)))
    expect_identical(schema$unit[schema$name == "yc"], "custom chill units")
    custom <- schema$component == "chill"
    expect_true(all(is.na(schema$unit[custom])))
    expect_true(all(startsWith(schema$description[custom], "Custom chill parameter:")))
    expect_true(all(!is.na(schema$description)))
  }
})
