prediction_list_example <- function() {
  first <- pheno_model("sequential", chill_dynamic("kinetic"),
                       parameters = c(yc = 0.1, zc = 5))
  second <- pheno_model("parallel", chill_dynamic("kinetic"), heat_gdh("scaled"),
                        parameters = c(yc = 0.2, zc = 150, Tb = 2))
  cold <- data.frame(Temp = rep(8, 480), Year = 2008,
                     JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  warm <- cold
  warm$Temp <- 10
  list(models = list(first = first, second = second), seasons = list(cold = cold, warm = warm))
}

test_that("single models predict season lists without changing single-season results", {
  x <- prediction_list_example()
  for (model in x$models) for (stop in c(TRUE, FALSE)) {
    expected <- vapply(x$seasons, function(w) predict_phenology(model, w, stopatzc = stop), numeric(1))
    expect_identical(predict_phenology(model, x$seasons, stopatzc = stop), expected)
    expect_true(all(is.finite(expected)))
    details <- lapply(x$seasons, function(w) predict_phenology(model, w,
                       stopatzc = stop, basic_output = FALSE))
    expect_identical(predict_phenology(model, x$seasons, stopatzc = stop, basic_output = FALSE), details)
  }
  model <- x$models[[1]]
  override <- model$parameters
  override["zc"] <- 1e12
  result <- predict_phenology(model, x$seasons, rev(override))
  expect_identical(result, c(cold = NA_real_, warm = NA_real_))
  expect_identical(predict_phenology(model, x$seasons[1]),
                   c(cold = predict_phenology(model, x$seasons[[1]])))
})

test_that("ordinary independent model lists support common or nested weather", {
  x <- prediction_list_example()
  before <- x$models
  for (basic in c(TRUE, FALSE)) {
    expected <- lapply(x$models, function(m) predict_phenology(m, x$seasons[[1]], basic_output = basic))
    expect_identical(predict_phenology(x$models, x$seasons[[1]], basic_output = basic), expected)
    expected <- lapply(x$models, function(m) predict_phenology(m, x$seasons, basic_output = basic))
    expect_identical(predict_phenology(x$models, x$seasons, basic_output = basic), expected)
    nested <- list(cultivar_a = x$seasons[2:1], cultivar_b = x$seasons[1])
    expected <- list(cultivar_a = predict_phenology(x$models[[1]], nested[[1]], basic_output = basic),
                      cultivar_b = predict_phenology(x$models[[2]], nested[[2]], basic_output = basic))
    expect_identical(predict_phenology(x$models, nested, basic_output = basic), expected)
    expect_identical(predict_phenology(x$models, unname(nested), basic_output = basic),
                     setNames(expected, names(x$models)))
  }
  expect_identical(x$models, before)
  one <- x$models[1]
  expect_identical(predict_phenology(one, x$seasons),
                   list(first = predict_phenology(one[[1]], x$seasons)))
})

test_that("ordinary model lists accept per-model or common complete parameter overrides", {
  x <- prediction_list_example()
  before <- x$models
  p <- lapply(x$models, function(m) {
    values <- m$parameters
    values["zc"] <- values[["zc"]] * 2
    rev(values)
  })
  expected <- lapply(seq_along(x$models), function(i) predict_phenology(x$models[[i]], x$seasons, p[[i]]))
  names(expected) <- names(x$models)
  expect_identical(predict_phenology(x$models, x$seasons, p), expected)
  expect_identical(x$models, before)
  models <- list(a = x$models[[1]], b = x$models[[1]])
  expect_identical(predict_phenology(models, x$seasons[[1]], p[[1]]),
                   lapply(models, function(m) predict_phenology(m, x$seasons[[1]], p[[1]])))
  expect_error(predict_phenology(x$models, x$seasons, p[1]), "one complete named vector per model")
  p[[2]] <- p[[2]][-1]
  expect_error(predict_phenology(x$models, x$seasons, p), "schema names")
})

test_that("calibration collections and combined models accept common single-season weather", {
  x <- prediction_list_example()
  stages <- stage_pheno_models(x$models[[1]], c(early = 5, late = 10))
  combined <- combined_pheno_model(x$models[[1]], 2)
  for (model in list(stages, combined)) {
    p <- model_parameters(model)
    p["zc2"] <- 15
    simple <- as_pheno_models(model, p)
    for (basic in c(TRUE, FALSE)) {
      expect_identical(predict_phenology(model, x$seasons[[1]], p, basic_output = basic),
                       predict_phenology(simple, x$seasons[[1]], basic_output = basic))
      expect_identical(predict_phenology(model, x$seasons, p, basic_output = basic),
                       predict_phenology(simple, x$seasons, basic_output = basic))
      nested <- list(one = x$seasons[1], two = x$seasons)
      expect_identical(predict_phenology(model, nested, p, basic_output = basic),
                       predict_phenology(simple, nested, basic_output = basic))
    }
  }
})

test_that("list predictions reject empty, mixed or mismatched input shapes", {
  x <- prediction_list_example()
  for (weather in list(list(), list(NULL), list(x$seasons[[1]], NULL),
                       list(x$seasons), list(x$seasons, list())))
    expect_error(predict_phenology(x$models, weather), "seasonlist")
  expect_error(predict_phenology(x$models[[1]], list()), "seasonlist")
  expect_error(predict_phenology(x$models[[1]], list(x$seasons)), "seasonlist")
  expect_error(predict_phenology(list(), x$seasons), "non-empty list")
  expect_error(predict_phenology(list(x$models[[1]], NULL), x$seasons), "pheno_model objects")
  bad <- x$seasons
  bad[[2]]$Temp[1] <- NA_real_
  expect_error(predict_phenology(x$models, bad), "finite numeric")
  expect_error(predict_phenology(x$models, x$seasons, stopatzc = NA), "Output flags")
  expect_error(predict_phenology(x$models[[1]], x$seasons, basic_output = 1), "Output flags")
})

test_that("population season lists retain each population prediction result", {
  x <- prediction_list_example()
  model <- population_pheno_model("sequential", chill_dynamic("kinetic"), n = 2)
  model$model <- x$models[[1]]
  expected <- lapply(x$seasons, function(w) predict_phenology(model, w))
  expect_identical(predict_phenology(model, x$seasons), expected)
})
