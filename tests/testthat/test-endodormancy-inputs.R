test_that("endodormancy evaluation requires six numeric parameters before conversion", {
  for (parameters in list(NULL, numeric(), 1:5, 1:7, rep("1", 6),
                          rep(TRUE, 6), as.list(1:6), matrix(1:6, nrow = 2))) {
    expect_error(eval_function_endodormancy_meigo(parameters),
                 "x must be a numeric vector of length 6", fixed = TRUE)
  }
})

test_that("endodormancy evaluation accepts numeric vectors with six parameters", {
  weather <- data.frame(Temp = rep(8, 480))
  seasons <- list(cultivar_a = list(weather, weather),
                  cultivar_b = list(weather, weather, weather))
  parameters <- c(theta_star = 279, theta_c = 286.1, tau = 47.7,
                  pie_c = 28, Tf = 4, slope = 1.6)
  result <- eval_function_endodormancy_meigo(parameters, seasons)
  expect_named(result, c("F", "g"))
  expect_equal(result$F, 0)
  expect_true(all(is.finite(result$g)))
  expect_no_error(eval_function_endodormancy_meigo(as.integer(parameters), seasons))
})

test_that("endodormancy evaluation requires an outer list before parameter conversion", {
  for (seasons in list(NULL, 1:2, "weather", data.frame(Temp = 8))) {
    expect_error(eval_function_endodormancy_meigo(rep(0, 6), seasons),
                 "SeasonList must be a list of lists", fixed = TRUE)
  }
})

test_that("endodormancy evaluation requires every cultivar entry to be a season list", {
  weather <- data.frame(Temp = 8)
  for (bad_entry in list(NULL, 1:2, "weather", weather)) {
    for (cultivar in 1:2) {
      seasons <- list(list(weather), list(weather))
      seasons[cultivar] <- list(bad_entry)
      expect_error(eval_function_endodormancy_meigo(rep(0, 6), seasons),
                   paste0("SeasonList[[", cultivar, "]] must be a list of seasons"),
                   fixed = TRUE)
    }
  }
})
