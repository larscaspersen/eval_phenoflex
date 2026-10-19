test_that("temperature responses use the selected native chill and heat kernels", {
  skip_if_not_installed("ggplot2")
  temperatures <- c(-5, 0, 4, 8, 20, 25, 36, 40)
  for (chill in list(chill_dynamic(), chill_hours(), chill_utah(), chill_positive_utah(), chill_none())) {
    for (scaling in c("scaled", "unscaled")) {
      model <- pheno_model(chill = chill, heat = heat_gdh(scaling))
      before <- model
      plot <- plot_temperature_response(model, temp_values = temperatures, chill_duration = 48)
      expect_s3_class(plot, "ggplot")
      data <- plot$data
      expected_chill <- vapply(temperatures, function(t) {
        p <- model$parameters
        times <- 0:48
        values <- switch(chill$name,
          dynamic = calculate_chill_dynamic(rep(t, 49), times, E0 = p[["E0"]], E1 = p[["E1"]],
            A0 = p[["A0"]], A1 = p[["A1"]], Tf = p[["Tf"]], slope = p[["slope"]])[, "y"],
          chill_hours = calculate_chill_hours(rep(t, 49), times, CH_lower = p[["CH_lower"]],
            CH_upper = p[["CH_upper"]])[, "y"],
          utah = calculate_chill_utah(rep(t, 49), times, unname(p[paste0("utah_w", 1:7)]))[, "y_raw"],
          positive_utah = calculate_chill_positive_utah(rep(t, 49), times,
            unname(p[paste0("pos_utah_w", 1:7)]))[, "y"], none = rep(0, 49))
        tail(values, 1) / 48
      }, numeric(1))
      expected_heat <- vapply(temperatures, function(t) {
        p <- model$parameters
        fn <- if (scaling == "scaled") calculate_heat_gdh else calculate_heat_gdh_unscaled
        fn(rep(t, 2), 0:1, Tb = p[["Tb"]], Tu = p[["Tu"]], Tc = p[["Tc"]])
      }, numeric(1))
      expect_equal(data$response[seq_along(temperatures)], expected_chill)
      expect_equal(tail(data$response, length(temperatures)), expected_heat)
      expect_identical(model, before)
      expect_silent(ggplot2::ggplot_build(plot))
    }
  }
})

test_that("transition plots agree with native heat effectiveness", {
  skip_if_not_installed("ggplot2")
  weather <- data.frame(Temp = 8, Year = 2001, JDay = rep(1:2, each = 24))
  for (structure in c("sequential", "parallel", "phenoflex")) {
    model <- pheno_model(structure)
    plot <- plot_model_structure(model)
    expect_s3_class(plot, "ggplot")
    data <- plot$data
    expected <- vapply(data$chill, function(y) {
      p <- model$parameters
      fn <- switch(structure, sequential = apply_sequential_structure,
        parallel = apply_parallel_structure, phenoflex = apply_phenoflex_structure)
      arguments <- list(chill = matrix(rep(y, 2), ncol = 1), heat = 1, yc = p[["yc"]],
        zc = 1e6, stopatzc = FALSE, basic_output = FALSE)
      if (structure == "parallel") arguments$kmin <- p[["kmin"]]
      if (structure == "phenoflex") arguments$s1 <- p[["s1"]]
      do.call(fn, arguments)$z[2]
    }, numeric(1))
    expect_equal(data$effectiveness, expected, tolerance = 1e-14)
    expect_silent(ggplot2::ggplot_build(plot))
    with_weather <- plot_model_structure(model, weather = weather)
    expect_equal(ggplot2::ggplot_build(with_weather)$data, ggplot2::ggplot_build(plot)$data)
    expect_identical(with_weather$labels$caption, plot$labels$caption)
  }
})

test_that("partial-overlap curve, defining points and overlap guide match the model", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model("partial_overlap", parameters = c(yc = 2, b1 = 10, b2 = 20, b3 = 0.5, ol = 0.75))
  plot <- plot_model_structure(model)
  expect_equal(plot$data$requirement,
    ifelse(plot$data$chill >= 2, 10 + 20 * exp(-0.5 * (plot$data$chill - 2)), NA_real_))
  points <- plot$layers[[3]]$data
  # Native classification immediately below and above the drawn requirement.
  for (additional in c(0, 2, 4)) {
    requirement <- 10 + 20 * exp(-0.5 * additional)
    for (delta in c(-0.001, 0.001)) {
      predicted <- apply_partial_overlap_structure(matrix(c(2, 2 + additional), ncol = 1),
        heat = requirement + delta, yc = 2, b1 = 10, b2 = 20, b3 = 0.5, ol = 0.75)
      expect_identical(predicted$bloomindex, if (delta > 0) 2L else 0L)
    }
  }
  # The point layer, after horizontal guides and the requirement curve.
  expect_equal(points$chill, 2)
  expect_equal(points$requirement, 30)
  expect_equal(plot$layers[[1]]$data$requirement, c(10, 7.5))
  expect_silent(ggplot2::ggplot_build(plot))
  for (values in list(c(ol = 0), c(b2 = 0), c(b3 = 0))) {
    p <- model$parameters
    p[names(values)] <- values
    fixed <- plot_model_structure(model, parameters = p)
    expected <- if ("b2" %in% names(values)) 10 else 30
    expect_equal(fixed$data$requirement, ifelse(fixed$data$chill >= 2, expected, NA_real_))
    expect_equal(nrow(fixed$layers[[3]]$data), 1L)
  }
})

test_that("weather overlays use the actual partial-overlap cutoff and exclude overshoot", {
  skip_if_not_installed("ggplot2")
  custom <- chill_custom(function(weather, parameters) {
    matrix(seq.int(0, nrow(weather) - 1L), ncol = 1)
  })
  model <- pheno_model("partial_overlap", chill = custom,
    parameters = c(yc = 2.5, b1 = 5, b2 = 10, b3 = 0.5, ol = 0.75))
  weather <- data.frame(Temp = 25, Year = 2001, JDay = rep(1:2, each = 24))
  before <- weather
  trajectory <- .partial_overlap_plot_trajectory(model, weather)
  # Chill first crosses at row 4 (3 units); the crossing interval contributes heat.
  expect_identical(trajectory$chill[4], 0)
  expect_identical(trajectory$chill[5], 1)
  expect_identical(trajectory$chill[6], 2)
  expect_identical(trajectory$chill[7], 3)
  expect_true(all(trajectory$chill[7:48] == 3))
  expect_identical(trajectory$event[7], "Overlap cutoff")
  expect_equal(trajectory$heat, predict_phenology(model, weather, stopatzc = FALSE, basic_output = FALSE)$z)
  plot <- plot_model_structure(model, weather = weather)
  expect_silent(ggplot2::ggplot_build(plot))
  season_path <- plot$layers[[7]]$data
  expect_equal(nrow(season_path), attr(trajectory, "bloomindex"))
  expect_true(any(grepl("Event|/ event", season_path$event)))
  expect_equal(season_path$chill[1:7], c(0, 1, 2, 3, 4, 5, 6))
  expect_identical(weather, before)
  p <- model$parameters
  p["ol"] <- 0
  expect_true(all(.partial_overlap_plot_trajectory(as_pheno_models(model, p)[[1]], weather)$chill == 0))
  p["yc"] <- 100
  unfulfilled <- as_pheno_models(model, p)[[1]]
  expect_true(all(.partial_overlap_plot_trajectory(unfulfilled, weather)$heat == 0))
  empty <- plot_model_structure(unfulfilled, weather = weather)
  expect_silent(ggplot2::ggplot_build(empty))
})

test_that("histogram and gradient overlays handle boundaries, sparse months and arbitrary grids", {
  skip_if_not_installed("ggplot2")
  weather <- data.frame(Temp = c(0, 1, 1, 1, 2, 50), Month = 1)
  for (type in c("histogram", "gradient")) {
    plot <- plot_temperature_response(pheno_model(chill = chill_hours()),
      temp_values = c(0, 0.3, 1, 1.7, 2), chill_duration = 24, hourtemps = weather,
      bins = 2, weather_freq_plot = type)
    expect_true(all(is.finite(plot$data$response)))
    expect_silent(ggplot2::ggplot_build(plot))
    expect_match(plot$labels$caption, "supplied weather time series", fixed = TRUE)
    labels <- unlist(plot$facet$params$labeller(
      data.frame(component = c("Chill response", "Heat response"))))
    expect_identical(unname(labels), c("Chill response\nMonths: Nov, Dec, Jan, Feb",
      "Heat response\nMonths: Jan, Feb, Mar, Apr, May"))
    if (type == "gradient") expect_equal(plot$data$frequency[1:5], c(1, 1, 1, 0.25, 0.25))
  }
  no_month <- data.frame(Temp = c(0, 1), Year = 2001, JDay = c(1, 32))
  plot <- plot_temperature_response(pheno_model(), temp_values = c(0, 1, 2), chill_duration = 24,
    hourtemps = no_month, chill_months = 12, heat_months = 1, weather_freq_plot = "gradient")
  expect_true(all(plot$data$frequency[1:3] == 0))
  expect_silent(ggplot2::ggplot_build(plot))
  labels <- unlist(plot$facet$params$labeller(
    data.frame(component = c("Chill response", "Heat response"))))
  expect_identical(unname(labels), c("Chill response\nMonths: Dec", "Heat response\nMonths: Jan"))
  for (hourtemps in list(NULL, weather)) {
    plain <- plot_temperature_response(pheno_model(), temp_values = c(0, 1, 2),
      chill_duration = 24, hourtemps = hourtemps, weather_freq_plot = "none")
    expect_null(plain$labels$caption)
    labels <- unlist(plain$facet$params$labeller(
      data.frame(component = c("Chill response", "Heat response"))))
    expect_identical(unname(labels), c("Chill response", "Heat response"))
  }
})

test_that("overview and component plots accept stored and characteristic fitted parameters", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model(parameters = c(theta_star = 6))
  initial <- model_parameters(model, "characteristic")
  fit <- fit_phenology(model, function(p, m) sum(p), numeric(), numeric(), parameters = initial)
  expect_identical(plot_model_structure(fit)$data, plot_model_structure(model)$data)
  expect_identical(plot_temperature_response(model, temp_values = c(0, 8, 25), chill_duration = 24)$data,
    plot_temperature_response(fit, temp_values = c(0, 8, 25), chill_duration = 24)$data)
  expect_identical(plot_temperature_response(model, parameters = initial,
    temp_values = c(0, 8, 25), chill_duration = 24)$data,
    plot_temperature_response(model, temp_values = c(0, 8, 25), chill_duration = 24)$data)
  expect_s3_class(plot_model_structure(population_pheno_model(n = 2)), "ggplot")
  skip_if_not_installed("patchwork")
  overview <- plot(model, temp_values = c(0, 8, 25), chill_duration = 24)
  expect_s3_class(overview, "patchwork")
  expect_silent(patchwork::patchworkGrob(overview))
})

test_that("heat-only plots show the interval start and malformed inputs fail clearly", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model(chill = chill_none(), parameters = c(heat_start_index = 3.5))
  plot <- plot_model_structure(model)
  expect_equal(plot$data$effectiveness, as.numeric(plot$data$index >= 3.5))
  expect_silent(ggplot2::ggplot_build(plot))
  for (values in list(c(0, NA_real_), c(0, Inf), c(-273, 0), c(8, 0), 8, matrix(0:1)))
    expect_error(plot_temperature_response(pheno_model(), temp_values = values), "temp_values")
  for (value in list(0, -1, 1.5, Inf, NA_real_))
    expect_error(plot_temperature_response(pheno_model(), chill_duration = value), "chill_duration")
  for (values in list(c(1, 0), c(-1, 5), c(0, Inf), c(0, 1, 2)))
    expect_error(plot_model_structure(pheno_model(), chill_range = values), "chill_range")
  expect_error(plot_model_structure(combined_pheno_model(pheno_model(), 2)), "select a collection")
  expect_error(plot_temperature_response(pheno_model(), hourtemps = data.frame(Temp = NA_real_)), "Temp")
  expect_error(plot_temperature_response(pheno_model(), hourtemps = data.frame(Temp = 8)), "Month")
  expect_error(plot_temperature_response(pheno_model(), hourtemps = data.frame(Temp = 8, Month = 13)), "Month")
  expect_error(plot_temperature_response(pheno_model(), hourtemps = data.frame(Temp = 8, Month = 1),
    chill_months = 0), "calendar months")
  expect_error(plot_temperature_response(pheno_model(), hourtemps = data.frame(Temp = 8, Year = 2001, JDay = 366)),
    "invalid Year or JDay")
  expect_equal(max(plot_model_structure(pheno_model("partial_overlap"))$data$chill), 120)
  expect_equal(max(plot_model_structure(pheno_model("partial_overlap"), chill_range = c(0, 150))$data$chill), 150)
  testthat::local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  expect_error(plot_model_structure(pheno_model()), "Install the 'ggplot2'")
})
