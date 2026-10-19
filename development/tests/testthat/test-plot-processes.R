process_weather <- function(temp = 6, year = 2001, days = 50:55) {
  data.frame(Temp = temp, Year = year, JDay = rep(days, each = 24),
    Hour = rep(0:23, length(days)))
}

test_that("process curves and bloom dates match native predictions for each structure", {
  skip_if_not_installed("ggplot2")
  weather <- process_weather(25)
  before <- weather
  for (structure in c("sequential", "phenoflex", "parallel", "partial_overlap")) {
    parameters <- if (structure == "partial_overlap")
      c(yc = 2, b1 = 1, b2 = 2, b3 = 0.2, ol = 3) else c(yc = 2, zc = 4)
    model <- pheno_model(structure, chill = chill_hours(),
      parameters = c(parameters, CH_upper = 30))
    native <- predict_phenology(model, weather, stopatzc = FALSE, basic_output = FALSE)
    p <- plot_processes(model, weather)
    expect_s3_class(p, "ggplot")
    panels <- levels(p$data$process)
    expect_equal(p$data$value[p$data$process == panels[1]], native$chill[, ncol(native$chill)])
    expect_equal(p$data$value[p$data$process == panels[3]], native$z)
    expect_equal(attr(p, "predictions")$day, predict_phenology(model, weather))
    bloom <- attr(p, "bloom_summary")
    expect_equal(bloom$day, predict_phenology(model, weather))
    expect_equal(bloom$value, native$z[native$bloomindex])
    transition <- p$data$value[p$data$process == panels[2]]
    if (structure == "partial_overlap") {
      expect_equal(transition[1:6], c(NA, NA, 3, 1 + 2 * exp(-0.2),
        1 + 2 * exp(-0.4), 1 + 2 * exp(-0.4)))
      expect_equal(tail(transition, 1), transition[5])
    } else {
      # At 25 C the unscaled heat increment is exactly one per hour.
      applied <- if (structure == "phenoflex") head(transition, -1) else tail(transition, -1)
      expect_equal(diff(native$z), applied, tolerance = 1e-12)
    }
    expect_silent(ggplot2::ggplot_build(p))
    expect_true("Feb" %in% ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels())
  }
  expect_identical(weather, before)
})

test_that("season lists have pointwise medians and bloom range whiskers", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model("sequential", chill = chill_hours(), parameters = c(yc = 24, zc = 1))
  seasons <- list(cool = process_weather(6), warm = process_weather(7, year = 2002))
  before <- seasons
  p <- plot_processes(model, seasons)
  medians <- attr(p, "median_processes")
  for (panel in levels(p$data$process)) {
    cool <- p$data$value[p$data$process == panel & p$data$season == "cool"]
    warm <- p$data$value[p$data$process == panel & p$data$season == "warm"]
    expect_equal(medians$value[medians$process == panel], (cool + warm) / 2)
  }
  predicted <- unname(predict_phenology(model, seasons))
  bloom <- attr(p, "bloom_summary")
  expect_equal(attr(p, "predictions")$day, predicted)
  expect_equal(bloom$day, median(predicted))
  expect_equal(c(bloom$bloom_min, bloom$bloom_max), range(predicted))
  expect_match(p$labels$caption, "earliest to latest prediction", fixed = TRUE)
  expect_identical(p$layers[[1]]$aes_params$colour, "grey70")
  expect_identical(p$layers[[2]]$aes_params$colour, "black")
  expect_silent(ggplot2::ggplot_build(p))
  expect_identical(seasons, before)

  # No extrapolation into the part of a season missing from the other one.
  shorter <- seasons
  shorter$warm <- shorter$warm[shorter$warm$JDay >= 51, ]
  p <- plot_processes(model, shorter)
  medians <- attr(p, "median_processes")
  early <- medians[medians$day < 50.5, ]
  original <- p$data[p$data$day < 50.5, ]
  expect_equal(early$value, original$value)
})

test_that("ribbons show the full seasonal range and preserve medians and bloom summaries", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model("sequential", chill = chill_hours(), parameters = c(yc = 24, zc = 1))
  seasons <- list(cold = process_weather(0), cool = process_weather(6), warm = process_weather(7))
  lines <- plot_processes(model, seasons)
  p <- plot_processes(model, seasons, ribbon = TRUE)
  expect_identical(attr(p, "median_processes"), attr(lines, "median_processes"))
  expect_identical(attr(p, "bloom_summary"), attr(lines, "bloom_summary"))
  expect_identical(p$data, lines$data)
  ranges <- attr(p, "median_processes")
  for (panel in levels(p$data$process)) {
    values <- lapply(names(seasons), function(season) {
      p$data$value[p$data$process == panel & p$data$season == season]
    })
    expect_equal(ranges$lower[ranges$process == panel], do.call(pmin, values))
    expect_equal(ranges$upper[ranges$process == panel], do.call(pmax, values))
  }
  expect_s3_class(p$layers[[1]]$geom, "GeomRibbon")
  expect_s3_class(p$layers[[2]]$geom, "GeomLine")
  expect_identical(p$layers[[2]]$aes_params$colour, "black")
  expect_match(p$labels$caption, "Grey ribbon: range", fixed = TRUE)
  expect_silent(ggplot2::ggplot_build(p))
  single <- plot_processes(model, seasons[[1]], ribbon = TRUE)
  expect_identical(single$data, plot_processes(model, seasons[[1]])$data)
  expect_false(any(vapply(single$layers, function(x) inherits(x$geom, "GeomRibbon"), logical(1))))
  overlap <- pheno_model("partial_overlap", chill = chill_hours())
  empty <- plot_processes(overlap, list(process_weather(25), process_weather(26)), ribbon = TRUE)
  ranges <- attr(empty, "median_processes")
  expect_true(all(is.na(ranges$lower[grepl("Transition", ranges$process)])))
  expect_true(all(is.na(ranges$upper[grepl("Transition", ranges$process)])))
  expect_silent(ggplot2::ggplot_build(empty))
  for (invalid in list(NA, 1, "TRUE", NULL, c(TRUE, FALSE)))
    expect_error(plot_processes(model, seasons, ribbon = invalid), "ribbon must be TRUE or FALSE")
})

test_that("dates use native noon and year-boundary conventions", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model(chill = chill_none(), parameters = c(zc = 1))
  dates <- as.Date("2020-12-30") + 0:3
  weather <- data.frame(Temp = 25, Year = rep(as.integer(format(dates, "%Y")), each = 24),
    JDay = rep(as.integer(format(dates, "%j")), each = 24))
  p <- plot_processes(model, weather)
  expect_equal(p$data$day[c(12, 36, 60, 84)], c(-1, 0, 1, 2))
  expect_true(all(diff(p$data$day) > 0))
  expect_s3_class(p$data$date, "Date")
  expect_equal(p$data$date[1], as.Date("2020-12-30"))
  expect_equal(p$data$date[25], as.Date("2020-12-31"))
  expect_equal(as.numeric(diff(p$data$date)), rep(1 / 24, nrow(weather) - 1), tolerance = 1e-10)
  expect_equal(nrow(attr(p, "median_processes")), nrow(weather))
  # expect_identical(p$scales$get_scales("x")$get_labels(as.Date(c("2020-12-01", "2021-01-01"))),
  #   c("Dec", "Jan"))
  expect_equal(attr(p, "predictions")$day, return_JDay(2, weather$JDay, weather$Year))
  expect_length(levels(p$data$process), 1)
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("day-of-year alignment retains all seasons and their range around leap day", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model(chill = chill_none(), parameters = c(zc = 1))
  seasons <- list(leap = process_weather(6, 2020, 58:62),
    ordinary = process_weather(18, 2021, 58:62))
  p <- plot_processes(model, seasons, ribbon = TRUE)
  medians <- attr(p, "median_processes")
  expect_equal(nrow(medians), 5 * 24)
  expect_true(all(table(as.numeric(p$data$date)) == 2L))
  for (day in 59:61) {
    rows <- p$data[abs(p$data$day - day) < 1e-8, ]
    summary <- medians[abs(medians$day - day) < 1e-8, ]
    expect_equal(nrow(rows), 2L)
    expect_equal(summary$value, median(rows$value))
    expect_equal(c(summary$lower, summary$upper), range(rows$value))
    expect_gt(summary$upper, summary$lower)
  }
  expect_true(all(diff(medians$value) >= 0))
  day60 <- p$data[abs(p$data$day - 60) < 1e-8, ]
  expect_equal(format(day60$date, "%m-%d"), c("03-01", "03-01"))
  expect_equal(attr(p, "predictions")$day, unname(predict_phenology(model, seasons)))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("missing bloom is retained in curves and excluded from date summaries", {
  skip_if_not_installed("ggplot2")
  model <- pheno_model("sequential", chill = chill_hours(), parameters = c(yc = 24, zc = 1))
  seasons <- list(bloom = process_weather(6), no_bloom = process_weather(0))
  p <- plot_processes(model, seasons)
  expect_true(is.na(attr(p, "predictions")$day[2]))
  expect_equal(attr(p, "bloom_summary")$day, predict_phenology(model, seasons[[1]]))
  expect_match(p$labels$caption, "Bloom not reached in 1 of 2 seasons", fixed = TRUE)
  expect_true(all(attr(p, "median_processes")$value >= 0))
  none <- plot_processes(model, seasons[2])
  expect_equal(nrow(attr(none, "bloom_summary")), 0)
  expect_silent(ggplot2::ggplot_build(none))
  overlap <- pheno_model("partial_overlap", chill = chill_hours())
  p <- plot_processes(overlap, process_weather(25))
  transition <- grepl("Transition", p$data$process)
  expect_true(all(is.na(p$data$value[transition])))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("process plots retain sequential onset and accept parameter overrides and fits", {
  skip_if_not_installed("ggplot2")
  custom <- chill_custom(function(weather, parameters) {
    matrix(rep(c(0, 2, 0.5, 0.25), length.out = nrow(weather)), ncol = 1)
  }, label = "Changing chill")
  model <- pheno_model("sequential", chill = custom, parameters = c(yc = 1, zc = 2))
  p <- plot_processes(model, process_weather(25))
  transition <- p$data$value[grepl("Transition", p$data$process)]
  expect_equal(transition, c(0, rep(1, length(transition) - 1)))
  model <- pheno_model()
  initial <- model_parameters(model, "characteristic")
  fit <- fit_phenology(model, function(p, m) sum(p), numeric(), numeric(), parameters = initial)
  expect_equal(plot_processes(fit, process_weather())$data,
    plot_processes(model, process_weather(), parameters = initial)$data)
  expect_error(plot_processes(model, list()), "non-empty seasonlist")
  expect_error(plot_processes(model, list(1)), "non-empty seasonlist")
  expect_error(plot_processes(model, data.frame(Temp = 6)), "complete hourly days")
})
