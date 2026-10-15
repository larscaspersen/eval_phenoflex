#' Plot seasonal chill, transition and heat processes
#'
#' Follow a model through one hourly weather season or compare its processes
#' across a list of seasons. Curves use the same native calculations as
#' [predict_phenology()], continued through the full supplied weather period.
#'
#' @param model A single [pheno_model()], a fitted single model, or a
#'   [population_pheno_model()] (its mean parameters are shown). Select a single
#'   member of a collection with [as_pheno_models()] before plotting.
#' @param weather One complete hourly weather data frame, or a non-empty list
#'   of such seasons, as accepted by [predict_phenology()]. Each season must
#'   span at most two calendar years.
#' @param parameters Complete named parameter vector, defaulting to the model's
#'   stored values. Dynamic Model parameters may use either representation.
#' @param ribbon Logical. With a season list, use a grey ribbon spanning the
#'   minimum and maximum available values at each hour instead of individual
#'   seasonal curves. The black median line is retained. Default FALSE.
#'
#' @details
#' Chill accumulation uses the coupling state (the last native chill column).
#' The transition panel shows heat effectiveness for PhenoFlex and parallel
#' models. Sequential heat remains enabled after the first chill-threshold
#' crossing, even if chill subsequently decreases. For partial overlap, the
#' transition is the heat requirement, using additional chill since initial
#' chilling completion and held fixed after overlap ends. It is undefined
#' before minimum chilling. Heat-only models show just the heat panel.
#'
#' The horizontal axis shows calendar dates with abbreviated month labels.
#' A single season retains its actual dates. Season lists align by native
#' fractional day of year, including its year-boundary convention, rather than
#' by calendar month and day. Month labels and bloom-date annotations use a
#' common non-leap reference calendar. Thus day 60 combines February 29 in a
#' leap year with March 1 in a non-leap year; after February, labels can differ
#' by one day from a leap season's actual calendar dates.
#' Grey curves show individual seasons and the black curve is the pointwise
#' median of available seasons at each hourly position; missing periods are
#' not extrapolated. With `ribbon = TRUE`, a grey minimum-to-maximum envelope
#' replaces the individual curves. A single season is drawn in black.
#'
#' The dot and vertical dashed line mark the median predicted bloom date.
#' Horizontal whiskers span the earliest to latest successful prediction;
#' they describe variation among seasons, not a confidence interval. The dot
#' sits on the median heat curve at that date. Seasons without bloom remain
#' in process medians, are excluded from date summaries, and are counted in
#' the caption. With no successful predictions, bloom markers are omitted.
#'
#' @return An editable ggplot with vertically stacked process panels. Its
#'   data contains season, season_id, date, day, process and value columns.
#'   The date column supplies the plotted calendar positions; day retains
#'   the native fractional day-of-year convention of [return_JDay()].
#'   Attributes `predictions`, `median_processes` and `bloom_summary` contain
#'   per-season bloom dates, pointwise medians, and the plotted date summary.
#'   The median_processes data also includes lower and upper range bounds.
#'   Plotting does not modify model parameters or weather.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   model <- pheno_model("sequential", chill = chill_hours(),
#'                        parameters = c(yc = 24, zc = 20))
#'   weather <- data.frame(Temp = rep(c(6, 18), each = 3 * 24),
#'                         Year = 2001, JDay = rep(60:65, each = 24))
#'   plot_processes(model, weather)
#'   warmer <- weather
#'   warmer$Temp <- warmer$Temp + 1
#'   plot_processes(model, list(original = weather, warmer = warmer))
#'   plot_processes(model, list(original = weather, warmer = warmer), ribbon = TRUE)
#' }
#' @md
#' @export
plot_processes <- function(model, weather, parameters = model_parameters(model), ribbon = FALSE) {
  .require_plot_package("ggplot2")
  if (!is.logical(ribbon) || length(ribbon) != 1L || is.na(ribbon))
    stop("ribbon must be TRUE or FALSE.", call. = FALSE)
  model <- .single_plot_model(model, parameters)
  seasons <- if (is.data.frame(weather)) list(weather) else weather
  if (!.is_pheno_seasonlist(seasons))
    stop("weather must be an hourly data frame or a non-empty seasonlist of data frames.", call. = FALSE)
  season_names <- names(seasons)
  if (is.null(season_names)) season_names <- as.character(seq_along(seasons))
  unnamed <- is.na(season_names) | !nzchar(season_names)
  season_names[unnamed] <- as.character(which(unnamed))
  reference_year <- if (length(seasons) > 1L) 2001L else NULL
  panels <- .process_panel_names(model)
  trajectories <- lapply(seq_along(seasons), function(i) {
    .season_processes(model, seasons[[i]], i, season_names[i], panels, reference_year)
  })
  predictions <- data.frame(season = season_names,
    day = vapply(trajectories, function(x) attr(x, "bloom_day"), numeric(1)),
    date = as.Date(vapply(trajectories, function(x) as.numeric(attr(x, "bloom_date")),
      numeric(1)), origin = "1970-01-01"))
  data <- do.call(rbind, trajectories)
  attr(data, "bloom_day") <- NULL
  attr(data, "bloom_date") <- NULL
  rownames(data) <- NULL
  median_input <- data.frame(season_hour = round(data$day * 24), process = data$process,
    value = data$value)
  medians <- stats::aggregate(value ~ season_hour + process, data = median_input, na.action = stats::na.pass,
    FUN = function(x) {
      x <- x[!is.na(x)]
      if (!length(x)) return(c(value = NA_real_, lower = NA_real_, upper = NA_real_))
      c(value = stats::median(x), lower = min(x), upper = max(x))
    })
  medians <- cbind(medians[c("season_hour", "process")], as.data.frame(medians$value))
  medians$day <- medians$season_hour / 24
  medians$season_hour <- NULL
  median_year <- if (is.null(reference_year)) max(seasons[[1]]$Year) else reference_year
  medians$date <- as.Date(sprintf("%04d-01-01", median_year)) + medians$day - 13 / 24
  medians <- medians[order(medians$process, medians$day), ]
  rownames(medians) <- NULL
  successful <- predictions$day[!is.na(predictions$day)]
  heat <- medians[medians$process == panels[["heat"]], ]
  bloom <- data.frame(day = numeric(), bloom_min = numeric(), bloom_max = numeric(),
    date = as.Date(character()), date_min = as.Date(character()), date_max = as.Date(character()),
    value = numeric(), process = factor(levels = unname(panels)), label = character())
  multiple <- length(seasons) > 1L
  if (length(successful)) {
    bloom_day <- stats::median(successful)
    dates <- predictions$date[!is.na(predictions$date)]
    bloom_date <- as.Date(stats::median(as.numeric(dates)), origin = "1970-01-01")
    bloom <- data.frame(day = bloom_day, bloom_min = min(successful), bloom_max = max(successful),
      date = bloom_date, date_min = min(dates), date_max = max(dates),
      value = stats::approx(as.numeric(heat$date), heat$value, xout = as.numeric(bloom_date))$y,
      process = factor(panels[["heat"]], levels = unname(panels)),
      label = paste0(if (multiple) "Median bloom: " else "Bloom: ",
        format(bloom_date, "%d %b")))
  }
  caption <- if (multiple)
    paste0(if (ribbon) "Grey ribbon: range of available seasons" else "Grey: individual seasons",
      "; black: median of available seasons.\n",
      "Dot: median bloom date; whiskers: earliest to latest prediction.") else "Dot: predicted bloom date."
  if (anyNA(predictions$day)) caption <- paste0(caption, "\nBloom not reached in ",
    sum(is.na(predictions$day)), " of ", length(seasons), " seasons.")
  result <- ggplot2::ggplot(data, ggplot2::aes(x = date, y = value))
  if (multiple) {
    if (ribbon) result <- result + ggplot2::geom_ribbon(data = medians,
      ggplot2::aes(x = date, ymin = lower, ymax = upper), inherit.aes = FALSE,
      fill = "grey70", alpha = 0.5, colour = NA, na.rm = TRUE) else
      result <- result + ggplot2::geom_line(ggplot2::aes(group = season_id),
        colour = "grey70", linewidth = 0.4, alpha = 0.7, na.rm = TRUE)
  }
  result <- result + ggplot2::geom_line(data = medians, colour = "black", linewidth = 0.9,
    na.rm = TRUE)
  if (nrow(bloom)) {
    result <- result + ggplot2::geom_vline(xintercept = bloom$date,
      colour = "grey45", linetype = "dashed", linewidth = 0.4)
    if (length(successful) > 1L) result <- result + ggplot2::geom_errorbar(data = bloom,
      ggplot2::aes(y = value, xmin = date_min, xmax = date_max), orientation = "y",
      inherit.aes = FALSE, width = 0.03 * max(diff(range(heat$value)), 1), linewidth = 0.7)
    result <- result + ggplot2::geom_point(data = bloom, size = 2.8) +
      ggplot2::geom_text(data = bloom, ggplot2::aes(label = label),
        nudge_y = -0.07 * max(diff(range(heat$value)), 1),
        nudge_x =  0.01 * max(diff(range(heat$day)), 1),
        hjust = 0, size = 3.4)
  }
  result <- result + ggplot2::facet_wrap(ggplot2::vars(process), ncol = 1, scales = "free_y") +
    ggplot2::scale_x_date(breaks = function(limits) {
      months <- seq(as.Date(format(limits[1], "%Y-%m-01")),
        as.Date(format(limits[2], "%Y-%m-01")), by = "month")
      months <- months[months >= limits[1] & months <= limits[2]]
      if (length(months)) months else as.Date(mean(as.numeric(limits)), origin = "1970-01-01")
    }, date_labels = "%b") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.13))) +
    ggplot2::labs(title = paste(.plot_structure_name(model), "seasonal processes"),
      x = "Date", y = NULL, caption = caption) +
    .pheno_plot_theme()
  attr(result, "predictions") <- predictions
  attr(result, "median_processes") <- medians
  attr(result, "bloom_summary") <- bloom
  result
}

.process_panel_names <- function(model) {
  heat <- paste0("Heat accumulation (", .plot_heat_units(model), ")")
  if (model$chill$name == "none") return(c(heat = heat))
  c(chill = paste0("Chill accumulation (", .plot_chill_units(model), ")"),
    transition = if (model$structure == "partial_overlap")
      paste0("Transition: heat requirement (", .plot_heat_units(model), ")") else
      "Transition: fraction of potential heat that is effective",
    heat = heat)
}

.season_processes <- function(model, weather, id, name, panels, reference_year = NULL) {
  native <- predict_phenology(model, weather, stopatzc = FALSE, basic_output = FALSE)
  years <- unique(weather$Year)
  if (length(years) > 2L)
    stop("Each weather season must span at most two calendar years.", call. = FALSE)
  day <- weather$JDay
  if (length(years) == 2L) {
    previous <- weather$Year == min(years)
    year <- min(years)
    leap <- year %% 4 == 0 & (year %% 100 != 0 | year %% 400 == 0)
    day[previous] <- day[previous] - 365L - as.integer(leap)
  }
  day <- day + (rep(1:24, length.out = nrow(weather)) - 12L) / 24
  dates <- as.Date(sprintf("%04d-01-01", weather$Year)) + weather$JDay - 1L
  dates <- dates + rep(0:23, length.out = nrow(weather)) / 24
  if (!is.null(reference_year))
    dates <- as.Date(sprintf("%04d-01-01", reference_year)) + day - 13 / 24
  values <- list(heat = native$z)
  if (model$chill$name != "none") {
    chill <- native$chill[, ncol(native$chill)]
    transition <- if (model$structure == "partial_overlap") {
      trajectory <- .partial_overlap_plot_trajectory(model, weather, native)
      requirement <- .partial_overlap_requirement(model$parameters, trajectory$chill)
      requirement[trajectory$accumulated_chill < model$parameters[["yc"]]] <- NA_real_
      requirement
    } else if (model$structure == "sequential") {
      cummax(as.numeric(chill >= model$parameters[["yc"]]))
    } else .chill_heat_effectiveness(model, chill)
    values <- list(chill = chill, transition = transition, heat = native$z)
  }
  data <- data.frame(season = name, season_id = id, day = rep(day, length(values)),
    date = rep(dates, length(values)),
    process = factor(rep(unname(panels), each = length(day)), levels = unname(panels)),
    value = unlist(values, use.names = FALSE))
  attr(data, "bloom_day") <- return_JDay(native$bloomindex, weather$JDay, weather$Year)
  attr(data, "bloom_date") <- if (native$bloomindex > 0) dates[native$bloomindex] else as.Date(NA)
  data
}

utils::globalVariables(c("day", "date", "value", "season_id", "process", "bloom_min", "bloom_max",
  "date_min", "date_max", "season_hour", "lower", "upper"))
