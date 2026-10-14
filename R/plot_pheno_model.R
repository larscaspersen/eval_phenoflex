#' Plot a model's temperature responses and chill-to-heat transition
#'
#' Visualize the stored parameters of a modular phenology model. The overview
#' combines chill and heat temperature responses with the structure's transition.
#' Individual plotting functions return editable ggplot objects; the overview
#' returns a patchwork composition. These functions require the optional packages
#' ggplot2 and, for the overview, patchwork.
#'
#' Dynamic Model chill response is accumulated chill divided by the duration of
#' constant-temperature exposure, calculated with the native model kernel. It is
#' a finite-exposure average, not an instantaneous or equilibrium response. Utah
#' shows signed chill units; positive Utah shows non-negative daily chill units.
#' Custom chill shows the average change in its final coupling state on synthetic
#' hourly weather beginning on 1 January 2001. Custom callbacks requiring other
#' weather fields cannot be evaluated this way. Heat uses the selected GDH scaling.
#'
#' Sequential coupling shows the threshold that first enables heat accumulation;
#' heat remains enabled after that crossing. Parallel and PhenoFlex show the
#' fraction of potential heat that is effective at a given chill state, including
#' PhenoFlex's native sigmoid saturation limits.
#' Sequential plots shade zero heat effectiveness before the first yc crossing.
#' PhenoFlex plots shade effectiveness below 1%, ending at
#' yc / (1 - qlogis(0.01) / (s1 * yc)). This shading is a visual reference;
#' it does not change the transition curve or prediction calculations.
#'
#' Partial overlap shows the heat requirement against accumulated chill, with
#' additional chill measured relative to the first yc crossing. Its single
#' defining point marks minimum chilling yc and the initial heat requirement
#' b1+b2. Below yc, heat cannot trigger the event. Above yc the requirement is
#' b1 + b2 * exp(-b3 * additional chill). Lines identify the lower heat asymptote
#' b1 and the overlap cutoff at accumulated heat b1*ol. The cutoff
#' chill value depends on the weather and cannot be inferred from ol alone.
#' With ol = 0, additional chill cannot affect the requirement, which is b1+b2.
#' Optional weather overlays the native heat trajectory against the additional
#' chill used for coupling, held constant after overlap ends. Hourly overshoot at
#' the first yc crossing is excluded from additional chill, as in prediction.
#' The displayed weather trajectory ends at the first event, or at the season's
#' endpoint when no event is reached.
#'
#' @param model A single [pheno_model()], a fitted single model, or a
#'   [population_pheno_model()] (its mean parameters are shown). Select a single
#'   member from collections using [as_pheno_models()] before plotting.
#' @param parameters Complete named vector, defaulting to stored values.
#'   Dynamic Model overrides accept kinetic or characteristic names, without mixing.
#' @param temp_values Finite, strictly increasing Celsius temperatures, above -273.
#' @param chill_duration Positive integer number of hours of constant-temperature
#'   exposure used for the chill response. Default 1200 hours, matching the intent
#'   of the original PhenoFlex temperature-response plot.
#' @param hourtemps Optional observed temperatures with Temp and Month columns.
#'   Month may instead be derived from Year and JDay. Temperatures outside the
#'   displayed range are omitted from the frequency overlay.
#' @param chill_months,heat_months Calendar months used for each observed-temperature
#'   overlay. Histogram height and gradient intensity are relative to the most
#'   frequent displayed temperature bin in that panel; they are not chill or heat.
#' @param weather_freq_plot "histogram", "gradient" or "none".
#' @param bins Number of equal-width bins for observed temperatures.
#' @param chill_range Optional two increasing non-negative endpoints for the
#'   transition x-axis, in accumulated chill units. NULL chooses a range from
#'   parameters, capped at 120 portions for Dynamic Model chill. An explicit
#'   range can extend beyond that display cap.
#' @param weather Optional complete hourly season passed to [predict_phenology()].
#'   For partial overlap, adds its effective chill/heat trajectory and overlap/event
#'   points. For other chill structures, marks observed chill states along the x-axis.
#' @param x A pheno_model for the S3 plot method.
#' @param y Unused; retained for the plot generic.
#' @param ... Arguments passed from plot() to plot_pheno_model().
#' @return plot_temperature_response() and plot_model_structure() return ggplot
#'   objects. Their data are accessible through $data; partial-overlap markers
#'   are in the point layer's data. plot_pheno_model() and plot() return a patchwork.
#'   Plotting does not modify model parameters or weather.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   model <- pheno_model("partial_overlap")
#'   plot_temperature_response(model)
#'   plot_model_structure(model)
#'   if (requireNamespace("patchwork", quietly = TRUE)) plot(model)
#' }
#' @md
#' @export
plot_pheno_model <- function(model, parameters = model_parameters(model),
                             temp_values = seq(-5, 40, by = 0.25),
                             chill_duration = 1200, hourtemps = weather,
                             chill_months = c(11:12, 1:2), heat_months = 1:5,
                             weather_freq_plot = c("histogram", "gradient", "none"),
                             bins = 30, chill_range = NULL, weather = NULL) {
  .require_plot_package("patchwork")
  model <- .single_plot_model(model, parameters)
  response <- plot_temperature_response(model, temp_values = temp_values,
    chill_duration = chill_duration, hourtemps = hourtemps,
    chill_months = chill_months, heat_months = heat_months,
    weather_freq_plot = weather_freq_plot, bins = bins)
  transition <- plot_model_structure(model, chill_range = chill_range, weather = weather)
  title <- paste(.plot_structure_name(model), "model")
  patchwork::wrap_plots(response, transition, ncol = 1, heights = c(1, 1.15)) +
    patchwork::plot_annotation(title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 17)))
}

#' @rdname plot_pheno_model
#' @export
plot.pheno_model <- function(x, y = NULL, ...) plot_pheno_model(x, ...)

#' @rdname plot_pheno_model
#' @export
plot_temperature_response <- function(model, parameters = model_parameters(model),
                                      temp_values = seq(-5, 40, by = 0.25),
                                      chill_duration = 1200, hourtemps = NULL,
                                      chill_months = c(11:12, 1:2), heat_months = 1:5,
                                      weather_freq_plot = c("histogram", "gradient", "none"),
                                      bins = 30) {
  .require_plot_package("ggplot2")
  model <- .single_plot_model(model, parameters)
  .plot_numeric_grid(temp_values, "temp_values", lower = -273, strict_lower = TRUE)
  .fit_positive_scalar(chill_duration, "chill_duration", integer = TRUE)
  .fit_positive_scalar(bins, "bins", integer = TRUE)
  weather_freq_plot <- match.arg(weather_freq_plot)
  response <- .temperature_response_data(model, temp_values, chill_duration)
  frequencies <- NULL
  if (!is.null(hourtemps) && weather_freq_plot != "none") {
    frequencies <- .temperature_frequency_data(hourtemps, chill_months, heat_months,
      temp_values, bins, levels(response$component))
    peaks <- vapply(split(response$response, response$component), function(x) max(abs(x)), numeric(1))
    frequencies$response <- frequencies$frequency * peaks[as.character(frequencies$component)]
    if (weather_freq_plot == "gradient") {
      response$frequency <- vapply(seq_len(nrow(response)), function(i) {
        rows <- frequencies[frequencies$component == response$component[i], ]
        j <- findInterval(response$temperature[i], c(rows$lower, tail(rows$upper, 1)),
          left.open = TRUE, all.inside = TRUE)
        rows$frequency[j]
      }, numeric(1))
    }
  }
  result <- ggplot2::ggplot(response, ggplot2::aes(x = temperature, y = response))
  if (!is.null(frequencies) && weather_freq_plot == "histogram")
    result <- result + ggplot2::geom_rect(data = frequencies,
      ggplot2::aes(xmin = lower, xmax = upper, ymin = 0, ymax = response),
      inherit.aes = FALSE, fill = "#cbd3db", alpha = 0.65)
  if (!is.null(frequencies) && weather_freq_plot == "gradient") {
    result <- result + ggplot2::geom_line(ggplot2::aes(colour = frequency), linewidth = 1.2) +
      ggplot2::scale_colour_gradient(low = "#c7dfe8", high = "#164b65",
        limits = c(0, 1), name = "Relative temperature\nfrequency")
  } else {
    result <- result + ggplot2::geom_line(ggplot2::aes(colour = component), linewidth = 1.1) +
      ggplot2::scale_colour_manual(values = c("#287c9b", "#cf6940"), guide = "none")
  }
  result + ggplot2::facet_wrap(ggplot2::vars(component), ncol = 2, scales = "free_y") +
    ggplot2::labs(x = "Temperature (\u00b0C)", y = "Accumulation per hour") +
    .pheno_plot_theme()
}

#' @rdname plot_pheno_model
#' @export
plot_model_structure <- function(model, parameters = model_parameters(model),
                                  chill_range = NULL, weather = NULL) {
  .require_plot_package("ggplot2")
  model <- .single_plot_model(model, parameters)
  if (identical(model$chill$name, "none")) return(.plot_heat_only_structure(model, weather))
  p <- model$parameters
  partial <- identical(model$structure, "partial_overlap")
  if (is.null(chill_range)) {
    upper <- if (partial && p[["b3"]] > 0 && p[["b2"]] > 0 && p[["ol"]] > 0)
      p[["yc"]] + 4 / p[["b3"]] else 2 * p[["yc"]]
    if (model$chill$name == "dynamic") upper <- min(upper, 120)
    chill_range <- c(0, upper)
  }
  .plot_numeric_grid(chill_range, "chill_range", lower = 0)
  if (length(chill_range) != 2L) stop("chill_range must have two endpoints.", call. = FALSE)
  if (partial) return(.plot_partial_overlap(model, chill_range, weather))
  zero_end <- switch(model$structure, sequential = p[["yc"]],
    phenoflex = p[["yc"]] / (1 - stats::qlogis(0.01) / p[["s1"]] / p[["yc"]]), NULL)
  grid <- sort(unique(c(seq(chill_range[1], chill_range[2], length.out = 401),
    zero_end[zero_end >= chill_range[1] & zero_end <= chill_range[2]],
    p[["yc"]][p[["yc"]] >= chill_range[1] & p[["yc"]] <= chill_range[2]])))
  data <- data.frame(chill = grid, effectiveness = .chill_heat_effectiveness(model, grid))
  result <- ggplot2::ggplot(data, ggplot2::aes(x = chill, y = effectiveness))
  if (!is.null(zero_end) && chill_range[1] < zero_end) {
    shade_end <- min(zero_end, chill_range[2])
    narrow <- (shade_end - chill_range[1]) / diff(chill_range) < 0.2
    result <- result + ggplot2::annotate("rect", xmin = chill_range[1], xmax = shade_end,
      ymin = 0, ymax = 1, fill = "#ccd4dc", alpha = 0.2) +
      ggplot2::annotate("text", x = mean(c(chill_range[1], shade_end)), y = 0.5,
        label = if (model$structure == "phenoflex") {
          if (narrow) "Less than 1% effective heat" else "Very low heat effectiveness\nLess than 1%"
        } else if (narrow) "Flowering not possible" else "Flowering not possible\nEffective heat = 0",
        angle = if (narrow) 90 else 0, size = 3.4, colour = "#657380")
  }
  result <- result +
    ggplot2::geom_vline(xintercept = p[["yc"]], linetype = "dashed", colour = "#94a3b0")
  if (model$structure == "sequential") {
    result <- result + ggplot2::geom_step(linewidth = 1.1, colour = "#6d5799")
    points <- data.frame(chill = p[["yc"]], effectiveness = 1,
      label = paste0("Heat starts at yc = ", .plot_number(p[["yc"]])))
    caption <- "Heat starts at the first yc crossing and stays enabled."
  } else if (model$structure == "parallel") {
    result <- result + ggplot2::geom_line(linewidth = 1.1, colour = "#6d5799")
    points <- data.frame(chill = c(0, p[["yc"]]), effectiveness = c(p[["kmin"]], 1),
      label = c(paste0("kmin = ", .plot_number(p[["kmin"]])), paste0("yc = ", .plot_number(p[["yc"]]))))
    caption <- "Heat effectiveness increases linearly from kmin to 1 at yc."
  } else {
    result <- result + ggplot2::geom_line(linewidth = 1.1, colour = "#6d5799")
    points <- data.frame(chill = p[["yc"]], effectiveness = 0.5,
      label = paste0("yc = ", .plot_number(p[["yc"]]), "\n50% effective heat"))
    caption <- paste0("PhenoFlex: logistic[s1 * yc * (y - yc) / y], s1 = ", .plot_number(p[["s1"]]),
      ".")
  }
  points <- points[points$chill >= chill_range[1] & points$chill <= chill_range[2], ]
  result <- result + ggplot2::geom_point(data = points, size = 2.8, colour = "#6d5799") +
    ggplot2::geom_text(data = points, ggplot2::aes(label = label),
      nudge_y = -0.075, nudge_x = 0.01 * diff(chill_range), hjust = 0, size = 3.5, colour = "#394956")
  if (!is.null(weather)) {
    native <- predict_phenology(model, weather, stopatzc = FALSE, basic_output = FALSE)
    chill <- native$chill[, ncol(native$chill)]
    result <- result + ggplot2::geom_rug(data = data.frame(chill = chill),
      ggplot2::aes(x = chill), inherit.aes = FALSE, sides = "b", alpha = 0.15, colour = "#cf6940")
    caption <- paste(caption, "Orange rug: chill states in the supplied season.")
  }
  result + ggplot2::coord_cartesian(xlim = chill_range, ylim = c(-0.04, 1.06), clip = "off") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(title = paste(.plot_structure_name(model), "transition"),
      x = paste0("Accumulated chill (", .plot_chill_units(model), ")"),
      y = "Fraction of potential heat that is effective", caption = caption) + .pheno_plot_theme()
}

.require_plot_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the '", package, "' package to use this plotting function.", call. = FALSE)
}

.single_plot_model <- function(model, parameters) {
  force(parameters)
  if (inherits(model, "phenology_fit")) model <- model$model
  if (inherits(model, "population_pheno_model")) model <- model$model
  if (!inherits(model, "pheno_model"))
    stop("Plot one pheno_model; select a collection member with as_pheno_models().", call. = FALSE)
  as_pheno_models(model, parameters)[[1]]
}

.plot_numeric_grid <- function(x, name, lower, strict_lower = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 2L || any(!is.finite(x)) ||
      any(diff(x) <= 0) || any(if (strict_lower) x <= lower else x < lower))
    stop(name, " must be a finite, strictly increasing numeric vector with values ",
      if (strict_lower) "above " else "at least ", lower, ".", call. = FALSE)
}

.plot_number <- function(x) format(signif(x, 4), trim = TRUE, scientific = FALSE)
.plot_chill_name <- function(model) switch(model$chill$name, dynamic = "Dynamic Model",
  chill_hours = "Chill hours", utah = "Utah", positive_utah = "Positive Utah",
  custom = model$chill$label, none = "None")
.plot_structure_name <- function(model) if (model$chill$name == "none") "Heat-only" else
  switch(model$structure, sequential = "Sequential", parallel = "Parallel",
    phenoflex = "PhenoFlex", partial_overlap = "Partial overlap")
.plot_chill_units <- function(model) switch(model$chill$name, dynamic = "chill portions",
  chill_hours = "chill hours", utah = "chill units", positive_utah = "positive chill units",
  custom = "custom chill units", none = "none")
.plot_heat_units <- function(model) if (model$heat$scaling == "scaled") "GDH" else "unscaled GDH"

.pheno_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "#edf0f3"),
    panel.border = ggplot2::element_rect(colour = "#b8c3cd"),
    strip.background = ggplot2::element_rect(fill = "#f3f6f8", colour = "#b8c3cd"),
    strip.text = ggplot2::element_text(face = "bold", colour = "#394956"),
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.caption = ggplot2::element_text(hjust = 0, colour = "#52606d", size = 9),
    plot.margin = ggplot2::margin(10, 18, 10, 10), legend.position = "bottom")
}

.temperature_response_data <- function(model, temperatures, duration) {
  p <- model$parameters
  i <- seq.int(0, duration)
  dates <- as.Date("2001-01-01") + i %/% 24L
  weather <- data.frame(Temp = 0, Year = as.integer(format(dates, "%Y")),
    JDay = as.integer(format(dates, "%j")), Month = as.integer(format(dates, "%m")), Hour = i %% 24L)
  chill <- vapply(temperatures, function(temp) {
    if (model$chill$name == "none") return(0)
    weather$Temp <- rep(temp, nrow(weather))
    values <- if (model$chill$name == "positive_utah")
      calculate_chill_positive_utah(weather$Temp, seq_len(nrow(weather)),
        weights = unname(p[paste0("pos_utah_w", 1:7)]), day_group = i %/% 24L) else
      .calculate_pheno_inputs(model, weather, p)$chill
    column <- if (model$chill$name == "utah") "y_raw" else ncol(values)
    (values[nrow(values), column] - values[1L, column]) / duration
  }, numeric(1))
  heat_fun <- if (model$heat$scaling == "scaled") calculate_heat_gdh else calculate_heat_gdh_unscaled
  heat <- heat_fun(c(temperatures, tail(temperatures, 1)), seq_len(length(temperatures) + 1L),
    Tb = p[["Tb"]], Tu = p[["Tu"]], Tc = p[["Tc"]])
  panels <- c("Chill response", "Heat response")
  if (any(!is.finite(c(chill, heat)))) stop("Temperature responses must be finite.", call. = FALSE)
  data.frame(temperature = rep(temperatures, 2), response = c(chill, heat),
    component = factor(rep(panels, each = length(temperatures)), levels = panels))
}

.temperature_frequency_data <- function(weather, chill_months, heat_months, temperatures, bins, panels) {
  for (months in list(chill_months, heat_months)) {
    if (!is.numeric(months) || !is.null(dim(months)) || !length(months) || any(!is.finite(months)) ||
        any(months != floor(months) | months < 1 | months > 12))
      stop("chill_months and heat_months must contain calendar months 1:12.", call. = FALSE)
  }
  if (!is.data.frame(weather) || !"Temp" %in% names(weather) ||
      !is.numeric(weather$Temp) || any(!is.finite(weather$Temp)))
    stop("hourtemps must contain finite numeric Temp values.", call. = FALSE)
  months <- weather$Month
  if (is.null(months)) {
    if (!all(c("Year", "JDay") %in% names(weather)))
      stop("hourtemps needs Month, or Year and JDay to derive months.", call. = FALSE)
    if (!is.numeric(weather$Year) || !is.numeric(weather$JDay) ||
        any(!is.finite(c(weather$Year, weather$JDay))) ||
        any(weather$Year != floor(weather$Year) | weather$Year < 1 | weather$Year > 9999) ||
        any(weather$JDay != floor(weather$JDay)))
      stop("hourtemps has invalid Year or JDay values.", call. = FALSE)
    leap <- weather$Year %% 4 == 0 & (weather$Year %% 100 != 0 | weather$Year %% 400 == 0)
    if (any(weather$JDay < 1 | weather$JDay > 365 + leap))
      stop("hourtemps has invalid Year or JDay values.", call. = FALSE)
    dates <- as.Date(sprintf("%04d-01-01", as.integer(weather$Year))) + weather$JDay - 1L
    months <- as.integer(format(dates, "%m"))
  }
  if (!is.numeric(months) || any(!is.finite(months)) || any(months != floor(months) | months < 1 | months > 12))
    stop("hourtemps Month must contain calendar months 1:12.", call. = FALSE)
  breaks <- seq(min(temperatures), max(temperatures), length.out = bins + 1L)
  result <- lapply(seq_along(panels), function(i) {
    use <- months %in% list(chill_months, heat_months)[[i]] &
      weather$Temp >= breaks[1] & weather$Temp <= tail(breaks, 1)
    counts <- graphics::hist(weather$Temp[use], breaks = breaks, plot = FALSE)$counts
    frequency <- if (max(counts) > 0) counts / max(counts) else rep(0, bins)
    data.frame(lower = head(breaks, -1), upper = tail(breaks, -1), frequency = frequency,
      component = factor(panels[i], levels = panels))
  })
  do.call(rbind, result)
}

.chill_heat_effectiveness <- function(model, chill) {
  p <- model$parameters
  switch(model$structure,
    sequential = as.numeric(chill >= p[["yc"]]),
    parallel = pmin(1, p[["kmin"]] + (1 - p[["kmin"]]) * (chill / p[["yc"]])),
    phenoflex = {
      result <- numeric(length(chill))
      positive <- chill > 0
      a <- p[["s1"]] * p[["yc"]] * (1 - p[["yc"]] / chill[positive])
      a[chill[positive] == p[["yc"]]] <- 0
      result[positive] <- ifelse(a >= 17, 1, ifelse(a <= -20, 0, stats::plogis(a)))
      result
    })
}

.partial_overlap_requirement <- function(p, additional_chill) {
  additional <- if (p[["ol"]] == 0) rep(0, length(additional_chill)) else additional_chill
  p[["b1"]] + p[["b2"]] * exp(-p[["b3"]] * additional)
}

.plot_partial_overlap <- function(model, chill_range, weather) {
  p <- model$parameters
  # The static curve starts at the exact minimum; native seasons can overshoot it.
  onset <- p[["yc"]]
  if (!is.null(weather)) {
    trajectory <- .partial_overlap_plot_trajectory(model, weather)
    if (!is.na(attr(trajectory, "chill_at_completion")))
      onset <- attr(trajectory, "chill_at_completion")
  }
  grid <- sort(unique(c(seq(chill_range[1], chill_range[2], length.out = 401),
    onset[onset >= chill_range[1] & onset <= chill_range[2]])))
  data <- data.frame(chill = grid, additional_chill = pmax(0, grid - onset))
  data$requirement <- .partial_overlap_requirement(p, data$additional_chill)
  data$requirement[data$chill < onset] <- NA_real_
  if (any(!is.finite(data$requirement[!is.na(data$requirement)])))
    stop("Partial-overlap heat requirement overflowed.", call. = FALSE)
  points <- data.frame(chill = onset, requirement = p[["b1"]] + p[["b2"]],
    label = paste0(if (onset == p[["yc"]]) "Minimum chilling" else "Minimum chilling reached",
      "\n", .plot_number(onset), " ", .plot_chill_units(model),
      "; heat required = ", .plot_number(p[["b1"]] + p[["b2"]])))
  points <- points[points$chill >= chill_range[1] & points$chill <= chill_range[2], ]
  guides <- data.frame(requirement = c(p[["b1"]], p[["b1"]] * p[["ol"]]),
    guide = c(paste0("Lower heat asymptote: ", .plot_number(p[["b1"]])),
      paste0("End of overlap at accumulated heat: ", .plot_number(p[["b1"]] * p[["ol"]]))))
  result <- ggplot2::ggplot(data, ggplot2::aes(x = chill, y = requirement)) +
    ggplot2::geom_hline(data = guides, ggplot2::aes(yintercept = requirement, linetype = guide),
      colour = "#94a3b0", linewidth = 0.6) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("dashed", "dotted"), guides$guide), name = NULL) +
    ggplot2::geom_line(linewidth = 1.2, colour = "#6d5799", na.rm = TRUE) +
    ggplot2::geom_point(data = points, size = 2.8, colour = "#6d5799") +
    ggplot2::geom_text(data = points, ggplot2::aes(label = label),
      nudge_y = 0.055 * (p[["b1"]] + p[["b2"]]), nudge_x = 0.01 * diff(chill_range),
      hjust = 0, vjust = 0, size = 3.4, colour = "#394956")
  if (chill_range[1] < p[["yc"]]) {
    shade_end <- min(p[["yc"]], chill_range[2])
    narrow <- (shade_end - chill_range[1]) / diff(chill_range) < 0.2
    result <- result + ggplot2::annotate("rect", xmin = chill_range[1],
      xmax = shade_end, ymin = 0, ymax = Inf, fill = "#ccd4dc", alpha = 0.2) +
      ggplot2::annotate("text", x = mean(c(chill_range[1], shade_end)),
        y = 0.5 * (p[["b1"]] + p[["b2"]]),
        label = if (narrow) "Insufficient chilling" else "Insufficient chilling\nHeat cannot trigger\nthe event",
        angle = if (narrow) 90 else 0, size = 3.4, colour = "#657380")
  }
  caption <- if (p[["ol"]] == 0) "No overlap: additional chilling does not reduce the heat requirement." else
    "Additional chilling reduces the heat requirement until the end of overlap."
  if (!is.null(weather)) {
    bloom <- attr(trajectory, "bloomindex")
    if (bloom > 0) trajectory <- trajectory[seq_len(bloom), ]
    trajectory$chill <- trajectory$accumulated_chill
    result <- result + ggplot2::geom_path(data = trajectory,
      ggplot2::aes(x = chill, y = heat, colour = "Supplied season"), inherit.aes = FALSE,
      linewidth = 0.8, alpha = 0.85) +
      ggplot2::scale_colour_manual(values = c("Supplied season" = "#cf6940"), name = NULL)
    marked <- trajectory[trajectory$event != "", ]
    if (nrow(marked)) result <- result + ggplot2::geom_point(data = marked,
      ggplot2::aes(x = chill, y = heat), inherit.aes = FALSE, size = 2.5, colour = "#cf6940") +
      ggplot2::geom_text(data = marked, ggplot2::aes(x = chill, y = heat, label = event),
        inherit.aes = FALSE, nudge_x = 0.015 * diff(chill_range), hjust = 0, vjust = -0.5,
        size = 3.3, colour = "#a34d2b")
    caption <- paste(caption, "Orange: accumulated heat versus chill used for the transition.", sep = "\n")
  }
  result + ggplot2::coord_cartesian(xlim = chill_range, ylim = c(0, NA), clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0), breaks = pretty(chill_range, n = 6)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.2))) +
    ggplot2::labs(title = "Partial-overlap transition", x = paste0("Accumulated chill (",
      .plot_chill_units(model), ")"),
      y = paste0(if (is.null(weather)) "Heat requirement (" else "Heat (", .plot_heat_units(model), ")"),
      caption = caption) + .pheno_plot_theme()
}

.partial_overlap_plot_trajectory <- function(model, weather) {
  native <- predict_phenology(model, weather, stopatzc = FALSE, basic_output = FALSE)
  p <- model$parameters
  y <- native$chill[, ncol(native$chill)]
  heat <- native$z
  additional <- numeric(length(y))
  event <- rep("", length(y))
  onset <- which(y >= p[["yc"]])[1]
  if (!is.na(onset)) {
    rows <- seq.int(onset, length(y))
    additional[rows] <- pmax(0, y[rows] - y[onset])
    cutoff <- rows[heat[rows] >= p[["b1"]] * p[["ol"]]][1]
    if (!is.na(cutoff)) {
      additional[seq.int(cutoff, length(y))] <- additional[cutoff]
      event[cutoff] <- "Overlap cutoff"
    }
  }
  if (native$bloomindex > 0) {
    i <- native$bloomindex
    event[i] <- if (nzchar(event[i])) paste(event[i], "/ event") else "Event"
  }
  accumulated <- y
  if (!is.na(onset)) accumulated[seq.int(onset, length(y))] <- y[onset] + additional[seq.int(onset, length(y))]
  result <- data.frame(chill = additional, accumulated_chill = accumulated, heat = heat, event = event)
  attr(result, "bloomindex") <- native$bloomindex
  attr(result, "chill_at_completion") <- if (is.na(onset)) NA_real_ else y[onset]
  result
}

.plot_heat_only_structure <- function(model, weather) {
  start <- model$parameters[["heat_start_index"]]
  active <- ceiling(start)
  finish <- if (is.null(weather)) max(2, active + 48) else nrow(weather) - 1L
  grid <- sort(unique(c(1, seq_len(min(finish, 1000)), active, finish)))
  data <- data.frame(index = grid, effectiveness = as.numeric(grid >= active))
  caption <- paste0("Chill is bypassed. heat_start_index = ", .plot_number(start),
    " enables heat from interval row ", .plot_number(active),
    "; event heat requirement zc = ", .plot_number(model$parameters[["zc"]]), " ", .plot_heat_units(model), ".")
  if (!is.null(weather)) .validate_sequential_weather(weather)
  ggplot2::ggplot(data, ggplot2::aes(x = index, y = effectiveness)) +
    ggplot2::geom_step(linewidth = 1.1, colour = "#6d5799") +
    ggplot2::geom_vline(xintercept = active, linetype = "dashed", colour = "#94a3b0") +
    ggplot2::labs(title = "Heat-only start", x = "Interval row (one-based)",
      y = "Fraction of potential heat that is effective", caption = caption) + .pheno_plot_theme()
}

utils::globalVariables(c("temperature", "response", "component", "lower", "upper", "frequency",
  "chill", "effectiveness", "label", "requirement", "guide", "heat", "event", "index"))
