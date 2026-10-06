#' Run the PhenoFlex bud population model
#'
#' Chill accumulation is shared by all buds; heat accumulation is computed
#' separately for each chill and heat requirement. Temperatures are in Celsius
#' and rows must represent consecutive hours. Forcing keeps chill fixed at
#' cutting and applies constant temperature, following the experimental model.
#' The original C++ `exp` output uses zero-based forcing step indices (9999 for
#' failure), whereas `bloomindex` is a one-based input row (0 for no bloom).
#' These conventions are retained for compatibility with experimental analyses.
#' @param temp_df Data frame with finite `Temp` and `JDay` columns and optional `Year`.
#' @param par Twelve PhenoFlex parameters in order: yc, zc, s1, Tu, E0, E1,
#' A0, A1, Tf, Tc, Tb, slope. A fully named vector may use any order.
#' @param yc_sd,zc_sd Standard deviations of bud requirements.
#' @param jday_cut Optional cutting days. Each day must identify exactly one date.
#' @param n Positive integer number of buds.
#' @param adjust_zc Positive multiplier for the heat requirement during forcing.
#' @param basic_output If TRUE, return bloom indices and sampled requirements only.
#' @param stop_at_zc Stop heat accumulation at bloom. With cutting days, complete
#' trajectories are always computed to avoid zero-filled values after bloom.
#' @param dist_chill,dist_heat `normal` or `normal_skewed`.
#' @param add_distpar Two skew-normal shape parameters; NULL means zero shapes.
#' @param seed Seed passed to [sample_bud_population()].
#' @param max_days_forcing Positive integer duration of a forcing experiment.
#' @param forcing_temperature Constant forcing temperature in Celsius.
#' @param population Optional list with positive `yc_pop` and `zc_pop` vectors
#' of length `n`. Supplying it bypasses sampling, useful for calibration.
#' @return A list with `bloomindex`, `yc_pop`, and `zc_pop`. Detailed output adds
#' shared chill vectors `x`, `y`, an hours-by-buds heat matrix `z`, cutting-days-by-buds
#' matrices `exp` and `z_effec`, and the potential forcing increment `exp_force_inc`.
#' @examples
#' weather <- data.frame(Temp = rep(15, 240), JDay = rep(1:10, each = 24))
#' par <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e13,
#'          4, 36, 4, 1.6)
#' result <- phenoflex_population(weather, par, n = 2)
#' @export
phenoflex_population <- function(temp_df, par, yc_sd = 0, zc_sd = 0,
                                 jday_cut = NULL, n = 100, adjust_zc = 1,
                                 basic_output = FALSE, stop_at_zc = TRUE,
                                 dist_chill = "normal", dist_heat = "normal",
                                 add_distpar = NULL, seed = 12345,
                                 max_days_forcing = 50, forcing_temperature = 23,
                                 population = NULL) {
  parameter_names <- c("yc", "zc", "s1", "Tu", "E0", "E1", "A0", "A1",
                       "Tf", "Tc", "Tb", "slope")
  if (!is.numeric(par) || length(par) != 12L || any(!is.finite(par))) {
    stop("par must contain twelve finite PhenoFlex parameters.", call. = FALSE)
  }
  if (!is.null(names(par))) {
    if (anyDuplicated(names(par)) || !setequal(names(par), parameter_names)) {
      stop("Named par must use the twelve documented parameter names.", call. = FALSE)
    }
    par <- par[parameter_names]
  }
  if (any(par[c(1:3, 5:8, 12)] <= 0) ||
      !(par[11] < par[4] && par[4] < par[10]) || any(par[c(4,9:11)] <= -273)) {
    stop("Requirements and rate parameters must be positive; Tb < Tu < Tc is required.",
         call. = FALSE)
  }
  if (!is.data.frame(temp_df) || !all(c("Temp", "JDay") %in% names(temp_df)) ||
      nrow(temp_df) < 2L || !is.numeric(temp_df$Temp) ||
      any(!is.finite(temp_df$Temp)) || any(temp_df$Temp <= -273) ||
      !is.numeric(temp_df$JDay) || any(!is.finite(temp_df$JDay))) {
    stop("temp_df must contain at least two valid hourly Temp and JDay rows.", call. = FALSE)
  }
  if (!is.null(temp_df$Year) &&
      (!is.numeric(temp_df$Year) || any(!is.finite(temp_df$Year)))) {
    stop("Year must contain finite years.", call. = FALSE)
  }
  .population_scalar(n, "n", positive = TRUE, integer = TRUE)
  .population_scalar(adjust_zc, "adjust_zc", positive = TRUE)
  .population_scalar(max_days_forcing, "max_days_forcing", positive = TRUE, integer = TRUE)
  .population_scalar(forcing_temperature, "forcing_temperature")
  if (forcing_temperature <= -273) stop("Invalid forcing temperature.", call. = FALSE)
  if (!is.logical(basic_output) || length(basic_output) != 1L || is.na(basic_output) ||
      !is.logical(stop_at_zc) || length(stop_at_zc) != 1L || is.na(stop_at_zc)) {
    stop("Output flags must be single non-missing logical values.", call. = FALSE)
  }
  if (is.null(population)) {
    # Preserve the experimental helper's deterministic single-bud convention.
    population <- sample_bud_population(n, par[1], par[2],
      if (n == 1L) 0 else yc_sd, if (n == 1L) 0 else zc_sd,
      dist_chill, dist_heat, if (is.null(add_distpar)) c(0, 0) else add_distpar, seed)
  }
  if (!is.list(population) || !all(c("yc_pop", "zc_pop") %in% names(population)) ||
      !all(vapply(population[c("yc_pop", "zc_pop")], function(x)
        is.numeric(x) && length(x) == n && all(is.finite(x)) && all(x > 0), logical(1)))) {
    stop("population must contain n positive yc_pop and zc_pop values.", call. = FALSE)
  }
  cuts <- NULL
  if (length(jday_cut)) {
    if (!is.numeric(jday_cut) || any(!is.finite(jday_cut))) {
      stop("jday_cut must contain finite days.", call. = FALSE)
    }
    cuts <- vapply(jday_cut, function(day) {
      rows <- which(temp_df$JDay == day)
      if (!length(rows) || any(diff(rows) != 1L) ||
          (!is.null(temp_df$Year) && length(unique(temp_df$Year[rows])) != 1L)) {
        stop("Each cutting day must identify one contiguous date in temp_df.", call. = FALSE)
      }
      floor(stats::median(rows)) - 1
    }, numeric(1))
  }
  out <- PhenoFlex_pop_slim(temp_df$Temp, seq_len(nrow(temp_df)),
    population$yc_pop, population$zc_pop, max_days_forcing = max_days_forcing,
    i_cut = cuts, forcing_temperature = forcing_temperature,
    s1 = par[3], E0 = par[5], E1 = par[6], A0 = par[7], A1 = par[8],
    Tf = par[9], slope = par[12], Tb = par[11], Tu = par[4], Tc = par[10],
    adjust_zc_forcing_exp = adjust_zc, basic_output = basic_output,
    stopatzc = stop_at_zc && !length(cuts))
  out$yc_pop <- population$yc_pop
  out$zc_pop <- population$zc_pop
  out
}

#' @rdname phenoflex_population
#' @export
helper_run_pop_model <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100,
                                 adjust_zc = 1, basic_output = FALSE, stop_at_zc = TRUE,
                                 dist_chill = "normal", dist_heat = "normal",
                                 add_distpar = NULL) {
  phenoflex_population(temp_df, par, yc_sd, zc_sd, jday_cut, n, adjust_zc,
                        basic_output, stop_at_zc, dist_chill, dist_heat, add_distpar)
}
