#' Sample bud-specific chill and heat requirements
#'
#' Independent normal or skew-normal requirements. Zero standard deviation
#' produces a constant population. Skew-normal parameters describe the actual
#' mean and standard deviation, with shape supplied in `skew`.
#' Negative draws are rejected with an error rather than silently truncated.
#' @param n Positive integer population size.
#' @param yc,zc Positive mean chill and heat requirements.
#' @param yc_sd,zc_sd Non-negative standard deviations.
#' @param dist_chill,dist_heat Either `normal` or `normal_skewed`.
#' @param skew Two finite shape parameters, for chill and heat respectively.
#' @param seed Optional non-negative integer seed. A supplied seed gives
#' reproducible draws without changing the caller's random-number state.
#' @return A list containing numeric vectors `yc_pop` and `zc_pop`, each of length `n`.
#' @export
sample_bud_population <- function(n = 100, yc = 40, zc = 190,
                                  yc_sd = 0, zc_sd = 0,
                                  dist_chill = "normal", dist_heat = "normal",
                                  skew = c(0, 0), seed = 12345) {
  .population_scalar(n, "n", positive = TRUE, integer = TRUE)
  .population_scalar(yc, "yc", positive = TRUE)
  .population_scalar(zc, "zc", positive = TRUE)
  .population_scalar(yc_sd, "yc_sd", nonnegative = TRUE)
  .population_scalar(zc_sd, "zc_sd", nonnegative = TRUE)
  dist_chill <- match.arg(dist_chill, c("normal", "normal_skewed"))
  dist_heat <- match.arg(dist_heat, c("normal", "normal_skewed"))
  if (!is.numeric(skew) || length(skew) != 2L || any(!is.finite(skew))) {
    stop("skew must contain two finite shape parameters.", call. = FALSE)
  }
  if (!is.null(seed)) {
    .population_scalar(seed, "seed", nonnegative = TRUE, integer = TRUE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  draw <- function(mean, sd, distribution, shape) {
    if (sd == 0) return(rep(mean, n))
    if (distribution == "normal") return(stats::rnorm(n, mean, sd))
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Install 'sn' to sample skew-normal requirements.", call. = FALSE)
    }
    delta <- shape / sqrt(1 + shape^2)
    omega <- sd / sqrt(1 - 2 * delta^2 / pi)
    xi <- mean - omega * delta * sqrt(2 / pi)
    as.vector(sn::rsn(n, xi = xi, omega = omega, alpha = shape))
  }
  out <- list(yc_pop = draw(yc, yc_sd, dist_chill, skew[1]),
              zc_pop = draw(zc, zc_sd, dist_heat, skew[2]))
  if (any(!is.finite(unlist(out))) || any(unlist(out) <= 0)) {
    stop("Sampled requirements must be positive; reduce dispersion or change the seed.",
         call. = FALSE)
  }
  out
}

.population_scalar <- function(x, name, positive = FALSE, nonnegative = FALSE,
                               integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (positive && x <= 0) || (nonnegative && x < 0) ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop("Invalid value for ", name, ".", call. = FALSE)
  }
  invisible(x)
}
