#' Stable residuals for the Dynamic Model parameter conversion
#'
#' Evaluate the two logarithmic residuals used to convert intermediate chill
#' parameters into E0 and E1. Uses expm1/log1p and log-scale rates to avoid
#' cancellation and intermediate overflow. Reference temperatures are 297 and
#' 279 K, with eta = 1/3. Equation 38 requires theta_c < 297 K for a real log.
#' These functions evaluate residuals; they do not themselves find a root.
#' @param x Numeric vector c(E0, E1), with finite 0 < E0 < E1.
#' @param params Numeric vector c(theta_star, theta_c, tau, pie_c), with
#' 0 < theta_star < theta_c < 297, tau > 0 and pie_c > 0; all finite.
#' @return Numeric vector of two logarithmic residuals, zero at a solution.
#' @seealso \code{\link{solve_nle_transformed}}, \code{\link{convert_parameters}}
#' @export
solve_nle <- function(x, params) {
  .validate_chill_params(params)
  if (!.valid_chill_energies(x))
    stop("x must contain finite energies satisfying 0 < E0 < E1.", call. = FALSE)
  .chill_residuals(x, params)
}

#' Dynamic Model residuals with positive, ordered energies
#'
#' Uses E0 = exp(z[1]) and E1 = E0 + exp(z[2]). Supply log-transformed starting
#' values and transform the solution back before using it as model parameters.
#' Extreme trials that overflow, underflow, or round E1 to E0 return non-finite
#' residuals, allowing nleqslv to backtrack rather than abort on an assertion.
#' @param z Finite numeric vector c(log(E0), log(E1 - E0)).
#' @inheritParams solve_nle
#' @return Two residuals, or c(Inf, Inf) for an unrepresentable trial point.
#' @examples
#' params <- c(279, 286.1, 47.7, 28)
#' fit <- nleqslv::nleqslv(log(c(500, 15000 - 500)),
#'                        solve_nle_transformed, params = params,
#'                        method = "Newton", xscalm = "auto")
#' E0 <- exp(fit$x[1])
#' E1 <- E0 + exp(fit$x[2])
#' solve_nle(c(E0, E1), params)
#' @export
solve_nle_transformed <- function(z, params) {
  .validate_chill_params(params)
  if (!is.numeric(z) || length(z) != 2L || any(!is.finite(z)))
    stop("z must contain two finite log-transformed energies.", call. = FALSE)
  x <- .decode_chill_energies(z)
  if (!.valid_chill_energies(x)) return(c(Inf, Inf))
  result <- .chill_residuals(x, params)
  if (any(!is.finite(result))) c(Inf, Inf) else result
}

.validate_chill_params <- function(params) {
  if (!is.numeric(params) || length(params) != 4L || any(!is.finite(params)) ||
      params[1] <= 0 || params[2] <= params[1] || params[2] >= 297 ||
      params[3] <= 0 || params[4] <= 0)
    stop("params must satisfy 0 < theta_star < theta_c < 297, tau > 0, pie_c > 0.",
         call. = FALSE)
}

.valid_chill_energies <- function(x) {
  is.numeric(x) && length(x) == 2L && all(is.finite(x)) && x[1] > 0 && x[2] > x[1]
}

.decode_chill_energies <- function(z) {
  E0 <- exp(z[1])
  c(E0, E0 + exp(z[2]))
}

# Stable log(1 - exp(-u)), for positive u.
.log1mexp_chill <- function(u) {
  if (u <= log(2)) log(-expm1(-u)) else log1p(-exp(-u))
}

# For u > 36, the relative correction to -log(1-exp(-u)) ~ exp(-u)
# is below double precision. Staying on the log scale avoids underflow.
.log_neg_log_chill <- function(u) {
  if (u > 36) -u else log(-.log1mexp_chill(u))
}

# log(1 - exp(-rate)), given log(rate), without constructing extreme rates.
.log_rate_probability <- function(log_rate) {
  if (log_rate < -36) return(log_rate)
  if (log_rate > log(36)) return(0)
  log(-expm1(-exp(log_rate)))
}

.chill_residuals <- function(x, params) {
  E0 <- x[1]
  E1 <- x[2]
  delta <- E1 - E0
  theta_star <- params[1]
  theta_c <- params[2]
  u <- (delta / theta_star) * ((theta_c - theta_star) / theta_c)
  if (!is.finite(u) || u <= 0) return(c(Inf, Inf))
  log_term <- .log_neg_log_chill(u)
  residual_1 <- if (u > 36) log1p(-E0 / E1) else
    log(delta) - log(expm1(u)) - log_term - log(E1)
  a <- (delta / theta_c) * ((297 - theta_c) / 297)
  b <- (delta / 279) * ((297 - 279) / 297)
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(c(Inf, Inf))
  log_lhs <- if (a > 36 && b > 36) {
    (delta / theta_c) * ((279 - theta_c) / 279) +
      .log1mexp_chill(a) - .log1mexp_chill(b)
  } else {
    log_a <- if (a > 36) a + .log1mexp_chill(a) else log(expm1(a))
    log_b <- if (b > 36) b + .log1mexp_chill(b) else log(expm1(b))
    log_a - log_b
  }
  log_k1 <- (E1 / theta_star) * ((297 - theta_star) / 297) - log(params[3]) + log_term
  log_k2 <- (E1 / theta_star) * ((279 - theta_star) / 279) - log(params[3]) + log_term
  log_num <- log_k2 + log(2 / 3) + log(params[4])
  log_other <- log_k1 + log(1 / 3) + log(params[4])
  m <- max(log_num, log_other)
  if (!is.finite(m)) return(c(Inf, Inf))
  log_den <- m + log1p(exp(min(log_num, log_other) - m))
  residual_2 <- log_lhs -
    (.log_rate_probability(log_num) - .log_rate_probability(log_den))
  c(residual_1, residual_2)
}

# Shared solver: physical energies in $x, legacy termination field for adapters.
.fit_chill_parameters <- function(params, start = c(500, 15000)) {
  failure <- function(message) list(x = c(NA_real_, NA_real_), A0 = NA_real_,
    A1 = NA_real_, fvec = c(Inf, Inf), termcd = 3L, message = message)
  valid <- tryCatch({ .validate_chill_params(params); TRUE }, error = conditionMessage)
  if (!isTRUE(valid)) return(failure(valid))
  if (!.valid_chill_energies(start)) return(failure("Invalid starting energies."))
  output <- tryCatch(nleqslv::nleqslv(
    log(c(start[1], start[2] - start[1])), solve_nle_transformed,
    params = params, xscalm = "auto", method = "Newton",
    control = list(trace = 0, allowSingular = TRUE, ftol = 1e-8)),
    error = function(e) failure(conditionMessage(e)))
  if (any(!is.finite(output$x))) return(output)
  output$z <- output$x
  output$x <- .decode_chill_energies(output$z)
  if (!.valid_chill_energies(output$x)) return(failure("Unrepresentable energies."))
  delta <- output$x[2] - output$x[1]
  u <- (delta / params[1]) * ((params[2] - params[1]) / params[2])
  log_A1 <- output$x[2] / params[1] - log(params[3]) + .log_neg_log_chill(u)
  output$A1 <- exp(log_A1)
  output$A0 <- exp(log_A1 - delta / params[2])
  output$fvec <- solve_nle(output$x, params)
  if (any(!is.finite(c(output$A0, output$A1, output$fvec))) ||
      output$A0 <= 0 || output$A1 <= 0 || max(abs(output$fvec)) > 1e-8) {
    output$termcd <- 3L
    output$message <- "Parameter conversion did not produce a usable root."
  }
  output
}
