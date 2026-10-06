#' Convert characteristic chill parameters to kinetic PhenoFlex parameters
#'
#' Uses stable residuals and the log-transformed energy solver. Positions outside
#' 5:8 are preserved. Failure policies retain the LarsChill calling convention.
#' @param par Numeric vector of length 12, ordered yc, zc, s1, Tu, theta_star,
#' theta_c, tau, pie_c, Tf, Tc, Tb, slope.
#' @param failure_return One of "MEIGO" (a penalty list), "NA" (NA at positions
#' 5:8), or "Ignore Error" (the last estimate, which may be invalid).
#' @return A vector with positions 5:8 replaced by E0, E1, A0, A1, or
#' list(F = 1e6, g = rep(1e6, 5)) on failure with policy "MEIGO".
#' @seealso \code{\link{solve_nle_transformed}}
#' @export
characteristic_to_kinetic <- function(par, failure_return = "MEIGO") {
  if (!is.numeric(par) || length(par) != 12L)
    stop("par must be a numeric PhenoFlex vector of length 12.", call. = FALSE)
  failure_return <- match.arg(failure_return, c("MEIGO", "NA", "Ignore Error"))
  output <- .fit_chill_parameters(unname(par[5:8]))
  if (output$termcd >= 3L && failure_return != "Ignore Error") {
    if (failure_return == "MEIGO") return(list(F = 1e6, g = rep(1e6, 5)))
    par[5:8] <- NA_real_
  } else {
    par[5:8] <- c(output$x, output$A0, output$A1)
  }
  if (!is.null(names(par))) names(par)[5:8] <- c("E0", "E1", "A0", "A1")
  par
}

#' @rdname characteristic_to_kinetic
#' @export
convert_parameters <- function(par, failure_return = "MEIGO") {
  characteristic_to_kinetic(par, failure_return = failure_return)
}
