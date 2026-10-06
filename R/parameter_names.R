#' PhenoFlex parameter names for the two chill parameterizations
#'
#' Kinetic parameters specify the reaction rates (E0, E1, A0, A1).
#' Characteristic parameters describe temperatures and times (theta_star,
#' theta_c, tau, pie_c). Only positions 5:8 differ; the other eight PhenoFlex
#' parameters are shared. The old/new names are compatibility aliases.
#' The same vectors can also be loaded with data().
#' @format Character vectors of length 12, in model input order.
#' @export
phenoflex_parnames_kinetic <- c("yc", "zc", "s1", "Tu", "E0", "E1", "A0", "A1",
                               "Tf", "Tc", "Tb", "slope")

#' @rdname phenoflex_parnames_kinetic
#' @export
phenoflex_parnames_characteristic <- c("yc", "zc", "s1", "Tu", "theta_star",
                                      "theta_c", "tau", "pie_c", "Tf", "Tc", "Tb", "slope")

#' @rdname phenoflex_parnames_kinetic
#' @export
phenoflex_parnames_old <- phenoflex_parnames_kinetic

#' @rdname phenoflex_parnames_kinetic
#' @export
phenoflex_parnames_new <- phenoflex_parnames_characteristic
