#' Evaluate phenology model parameters for GenSA or DEoptim
#'
#' Evaluate parameters for a PhenoFlex, sequential, parallel or partial overlap
#' model during calibration. Return the residual sum of squares, or Inf when
#' parameter constraints are violated.
#'
#' Optional controls convert characteristic Dynamic Model parameters to kinetic
#' parameters and check additional parameter constraints before prediction.
#'
#' @param x Numeric model parameter vector in the order required by modelfn.
#'
#' @param modelfn Function that predicts a bloom date from a season's weather
#' and model parameters, such as \code{\link{custom_PhenoFlex_GDHwrapper}}.
#'
#' @param bloomJDays Numeric vector of observed bloom dates, expressed as days of
#' the year and ordered to match SeasonList.
#'
#' @param SeasonList List of hourly weather data frames, one per phenological season.
#' Each data frame contains Temp (hourly temperature) and JDay (day of the year).
#' Season lists can be generated with \code{\link[chillR]{genSeasonList}}.
#'
#' @param na_penalty Numeric replacement for a missing bloom prediction. Defaults
#' to 365; this value is used in place of the missing date when calculating loss.
#'
#' @param intermed_chill Logical. If TRUE, convert characteristic chill parameters
#' (theta_star, theta_c, tau and pie_c) to kinetic parameters (E0, E1, A0 and A1).
#' Defaults to FALSE. Characteristic parameters provide a more narrowly defined
#' search space.
#'
#' @param pos_int Numeric vector of four positions in x to convert when
#' intermed_chill is TRUE, in the order theta_star, theta_c, tau, pie_c.
#'
#' @param check_constrain Logical. If TRUE, check Q10 values and GDH temperature
#' constraints before prediction. Defaults to FALSE. Q10 values for E0 and E1
#' must lie within the supplied tolerance limits, and GDH temperatures must be
#' ordered. Return Inf if a constraint is violated.
#'
#' @param pos_E0_E1 Numeric vector giving the positions of E0 and E1 in x.
#' Used when checking constraints.
#'
#' @param pos_heat Numeric vector giving the positions of Tb, Tu and Tc in x.
#' Used when checking constraints.
#'
#' @param c_L Numeric vector of five lower constraint limits.
#'
#' @param c_U Numeric vector of five upper constraint limits.
#'
#' @return Numeric residual sum of squares between predicted and observed dates. Return Inf if parameter constraints are violated.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#'
#' @section Deprecated:
#' \code{eval_all_daoptim()} is deprecated since evalpheno 0.3.1.
#' Use fit_phenology() with phenology_rss(); choose fixed parameters through named bounds.
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export eval_all_daoptim
eval_all_daoptim <- function(x, modelfn, bloomJDays, SeasonList, na_penalty = 365, intermed_chill = FALSE,
                     pos_int = c(5:8), check_constrain = FALSE,
                     pos_E0_E1 = c(5,6),
                     pos_heat = c(11, 4, 10),
                     c_L = c(  0,   0,   0,     1.5, 1.5), 
                     c_U = c(Inf, Inf, Inf,     3.5, 3.5)){
  .Deprecated(new = "phenology_rss", package = "evalpheno",
              msg = "eval_all_daoptim() is deprecated. Use fit_phenology() with phenology_rss(); choose fixed parameters through named bounds. See help('evalpheno-deprecated') for migration differences.")
  if(intermed_chill){
    int <- characteristic_to_kinetic(c(rep(0, 4), x[pos_int], rep(0,4)))
    
    if(is.list(int)) return(Inf)
    
    x[pos_int] <- int[5:8]
  }
  
  #wrapper_seq_model(x = SeasonList[[1]], par = par)
  
  par <- x
  if (!check_constrain) {
    return(.phenology_rss(par, modelfn, bloomJDays, SeasonList, na_penalty))
  }
  g <- rep(0, 5)
  g[1] <- par[pos_heat[2]] - par[pos_heat[1]]
  g[2] <- par[pos_heat[3]] - par[pos_heat[1]]
  g[3] <- par[pos_heat[3]] - par[pos_heat[2]]
  g[4] <- exp((10 * par[pos_E0_E1[1]])/(297 * 279)) 
  g[5] <- exp((10 * par[pos_E0_E1[2]])/(297 * 279))
  
  if(any(!is.finite(g)) || any(g <= c_L) || any(g >= c_U)){
    return(Inf)
  }
  .phenology_rss(par, modelfn, bloomJDays, SeasonList, na_penalty)
}
