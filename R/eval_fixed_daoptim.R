#' Evaluate structure parameters with fixed chill and heat submodels
#'
#' Evaluate the structure parameters of a PhenoFlex, sequential, parallel or
#' partial overlap model while keeping chill and heat submodel parameters fixed.
#' Return the residual sum of squares in a format compatible with GenSA and
#' DEoptim.
#'
#' @param x Numeric vector containing only the structure parameters required
#' by modelfn. Chill and heat submodel parameters are supplied separately.
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
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature (in Kelvin) for the sigmoidal function
#' @param Tb numeric. GDH base temperature (lower threshold)
#' @param Tu numeric. GDH optimal temperature
#' @param Tc Numeric upper temperature threshold of the GDH model, in degrees Celsius.
#'
#' @return Numeric residual sum of squares between predicted and observed dates.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#'
#' @export eval_fixed_daoptim
eval_fixed_daoptim <- function(x, modelfn, bloomJDays, SeasonList, na_penalty = 365,
                       E0 = 4153.5, E1 = 12888.8, A0 = 139500, A1 = 2567000000000000000, 
                       Tf=4, slope=1.6, Tb = 4, Tu = 25, Tc = 36){
  par <- c(x, E0, E1, A0, A1, Tf, slope, Tb, Tu, Tc)
  
  .phenology_rss(par, modelfn, bloomJDays, SeasonList, na_penalty)
  
}
