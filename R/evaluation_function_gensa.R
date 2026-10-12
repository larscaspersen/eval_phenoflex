#' Evaluate PhenoFlex parameters for the GenSA optimizer
#'
#' Predict bloom dates for the supplied seasons, compare them with observed
#' dates and return a residual sum of squares for GenSA to minimize.
#'
#' Parameters may use the kinetic chill representation (E0, E1, A0 and A1)
#' or the characteristic representation (theta_star, theta_c, tau and pie_c).
#' Characteristic parameters are converted before bloom dates are predicted.
#'
#' @param x Numeric PhenoFlex parameter vector. With intermed_par = TRUE,
#' use yc, zc, s1, Tu, theta_star, theta_c, tau, pie_c, Tf, Tc, Tb, slope.
#' With intermed_par = FALSE, use yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope.
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
#' @param intermed_par Logical. If TRUE (the default), convert characteristic
#' chill parameters before prediction. If FALSE, use kinetic parameters directly.
#'
#' @param q10_lower Lower allowed Q10 value for E0 and E1. Defaults to 1.5.
#' Candidates outside the allowed range return q10_penalty regardless of prediction accuracy.
#'
#' @param q10_upper Upper allowed Q10 value for E0 and E1. Defaults to 3.5.
#' Candidates outside the allowed range return q10_penalty regardless of prediction accuracy.
#'
#' @param q10_penalty Objective value returned when chill parameter conversion
#' fails or a Q10 value is outside the allowed range. Defaults to 10^10.
#'
#' @param na_penalty Numeric replacement for a missing bloom prediction. Defaults
#' to 365; this value is used in place of the missing date when calculating loss.
#'
#' @return Numeric residual sum of squares between predicted and observed bloom
#' dates, or q10_penalty if conversion fails or a Q10 constraint is violated.
#'
#' @author Lars Caspersen
#' @keywords utility
#' @importFrom nleqslv nleqslv
#'
#' @section Deprecated:
#' \code{evaluation_function_gensa()} is deprecated since evalpheno 0.0.3.1.
#' Use fit_phenology(optimizer = 'GenSA') with phenology_rss().
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export evaluation_function_gensa

evaluation_function_gensa <- function (x, 
                               modelfn, 
                               bloomJDays, 
                               SeasonList, 
                               intermed_par = TRUE,
                               q10_lower = 1.5, 
                               q10_upper = 3.5,
                               q10_penalty = 10^10,
                               na_penalty = 365) 
{
  .Deprecated(new = "phenology_rss", package = "evalpheno",
              msg = "evaluation_function_gensa() is deprecated. Use fit_phenology(optimizer = 'GenSA') with phenology_rss(). See help('evalpheno-deprecated') for migration differences.")
  par <- x
  
  if(intermed_par){
    params <- numeric(4)
    params[1] <- x[5]
    params[2] <- x[6]
    params[3] <- x[7]
    params[4] <- x[8]
    
    #convert parameters
    output <- .fit_chill_parameters(params)
    if (output$termcd >= 3) {
      E0 <- NA
      E1 <- NA
      A0 <- NA
      A1 <- NA
      return(q10_penalty)
    }
    else {
      E0 <- output$x[1]
      E1 <- output$x[2]
      q = 1/params[1] - 1/params[2]
      A1 <- output$A1
      A0 <- output$A0
    }
    par[5:8] <- c(E0, E1, A0, A1)
  }

  
  #check for q10 criterion, if failed then penalize
  q10_e0 <- exp((10 * par[5])/(297 * 279))
  q10_e1 <- exp((10 * par[6])/(297 * 279))
  if(q10_e0 < q10_lower | q10_e0 > q10_upper | q10_e1 < q10_lower | q10_e1 > q10_upper){
    return(q10_penalty)
  }
  
  #predict bloom
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, 
                              par = par))
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, 
                       no = pred_bloom)
  return(sum((pred_bloom - bloomJDays)^2))
}
