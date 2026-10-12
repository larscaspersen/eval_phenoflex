#' Predict a bloom date with PhenoFlex
#'
#' Predict the day of the year of bloom from hourly weather and PhenoFlex
#' parameters for one season.
#'
#' This wrapper is used by calibration objective functions such as
#' \code{\link{evaluation_function_meigo}},
#' \code{\link{evaluation_function_meigo_nonlinear}} and
#' \code{\link{evaluation_function_meigo_vns}}.
#'
#' @param x Hourly weather data frame with Temp (hourly temperature)
#' and JDay (day of the year) columns.
#'
#' @param par Numeric PhenoFlex parameter vector of length 12, ordered
#' yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope. Convert characteristic chill
#' parameters (theta_star, theta_c, tau and pie_c) with
#' \code{\link{convert_parameters}} before supplying them.
#'
#' @param constraints Logical. If TRUE, check selected model parameter
#' constraints before predicting bloom. Defaults to FALSE. See Details.
#'
#' @return A numeric predicted bloom date, expressed as a day of the year.
#'
#' @details If constraints is TRUE, check selected Dynamic Model and GDH parameters.
#' If a check fails, return NA instead of predicting bloom. The checks require:
#' \itemize{
#'   \item GDH thresholds satisfying \code{Tu > Tb} and \code{Tc > Tb}.
#'   \item \code{exp(10 * E0 / (297 * 279))} between 1.5 and 3.5, inclusive.
#'   \item \code{exp(10 * E1 / (297 * 279))} at least 1.5.
#' }
#'
#' The Dynamic Model parameter checks follow Egea et al. (2021), who used
#' Q10 to constrain E0 and E1. Q10 describes the factor by which a reaction rate
#' increases with a 10-degree Celsius temperature increase. Biological reactions
#' typically have Q10 values between 2 and 3. These legacy wrappers check both
#' limits for E0 and only the lower limit for E1.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @import chillR
#' @export custom_PhenoFlex_GDHwrapper
custom_PhenoFlex_GDHwrapper <- function (x, par, constraints = FALSE){
  #x is one of the elements in season List
  #par are the parameters of the model
  
  #make explicit what which parameters is
  yc = par[1]
  zc = par[2]
  s1 = par[3]
  Tu = par[4]
  E0 = par[5]
  E1 = par[6]
  A0 = par[7]
  A1 = par[8]
  Tf = par[9]
  Tc = par[10] 
  Tb = par[11]
  slope = par[12]
  
  
  #in case the parameters do not make sense, return NA
  if(constraints){
    
    t1 <- Tu <= Tb
    t2 <- Tc <= Tb
    t3 <- exp((10 * E0)/(297 * 279)) < 1.5 | exp((10 * E0)/(297 * 279)) > 3.5
    t4 <- exp((10 * E1)/(297 * 279)) < 1.5 | exp((10 * E0)/(297 * 279)) > 3.5
    
    if(any(c(t1, t2, t3, t4))){
      return(NA)
    }
    
  }
  
  #calculat the bloom day
  bloomindex <- chillR::PhenoFlex(temp = x$Temp, 
                          times = seq_along(x$Temp), 
                          yc = yc, 
                          zc = zc, 
                          s1 = s1, 
                          Tu = Tu, 
                          E0 = E0, 
                          E1 = E1, 
                          A0 = A0, 
                          A1 = A1, 
                          Tf = Tf, 
                          Tc = Tc, 
                          Tb = Tb, 
                          slope = slope, 
                          Imodel = 0L, 
                          basic_output = TRUE)$bloomindex
  
  

  #return values
  return_JDay(bloomindex, x$JDay, x$Year)
}
