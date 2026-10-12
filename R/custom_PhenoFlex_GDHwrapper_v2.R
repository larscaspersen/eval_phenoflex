#' Predict a PhenoFlex bloom date and chill and heat trajectories
#'
#' Predict the day of the year of bloom and return accumulated chill and heat
#' from hourly weather and PhenoFlex parameters for one season.
#'
#' This wrapper supports calibration functions that need the full trajectories,
#' including \code{\link{eval_phenoflex_three_stages}}.
#'
#' @param x Hourly weather data frame with Temp (hourly temperature)
#' and JDay (day of the year) columns. A Year column is required for the returned trajectory matrix.
#'
#' @param par Numeric PhenoFlex parameter vector of length 12, ordered
#' yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope. Convert characteristic chill
#' parameters (theta_star, theta_c, tau and pie_c) with
#' \code{\link{convert_parameters}} before supplying them.
#'
#' @param constraints Logical. If TRUE, check selected model parameter
#' constraints before predicting bloom. Defaults to FALSE. See Details.
#'
#' @return A list with JDay (the predicted bloom date) and chill_heat (a matrix
#' with columns Year, JDay, chill_accumulated and heat_accumulated).
#'
#' @details If constraints is TRUE, check selected Dynamic Model and GDH parameters.
#' If a check fails, return NA instead of predicting bloom. The checks require:
#' \itemize{
#'   \item GDH thresholds satisfying \code{Tu > Tb} and \code{Tc > Tb}.
#'   \item \code{exp(10 * E0 / (297 * 279))} between 1.5 and 3.5, inclusive.
#'   \item \code{exp(10 * E1 / (297 * 279))} between 1.5 and 3.5, inclusive.
#' }
#'
#' The Dynamic Model parameter checks follow Egea et al. (2021), who used
#' Q10 to constrain E0 and E1. Q10 describes the factor by which a reaction rate
#' increases with a 10-degree Celsius temperature increase. Biological reactions
#' typically have Q10 values between 2 and 3. These wrappers check both
#' limits, 1.5 and 3.5, for both E0 and E1.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @import chillR
#' @section Deprecated:
#' \code{custom_PhenoFlex_GDHwrapper_v2()} is deprecated since evalpheno 0.0.3.1.
#' Use predict_phenology() with a PhenoFlex model and basic_output = FALSE for trajectories.
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export custom_PhenoFlex_GDHwrapper_v2

custom_PhenoFlex_GDHwrapper_v2 <- function (x, par, constraints = FALSE){
  .Deprecated(new = "predict_phenology", package = "evalpheno",
              msg = "custom_PhenoFlex_GDHwrapper_v2() is deprecated. Use predict_phenology() with a PhenoFlex model and basic_output = FALSE for trajectories. See help('evalpheno-deprecated') for migration differences.")
  #x is one of the elements in season List
  #par are the parameters of the model
  
  
  #x <- SeasonList[[1]]
  
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
    t4 <- exp((10 * E1)/(297 * 279)) < 1.5 | exp((10 * E1)/(297 * 279)) > 3.5
    
    if(any(c(t1, t2, t3, t4))){
      return(NA)
    }
    
  }
  
  #calculat the bloom day
  pheno_out <- chillR::PhenoFlex(temp = x$Temp, 
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
                                 stopatzc = FALSE,
                                 basic_output = FALSE)
  
  #when I do not stop the accumulation at zc, then the returned bloomdate is not when the flowering is predicted
  #(hopefully)
  if(any(pheno_out$z >= zc)){
    bloomindex <- min(which(pheno_out$z >= zc))
  } else {
    bloomindex <- 0
  }
  
  chill_heat_out <- matrix(c(x$Year, x$JDay, pheno_out$y, pheno_out$z),ncol = 4)
  colnames(chill_heat_out) <- c('Year', 'JDay', 'chill_accumulated', 'heat_accumulated')
  
  #make the bloomindex a fraction of the Jday, depending on when exactly flower was predicted
  
  #return values
  list(JDay = return_JDay(bloomindex, x$JDay, x$Year),
       chill_heat = chill_heat_out)
}
