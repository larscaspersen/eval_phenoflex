#' Wrapper for Parallel Model
#' 
#' Takes hourly weather time series and parameters from Parallel model and returns day of the year with bloom. 
#' 
#' This function is mainly used by function called by the optimizer during calibration
#' 
#' @param x data.frame with at least the column "Temp" (hourly temperature data) and "JDay" (day of the year)
#' @param par vector of length 12 with the parameters of the Partial Overlap model (3), chill submodel (6) and heat submodel (3):
#' parameters from parallel model include: yc (chill requirement), zc (heat requirement), kmin (share of buds that can flower without any chill)
#' parameters of chill submodel (Dynamic Model) include: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' parameters of heat submodel (Growing Degree Hour model): Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return single numeric value with the day of the year, for which the model predicts bloom
#' with given temperature data and model parameters. 
#' 
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @references Landsberg, J. J. (1974). Apple fruit bud development and growth; analysis and an empirical model. Annals of Botany, 38(5), 1013-1023.
#' @import chillR
#' @export wrapper_parallel_model

wrapper_parallel_model <- function(x, par){
  yc = par[1]
  zc = par[2]
  kmin = par[3]
  E0 = par[4]
  E1 = par[5]
  A0 = par[6]
  A1 = par[7]
  Tf = par[8]
  slope = par[9]
  Tb = par[10]
  Tu = par[11]
  Tc = par[12]
  
  bloomindex <- parallel_model(temp = x$Temp, times = seq_along(x$Temp), 
                               yc = yc, zc = zc, kmin = kmin,
                               E0 = E0, E1 = E1, A0 = A0, A1 = A1, Tf = Tf, slope = slope,
                               Tb = Tb, Tu = Tu, Tc = Tc, 
                               basic_output = TRUE)$bloomindex
  if (bloomindex == 0) {
    return(NA)
  }
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  if (length(unique(x$Year)) == 2 & x$Year[bloomindex] == min(x$Year)) {
    JDay <- JDay - 365
  }
  n <- length(JDaylist)
  if (n == 1) {
    return(JDay)
  }
  return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
}