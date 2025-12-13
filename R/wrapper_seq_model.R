#' Wrapper for Sequential Model
#' 
#' Takes hourly weather time series and parameters from Sequential model and returns day of the year with bloom. 
#' 
#' This function is mainly used by function called by the optimizer during calibration
#' 
#' @param x data.frame with at least the column "Temp" (hourly temperature data) and "JDay" (day of the year)
#' @param par vector of length 11 with the parameters of the Sequential model (2), chill submodel (6) and heat submodel (3):
#' parameter of sequenital model include: yc (chill requirement), zc (heat requirement).
#' parameters of chill submodel (Dynamic Model) include: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' parameters of heat submodel (Growing Degree Hour model): Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return single numeric value with the day of the year, for which the model predicts bloom
#' with given temperature data and model parameters. 
#' 
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @references Ashcroft, G. L., Richardson, E. A., & Seeley, S. D. (1977). A statistical method of determining chill unit and growing degree hour requirements for deciduous fruit Trees1. HortScience, 12(4), 347-348.
#' @import chillR
#' @export wrapper_seq_model

wrapper_seq_model <- function(x, par){
  yc = par[1]
  zc = par[2]
  E0 = par[3]
  E1 = par[4]
  A0 = par[5]
  A1 = par[6]
  Tf = par[7]
  slope = par[8]
  Tb = par[9]
  Tu = par[10]
  Tc = par[11]
  
  bloomindex <- seq_model(temp = x$Temp, times = seq_along(x$Temp), 
                          yc = yc, zc = zc, 
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