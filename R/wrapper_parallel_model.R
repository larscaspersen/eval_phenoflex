#' Wrapper for Parallel Model
#' 
#' Takes hourly weather time series and parameters from Parallel model and returns day of the year with bloom. 
#' 
#' This function is mainly used by function called by the optimizer during calibration
#' 
#' @param x data.frame with at least the column "Temp" (hourly temperature data) and "JDay" (day of the year)
#' @param par vector of length 12 with the parameters of the Parallel model (3), chill submodel (6) and heat submodel (3):
#' parameters from parallel model include: yc (chill requirement), zc (heat requirement), kmin (minimum relative heat effectiveness at zero chill)
#' parameters of chill submodel (Dynamic Model) include: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' parameters of heat submodel (Growing Degree Hour model): Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return single numeric value with the day of the year, for which the model predicts bloom
#' with given temperature data and model parameters. 
#' 
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @references Hänninen, H. and Kramer, K. (2007). A framework for modelling the annual cycle of trees in boreal and temperate regions. Silva Fennica 41(1), 167-205. Equation B4b. doi:10.14214/sf.313.
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
  return_JDay(bloomindex, x$JDay, x$Year)
}
