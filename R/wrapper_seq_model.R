#' Predict a bloom date with the sequential model
#'
#' Predict the day of the year of bloom from hourly weather and sequential
#' model parameters. This wrapper is used by objective functions during
#' model calibration.
#'
#' @param x Hourly weather data frame with Temp (hourly temperature)
#' and JDay (day of the year) columns.
#'
#' @param par Numeric vector of length 11 with the parameters of the sequential model (2), chill submodel (6) and heat submodel (3):
#' Structure parameters: yc (chill requirement), zc (heat requirement).
#' Dynamic Model chill parameters: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' GDH heat parameters: Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return A numeric predicted bloom date, expressed as a day of the year,
#' for the supplied weather and model parameters.
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
  return_JDay(bloomindex, x$JDay, x$Year)
}
