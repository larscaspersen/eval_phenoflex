#' Predict a bloom date with the partial overlap model
#'
#' Predict the day of the year of bloom from hourly weather and partial overlap
#' model parameters. This wrapper is used by objective functions during
#' model calibration.
#'
#' @param x Hourly weather data frame with Temp (hourly temperature)
#' and JDay (day of the year) columns.
#'
#' @param par Numeric vector of length 14 with the parameters of the partial overlap model (5), chill submodel (6) and heat submodel (3):
#' yc (initial chill requirement), b1 (minimum heat requirement), b2 (additional
#' heat requirement at zero additional chill), b3 (exponential decay coefficient),
#' ol (overlap multiplier; additional chill stops affecting the heat
#' requirement at b1*ol accumulated heat).
#' Additional chill is counted from the first row reaching yc; total chill and
#' internal pools are not reset.
#' Dynamic Model chill parameters: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' GDH heat parameters: Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return A numeric predicted bloom date, expressed as a day of the year,
#' for the supplied weather and model parameters.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @references Pope, K.S., Da Silva, D., Brown, P.H., DeJong, T.M. (2014). A biologically based approach to modeling spring phenology in temperate deciduous trees. Agricultural and Forest Meteorology, 198, 15-23. https://doi.org/10.1016/j.agrformet.2014.07.009
#' @import chillR
#' @export wrapper_po_model

#helper function partial overlap model
wrapper_po_model <- function(x, par){
  yc = par[1]
  b1 = par[2]
  b2 = par[3]
  b3 = par[4]
  ol = par[5]
  E0 = par[6]
  E1 = par[7]
  A0 = par[8]
  A1 = par[9]
  Tf = par[10]
  slope = par[11]
  Tb = par[12]
  Tu = par[13]
  Tc = par[14]
  
  bloomindex <- po_model(temp = x$Temp, times = seq_along(x$Temp), 
                         yc = yc, b1 = b1, b2 = b2, b3 = b3, ol = ol,
                         E0 = E0, E1 = E1, A0 = A0, A1 = A1, Tf = Tf, slope = slope,
                         Tb = Tb, Tu = Tu, Tc = Tc, 
                         basic_output = TRUE)$bloomindex
  return_JDay(bloomindex, x$JDay, x$Year)
}
