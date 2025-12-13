#' Wrapper for Partial Overlap Model
#' 
#' Takes hourly weather time series and parameters from Partial Overlap model and returns day of the year with bloom. 
#' 
#' This function is mainly used by function called by the optimizer during calibration
#' 
#' @param x data.frame with at least the column "Temp" (hourly temperature data) and "JDay" (day of the year)
#' @param par vector of length 14 with the parameters of the Partial Overlap model (5), chill submodel (6) and heat submodel (3):
#' yc (chill requirement), b1 (heat requirement at marginal chilling), b2 (heat requirement at optimal chilling), b3 (steepness of chill and heat compensation (low b3 = linear compensation, high b3 = much more heat needed at marginal chill), ol = extra chill accumulation after reachining minimum chill requirement, up to overlap * heat requirement.
#' parameters of chill submodel (Dynamic Model) include: E0, E1, A0, A1, Tf, slope. See \link[chillR]{Dynamic_Model} for more information.
#' parameters of heat submodel (Growing Degree Hour model): Tb, Tu, Tc. See \link[chillR]{GDH} for more information.
#' @return single numeric value with the day of the year, for which the model predicts bloom
#' with given temperature data and model parameters. 
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