#' Returns Julian Day (with fraction)
#' 
#' Function needed to handle the calibration function for several phenology stages
#' 
#' @param index integer, when heat requirement is met 
#' @param Jday_vec vector of julian days
#' @param year_vec vector with the year
#' @return day of the year, with a fraction to account for when (at which hour)
#' the requirement is met. 12:00 midday corresponds to .0, in the afternoon up to +0.5 is added to the julian day, 
#' if it is before midday up to -0.5 is added to the Julian Day
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom assertthat are_equal
#' @importFrom nleqslv nleqslv
#' @examples 
#' \dontrun{
#'  
#'  i <- 56
#'  Jday_vec <- rep(1:30, each = 24)
#'  year_vec <- rep(2024, length(Jday_vec))
#'  return_JDay(i, return_JDay, year_vec)
#' }
#' @export return_JDay
return_JDay <- function(index, Jday_vec, year_vec){
  
  # index <- index_budburst
  # Jday_vec <- chill_heat[,2]
  # year_vec <- chill_heat[,1]
  
  if (index == 0){
    JDay_out <- NA
  } else {
    JDay <- Jday_vec[index]
    JDaylist <- which(Jday_vec == JDay)
    
    #if we are in the norhtern hemisphere and the year corresponding to the index is the smaller one, return negative numbers relative to Jan-1 being 1
    if(length(unique(year_vec)) == 2 & year_vec[index] == min(year_vec)){
      
      JDay <- JDay - 365
      
    } 
    
    n <- length(JDaylist)
    if (n == 1){
      JDay_out <- JDay
      return(JDay)
    } else {
      JDay_out <- JDay + which(JDaylist == index)/n - 1/(n/ceiling(n/2))
    }
  }
  return(JDay_out)
}