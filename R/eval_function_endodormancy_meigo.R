#' Evaluate Dynamic Model parameters using endodormancy observations
#'
#' Calculate accumulated chill portions for each supplied temperature series,
#' then measure the dispersion of these totals across seasons within each
#' cultivar. Calibration minimizes the sum of the cultivar-specific dispersion
#' values. Characteristic parameters are converted to kinetic Dynamic Model
#' parameters before chill accumulation is calculated.
#'
#' Supply a nested list of seasons, with one outer entry per cultivar.
#'
#' @param x Numeric characteristic Dynamic Model parameter vector in the order
#' theta_star, theta_c, tau, pie_c, Tf, slope. Characteristic parameters are
#' converted to E0, E1, A0 and A1 before calculating chill accumulation.
#'
#' @param SeasonList Nested list with one season list per cultivar. Each season
#' contains an hourly temperature series in a Temp column for the observation
#' period used to assess endodormancy release.
#'
#' @param dispersion_fun Name of a function that measures dispersion of accumulated
#' chill portions across seasons within a cultivar. The function must be available
#' by name. Defaults to "calc_cv", which calculates the coefficient of variation.
#'
#' @return A list with F (the sum of cultivar-specific dispersion values)
#' and g (the Q10 constraint values).
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @importFrom chillR Dynamic_Model
#' @importFrom utils tail
#' @examples
#' \dontrun{
#' #          theta_star, theta_c, pie_c, tau, Tf, slope
#' par <-   c(279,        281,      ?      ?    4,  1.6)
#'
#' #prepare weather data
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)
#'
#' #get endodormancy / ecodormancy release date
#' #--> prepare SeasonList for different cultivars
#'
#' eval_function_endodormancy_meigo(x = par, par = SeasonList)
#'
#'
#' }
#' @export eval_function_endodormancy_meigo

eval_function_endodormancy_meigo <- function(x, SeasonList, dispersion_fun = 'calc_cv'){
  
  #convert the parameters from new to old
  x_new <- convert_parameters(c(0,0,0,0, x[1:5], 0, 0, x[6]))
  
  #inequality constrains
  g <- rep(0,2)
  
  #in case the conversion fails, return simply a large error for the evaluation function
  if(is.list(x_new)){
    return(list(F=x_new$F, g=rep(x_new$F, length(g))))
  }
  
  x_new <- x_new[c(5:9, 12)]
  
  #iterate over seasonlist of cultivars
  coefficinets_of_variation <- purrr::map_dbl(SeasonList, function(pheno){
    #iterate over seasons of one cultivar and sum gdh
    gdh <- purrr::map_dbl(pheno, function(p){
      dynamic_model_wrapper(p$Temp, x_new)
    })
    #calculate the dispersion function
      return( get(dispersion_fun)(gdh))
  })
  
  F <- sum(coefficinets_of_variation)
  
  
  #E0 and E1
  g[1] <- exp((10 * x_new[1]) / (297 * 279))
  g[2] <- exp((10 * x_new[2]) / (297 * 279))
  
  
  
  #check the paraemeters
  
  return(list(F=F, g=g))
}

dynamic_model_wrapper <- function(Temp, par, summ = TRUE){
  
  # par <- c(4153.5, 12888.8, 139500, 2.567e+18, 1.6, 4)
  # Temp <- seasonlist_in[[1]][[1]]$Temp
  # summ <-  TRUE
  out <- chillR::Dynamic_Model(HourTemp = Temp, 
                        summ = summ, 
                        E0 = par[1], 
                        E1 = par[2], 
                        A0 = par[3], 
                        A1 = par[4], 
                        slope = par[5], 
                        Tf = par[6] + 273) 

    return(utils::tail(out, n= 1))
}

