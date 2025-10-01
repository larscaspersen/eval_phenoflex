#' Evaluation function for endodormancy 
#' 
#' Takes model parameters for Dynamic Model and list of seasons with daily temperature time series. 
#' 
#' The function evaluates a set of parameters and calculates accumulated chill portions (CP) for each of the supplied temperature time series.
#' Then it returns a dispersion metric measuring how different accumulated CP are among the seasons.
#' The function is usually used in an optimization problem, where the goal is to estimate
#' Dynamic Model parameters most suitable to minimize the dispersion of accumulated CP.
#' 
#' The function has very specific requirements on the structure of the Seasonlist.
#' In fact it is a list of Seasonlists created for several cultivars.
#' For each cultivar, the seasonlist contains exactly the period from endodormancy 
#' release to observed ecodormancy release.  
#' 
#' 
#' @param x model parameters of Dynamic Model, new format: theta_star, theta_c, pie_c, tau, Tf, slope
#' will get converted to the format: E0, E1, A0, A1, Tf, slope to run the Dynamic Model 
#' @param SeasonList List of SeasonLists. For each cultivar one entry. Individual seasonlists 
#' per cultivar contain exactly period of: measured end of endodormancy release to measured
#' ecodormancy release.
#' @param dispersion_fun name of function to measure dispersion of accumulated GDH per
#' cultivar. Needs to be a function loaded in the environment. Can also be a custom function.
#' By default "calc_cv", which calculates the coefficient of variation.
#' @return sum of the dispersion function. This will be subject to global optimization. 
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @details
#' Used this function for the Olive study of Efe Deger
#' 
#' @importFrom purrr map_dbl
#' @importFrom chillR Dynamic_Model
#' @importFrom LarsChill convert_parameters
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
  x_new <- LarsChill::convert_parameters(c(0,0,0,0, x[1:5], 0, 0, x[6]))
  
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

