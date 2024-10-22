#' Evaluation function for ecodormancy 
#' 
#' Takes model parameters for GDH model and list of seasons with daily temperature time series. 
#' 
#' The function evaluates a set of parameters and calculates accumulated GDH for each of the supplied temperature time series.
#' Then it returns a dispersion metric measuring how different accumulated GDH are among the seasons.
#' The function is usually used in an optimization probelm, where the goal is to estimate
#' GDH model parameters most suitable to minimize the dispersion of accumulated GDH.
#' 
#' The function has very specific requirements on the structure of the Seasonlist.
#' In fact it is a list of Seasonlists created for several cultivars.
#' For each cultivar, the seasonlist contains exactly the period from endodormancy 
#' release to observed ecodormancy release.  
#' 
#' 
#' @param x model parameters of GDH model, order: Tb, Tu, Tc
#' @param SeasonList List of SeasonLists. For each cultivar one entry. Individual seasonlists 
#' per cultivar contain exactly period of: measured end of endodormancy release to measured
#' ecodormancy release.
#' @param dispersion_fun name of function to measure dispersion of accumulated GDH per
#' cultivar. Needs to be a function loaded in the environment. Can also be a custom function.
#' By default "calc_cv", which calculates the coefficient of variation.
#' @return sum of the dispersion function. This will be subject to global optimization. 
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @examples 
#' \dontrun{
#' #          yc,  zc,  s1, Tu,    E0,      E1,     A0,         A1,   Tf, Tc, Tb,  slope
#' par <-   c(40, 190, 0.5, 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
#' 
#' #prepare weather data
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)
#' 
#' custom_PhenoFlex_GDHwrapper(x = seasonList[[1]], par = par)
#' 
#' 
#' }
#' @export eval_function_gdh

eval_function_gdh <- function(x, SeasonList, dispersion_fun = 'calc_cv'){
  #iterate over seasonlist of cultivars
  coefficinets_of_variation <- purrr::map_dbl(SeasonList, function(pheno){
    #iterate over seasons of one cultivar and sum gdh
    gdh <- purrr::map_dbl(pheno, function(p){
      GDH_wrapper(p$Temp, x)
    })
    #calculate the dispersion function
    get(dispersion_fun)(gdh) %>%
      return()
  })
  return(sum(coefficinets_of_variation))
}


#wrapper to optimize the GDH in Efes Olive paper
GDH_wrapper <- function(Temp, par, summ = TRUE){
  Stress <- 1
  Tb <- par[1]
  Tu <- par[2]
  Tc <- par[3]
  GDH_weight <- rep(0, length(Temp))
  GDH_weight[which(Temp >= Tb & Temp <= Tu)] <- Stress *
    (Tu - Tb)/2 * (1 + cos(pi + pi * (Temp[which(Temp >=
                                                   Tb & Temp <= Tu)] - Tb)/(Tu - Tb)))
  GDH_weight[which(Temp > Tu & Temp <= Tc)] <- Stress *
    (Tu - Tb) * (1 + cos(pi/2 + pi/2 * (Temp[which(Temp >
                                                     Tu & Temp <= Tc)] - Tu)/(Tc - Tu)))
  if (summ)
    return(sum(GDH_weight))
  else return(GDH_weight)
}

calc_cv <- function(x){
  sd(x) / mean(x) * 100 %>%
    return()
}