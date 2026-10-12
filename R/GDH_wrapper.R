#' Evaluate GDH parameters using ecodormancy release observations
#'
#' Calculate accumulated growing degree hours (GDH) for each supplied
#' temperature series, then measure the dispersion of these totals across
#' seasons within each cultivar. Calibration minimizes the sum of the
#' cultivar-specific dispersion values.
#'
#' Supply a nested list of seasons, with one outer entry per cultivar.
#' Each temperature series should cover the period from observed endodormancy
#' release to observed ecodormancy release.
#'
#' @param x Numeric GDH parameter vector in the order Tb, Tu, Tc.
#'
#' @param SeasonList Nested list with one season list per cultivar. Each hourly
#' temperature series covers the period from observed endodormancy release to
#' observed ecodormancy release.
#'
#' @param dispersion_fun Name of a function that measures dispersion of accumulated
#' GDH across seasons within a cultivar. The function must be available by name.
#' Defaults to "calc_cv", which calculates the coefficient of variation.
#'
#' @return Numeric sum of the cultivar-specific dispersion values to minimize.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @examples
#' \dontrun{
#' #          Tb, Tu, Tc
#' par <-   c(4,  26, 36)
#'
#' #prepare weather data
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)
#'
#' #get endodormancy / ecodormancy release date
#' #--> prepare SeasonList for different cultivars
#'
#' eval_function_gdh(x = par, par = SeasonList)
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
      return(get(dispersion_fun)(gdh))
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
  stats::sd(x) / mean(x) * 100 %>%
    return()
}