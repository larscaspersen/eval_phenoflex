#' Evaluate PhenoFlex parameters for the MEIGO optimizer
#'
#' Predict bloom dates from hourly weather and model parameters, then return
#' the residual sum of squares and parameter constraint values for calibration.
#'
#' @param x traditional model parameters of PhenoFlex in the order yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope.
#' @param modelfn Function that predicts a bloom date from a season's weather
#' and model parameters, such as \code{\link{custom_PhenoFlex_GDHwrapper}}.
#'
#' @param bloomJDays Numeric vector of observed bloom dates, expressed as days of
#' the year and ordered to match SeasonList.
#'
#' @param SeasonList List of hourly weather data frames, one per phenological season.
#' Each data frame contains Temp (hourly temperature) and JDay (day of the year).
#' Season lists can be generated with \code{\link[chillR]{genSeasonList}}.
#'
#' @param na_penalty Numeric replacement for a missing bloom prediction. Defaults
#' to 365; this value is used in place of the missing date when calculating loss.
#'
#' @return A list with F (the residual sum of squares) and g (the model's
#' constraint values).
#'
#' @author Lars Caspersen
#' @keywords utility
#'
#' @section Deprecated:
#' \code{evaluation_function_meigo()} is deprecated since evalpheno 0.0.3.1.
#' Use pheno_model(), phenology_rss() and fit_phenology().
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export evaluation_function_meigo

evaluation_function_meigo <- function(x, 
                                      modelfn,
                                      bloomJDays,
                                      SeasonList,
                                      na_penalty = 365){
  .Deprecated(new = "phenology_rss", package = "evalpheno",
              msg = "evaluation_function_meigo() is deprecated. Use pheno_model(), phenology_rss() and fit_phenology(). See help('evalpheno-deprecated') for migration differences.")
  
  #innput:
  #         x is the parameters in meigo
  #         modelfn is the function used to calculate the bloomdays
  #         SeasonList contains the weather data for the different years
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  #output: inequality constraints g
  #        model performance value F
  
  #change the name of the parameters, dont know if necessary
  par <- x
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  #calculate the model performance value
  F <- sum((pred_bloom - bloomJDays)^2)
  
  
  
  #####
  #inequality constraints
  #####
  
  #this is the vector containing the values for the inequality constraints
  #at first initialize the vector
  g <- rep(0,5)
  
  
  #equality constraints should be always stated before inequality constraints, according to meigo vignette
  #(but we dont have any in this case)
  
  
  #inequality constraints are reformulated as differences
  
  #Tu >= Tb
  g[1] <- x[4] - x[11]
  #Tx >= Tb
  g[2] <- x[10] - x[11]
  #Tc >= Tu
  g[3] <- x[10] - x[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * x[5]) / (297 * 279))
  
  g[5] <- exp((10 * x[6]) / (297 * 279))
  
  
  return(list(F=F, g=g))

}
