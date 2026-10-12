#' Evaluate PhenoFlex structure parameters with fixed submodels
#'
#' Evaluate the chill requirement (yc), heat requirement (zc) and transition
#' parameter (s1) while keeping all chill and heat submodel parameters fixed.
#' The function returns an objective value and parameter constraint values
#' for use during calibration.
#'
#' @param x Numeric PhenoFlex structure parameter vector in the order
#' yc (chill requirement), zc (heat requirement), s1 (transition parameter).
#'
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
#' @param par_fixed Fixed chill and heat submodel parameters in the order
#' Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope.
#'
#' @param constrains Logical compatibility argument, defaulting to TRUE.
#' The current implementation returns constraint values regardless of this setting.
#'
#' @param return_pred Logical. If TRUE, return predicted bloom dates instead of
#' the objective value and constraint values. Defaults to FALSE.
#'
#' @param na_penalty Numeric replacement for a missing bloom prediction. Defaults
#' to 365; this value is used in place of the missing date when calculating loss.
#'
#' @return A list with F (the residual sum of squares) and g (the model's
#' constraint values). If return_pred is TRUE, return predicted bloom dates instead.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#'
#' @section Deprecated:
#' \code{eval_phenoflex_onlyreq()} is deprecated since evalpheno 0.3.1.
#' Use pheno_model() and phenology_rss(); bound only the structure parameters in fit_phenology().
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export eval_phenoflex_onlyreq
eval_phenoflex_onlyreq  <- function(x,
                                    modelfn,
                                    bloomJDays,
                                    SeasonList,
                                    par_fixed,
                                    constrains = TRUE,
                                    return_pred = FALSE,
                                    na_penalty = 365){
  .Deprecated(new = "phenology_rss", package = "evalpheno",
              msg = "eval_phenoflex_onlyreq() is deprecated. Use pheno_model() and phenology_rss(); bound only the structure parameters in fit_phenology(). See help('evalpheno-deprecated') for migration differences.")
  
  
  par <- c(x, par_fixed)
  pred_bloom <- NULL
  
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
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
  g[1] <- par[length(par)-8] - par[length(par)-1]
  #Tx >= Tb
  g[2] <- par[length(par)-2] - par[length(par)-1]
  #Tc >= Tu
  g[3] <- par[length(par)-2] - par[length(par)-8]
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * par[length(par)-7]) / (297 * 279))
  g[5] <- exp((10 * par[length(par)-6]) / (297 * 279))
  if(return_pred){
    pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    return(pred_bloom)
  } else{
    return(list(F=F, g=g))
  }
}
