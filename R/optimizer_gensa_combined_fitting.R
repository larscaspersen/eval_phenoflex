#' Calibrate a combined cultivar model with GenSA
#'
#' Fit cultivar-specific chill requirements (yc), heat requirements (zc) and
#' transition parameters (s1), with shared chill and heat submodel parameters.
#' The GenSA optimizer minimizes the total residual sum of squares across
#' cultivars and seasons.
#'
#' @param par.guess Numeric initial parameter vector containing one yc value per
#' cultivar, followed by one zc and one s1 value per cultivar, then shared parameters
#' Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope for par_type = "old".
#' For par_type = "new", replace E0, E1, A0 and A1 with theta_star, theta_c,
#' tau and pie_c.
#'
#' @param modelfn Function that predicts a bloom date from weather and parameters.
#' Defaults to \code{\link[chillR]{PhenoFlex_GDHwrapper}}.
#'
#' @param bloomJDays List of numeric observed bloom-date vectors, one per cultivar.
#' Dates must match the order of seasons in SeasonList.
#'
#' @param SeasonList Nested list of hourly weather data frames, with one season
#' list per cultivar in the same order as bloomJDays.
#'
#' @param n_cult Number of cultivars to calibrate together.
#'
#' @param par_type Character string selecting the chill parameter representation:
#' "old" for kinetic parameters or "new" for characteristic parameters.
#'
#' @param control List of GenSA optimizer controls. See \code{\link[GenSA]{GenSA}}
#' for available settings.
#'
#' @param lower Numeric vector of lower parameter bounds.
#'
#' @param upper Numeric vector of upper parameter bounds.
#'
#' @param seed Numeric seed for reproducible optimization.
#'
#' @param ... Additional arguments accepted by the fitting interface.
#'
#' @return A phenologyFit result containing model_fit (the native GenSA result),
#' par (the fitted parameters), initial settings, weather and observed bloom dates.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom GenSA GenSA
#' @importFrom purrr map_dbl
#' @importFrom chillR PhenoFlex_GDHwrapper
#'
#' @section Deprecated:
#' \code{phenologyFitter_combined_fitting()} is deprecated since evalpheno 0.0.3.1.
#' Use a model collection with fit_phenology(optimizer = 'GenSA') and phenology_rss_combined().
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export phenologyFitter_combined_fitting
phenologyFitter_combined_fitting <- function (par.guess = NULL, modelfn = chillR::PhenoFlex_GDHwrapper, bloomJDays, 
                                              SeasonList, 
                                              n_cult, 
                                              par_type = 'old',
                                              control = list(smooth = FALSE, verbose = TRUE, 
                                                             maxit = 1000, nb.stop.improvement = 250), lower, upper, 
                                              seed = 1235433, ...) 
{
  .Deprecated(new = "fit_phenology", package = "evalpheno",
              msg = "phenologyFitter_combined_fitting() is deprecated. Use a model collection with fit_phenology(optimizer = 'GenSA') and phenology_rss_combined(). See help('evalpheno-deprecated') for migration differences.")
  control$seed <- seed
  stopifnot(is.list(SeasonList))
  stopifnot(length(SeasonList) == length(bloomJDays))
  res <- phenologyFit()
  res$par.guess <- par.guess
  res$modelfn <- modelfn
  #bloomJdays is now a list of vectors, one for each cultivar
  res$bloomJDays <- bloomJDays
  res$SeasonList <- SeasonList
  res$lower <- lower
  res$upper <- upper
  res$control <- control
  
  if(par_type == 'old'){
    convert_par <- FALSE
  } else if(par_type == 'new'){
    convert_par <- TRUE
  }
  
  #seasonlist needs to be a list of season listst
  for(j in 1:length(SeasonList)){
    for (i in c(1:length(SeasonList[[j]]))) {
      minJDay <- SeasonList[[j]][[i]]$JDay[1]
      maxJDay <- SeasonList[[j]][[i]]$JDay[length(SeasonList[[j]][[i]]$JDay)]
      if (maxJDay > minJDay) {
        stop(paste0("Season ", i, " is overlapping with the previous or following one. Aborting!"))
      }
      if (bloomJDays[[j]][i] > maxJDay && bloomJDays[[j]][i] < minJDay) {
        stop(paste0("In season ", i, " the bloomJDay is outside the provided JDay vector. Aborting!"))
      }
      dx <- diff(SeasonList[[j]][[i]]$JDay)
      mx <- min(dx)
      kmx <- which(dx == mx)
      SeasonList[[j]][[i]]$JDay[1:kmx] <- SeasonList[[j]][[i]]$JDay[1:kmx] + 
        mx - 1
      if (bloomJDays[[j]][i] > minJDay) 
        bloomJDays[[j]][i] <- bloomJDays[[j]][i] + mx - 1
      res$SeasonList[[j]][[i]]$JDayunwrapped <- SeasonList[[j]][[i]]$JDay
      res$bloomJDaysunwrapped <- bloomJDays
    }
  }
  
  
  res$model_fit <- GenSA::GenSA(par = par.guess, fn = chifull_combined_fitting, 
                                bloomJDays = bloomJDays, SeasonList = SeasonList, modelfn = modelfn, 
                                control = control, lower = lower, upper = upper, n_cult = n_cult,
                                convert_par = convert_par, ...)
  res$par <- res$model_fit$par
  
  
  #res$pbloomJDays <- predictBloomDays(par = res$par, SeasonList = res$SeasonList, 
  #                                    modelfn = res$modelfn)
  return(res)
}

chifull_combined_fitting <- function (par,
                                      modelfn, 
                                      bloomJDays, 
                                      SeasonList, 
                                      na_penalty = 365, 
                                      n_cult,
                                      convert_par = FALSE,
                                      ...) 
{
  rss <- purrr::map_dbl(1:length(SeasonList), function(i){
    
    #extract the parameters
    par_cult <- par[c(i, 1+((n_cult -1))+i, 2+((n_cult -1)*2)+i, (length(par)-8):length(par))]
    
    #convert parameters new to old if needed
    if(convert_par){
      par_cult <- characteristic_to_kinetic(par_cult)
      if (is.list(par_cult)) return(Inf)
    }
    
    #check if parameters violate any 
    
    #predict the bloom
    pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    #calculate for each cultivar the rss
    rss <- sum((pred_bloom - bloomJDays[[i]])^2)
    return(rss)
    
  })
  
  return(sum(rss))
}
