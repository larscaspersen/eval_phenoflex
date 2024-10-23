#' Evaluation function for combined fitting of cultivars, using GenSA optimizer
#'  
#' Takes model parameters for PhenoFlex, parameters of yc, zc, s1 for individual cultivars and the remaining ones are shared. 
#' Used the function for the apple study in Bonn and Asturias (Hajar Mojahid)  
#' 
#' The function evaluates a set of parameters, it assumes the chill requirement (yc), heat
#' requirement (zc) and transition parameter (s1) cultivar-specific and the remaining
#' parameters for chill and heat submodels are shared. The function is mainly used
#' in model calibration and gets called by the global optimization algorithm. 
#' 
#' @param par.guess vector, numeric. Initial guess of the model parameters 
#' yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope
#' @param modelfn function used within the evaluation function to calculate the actual bloomday, by default we
#' are using the PhenoFlex_GDHwrapper f chillR package
#' @param bloomJDays numeric containing the days of the year with the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param n_cult numeric, indicates the number of cultivars combined in the study
#' @param par_type character, controls if we use the old type (with E0, E1, A0, A1) or 
#' when equal to "new" for the format with (theta_star, theta_c, pie_c, tau)
#' @param control list, contains arguments controling the GenSA optimization algorithm.
#' For more details see in the documentation of the GenSA package.
#' @param lower vector, contains lower bound of the search space for the parameters
#' @param upper vector, contains upper bound of the search space for the parameters
#' @param seed numeric, sets seed to make the random process of optimization reproducible
#' @param ... further arguments for the modelfn function
#' @return list with two elements. First is called 'f' and contains the residual sum of squares of the model. The 
#' second is 'g' which is the values of the additional model constraints defined in the function.
#' If the flag for return_pred set TRUE, then it returns bloom dates
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom GenSA GenSA
#' @importFrom purrr map_dbl
#' @importFrom chillR PhenoFlex_GDHwrapper
#' 
#' @export phenologyFitter_combined_fitting
phenologyFitter_combined_fitting <- function (par.guess = NULL, modelfn = chillR::PhenoFlex_GDHwrapper, bloomJDays, 
                                              SeasonList, 
                                              n_cult, 
                                              par_type = 'old',
                                              control = list(smooth = FALSE, verbose = TRUE, 
                                                             maxit = 1000, nb.stop.improvement = 250), lower, upper, 
                                              seed = 1235433, ...) 
{
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
                                convert_par = convert_par)
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
      par_cult <- convert_parameters(par_cult)
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