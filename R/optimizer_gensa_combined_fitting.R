#wrapper for the GenSA optimization of combined fitting

phenologyFitter_combined_fitting <- function (par.guess = NULL, modelfn = PhenoFlex_GDHwrapper, bloomJDays, 
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