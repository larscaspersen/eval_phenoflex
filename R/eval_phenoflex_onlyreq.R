#used in the lake constance study and the almond study

eval_phenoflex_onlyreq  <- function(x,
                                    modelfn,
                                    bloomJDays,
                                    SeasonList,
                                    par_fixed,
                                    constrains = TRUE,
                                    return_pred = FALSE,
                                    na_penalty = 365){
  
  
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