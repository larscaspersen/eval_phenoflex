#' Evaluation function for combined fitting of cultivars
#'  
#' Takes model parameters for PhenoFlex, parameters of yc, zc, s1. Other parameters are fixed.
#' Used this function 
#' 
#' The function evaluates a set of parameters, it assumes the chill requirement (yc), heat
#' requirement (zc) and transition parameter (s1) cultivar-specific and the remaining
#' parameters for chill and heat submodels are shared. The function is mainly used
#' in model calibration and gets called by the global optimization algorithm. 
#' 
#' @param x model parameters of PhenoFlex, with chill requirement (yc), 
#' heat requirement (zc), and transition parameter (s1) new format: 
#' @param modelfn function used within the evaluation function to calculate the actual bloomday, often we use
#' the 'custom_GDH_wrapper' function for that
#' @param bloomJDays numeric containing the days of the year with the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param par_fixed parameters of chill and heat submodel, order: Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope 
#' @param constraints logical, by default set TRUE. Indicates if there are further
#' checks on the parameters except the boundaries of the optimization algorithm
#' @param return_pred logical, by default FALSE. If set TRUE the function will return
#' bloom dates instead of the output required by the global optimzation algorithm.
#' @param na_penalty numeric, by default 365. Penalty for the phenology
#' prediction function when it fails return a bloom prediction
#' @return list with two elements. First is called 'f' and contains the residual sum of squares of the model. The 
#' second is 'g' which is the values of the additional model constraints defined in the function.
#' If the flag for return_pred set TRUE, then it returns bloom dates
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' 
#' @export eval_phenoflex_onlyreq
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