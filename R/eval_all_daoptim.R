#' Evaluation funciton full model, general format
#'  
#' Takes model parameters from PhenoFlex, Sequential Model, Parallel Model or Partial Overlap Model
#' 
#' The function is parsed to the optimization algorithm during calibration. Format is compatible with
#' GenSA and DEoptim optimizers. Allows for additional parameter constraints. Also can allow
#' parameter conversion for chill submodel (Dynamic Model). Returns the squared sum of prediction errors or
#' Inf if constraints are violated.
#' 
#' @param x model parameters needs to fit the format of PhenoFlex, Parallel Model, Sequential Model or Partial Overlap Model
#' @param modelfn function used within the evaluation function to calculate the actual bloomday
#' @param bloomJDays numeric containing the days of the year with the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param na_penalty numeric, by default 365. Penalty for the phenology
#' prediction function when it fails return a bloom prediction
#' @param intermed_chill boolean, by default `FALSE`. If `TRUE` converts intermediate parameters of chill submodel (theta_star, theta_c, pie_c, tau) to E0, E1, A0, A1. 
#' Intermediate parameters have tighter search space and optimization algorithms find faster combinations.
#' @param pos_int numeric vector of length 4. Tells the position of intermediate parameters in vector `x` that should be converted 
#' if `intermed_chill=TRUE`. Assumes order: theta_star, theta_c, piec and tau. 
#' @param check_constrain booloean, by default `FALSE`. If set `TRUE` additional constraints will be evaluated. This includes
#' Q10 metric for parameters E0 and E1 should be in plausible range for biological processes (between 1.5 and 3.5) and 
#' logical constraints in GDH model (Tb <= Tu <= Tc). If constraints are violated, function does not compute prdiction error
#' but Inf
#' @param pos_E0_E1 numeric vector of length 2, indicates position of E0 and E1 in the parameter vector `x`. Only relevant when checking for constraints.
#' @param pos_heat numeric vector of length 3, indicates position of Tb, Tu and Tc in the parameter vector `x`. Only relevant when checking for constraints.
#' @param c_L numeric vector of length 5, lower tolerance level of constraints.
#' @param c_U numeric vector of length 5, upper tolerance level of constraints.
#' @return numeric, sum of squared difference between predicted and observed phenology observations.
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' 
#' @export eval_all_daoptim
eval_all_daoptim <- function(x, modelfn, bloomJDays, SeasonList, na_penalty = 365, intermed_chill = FALSE,
                     pos_int = c(5:8), check_constrain = FALSE,
                     pos_E0_E1 = c(5,6),
                     pos_heat = c(11, 4, 10),
                     c_L = c(  0,   0,   0,     1.5, 1.5), 
                     c_U = c(Inf, Inf, Inf,     3.5, 3.5)){
  if(intermed_chill){
    int <- LarsChill::convert_parameters(c(rep(0, 4), x[pos_int], rep(0,4)))
    
    if(is.list(int)) return(Inf)
    
    x[pos_int] <- int[5:8]
  }
  
  #wrapper_seq_model(x = SeasonList[[1]], par = par)
  
  par <- x
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, 
                              par = par))
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, 
                       no = pred_bloom)
  F <- sum((pred_bloom - bloomJDays)^2)
  g <- rep(0, 5)
  g[1] <- par[pos_heat[2]] - par[pos_heat[1]]
  g[2] <- par[pos_heat[3]] - par[pos_heat[1]]
  g[3] <- par[pos_heat[3]] - par[pos_heat[2]]
  g[4] <- exp((10 * par[pos_E0_E1[1]])/(297 * 279)) 
  g[5] <- exp((10 * par[pos_E0_E1[2]])/(297 * 279))
  
  if((any(g <= c_L) | any(g >= c_U)) & check_constrain){
    return(Inf)
  }
  return(F)
}