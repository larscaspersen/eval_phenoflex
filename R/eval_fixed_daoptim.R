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
#' Assumes only model specific parameters, excluding parameters of the submodel for chill and heat accumulation.
#' @param modelfn function used within the evaluation function to calculate the actual bloomday
#' @param bloomJDays numeric containing the days of the year with the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param na_penalty numeric, by default 365. Penalty for the phenology
#' prediction function when it fails return a bloom prediction
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature (in degree Kelvin) for the sigmoidal function
#' @param Tb numeric. GDH base temperature (lower threshold) 
#' @param Tu numeric. GDH optimal temperature 
#' @param Tc numeric. GDH upper temperature (upper threshold)
#' @return numeric, sum of squared difference between predicted and observed phenology observations.
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' 
#' @export eval_fixed_daoptim
eval_fixed_daoptim <- function(x, modelfn, bloomJDays, SeasonList, na_penalty = 365,
                       E0 = 4153.5, E1 = 12888.8, A0 = 139500, A1 = 2567000000000000000, 
                       Tf=4, slope=1.6, Tb = 4, Tu = 25, Tc = 36){
  par <- c(x, E0, E1, A0, A1, Tf, slope, Tb, Tu, Tc)
  
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, 
                              par = par))
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, 
                       no = pred_bloom)
  
  F <- sum((pred_bloom - bloomJDays)^2)
  return(F)
  
}