#' Evaluation function for phenoflex formatted for GenSA optimizer
#' 
#' This function communicates with the GenSA optimizer during the model calibration. It
#' takes model parameters and seasonlist, predicts bloom dates, compares it to the 
#' observed bloom dates and returns an error score. The optimizer
#' tries to reduce the error score.
#' 
#' Parameters can follow the traditional format (with the E0 to A1 parameters),
#' or include the intermediate parameter format (with theta_star to tau). In case of
#' intermediate parameters, they get converted to the traditional format before bloom dates are calculated.
##' 
#' @param x model parameters of PhenoFlex. Expected order when intermed_par = TRUE (default) yc, zc, s1, Tu, theta_star, theta_c, tau, pie_c, Tf, Tc, Tb, slope.
#' Expected order when intermed_par = FALSE yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope.
#' @param modelfn function used within the evaluation function to calculate the actual bloomday, often we use
#' the 'custom_GDH_wrapper' function for that
#' @param bloomJDays numeric containing the days of the year witht the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phenological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param intermed_par boolean, by default TRUE. Decides if intermediate parameter format is used, when set equal FALSE, then old format without conversion step is used
#' @param q10_lower lower boundary of Q10 metric, affects the E0 to A1 parameters. 
#' By default set to 1.5. If parameters lead to lower Q10 score, the evalfunction punishes the parameters, despite their 
#' predictive performance
#' @param q10_upper upper boundary of Q10 metric, affects the E0 to A1 parameters. 
#' By default set to 3.5. If parameters lead to higher Q10 score, the evalfunction punishes the parameters, despite their 
#' predictive performance
#' @param q10_penalty by default 10^10, the penalty when the calculated Q10 metric is outside the upper and lower boundary
#' @param na_penalty numeric, value which is used when the model fails to generate a prediction
#' for the bloom date. By default 365
#' @return residual sum of predicted and observed bloom dates. This is the score that optimizer tries to reduce
#' 
#' @author Lars Caspersen
#' @keywords utility
#' @importFrom nleqslv nleqslv
#'  
#' @export evaluation_function_gensa

evaluation_function_gensa <- function (x, 
                               modelfn, 
                               bloomJDays, 
                               SeasonList, 
                               intermed_par = TRUE,
                               q10_lower = 1.5, 
                               q10_upper = 3.5,
                               q10_penalty = 10^10,
                               na_penalty = 365) 
{
  par <- x
  
  if(intermed_par){
    params <- numeric(4)
    params[1] <- x[5]
    params[2] <- x[6]
    params[3] <- x[7]
    params[4] <- x[8]
    
    #convert parameters
    output <- nleqslv::nleqslv(c(500, 15000), solve_nle, 
                               jac = NULL, params, xscalm = "auto", method = "Newton", 
                               control = list(trace = 0, allowSingular = TRUE))
    if (output$termcd >= 3) {
      E0 <- NA
      E1 <- NA
      A0 <- NA
      A1 <- NA
      return(q10_penalty)
    }
    else {
      E0 <- output$x[1]
      E1 <- output$x[2]
      q = 1/params[1] - 1/params[2]
      A1 <- -exp(E1/params[1])/params[3] * log(1 - exp((E0 - 
                                                          E1) * q))
      A0 <- A1 * exp((E0 - E1)/params[2])
    }
    par[5:8] <- c(E0, E1, A0, A1)
  }

  
  #check for q10 criterion, if failed then penalize
  q10_e0 <- exp((10 * par[5])/(297 * 279))
  q10_e1 <- exp((10 * par[6])/(297 * 279))
  if(q10_e0 < q10_lower | q10_e1 > q10_upper | q10_e1 < q10_lower | q10_e1 > q10_upper){
    return(q10_penalty)
  }
  
  #predict bloom
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, 
                              par = par))
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, 
                       no = pred_bloom)
  return(sum((pred_bloom - bloomJDays)^2))
}