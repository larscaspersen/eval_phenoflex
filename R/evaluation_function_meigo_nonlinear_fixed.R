#' Evaluate multiple cultivars with selected PhenoFlex parameters fixed
#'
#' Evaluate cultivar-specific yc, zc and s1 parameters, with shared chill and
#' heat submodel parameters. The parameters theta_star, pie_c, Tc, Tb and slope
#' are supplied separately and held fixed. Characteristic chill parameters
#' are converted to kinetic parameters before prediction.
#'
#' Use this objective function with \link[MEIGOR]{MEIGO}. Cultivars in
#' bloomJDays and SeasonList must appear in the same order as their parameter
#' values; seasons must also appear in matching order within each cultivar.
#'
#' @param x Numeric parameter vector containing one yc value per cultivar,
#' followed by one zc and one s1 value per cultivar, then the shared parameters
#' Tu, theta_c, tau and Tf. Other parameters are supplied separately.
#'
#' @param modelfn Function that predicts a bloom date from a season's weather
#' and model parameters, such as \code{\link{custom_PhenoFlex_GDHwrapper}}.
#'
#' @param bloomJDays List of data frames with observed bloom dates in a pheno
#' column, one data frame per cultivar. Dates are expressed as days of the year.
#'
#' @param SeasonList Nested list with one season list per cultivar.
#' Each season is a weather data frame with Temp and JDay columns, ordered to
#' match the observations for that cultivar.
#'
#' @param n_cult Number of cultivars to calibrate together.
#'
#' @param na_penalty Numeric replacement for a missing bloom prediction. Defaults
#' to 365; this value is used in place of the missing date when calculating loss.
#'
#' @param theta_star Numeric optimal temperature for chill accumulation, in Kelvin.
#'
#' @param pie_c Numeric duration, in hours, of the alternating-temperature
#' experiment that produces complete chill negation. Defaults to 28.
#'
#' @param Tc Numeric upper temperature threshold of the GDH model, in degrees Celsius.
#'
#' @param Tb Numeric lower temperature threshold of the GDH model, in degrees Celsius.
#'
#' @param slope Numeric slope of the transition function that converts labile
#' chill to stable chill portions. Defaults to 1.6.
#'
#' @return A list with F (the residual sum of squares) and g (the model's
#' constraint values).
#'
#' @author Lars Caspersen
#' @keywords utility
#' @importFrom nleqslv nleqslv
#'
#' @section Deprecated:
#' \code{evaluation_function_meigo_nonlinear_fixed()} is deprecated since evalpheno 0.3.1.
#' Use a characteristic model collection with phenology_rss_combined(); choose fixed parameters through named bounds in fit_phenology().
#' The legacy implementation remains available and warns on use.
#' See \code{\link{evalpheno-deprecated}} for migration differences.
#' @export evaluation_function_meigo_nonlinear_fixed

evaluation_function_meigo_nonlinear_fixed <- function(x, 
                                                      modelfn,
                                                      bloomJDays,
                                                      SeasonList,
                                                      na_penalty = 365,
                                                      n_cult,
                                                      theta_star = 286,
                                                      pie_c = 28,
                                                      Tc = 36,
                                                      Tb  = 2,
                                                      slope = 1.6){
  .Deprecated(new = "phenology_rss_combined", package = "evalpheno",
              msg = "evaluation_function_meigo_nonlinear_fixed() is deprecated. Use a characteristic model collection with phenology_rss_combined(); choose fixed parameters through named bounds in fit_phenology(). See help('evalpheno-deprecated') for migration differences.")
  
  
  #bring x to another variable
  x_supplied <- x
  
  x <- c(x_supplied[1:((n_cult*3)+1)], #repeated yv, zc, s1 and lastly single Tu
         theta_star, #fixed
         x_supplied[c(length(x_supplied)-2, length(x_supplied)-1)], #theta_c and tau
         pie_c, #fixed
         x_supplied[length(x_supplied)], #Tf
         Tc,
         Tb,
         slope) 
  
  
  params<-numeric(4)
  
  params[1] <- theta_star  #theta*
  params[2] <- x[6+((n_cult -1) * 3)]    #theta_c
  params[3] <- x[7+((n_cult -1) * 3)]    #Tau(thetha*)
  params[4] <- pie_c     #pi_c
  
  
  output <- .fit_chill_parameters(params)
  
  
  #This is a numerical method which can produce non-convergence. Check this
  if (output$termcd >= 3){
    #if the nle algorithm has stalled just discard this solution
    E0<-NA; E1<-NA; A0<-NA; A1<-NA
    return(list(F=10^6, g=rep(10^6,5)))
    
    #You would add here a flag to let your optimization procedure know
    #That this solution should be ignored by lack of convergence
    
  } else {
    
    E0 <- output$x[1]
    E1 <- output$x[2]
    
    #A1 and A0 can be calculated through Equations 36 and 37
    
    q=1/params[1]-1/params[2]
    
    A1 <- output$A1
    A0 <- output$A0
  }
  
  
  #change the name of the parameters, dont know if necessary
  par <- x
  par[(5:8)+((n_cult -1) * 3)] <- c(E0, E1, A0, A1)
  
  
  
  #loop over the cultivars, calculate predicted days for each cultivar
  rss <- purrr::map_dbl(1:length(SeasonList), function(i){
    
    #extract the parameters
    par_cult <- par[c(i, 1+((n_cult -1))+i, 2+((n_cult -1)*2)+i, (length(par)-8):length(par))]
    #predict the bloom
    pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    #calculate for each cultivar the rss
    rss <- sum((pred_bloom - bloomJDays[[i]]$pheno)^2)
    return(rss)
    
  })
  
  
  #calculate the model performance value
  F <- sum(rss)
  
  
  
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
  
  
  
  return(list(F=F, g=g))
  
}
