#' Evaluate characteristic PhenoFlex parameters for multiple cultivars
#'
#' Evaluate cultivar-specific yc, zc and s1 parameters, with shared chill and
#' heat submodel parameters. Characteristic chill parameters (theta_star,
#' theta_c, tau and pie_c) are converted to kinetic parameters before prediction.
#'
#' Use this objective function with \link[MEIGOR]{MEIGO}. Cultivars in
#' bloomJDays and SeasonList must appear in the same order as their parameter
#' values; seasons must also appear in matching order within each cultivar.
#'
#' @param x Numeric parameter vector containing one yc value per cultivar,
#' followed by one zc and one s1 value per cultivar, then shared parameters
#' Tu, theta_star, theta_c, tau, pie_c, Tf, Tc, Tb, slope.
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
#' @param return_bloom_days Logical. If TRUE, return predicted bloom dates instead
#' of the objective value and constraint values. Defaults to FALSE.
#'
#' @return A list with F (the residual sum of squares) and g (the model's
#' constraint values) when return_bloom_days is FALSE. If parameter conversion
#' fails, F is 1e6 and g contains five values of 1e6.
#' If return_bloom_days is TRUE, return a list of numeric bloom-date vectors,
#' one per cultivar in SeasonList order. Each vector follows that cultivar's
#' season order. If parameter conversion fails, return the same list structure
#' with NA_real_ for every season. After successful conversion, missing bloom
#' predictions are replaced by na_penalty.
#'
#' @author Lars Caspersen
#' @keywords utility
#' @importFrom nleqslv nleqslv
#'
#' @export evaluation_function_meigo_nonlinear_combined

#bloomJDays is now a dataframe with the column pheno
#seasonlist elements must be in the same order as bloomJDays elements (also within the elements)
evaluation_function_meigo_nonlinear_combined <- function(x, 
                                                         modelfn,
                                                         bloomJDays,
                                                         SeasonList,
                                                         n_cult,
                                                         na_penalty = 365,
                                                         return_bloom_days = FALSE){
  
  #innput:
  #         x is the parameters in meigo
  #         modelfn is the function used to calculate the bloomdays
  #         SeasonList contains the weather data for the different years
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  #output: inequality constraints g
  #        model performance value F
  
  
  
  #instead of A0, A1, E0 and E1 we have now theta*, theta_c, Tau(thetha*) and pie_c
  #--> we need to solve now a non-linear system to calculate A0, A1, E0 and E1 based on the parameters
  #    ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c', 'Tf', 'Tc', 'Tb',  'slope')
  #x<- c(40,   190,   0.5,  25,   279,      287,       28,             26,       4,   36,    4,    1.60)
  
  params<-numeric(4)
  
  params[1] <- x[5+((n_cult -1) * 3)]   #theta*
  params[2] <- x[6+((n_cult -1) * 3)]    #theta_c
  params[3] <- x[7+((n_cult -1) * 3)]    #Tau(thetha*)
  params[4] <- x[8+((n_cult -1) * 3)]     #pi_c
  
  
  output <- .fit_chill_parameters(params)
  
  
  #This is a numerical method which can produce non-convergence. Check this
  if (output$termcd >= 3) {
    if (return_bloom_days) {
      return(lapply(seq_along(SeasonList), function(i) {
        rep(NA_real_, length(SeasonList[[i]]))
      }))
    }
    return(list(F = 1e6, g = rep(1e6, 5)))
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
  
  #split the parameters and reconstruct them for the different cultivar
  
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
  
  
  if(return_bloom_days == FALSE){
    #output
    return(list(F=F, g=g))
  } else{
    
    pred_bloom <- purrr::map(1:length(SeasonList), function(i){
      
      #extract the parameters
      par_cult <- par[c(i, 1+((n_cult -1))+i, 2+((n_cult -1)*2)+i, (length(par)-8):length(par))]
      #predict the bloom
      pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
      pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
      #calculate for each cultivar the rss
      return(pred_bloom)
      
    })
    
    
    
    return(pred_bloom)
  }
}
