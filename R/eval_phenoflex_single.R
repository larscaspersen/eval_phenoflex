#' Evaluation function for single fitting of cultivars
#'  
#' Takes model parameters for PhenoFlex, parameters of yc, zc, s1, Tu, theta_c, pie_c, tau, Tf, Tb, slope. 
#' 
#' The function is used to calibrate the PhenoFlex model. We assumed that the Tc and theta_star parameter
#' do not contribute much to the calibrated model performance and fixed them
#' to their default values. 
#' 
#' I used the function for the analysis on apple cultivars at lake Constance, though the 
#' results will probaly be not part of the final paper.
#' 
#' @param x model parameters of PhenoFlex, new format. yc, zc, s1, Tu, theta_c, pie_c, tau, Tf, Tb, slope. 
#' @param modelfn function used within the evaluation function to calculate the actual bloomday, often we use
#' the 'custom_GDH_wrapper' function for that
#' @param bloomJDays numeric containing the days of the year with the observed bloom
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param Tc numeric, by default 36. Critical temperature (°C) parameter of the GDH model.
#' @param theta_star numeric, by default 279. Optimal temperature (K) for chill accumulation of the
#' Dynamic Model. 
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
eval_phenoflex_single  <- function(x, 
                                   modelfn,
                                   bloomJDays,
                                   SeasonList,
                                   Tc = 36,
                                   theta_star = 279,
                                   return_pred = FALSE,
                                   na_penalty = 365){
  
  x <-  c(x[1:4], theta_star, x[5:8], Tc, x[9:10])
  
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
  
  params[1] <- x[5]   #theta*
  params[2] <- x[6]    #theta_c
  params[3] <- x[7]    #Tau(thetha*)
  params[4] <- x[8]     #pi_c
  
  
  output<-nleqslv::nleqslv(c(500, 15000), LarsChill::solve_nle, jac=NULL, params, xscalm="auto", method="Newton",
                           control=list(trace=0,allowSingular=TRUE))
  
  
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
    
    A1 <- -exp(E1/params[1])/params[3]*log(1-exp((E0-E1)*q))
    A0 <- A1*exp((E0-E1)/params[2])
  }
  
  
  #change the name of the parameters, dont know if necessary
  par <- x
  par[5:8] <- c(E0, E1, A0, A1)
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  #calculate the model performance value
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
  g[1] <- par[4] - par[11]
  #Tx >= Tb
  g[2] <- par[10] - par[11]
  #Tc >= Tu
  g[3] <- par[10] - par[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * par[5]) / (297 * 279))
  
  g[5] <- exp((10 * par[6]) / (297 * 279))
  
  
  
  if(return_pred){
    return(pred_bloom)
  } else{
    return(list(F=F, g=g))
  }
  
}