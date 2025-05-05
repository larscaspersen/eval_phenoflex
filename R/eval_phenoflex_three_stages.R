#' Evaluation function for three phenological stages at once
#'  
#' Takes model parameters for PhenoFlex, however, three individual heat requirements for 
#' each modelled stage individually.
#' 
#' The function returns for one season the timing of the three phenological stages. 
#' E0 to A1 parameters of the chill submodel are replaced by intermediate parameters. 
#' Parameters Tc and theta_star are fixed, as they contribute little to
#' model performance. The function is usually called by an optimization algorithm
#' during model calibration.
#' 
#' @param x model parameters of PhenoFlex, new format: 
#' c(yc, zc1, zc2, zc3 s1,Tu, theta_c, pie_c, tau, Tf, Tb, slope)
#' will get converted to the format: E0, E1, A0, A1, Tf, slope to run the Dynamic Model 
#' @param modelfn function used within the evaluation function to calculate the actual bloomday, often we use
#' the 'custom_GDH_wrapper' function for that
#' @param bloomdates_df data.frame, with four columns: year, firstbloom, fullbloom and budbreak.
#' The phenological stages have the julian day of the occurence as values. 
#' @param SeasonList list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
#' generated using \link[chillR]{genSeasonList}
#' @param Tc numeric, by default 36. Critical temperature (°C) parameter of the GDH model.
#' @param theta_star numeric, by default 279. Optimal temperature (K) for chill accumulation of the
#' Dynamic Model. 
#' @param na_penalty numeric, by default 365. Penalty for the phenology
#' prediction function when it fails return a bloom prediction
#' @return list with two elements. First is called 'f' and contains the residual sum of squares of the model. The 
#' second is 'g' which is the values of the additional model constraints defined in the function.
#' If the flag for return_pred set TRUE, then it returns bloom dates
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom assertthat are_equal
#' @importFrom nleqslv nleqslv
#' @examples 
#' \dontrun{
#' 
#' ncult <- 3
#' #          theta_star, theta_c, pie_c, tau, Tf, slope
#' par <-   c(40, 50, 150, 200, 0.5, 26, 281, 24, 45, 4, 1.6)
#' 
#' #prepare weather data
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)
#' 
#' eval_phenoflex_three_stages(x = par, par = SeasonList)
#' }
#' @export eval_phenoflex_three_stages
eval_phenoflex_three_stages  <- function(x, 
                                          modelfn,
                                          bloomdates_df,
                                          SeasonList,
                                          Tc = 36,
                                          theta_star = 279,
                                          na_penalty = 365){
  
  #x <- c(64.6278698, 344.1353029, 351.9283211, 487.7301338,   0.6457447,  25.7247364, 286.2740348,  35.3590167,  43.9834312,   7.3949351,   7.2267818,   3.7180116)
  
  x <-  c(x[1:6], theta_star, x[7:10], Tc, x[11:12])
  
  #x <- c(rep(24.79, ncult),	rep(337.04, ncult),	rep(0.2529,ncult),	17.72,	285.54,	   45.67,	  29.49,	2.97,	1.87,	2.69)
  
  #assume that I have three cultivars with their own yc, zc, s1
  params<-numeric(4)
  
  
  #findout position of theta c
  
  params[1] <- x[5+2]   #theta*
  params[2] <- x[6+2]    #theta_c
  params[3] <- x[7+2]    #Tau(thetha*)
  params[4] <- x[8+2]     #pi_c
  
  
  output<-nleqslv::nleqslv(c(500, 15000), LarsChill::solve_nle, jac=NULL, params, xscalm="auto", method="Newton",
                           control=list(trace=0,allowSingular=TRUE))
  
  
  #This is a numerical method which can produce non-convergence. Check this
  if (output$termcd >= 3){
    #if the nle algorithm has stalled just discard this solution
    E0<-NA; E1<-NA; A0<-NA; A1<-NA
    return(list(F=10^6, g=rep(10^6,8)))
    
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
  
  
  par <- x
  
  #add the converted parameters
  par[(5+2):(8+2)] <- c(E0, E1, A0, A1)
  
  #make the predictions for only the first bloom, and make checks for the others
  par_run_pheno <- par[c(1,3,5:14)]
  #par_run_pheno <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e+13, 4, 36, 4, 1.6)
  
  #calculate the predicted flower dates
  model_out <- lapply(X = SeasonList, FUN = modelfn, par = par_run_pheno)
  
  
  #heat requirement for full bloom
  zc1 <- par[2]
  zc2 <- par[3]
  zc3 <- par[4]
  yc <- par[1]
  
  #calculate the predicted full bloom by comparing the accumulated heat with heat requirement zc3
  pred_rest <- purrr::map(model_out, 'chill_heat') %>% 
    purrr::map(function(x){
      
      #x <- model_out[[1]]$chill_heat
      
      #julian day of fullbloom
      fullbloom_index <- which(x[,'heat_accumulated'] >= zc3) 
      
      if(length(fullbloom_index)){
        JDay_fullbloom <- fullbloom_index %>% 
          min() %>% 
          return_JDay(Jday_vec =  x[,'JDay'], year_vec = x[, 'Year'])
      } else {
        JDay_fullbloom <- na_penalty
      }
      
      
      budburst_index <- which(x[,'heat_accumulated'] >= zc1)
      
      if(length(budburst_index)){
        JDay_budburst <- budburst_index %>% 
          min() %>% 
          return_JDay(Jday_vec =  x[,'JDay'], year_vec = x[, 'Year'])
      } else {
        JDay_budburst <- na_penalty
      }
      
      
      chill_met_at_budburst <- x[JDay_budburst, 3] >= yc
      
      return(list(JDay_fullbloom = JDay_fullbloom,
                  JDay_budburst = JDay_budburst,
                  chill_met_at_budburst = chill_met_at_budburst))
    })
  
  #extract the individual obtained julian days
  pred_fullbloom <- purrr::map_dbl(pred_rest, 'JDay_fullbloom') %>% 
    replace_na(na_penalty)
  
  pred_budburst <- purrr::map_dbl(pred_rest, 'JDay_budburst') %>% 
    replace_na(na_penalty)
  
  #extract the predicted start of the bloom
  pred_firstbloom <- purrr::map_dbl(model_out, 'JDay') %>% 
    replace_na(na_penalty)
  
  
  #in cases where the chill requirement was not met at budburst, the prediction becomes 365 instead
  chill_at_budburst <- purrr::map_lgl(pred_rest, 'chill_met_at_budburst')
  
  if(all(chill_at_budburst) == FALSE){
    pred_fullbloom[chill_at_budburst == FALSE] <- na_penalty
    pred_budburst[chill_at_budburst == FALSE] <- na_penalty
    pred_firstbloom[chill_at_budburst == FALSE] <- na_penalty
  }
  

  #calculate the combined rss
  rss_firstbloom <- sum((pred_firstbloom - bloomdates_df$firstbloom)^2)
  rss_fullbloom <- sum((pred_fullbloom - bloomdates_df$fullbloom)^2)
  rss_budbreak <- sum((pred_budburst - bloomdates_df$budbreak)^2)
  
  
  #calculate the model performance value
  F <- rss_firstbloom + rss_fullbloom + rss_budbreak
  
  
  
  #####
  #inequality constraints
  #####
  
  #this is the vector containing the values for the constraints
  #at first initialize the vector
  g <- rep(0,8)
  
  #inequality constraints are reformulated as differences
  
  #heat requirements must not be smaller than from the previous stage
  g[1] <- zc2 - zc1
  g[2] <- zc3 - zc1
  g[3] <- zc3 - zc2
  
  #Tu >= Tb
  g[4] <- par[4+2] - par[11+2]
  #Tx >= Tb
  g[5] <- par[10+2] - par[11+2]
  #Tc >= Tu
  g[6] <- par[10+2] - par[4+2]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[7] <- exp((10 * par[5+2]) / (297 * 279))
  
  g[8] <- exp((10 * par[6+2]) / (297 * 279))
  
  
  
  return(list(F=F, g=g))
  
}



