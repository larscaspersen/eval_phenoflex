#used in the almond study, and the lake constance study

eval_phenoflex_combined  <- function(x,
                                     modelfn,
                                     bloomJDays,
                                     SeasonList,
                                     ncult = 3,
                                     Tc = 36,
                                     theta_star = 279,
                                     return_pred = FALSE,
                                     na_penalty = 365){
  #calculate how long x should be
  assertthat::are_equal(length(x), 7 + (3*ncult))
  #x <- c(rep(24.79, ncult),	rep(337.04, ncult),	rep(0.2529,ncult),	17.72,	285.54,	   45.67,	  29.49,	2.97,	1.87,	2.69)
  #assume that I have three cultivars with their own yc, zc, s1
  params<-numeric(4)
  #findout position of theta c
  params[1] <- theta_star   #theta*
  params[2] <- x[5+((ncult -1) * 3)]    #theta_c
  params[3] <- x[6+((ncult -1) * 3)]    #Tau(thetha*)
  params[4] <- x[7+((ncult -1) * 3)]     #pi_c
  output<-nleqslv::nleqslv(c(500, 15000), solve_nle, jac=NULL, params, xscalm="auto", method="Newton",
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
  #       yc, zc, s1, Tu                            Tf                    Tb                   slope
  par <- c(x[1:((ncult * 3)+1)], E0, E1, A0, A1, x[length(x) - 2], Tc,  x[length(x) - 1], x[length(x)])
  pred_bloom <- NULL
  #loop over the cultivars, calculate predicted days for each cultivar
  rss <- purrr::map_dbl(1:length(SeasonList), function(i){
    #extract the parameters
    par_cult <- par[c(i, 1+((ncult -1))+i, 2+((ncult -1)*2)+i, (length(par)-8):length(par))]
    #predict the bloom
    pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    #calculate for each cultivar the rss
    rss <- sum((pred_bloom - bloomJDays[[i]])^2)
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
  if(return_pred){
    pred_bloom <- purrr::map(1:length(SeasonList), function(i){
      #extract the parameters
      par_cult <- par[c(i, 1+((ncult -1))+i, 2+((ncult -1)*2)+i, (length(par)-8):length(par))]
      #predict the bloom
      pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
      pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
      #calculate for each cultivar the rss
      return(pred_bloom)
    })
    pred_bloom <- do.call('c', pred_bloom)
    return(pred_bloom)
  } else{
    return(list(F=F, g=g))
  }
}

