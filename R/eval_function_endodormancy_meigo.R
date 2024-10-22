dynamic_model_wrapper <- function(Temp, par, summ = TRUE){
  
  # par <- c(4153.5, 12888.8, 139500, 2.567e+18, 1.6, 4)
  # Temp <- seasonlist_in[[1]][[1]]$Temp
  # summ <-  TRUE
  chillR::Dynamic_Model(HourTemp = Temp, 
                        summ = summ, 
                        E0 = par[1], 
                        E1 = par[2], 
                        A0 = par[3], 
                        A1 = par[4], 
                        slope = par[5], 
                        Tf = par[6] + 273) %>% 
    tail(n= 1) %>% 
    return()
}

eval_function_endodormancy_meigo <- function(x, SeasonList, dispersion_fun = 'calc_cv'){
  
  #convert the parameters from new to old
  x_new <- convert_parameters(c(0,0,0,0, x[1:5], 0, 0, x[6]))
  
  #inequality constrains
  g <- rep(0,2)
  
  #in case the conversion fails, return simply a large error for the evaluation function
  if(is.list(x_new)){
    return(list(F=x_new$F, g=rep(x_new$F, length(g))))
  }
  
  x_new <- x_new[c(5:9, 12)]
  
  #iterate over seasonlist of cultivars
  coefficinets_of_variation <- purrr::map_dbl(SeasonList, function(pheno){
    #iterate over seasons of one cultivar and sum gdh
    gdh <- purrr::map_dbl(pheno, function(p){
      dynamic_model_wrapper(p$Temp, x_new)
    })
    #calculate the dispersion function
    get(dispersion_fun)(gdh) %>%
      return()
  })
  
  F <- sum(coefficinets_of_variation)
  
  
  #E0 and E1
  g[1] <- exp((10 * x_new[1]) / (297 * 279))
  g[2] <- exp((10 * x_new[2]) / (297 * 279))
  
  
  
  #check the paraemeters
  
  return(list(F=F, g=g))
}