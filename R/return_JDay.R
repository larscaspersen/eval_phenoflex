return_JDay <- function(index, Jday_vec, year_vec){
  
  # index <- index_budburst
  # Jday_vec <- chill_heat[,2]
  # year_vec <- chill_heat[,1]
  
  if (index == 0){
    JDay_out <- NA
  } else {
    JDay <- Jday_vec[index]
    JDaylist <- which(Jday_vec == JDay)
    
    #if we are in the norhtern hemisphere and the year corresponding to the index is the smaller one, return negative numbers relative to Jan-1 being 1
    if(length(unique(year_vec)) == 2 & year_vec[index] == min(year_vec)){
      
      JDay <- JDay - 365
      
    } 
    
    n <- length(JDaylist)
    if (n == 1){
      JDay_out <- JDay
      return(JDay)
    } else {
      JDay_out <- JDay + which(JDaylist == index)/n - 1/(n/ceiling(n/2))
    }
  }
  return(JDay_out)
}