library(tidyverse)
library(chillR)

#Rcpp::sourceCpp("src/testfun.cpp")
Rcpp::sourceCpp("src/phenoflex_pop.cpp")

seasonlist <- KA_weather[which(KA_weather$Year>2004),] %>% 
  chillR::fix_weather() %>% 
  chillR::stack_hourly_temps(latitude=50.4) %>% 
  purrr::pluck('hourtemps') %>% 
  chillR::genSeasonList(years = 2006)
s <- seasonlist[[1]]

PhenoFlex(temp = s$Temp, times = seq_along(s$Temp), 
          yc = 40,
          zc = 200,
          s1 = 0.5,
          E0 = 4153.5,
          E1 = 12888.8,
          A0 = 139500,
          A1 = 2567000000000000000,
          slope = 1.6,
          Tf = 277-273)

yc_mean <- 20
yc_sd <- 5
zc_mean = 200
zc_sd = 10
n <- 10

set.seed(12345)
yc_pop <- rnorm(n = n, mean = yc_mean, sd = yc_sd)
zc_pop <- rnorm(n = n, mean = zc_mean, sd = zc_sd)

#the function expects to get already the position at which the forcing experiment should take place
jday_cut <- c(320)
i_cut <-purrr::map_int(jday_cut, function(x){
  floor(median(which(x == s$JDay)))
})  
#c++ starts counting at zero, correct for that
i_cut <- i_cut -1

pop_out <- PhenoFlex_pop(temp = s$Temp, 
                         times = seq_along(s$Temp), 
                         yc = yc_pop,
                         zc = zc_pop,
                         i_cut = i_cut,
                         max_days_forcing = 40,
                         force_temp = 23, 
                         s1 = 0.5,
                         E0 = 4153.5,
                         E1 = 12888.8,
                         A0 = 139500,
                         A1 = 2567000000000000000,
                         slope = 1.6,
                         Tf = 277-273, 
                         basic_output = FALSE)



i_cut_debug <- i_cut +1


PFcn(pop_out$y[[10]][i_cut_debug], yc_pop[10], 0.5)

yc_met_i <- min(which(pop_out$y[[1]] >= yc_pop[1]))
s$JDay[yc_met_i]

#variation in chill accumulation the same, just the requirement different
#that is why we always get the same curve



z_df <- do.call('cbind', pop_out$z) %>% 
  as.data.frame() %>% 
  mutate(i = 1:nrow(.)) %>% 
  pivot_longer(cols = -i) 

z_df %>% 
  group_by(i) %>% 
  summarise(min = min(value), max  = max(value)) %>% 
  ggplot(aes(x = i)) +
  geom_ribbon(aes(ymin = min, ymax = max)) +
  geom_line(data = z_df[z_df$name == 'V1',], aes(y = value), col = 'red')

y_df <- data.frame(y = pop_out$y[[1]],
               i = 1:length(pop_out$y[[1]]))

y_df %>% 
  ggplot(aes(x = i)) +
  geom_ribbon(aes(ymin = min(yc_pop), ymax = max(yc_pop))) +
  geom_line(aes(y = y))


out <- testfun(temp = 1:10, 
        times = 1:10,
        jday = 1:10)

hist(out$yc)
hist(out$zc)

#random number gen works

#next: identify the jday if I supply a vector of Jdays and
#cutting days



