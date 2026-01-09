library(tidyverse)
library(chillR)

#Rcpp::sourceCpp("src/testfun.cpp")
Rcpp::sourceCpp("src/phenoflex_pop.cpp")

# s <- read.csv('experimental/Ravensburg_hourly_temp.csv') %>%
#   genSeasonList(years = 2022) %>%
#   purrr::pluck(1)
s <- read.csv('experimental/hohenheim_aug21-jul22.csv', sep = ';', dec = ',') %>% 
  mutate(Date = lubridate::dmy(Tag),
         Hour = lubridate::hm(Stunde) %>% hour(),
         JDay = lubridate::yday(Date),
         Year = lubridate::year(Date),
         Temp = AVG_TA200,
         ) %>% 
  select(Temp, JDay, Year)

#               yc          zc            s1      Tu          theta_star  theta_c
par_topaz <- c(40.0336321, 181.2843981, 0.1473177, 21.0231964, 279, 285.6807267,
#               tau         pie_c       Tf        Tc  Tb        slope        
               34.2354385, 32.6464321, 8.0650228, 36, 5.9122156, 1.9426370)
par_topaz <- LarsChill::convert_parameters(par_topaz)


yc_mean <- par_topaz[1] 
yc_sd <- 5
zc_mean = par_topaz[2]
zc_sd = 10
n <- 100

set.seed(12345)
yc_pop <- rnorm(n = n, mean = yc_mean, sd = yc_sd)
zc_pop <- rnorm(n = n, mean = zc_mean, sd = zc_sd)

#---------------------------------#
#OBSERVED

exp_obs <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                              sheet = 'T_2022_term+spur')

exp_obs <- exp_obs %>% 
  group_by(Data) %>% 
  mutate(cumsum = cumsum(Bubble_size))

exp_obs$Date <- lubridate::parse_date_time(exp_obs$Data, orders = 'dmy')
exp_obs$yday <- lubridate::yday(exp_obs$Date)
exp_obs$h_after_cut <- exp_obs$Days*24

exp_obs <- exp_obs %>% 
  filter(yday >= 300 | yday < 55)

#-----------------------------------#


#the function expects to get already the position at which the forcing experiment should take place
jday_cut <- unique(exp_obs$yday)
jday_name <- lubridate::stamp("Nov 03", orders = '%b %d')(unique(exp_obs$Date))
#jday_name <- c('Nov 03', 'Nov 17', 'Dec 01', 'Dec 15', 'Dec 29', 'Jan 12', 'Jan 26', 'Feb 09', 'Feb23')
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
                         max_days_forcing = 50,
                         forcing_temperature = 23, 
                         s1 = par_topaz[3],
                         E0 = par_topaz[5],
                         E1 = par_topaz[6],
                         A0 = par_topaz[7],
                         A1 = par_topaz[8],
                         slope = par_topaz[12],
                         Tf = par_topaz[9],
                         Tb = par_topaz[11],
                         Tu = par_topaz[4],
                         Tc = par_topaz[10],
                         placeholder_fail = 9999,
                         basic_output = FALSE)


pop_out$exp %>% 
  as.data.frame() %>% 
  mutate(jday = jday_cut) %>% 
  pivot_longer(cols = -jday) %>%
  mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
         value = ifelse(value == 9999, yes = 24*51, no = value),
         value_mod = value / 24,
         jday_fact = factor(jday, levels = jday_cut, 
                            labels = jday_name)) %>% 
  ggplot(aes(x = value_mod)) +
  stat_ecdf(aes(color = jday_fact),
            geom = 'step', size = 1.5) +
  coord_cartesian(xlim = c(0,51)) +
  xlab('Days after Cutting') + 
  ylab('Share of buds that reach Flowering') +
  facet_wrap(~jday_fact)


#add missing days to antons table
miss_df <- exp_obs %>% 
  group_by(Date, yday) %>% 
  summarise(max = max(cumsum),
            miss = 1 - max) %>% 
  ungroup() %>% 
  mutate(h_after_cut = 51 *24,
         Bubble_size = miss,
         Days = 51,
         cumsum = 1) %>% 
  select(Days, Date, yday, h_after_cut, Bubble_size, cumsum)

exp_obs_plot <- exp_obs %>% 
  select(Days, Date, yday, h_after_cut, Bubble_size, cumsum) %>% 
  rbind(miss_df)

exp_obs_plot_sub <- exp_obs_plot %>% 
  mutate(jday_fact = factor(yday, levels = jday_cut, 
                            labels = jday_name)) %>% 
  filter(Bubble_size != 0) 

exp_obs_plot_sub %>% 
  ggplot(aes(x = Days)) +
  stat_ecdf(aes(color = jday_fact),
            geom = 'step', size = 1.5) +
  #geom_point(col = 'black', shape = 4) +
  coord_cartesian(xlim = c(0,51)) +
  facet_wrap(~jday_fact)

test <- exp_obs_plot %>% 
  mutate(jday_fact = factor(yday, levels = jday_cut, 
                            labels = jday_name)) 

test %>%  
  ggplot(aes(x = Days, y = cumsum))+
  geom_step() +
  facet_wrap(~jday_fact)

ribbon_df <- data.frame(x = 50:60,
                        ymin = -Inf,
                        ymax = Inf)


pop_out$exp %>% 
  as.data.frame() %>% 
  mutate(jday = jday_cut) %>% 
  pivot_longer(cols = -jday) %>%
  mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
         value = ifelse(value == 9999, yes = 24*51, no = value),
         value_mod = value / 24,
         jday_fact = factor(jday, levels = jday_cut, 
                            labels = jday_name)) %>% 
  ggplot() +
  stat_ecdf(aes(x = value_mod,
                color = 'Modelled',
                linetype = 'Modelled'),
            geom = 'step', size = 1.5) +
  geom_step(data = test, aes(y = cumsum, x = Days, col = 'Observed', linetype = 'Observed'),
            size = 1.5) +
  geom_ribbon(data = ribbon_df, aes(ymin = ymin, ymax = ymax, x = x),
              fill = 'grey') +
  scale_color_discrete(name = 'Data Source') +
  scale_linetype_discrete(name = 'Data Source') +
  coord_cartesian(xlim = c(0,51)) +
  facet_wrap(~jday_fact) +
  theme_bw() +
  theme(legend.position = 'bottom')

helper_run_pop_model <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100){
  
  #--------------#
  #draw population of yc, zc
  
  set.seed(12345)
  yc_pop <- rnorm(n = n, mean = par[1], sd =  yc_sd)
  zc_pop <- rnorm(n = n, mean = par[2], sd = zc_sd)
  
  #---------------#
  #identify timepoints of cutting in temperature data
  
  i_cut <-purrr::map_int(jday_cut, function(x){
    floor(median(which(x == temp_df$JDay)))
  })  
  #c++ starts counting at zero, correct for that
  i_cut <- i_cut -1
  
  #----------------#
  #run model
  
  PhenoFlex_pop(temp = temp_df$Temp, 
                times = seq_along(temp_df$Temp), 
                yc = yc_pop,
                zc = zc_pop,
                i_cut = i_cut,
                max_days_forcing = 50,
                forcing_temperature = 23, 
                s1 = par[3],
                E0 = par[5],
                E1 = par[6],
                A0 = par[7],
                A1 = par[8],
                slope = par[12],
                Tf = par[9],
                Tb = par[11],
                Tu = par[4],
                Tc = par[10],
                placeholder_fail = 9999,
                basic_output = FALSE) %>% 
    return()
  
}

helper_plot_forcing_exp <- function(model_res, obs, jday_cut){
  
  ribbon_df <- data.frame(x = 50:60,
                          ymin = -Inf,
                          ymax = Inf)
  
  mod_df <- model_res$exp %>% 
    as.data.frame() %>% 
    mutate(jday = jday_cut) %>% 
    pivot_longer(cols = -jday) %>%
    mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
           value = ifelse(value == 9999, yes = 24*51, no = value),
           value_mod = value / 24,
           jday_fact = factor(jday, levels = jday_cut, 
                              labels = jday_name)) 
  
  ggplot(mod_df) +
    stat_ecdf(aes(x = value_mod,
                  color = 'Modelled',
                  linetype = 'Modelled'),
              geom = 'step', size = 1.5) +
    geom_step(data = obs, aes(y = cumsum, x = Days, col = 'Observed', linetype = 'Observed'),
              size = 1.5) +
    geom_ribbon(data = ribbon_df, aes(ymin = ymin, ymax = ymax, x = x),
                fill = 'grey') +
    scale_color_discrete(name = 'Data Source') +
    scale_linetype_discrete(name = 'Data Source') +
    coord_cartesian(xlim = c(0,51)) +
    facet_wrap(~jday_fact) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

#               yc          zc            s1      Tu          theta_star  theta_c
par_topaz <- c(40.0336321, 181.2843981, 0.1473177, 21.0231964, 279, 285.6807267,
               #               tau         pie_c       Tf        Tc  Tb        slope        
               34.2354385, 32.6464321, 8.0650228, 36, 5.9122156, 1.9426370) %>% 
  LarsChill::convert_parameters()

pop_out <- helper_run_pop_model(par = par_topaz, yc_sd = 5, zc_sd = 10, jday_cut = jday_cut, temp_df = s, n = 100) 

pop_out %>% helper_plot_forcing_exp(obs = test, jday_cut = jday_cut)


forcing_days <- 50
placeholder_fail <- 9999

exp_obs_plot_sub_list <- list()

#anton said every two to three days it got checked. 
#I choose every three days
#add that info to the figure, even if measurement did not change
for(jd_f in levels(exp_obs_plot_sub$jday_fact)){
  
  #jd_f <- levels(exp_obs_plot_sub$jday_fact)[1]
  
  #subset
  sub <- exp_obs_plot_sub[exp_obs_plot_sub$jday_fact == jd_f,]
  sub_merged <- data.frame()
  
  #assume every second day check
  example_days <- seq(0, 50, by = 3)
  
  #add empty df to list, will get filled during loop
  exp_obs_plot_sub_list[[jd_f]] <- data.frame()
  
  #merge example days with actual noted days
  for(i in 1:nrow(sub)){
    #at first entry, check if there are previous possible sample dates, if so give them the value zero
    if(i == 1){
      if(any(example_days < sub$Days[i])){
        
        cumsum_fill <- 0
        example_days_i <- which(example_days < sub$Days[i])
        
      }
    } else {
      cumsum_fill <- exp_obs_plot_sub_list[[jd_f]]$cumsum %>% tail(n = 1)
      example_days_i <- which(example_days < sub$Days[i] & example_days > max(exp_obs_plot_sub_list[[jd_f]]$Days))
      
      #if there are no days to add / fill, go to next
      if(length(example_days_i) == 0) next
    }
    
    #add data to the list of data.frames
    df_add <- data.frame(Days = example_days[example_days_i],
                                 h_after_cut = example_days[example_days_i] * 24,
                                 cumsum = cumsum_fill,
                                 jday_fact = jd_f,
                         yday = sub$yday[1]) %>% 
      rbind(sub[i, c('Days', 'h_after_cut', 'cumsum', 'jday_fact', 'yday')]) %>% 
      mutate( jday_fact = factor(jday_fact, levels = jday_name))
    
    exp_obs_plot_sub_list[[jd_f]] <- rbind(exp_obs_plot_sub_list[[jd_f]],
                                           df_add)
  }
}
exp_obs_plot_filled_df <- bind_rows(exp_obs_plot_sub_list) %>% 
  mutate(           jday_fact = factor(jday_fact, levels = jday_name))


#prepare forcing observations
forcing_obs_s1 <- purrr::map(levels(exp_obs_plot_filled_df$jday_fact), function(x){
  hours <- exp_obs_plot_filled_df$h_after_cut[exp_obs_plot_filled_df$jday_fact == x]
  cumsum <- exp_obs_plot_filled_df$cumsum[exp_obs_plot_filled_df$jday_fact == x]
  
  hours <- ifelse(hours > forcing_days * 24, yes = placeholder_fail, no = hours)
  return(list(yday_cut = exp_obs_plot_filled_df$yday[exp_obs_plot_filled_df$jday_fact == x][1],
              hours = hours,
              cumsum = cumsum))
})

#caluclate sum of squared difference, calculate dynamic time warping

exp_perf <- purrr::map(1:nrow(pop_out$exp), function(j){
  
  #replace placeholder with one day after max forcing
  exp_in <- pop_out$exp[j,]
  exp_in <- ifelse(exp_in == placeholder_fail, yes = (forcing_days + 1)*24, no = exp_in)
  
  fun <- ecdf(exp_in)
  exp_mod <- fun(exp_obs_plot_sub_list[[j]]$h_after_cut)
  return(exp_mod)
})


performance_df <- purrr::map(seq_along(exp_perf), function(i){
  mae <- mean(abs(exp_perf[[i]] - exp_obs_plot_sub_list[[i]]$cumsum)) %>% round(digits = 1)
  rmse <- chillR::RMSEP(exp_perf[[i]], exp_obs_plot_sub_list[[i]]$cumsum) %>% round(digits = 1)
  ssd <- sum((exp_perf[[i]] - exp_obs_plot_sub_list[[i]]$cumsum)^2) %>% round(digits = 1)
  dtw <- dtw::dtw(exp_perf[[i]], exp_obs_plot_sub_list[[i]]$cumsum)$distance %>% round(digits = 1)
  return(data.frame(mae = mae, rmse = rmse, ssd = ssd, dtw = dtw))
}) %>% 
  bind_rows() %>% 
  mutate(jday_fact = factor(jday_name, levels = jday_name))



helper_plot_forcing_exp(pop_out, obs = test, jday_cut = jday_cut) +
  geom_text(data = performance_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, nsmall = 1))), vjust = 1, hjust = 0) +
  #geom_text(data = performance_df, aes(x = 0, y = 0.8, label = paste('RMSE:' ,format(rmse, digits = 2))), vjust = 1, hjust = 0) +
  ggtitle(paste('Forcing Experiment. Total Sum of Squared Differences (SSD):', format(sum(performance_df$ssd), nsmall = 1)))


#make plot for all seasons in KOB





helper_run_pop_model(par = par_topaz, yc_sd = 10, zc_sd = 0, jday_cut = jday_cut, temp_df = s, n = 100) %>% 
  helper_plot_forcing_exp(obs = test, jday_cut = jday_cut)


#no population effect
helper_run_pop_model(par = par_topaz, yc_sd = 0, zc_sd = 0, jday_cut = jday_cut, temp_df = s, n = 100) %>% 
  helper_plot_forcing_exp(obs = test, jday_cut = jday_cut)

#--> need to move the curve to the left

start <- Sys.time()
par_topaz1 <- c(34, 165, 0.35, 21.0231964, 279, 285.6807267,
                #tau         pie_c       Tf        Tc  Tb        slope        
                34.2354385, 32.6464321, 8.0650228, 36, 5.9122156, 1.9426370) %>% 
  LarsChill::convert_parameters()

helper_run_pop_model(par = par_topaz1, yc_sd = 4, zc_sd = 0, jday_cut = jday_cut, temp_df = s, n = 100) %>% 
  helper_plot_forcing_exp(obs = test, jday_cut = jday_cut)
end <-  Sys.time()

end-start

shiny::runApp('experimental/app_pop-model/')







helper_bloomint_to_jday <- function(bloomindex, x, return_fail = 9999){
  if (bloomindex == 0) {
    return(return_fail)
  }
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  if (length(unique(x$Year)) == 2 & x$Year[bloomindex] == min(x$Year)) {
    JDay <- JDay - 365
  }
  n <- length(JDaylist)
  if (n == 1) {
    return(JDay)
  }
  return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
}

bloom_df <- data.frame(pheno = purrr::map_dbl(pop_out$bloomindex, helper_bloomint_to_jday, x = s))

yday_obs_10 <- lubridate::ymd('2022-04-14') %>% lubridate::yday()
yday_obs_50 <- lubridate::ymd('2022-04-18') %>% lubridate::yday()
yday_x_scale <- (floor(min(bloom_df$pheno))):(ceiling(max(bloom_df$pheno)))
yday_x_scale_label <- as.Date(yday_x_scale,
                              origin = '2021-12-31') %>% 
  format("%b %d")
annotate_obs <- data.frame(x = c(yday_obs_10, yday_obs_50),
                          y = c(0.1, 0.5),
                          label = c('Observed\nFirst Flowering', 'Observed\nFull Flowering'))
annotate_pred <- data.frame(x = quantile(bloom_df$pheno, probs = c(0.1, 0.5)),
                           y = c(0.1, 0.5),
                           label = c('Modelled\nFirst Flowering', 'Modelled\nFull Flowering'))

ggplot(bloom_df) +
  stat_ecdf(aes(x = pheno,
                color = 'Modelled'),
            geom = 'step', size = 1.5) +
  geom_point(data = annotate_obs, aes(y = y, x = x, col = 'Observed'),
             size = 3, shape = 25, stroke = 2, fill = 'white') +
  geom_point(data = annotate_pred, aes(y = y, x = x, col = 'Modelled'),
             size = 3, shape = 25, stroke = 2, fill = 'white') +
  ggrepel::geom_label_repel(data = annotate_obs, aes(x = x, y = y, label = label),
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.3, "lines"),
                            nudge_y = 0.1,
                            nudge_x = -0.1) +
  ggrepel::geom_label_repel(data = annotate_pred, aes(x = x, y = y, label = label),
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.3, "lines"),
                            nudge_y = -0.1,
                            nudge_x = 0.1) +
  scale_color_discrete(name = 'Data Source') +
  scale_x_continuous(breaks = yday_x_scale, labels = yday_x_scale_label) +
  theme_bw() +
  theme(legend.position = 'bottom')


#max is roughly equal to 75%

max(exp_obs$Bubble_size)

exp_obs %>% 
  group_by(Data) %>% 
  summarise(sum = sum(Bubble_size))


wrapper_phenoflex_pop <- function(yc_pop,
                                  zc_pop,
                                  par_rest,
                                  s, 
                                  cutting_days, 
                                  forcing_days = 50, forcing_temp = 23,
                                  placeholder_fail = 9999){
  
  pop_out <- PhenoFlex_pop(temp = s$Temp, 
                           times = seq_along(s$Temp), 
                           yc = yc_pop,
                           zc = zc_pop,
                           i_cut = cutting_days,
                           max_days_forcing = forcing_days,
                           forcing_temperature = forcing_temp, 
                           s1 = par_rest[1],
                           E0 = par_rest[3],
                           E1 = par_rest[4],
                           A0 = par_rest[5],
                           A1 = par_rest[6],
                           slope = par_rest[10],
                           Tf = par_rest[7],
                           Tb = par_rest[9],
                           Tu = par_rest[2],
                           Tc = par_rest[8],
                           placeholder_fail = placeholder_fail,
                           basic_output = FALSE)
  
  return(pop_out[c('bloomindex', 'exp')])
}

#par: yc_mean, yc_sd, zc_mean, zc_sd, s1
par_topaz <- c(40.0336321, 181.2843981, 0.1473177, 21.0231964, 279, 285.6807267,
               #               tau         pie_c       Tf        Tc  Tb        slope        
               34.2354385, 32.6464321, 8.0650228, 36, 5.9122156, 1.9426370)
par_topaz <- LarsChill::convert_parameters(par_topaz)


par <- c(40.0336321, 5, 181.2843981, 10, 0.1473177)
par_fixed <- c(2.102320e+01, 3.855167e+03, 1.025983e+04, 4.399565e+04, 2.397905e+14, 8.065023e+00,
               3.600000e+01, 5.912216e+00, 1.942637e+00)
forcing_days <- 50
placeholder_fail <- 9999

forcing_obs_s1 <- purrr::map(levels(exp_obs_plot_sub$jday_fact), function(x){
  hours <- exp_obs_plot_sub$h_after_cut[exp_obs_plot_sub$jday_fact == x]
  cumsum <- exp_obs_plot_sub$cumsum[exp_obs_plot_sub$jday_fact == x]
  
  hours <- ifelse(hours > forcing_days * 24, yes = placeholder_fail, no = hours)
  return(list(yday_cut = exp_obs_plot_sub$yday[exp_obs_plot_sub$jday_fact == x][1],
              hours = hours,
              cumsum = cumsum))
})

forcing_obs <- list(forcing_obs_s1)
seasonlist <- list(s)
yday <- c('14-04-2022', '18-04-2022') %>% lubridate::dmy() %>% yday()
#         first flower  full flower
share <- c(0.1,         0.5)
bloom_list <- list(list(yday = yday,
                        share = share))

eval_fun_phenopop_fixed <- function(par, seasonlist, bloom_obs, forcing_obs,
                                    modelfn = wrapper_phenoflex_pop,
                                    weight_exp = 1, #scales the weight of experiment performance, values > 1 place more weight on experiment
                                    par_fixed = c(21.0231964, 279, 285.6807267,
                                                  34.2354385, 32.6464321, 8.0650228, 36, 5.9122156, 1.9426370),
                                    forcing_days = 50, forcing_temp = 23, placeholder_fail = 9999,
                                    population_size = 100){
  
  ##generate population of chill and forcing requirements
  set.seed(12345)
  yc_pop <- rnorm(n = population_size, mean = par[1], sd = par[2])
  zc_pop <- rnorm(n = population_size, mean = par[3], sd = par[4])
  
  #combine parameters
  par_rest <- c(par[5], par_fixed)
  
  #iterate over seasons
  f <- purrr::map_dbl(1:seq_along(seasonlist), function(i){
    
    #return predicted results of bloom and experiments
    out <- modelfn(yc_pop = yc_pop,
                   zc_pop = zc_pop, 
                   par_rest = par_rest, 
            s = seasonlist[[i]], 
            cutting_days = purrr::map_int(forcing_obs[[i]], 'yday_cut'),
            forcing_days = forcing_days, 
            forcing_temp = forcing_temp,
            placeholder_fail = placeholder_fail)
    
    #compare observed and predicted results of forcing experiments.
    #return average squared distance of obs and predicted budbreaks 
    #number of observation vary by forcing experimen, because only contains entries when something changed
    
    exp_perf <- purrr::map_dbl(1:nrow(out$exp), function(j){
      fun <- ecdf(out$exp[j,])
      exp_mod <- fun(forcing_obs[[i]][[j]]$hours)
      exp_obs <- forcing_obs[[i]][[j]]$cumsum
      return(mean(abs((exp_mod - exp_obs)*100)^2))
           
    }) %>% sum()
    
    #convert bloomint (position in table when bloom reached) to julian days
    bloom_ydays_fun <- purrr::map_dbl(out$bloomindex, helper_bloomint_to_jday, x = s) %>% 
      ecdf()
    
    #return predicted bloomydays for time points of observations
    bloom_pred <- bloom_ydays_fun(bloom_list[[i]]$yday)
    bloom_perf <- mean(abs((bloom_pred - bloom_list[[i]]$share)*100)^2) 
    
    #calculate avergae experiment performance
    exp_perf_mean <- (exp_perf / nrow(out$exp))
    return(bloom_perf + (exp_perf_mean * weight_exp))
    
  })
  
  return(sum(f))
  
}

par <- c(40.0336321, 5, 181.2843981, 10, 0.1473177)
par_fixed <- c(2.102320e+01, 3.855167e+03, 1.025983e+04, 4.399565e+04, 2.397905e+14, 8.065023e+00,
               3.600000e+01, 5.912216e+00, 1.942637e+00)
forcing_days <- 50
placeholder_fail <- 9999

eval_fun_phenopop_fixed(par = par, 
                        seasonlist = seasonlist, 
                        bloom_obs = bloom_list, 
                        forcing_obs = forcing_obs, 
                        modelfn = wrapper_phenoflex_pop, 
                        weight_exp = 1, 
                        par_fixed =  par_fixed, 
                        forcing_days = 50, 
                        forcing_temp = 23, 
                        placeholder_fail = 9999, 
                        population_size = 100)


library(DEoptim)
set.seed(123456)

res <- DEoptim(fn = eval_fun_phenopop_fixed,
               lower = c(30, 0, 100, 0, 0),
               upper = c(70, 10, 400, 20, 1.5),
               control = DEoptim.control(itermax = 200),
               seasonlist = seasonlist, 
               bloom_obs = bloom_list, 
               forcing_obs = forcing_obs, 
               modelfn = wrapper_phenoflex_pop, 
               weight_exp = 1, 
               par_fixed =  par_fixed, 
               forcing_days = 50, 
               forcing_temp = 23, 
               placeholder_fail = 9999, 
               population_size = 100)

par <- c(51.394318,    8.344372,  178.120718,   11.550661,    1.079505)
n <- 100
set.seed(12345)
yc_pop <- rnorm(n = n, mean = par[1], sd = par[2])
zc_pop <- rnorm(n = n, mean = par[3], sd = par[4])

pop_out <- PhenoFlex_pop(temp = s$Temp, 
                         times = seq_along(s$Temp), 
                         yc = yc_pop,
                         zc = zc_pop,
                         i_cut = i_cut,
                         max_days_forcing = 50,
                         forcing_temperature = 23, 
                         s1 = par[5],
                         E0 = par_fixed[2],
                         E1 = par_fixed[3],
                         A0 = par_fixed[4],
                         A1 = par_fixed[5],
                         slope = par_fixed[9],
                         Tf = par_fixed[6],
                         Tb = par_fixed[8],
                         Tu = par_fixed[1],
                         Tc = par_fixed[7],
                         placeholder_fail = 9999,
                         basic_output = FALSE)

pop_out$exp %>% 
  as.data.frame() %>% 
  mutate(jday = jday_cut) %>% 
  pivot_longer(cols = -jday) %>%
  mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
         value = ifelse(value == 9999, yes = 24*51, no = value),
         value_mod = value / 24,
         jday_fact = factor(jday, levels = jday_cut, 
                            labels = jday_name)) %>% 
  ggplot() +
  stat_ecdf(aes(x = value_mod,
                color = 'modelled',
                linetype = 'modelled'),
            geom = 'step', size = 1.5) +
  geom_step(data = test, aes(y = cumsum, x = Days, col = 'observed', linetype = 'observed'),
            size = 1.5) +
  coord_cartesian(xlim = c(0,51)) +
  facet_wrap(~jday_fact)

#LarsChill::convert_parameters_old_to_new(c(0,0,0,0, 0.4153e4, 0.1289e5, 0.1395e6, 0.2567e19, 0, 0,0,0))
