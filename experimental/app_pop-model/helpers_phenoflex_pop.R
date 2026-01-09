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

kob_season <- read.csv('experimental/Ravensburg_hourly_temp.csv') %>% 
  genSeasonList(years = 2004:2022) %>% 
  setNames(2004:2022)

kob_bloom <- read.csv('experimental/Ravensburg_bloom_dates.csv') %>% 
  filter(variety == 'Topaz') %>% 
  mutate(firstbloom = lubridate::yday(firstbloom),
         fullbloom = lubridate::yday(fullbloom))

helper_run_pop_model <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100,
                                 basic_output = FALSE){
  
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
                basic_output = basic_output) %>% 
    return()
  
}

helper_plot_forcing_exp <- function(model_res, obs, jday_cut, jday_name, forcing_days = 50, placeholder_fail = 9999){
  
  # #hide the added fake observation, to force cumulative curve to reach 1.0 at day after maximum forcing days
  # ribbon_df <- data.frame(x = (forcing_days):(forcing_days + 10),
  #                         ymin = -Inf,
  #                         ymax = Inf)
  
  #prepare modeled forcing data
  mod_df <- model_res$exp %>% 
    as.data.frame() %>% 
    mutate(jday = jday_cut) %>% 
    pivot_longer(cols = -jday) %>%
    mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
           value = ifelse(value == placeholder_fail, yes = 24* (forcing_days+1), no = value),
           value_mod = value / 24,
           jday_fact = factor(jday, levels = jday_cut, 
                              labels = jday_name)) 
  
  #calculate ecdf based on the returned days when budbreak is reached
  mod_ecdf <- purrr::map(levels(mod_df$jday_fact), function(jdf){
    #set ecdf function
    ecdf_fun <- mod_df %>% 
      filter(jday_fact == jdf) %>% 
      pull(value_mod) %>% 
      ecdf()
    
    return(data.frame(cumsum_mod = ecdf_fun(seq(0, forcing_days, by = 0.1)),
                      Days = seq(0, forcing_days, by = 0.1),
                      jday_fact = jdf))
  }) %>% 
    bind_rows()
  
  #bring observed and modelled ecdf together
  ssd_df <- obs[,c('Days', 'cumsum', 'jday_fact')] %>% 
    merge(mod_ecdf, by = c('Days', 'jday_fact')) %>% 
    group_by(jday_fact) %>% 
    summarise(ssd = sum((cumsum - cumsum_mod)^2) %>% round(digits = 1))
  
  mod_ecdf %>% 
    rename(cumsum = 'cumsum_mod') %>% 
    mutate(source = 'Modelled') %>% 
    rbind(cbind(obs[,c('Days', 'cumsum', 'jday_fact')], source = 'Observed')) %>% 
    mutate(jday_fact = factor(jday_fact, levels = jday_name)) %>% 
    filter(Days <= forcing_days) %>% 
    ggplot() +
    geom_line(aes(x = Days, y = cumsum, col = source, linetype = source),
              size = 1.5) +
    # geom_ribbon(data = ribbon_df, aes(ymin = ymin, ymax = ymax, x = x),
    #             fill = 'grey') +
    #annotate(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, nsmall = 1))), vjust = 1, hjust = 0) +
    geom_text(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, digits = 2))), vjust = 1, hjust = 0) +
    ggtitle(paste('Forcing Experiment. Total Sum of Squared Differences (SSD):', format(sum(ssd_df$ssd), nsmall = 1))) + 
    scale_color_discrete(name = 'Data Source') +
    scale_linetype_discrete(name = 'Data Source') +
    ylab('Share of buds flowering') + 
    xlab('Days after cutting under forcing conditions') +
    coord_cartesian(xlim = c(0,forcing_days+1)) +
    facet_wrap(~jday_fact) +
    theme_bw() +
    theme(legend.position = 'bottom')
    
  
  # ggplot(mod_df) +
  #   stat_ecdf(aes(x = value_mod,
  #                 color = 'Modelled',
  #                 linetype = 'Modelled'),
  #             geom = 'step', size = 1.5) +
  #   geom_step(data = obs, aes(y = cumsum, x = Days, col = 'Observed', linetype = 'Observed'),
  #             size = 1.5) +
  #   geom_ribbon(data = ribbon_df, aes(ymin = ymin, ymax = ymax, x = x),
  #               fill = 'grey') +
  #   scale_color_discrete(name = 'Data Source') +
  #   scale_linetype_discrete(name = 'Data Source') +
  #   ylab('Share of buds flowering') + 
  #   xlab('Days after cutting under forcing conditions') +
  #   coord_cartesian(xlim = c(0,51)) +
  #   facet_wrap(~jday_fact) +
  #   theme_bw() +
  #   theme(legend.position = 'bottom')
}

#---------------------------------#
#OBSERVED

exp_obs <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                              sheet = 'T_2022_term+spur') %>% 
  group_by(Data) %>% 
  mutate(cumsum = cumsum(Bubble_size)) %>% 
  mutate(Date = lubridate::parse_date_time(Data, orders = 'dmy'),
         yday = lubridate::yday(Date),
         h_after_cut = Days*24) %>% 
  filter(yday >= 300 | yday < 55)

#-----------------------------------#
#cutting experiments

#the function expects to get already the position at which the forcing experiment should take place
jday_cut <- unique(exp_obs$yday)
jday_name <- lubridate::stamp("Nov 03", orders = '%b %d')(unique(exp_obs$Date))
#jday_name <- c('Nov 03', 'Nov 17', 'Dec 01', 'Dec 15', 'Dec 29', 'Jan 12', 'Jan 26', 'Feb 09', 'Feb23')
i_cut <-purrr::map_int(jday_cut, function(x){
  floor(median(which(x == s$JDay)))
})  
#c++ starts counting at zero, correct for that
i_cut <- i_cut -1

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

exp_obs <- exp_obs %>% 
  select(Days, Date, yday, h_after_cut, Bubble_size, cumsum) %>% 
  rbind(miss_df) %>%  
  mutate(jday_fact = factor(yday, levels = jday_cut, 
                            labels = jday_name)) 


#-----------------------------------#



#----------------------------------#
#prepare flower plot

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

helper_plot_flowering_old <- function(model_res, temp_df, obs = c('2022-04-14', '2022-04-18')){
  bloom_df <- data.frame(pheno = purrr::map_dbl(model_res$bloomindex, helper_bloomint_to_jday, x = temp_df))
  
  yday_obs_10 <- lubridate::ymd(obs[1]) %>% lubridate::yday()
  yday_obs_50 <- lubridate::ymd(obs[2]) %>% lubridate::yday()
  yday_x_scale <- (floor(min(c(bloom_df$pheno, yday_obs_10, yday_obs_50)))):(ceiling(max(c(bloom_df$pheno, yday_obs_10, yday_obs_50))))
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
    ylab('Share of buds flowering') +
    xlab('Date') +
    scale_x_continuous(breaks = yday_x_scale, labels = yday_x_scale_label) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

helper_plot_flowering <- function(bloom_df, obs_df, year_select = 'all', placeholder_fail = 9999){
  
  years <- as.numeric(year_select)
  if(year_select == 'all') years <- unique(bloom_df$year)
  
  obs_df <- obs_df %>% 
    filter(year %in% years) %>% 
    select(year, firstbloom, fullbloom) %>% 
    pivot_longer(cols = -year) %>% 
    mutate(share = ifelse(name == 'firstbloom', yes = 0.1, no = 0.5))
  
  bloom_df$value_mod <- ifelse(bloom_df$value == placeholder_fail, yes = NA, no = bloom_df$value)
  
  yday_x_scale <- (floor(min(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE))):(ceiling(max(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE)))
  yday_x_scale_label <- as.Date(yday_x_scale,
                                origin = '2021-12-31') %>% 
    format("%b %d")
  
  yday_x_scale_sub <- yday_x_scale[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  yday_x_scale_label_sub <- yday_x_scale_label[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  
  #merge predicted and observed
  rmse_df <- bloom_df %>% 
    filter(year %in% years) %>% 
    group_by(year) %>% 
    summarise(firstbloom = quantile(value, 0.1),
              fullbloom = quantile(value, 0.5)) %>% 
    pivot_longer(cols = -year, values_to = 'pred') %>% 
    merge(obs_df, by = c('year', 'name')) %>% 
    group_by(name) %>% 
    summarise(rmse = RMSEP(predicted = pred, observed = value) %>% round(digits = 1)) %>% 
    ungroup() %>% 
    pivot_wider(values_from = rmse)
  
  #calculate difference of predicted and observed
  performance_df <- bloom_df %>% 
    filter(year %in% years) %>% 
    group_by(year) %>% 
    summarise(firstbloom = quantile(value, 0.1),
              fullbloom = quantile(value, 0.5)) %>% 
    pivot_longer(cols = -year, values_to = 'pred') %>% 
    merge(obs_df, by = c('year', 'name')) %>% 
    mutate(diff = pred - value) %>% 
    select(year, name, diff) %>% 
    pivot_wider(names_from = name, values_from = diff)
  
  bloom_df %>% 
    filter(year %in% years) %>% 
    ggplot() +
    stat_ecdf(aes(x = value,
                  color = 'Modelled'),
              geom = 'step', size = 1.5) +
    scale_color_discrete(name = 'Data Source') +
    ylab('Share of buds flowering') +
    xlab('Date') +
    scale_x_continuous(breaks = yday_x_scale_sub, labels = yday_x_scale_label_sub,
                       limits = c(min(yday_x_scale), max(yday_x_scale))) +
    geom_point(data = obs_df, aes(x = value, y = share, col = 'Observed'),
               size = 1.5, shape = 25, stroke = 2, fill = 'white') +
    annotate(geom = 'text', x = min(yday_x_scale), y = 1, label = 'Error',
             vjust = 1, hjust = 0) +
    geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.8, 
                                         label = paste('F1:', format(round(firstbloom, digits = 1), nsmall = 1))),
              vjust = 1, hjust = 0) +
    geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.6, 
                                         label = paste('F2:', format(round(fullbloom, digits = 1), nsmall = 1))),
              vjust = 1, hjust = 0)+
    facet_wrap(~year) +
    ggtitle(paste('Predicted Bloom KOB. RMSE Firstbloom (F1):', format(rmse_df$firstbloom, nsmall = 1),
            'RMSE Fullbloom (F2):', format(rmse_df$fullbloom, nsmall = 1) ))+
    theme_bw() +
    theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
}

helper_combined_response_plot <- function(par, temp_values){
  
  response_df <- LarsChill::get_temp_response_df(par, temp_values)
  
  max_chill <- max(response_df$Chill_response)
  
  response_df %>% 
    mutate(Heat_response = Heat_response * max_chill) %>% 
    pivot_longer(cols = c('Chill_response', 'Heat_response')) %>% 
    mutate(name_plot = factor(name, 
                              levels = c('Chill_response', 'Heat_response'),
                              labels = c('Chill Response', 'Heat Response'))) %>% 
    ggplot(aes(x = Temperature, y = value, group = name_plot, col = name_plot)) +
    geom_line(aes(linetype = name_plot), show.legend = FALSE,
              size = 1.5) +
    scale_y_continuous(
      
      # Features of the first axis
      name = "Chill Response",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis( transform=~./40, name="Heat Response")
    ) +
    theme_bw(base_size = 15) +
    scale_color_manual(values = c('#377eb8',  '#e41a1c')) +
    scale_linetype_manual(values = c('dashed', 'solid')) +
    theme(
      # Primary Y-axis (left)
      axis.text.y.left = element_text(color = '#377eb8'),
      axis.title.y.left = element_text(color = '#377eb8'),
      axis.line.y.left = element_line(color = '#377eb8'),
      
      # Secondary Y-axis (right)
      axis.text.y.right = element_text(color = '#e41a1c'),
      axis.title.y.right = element_text(color = '#e41a1c'),
      axis.line.y.right = element_line(color = '#e41a1c'),
      legend.position = 'none'
    ) 
  
  #rescale heat
}
