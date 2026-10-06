library(tidyverse)
library(chillR)
library(ggrepel)

#Rcpp::sourceCpp("src/testfun.cpp")
Rcpp::sourceCpp("src/phenoflex_pop.cpp")

# s <- read.csv('experimental/Ravensburg_hourly_temp.csv') %>% 
#   genSeasonList(years = 2022) %>% 
#   purrr::pluck(1)


kob_season <- read.csv('experimental/Ravensburg_hourly_temp.csv') %>% 
  genSeasonList(years = 2004:2022) %>% 
  setNames(2004:2022)

kob_bloom <- read.csv('experimental/Ravensburg_bloom_dates.csv') %>% 
  filter(variety == 'Topaz') %>% 
  mutate(firstbloom = lubridate::yday(firstbloom),
         fullbloom = lubridate::yday(fullbloom),
         budbreak = lubridate::yday(budbreak))

get_skewed_dist <- function(mean, sd, skew, n=100){
  
  delta <- skew / sqrt(1 + skew^2)
  
  omega <- sd /
    sqrt(1 - 2 * delta^2 / pi)
  
  xi <- mean -
    omega * delta * sqrt(2 / pi)
  
  set.seed(123)
  
  return(sn::rsn(
    n = n,
    xi = xi,
    omega = omega,
    alpha = skew
  ))
}

helper_sample_dist <- function(n, par, yc_sd, zc_sd, dist_chill, dist_heat, add_par, seed = 12345){
  if(n > 1){
    set.seed(seed)
    if(yc_sd > 0){
      if(dist_chill == 'normal'){
        yc_pop <- rnorm(n = n, mean = par[1], sd =  yc_sd)
      } else if(dist_chill == 'normal_skewed'){
        yc_pop <- get_skewed_dist(mean = par[1], sd = yc_sd, skew = add_par[1],  n = n)
      }
    }
    if(zc_sd > 0){
      if(dist_heat == 'normal'){
        zc_pop <- rnorm(n = n, mean = par[2], sd = zc_sd)
      } else if(dist_heat == 'normal_skewed'){
        zc_pop <- get_skewed_dist(mean = par[2], sd = zc_sd, skew = add_par[2],  n = n)
      }
    }
  }
  
  return(list('yc_pop' = yc_pop |> as.vector(),
              'zc_pop' = zc_pop |> as.vector()))
}

helper_run_pop_model <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100,
                                 adjust_zc = 1, basic_output = FALSE, stop_at_zc = TRUE,
                                 dist_chill = 'normal', dist_heat = 'normal', 
                                 add_distpar = NULL){
  
  #--------------#
  #draw population of yc, zc
  
  yc_pop <- par[1]
  zc_pop <- par[2]
  
  #in case there is a distribution to sample from, sample yc_pop and zc_pop
  if(n > 1){
    dist_out <- helper_sample_dist(n = n, par = par,
                       yc_sd = yc_sd, zc_sd = zc_sd, dist_chill = dist_chill,
                       dist_heat = dist_heat, add_par = add_distpar)
    
    yc_pop <- dist_out$yc_pop
    zc_pop <- dist_out$zc_pop
  }
  
  #---------------#
  #identify timepoints of cutting in temperature data
  
  i_cut <-purrr::map_int(jday_cut, function(x){
    floor(median(which(x == temp_df$JDay)))
  })  
  #c++ starts counting at zero, correct for that
  i_cut <- i_cut -1
  
  #----------------#
  #run model
  
  mod_out <- PhenoFlex_pop_slim(temp = temp_df$Temp, 
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
                basic_output = basic_output, 
                adjust_zc_forcing_exp = adjust_zc, 
                stopatzc = stop_at_zc) 
  
  #attach yc_pop, zc_pop
  mod_out$yc_pop <- yc_pop
  mod_out$zc_pop <- zc_pop
  
  return(mod_out)
  
}

helper_run_pop_model_old <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100,
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
#modify function so that it can take two model outputs

#function to prepare data ready for plotting
helper_prepare_forcing_data <- function(model_res_list, obs, 
                                        jday_cut, 
                                        jday_name, 
                                        forcing_days = 50, placeholder_fail = 9999,
                                        legend_names=NULL, annotate_performance = TRUE){
  
  # Validate inputs
  if(any(names(model_res_list) %in% c("bloomindex","x","y","z"))) model_res_list <- list(model_res_list)
  
  # Handle legend_names
  n_models <- length(model_res_list)
  if (!is.null(legend_names)) {
    if (length(legend_names) != n_models) {
      stop("'legend_names' must have the same length as 'model_res_list'")
    }
  } else {
    legend_names <- if (n_models == 1) "Modelled" else paste("Model", seq_len(n_models))
  }
  
  # Process each model run and bind together
  mod_ecdf <- purrr::map2(model_res_list, legend_names, function(model_res, leg_name) {
    
    mod_df <- model_res$exp %>%
      as.data.frame() %>%
      mutate(jday = jday_cut) %>%
      pivot_longer(cols = -jday) %>%
      mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
             value = ifelse(value == placeholder_fail, yes = 24 * (forcing_days + 1), no = value),
             value_mod = value / 24,
             jday_fact = factor(jday, levels = jday_cut,
                                labels = jday_name))
    
    purrr::map(levels(mod_df$jday_fact), function(jdf) {
      ecdf_fun <- mod_df %>%
        filter(jday_fact == jdf) %>%
        pull(value_mod) %>%
        ecdf()
      
      data.frame(cumsum_mod = ecdf_fun(seq(0, forcing_days, by = 0.1)),
                 Days        = seq(0, forcing_days, by = 0.1),
                 jday_fact   = jdf,
                 model_run   = leg_name)
    }) %>%
      bind_rows()
    
  }) %>%
    bind_rows()
  
  ssd_df <- NULL
  # SSD per model run x chilling date combination
  if(annotate_performance){
    ssd_df <- purrr::map(legend_names, function(leg_name) {
      obs[, c('Days', 'cumsum', 'jday_fact')] %>%
        merge(mod_ecdf %>% filter(model_run == leg_name),
              by = c('Days', 'jday_fact')) %>%
        group_by(jday_fact) %>%
        summarise(ssd = sum(((cumsum - cumsum_mod) * 100)^2) %>% round(digits = 0),
                  .groups = "drop") %>%
        mutate(model_run = leg_name)
    }) %>%
      bind_rows()
  }

  
  # Build curve dataset: one row per Days x jday_fact x source
  curve_df <- mod_ecdf %>%
    rename(cumsum = cumsum_mod) %>%
    mutate(source = model_run) %>%
    select(-model_run) %>%
    rbind(
      cbind(obs[, c("Days", "cumsum", "jday_fact")],
            source = "Observed")
    ) %>%
    mutate(jday_fact = factor(jday_fact, levels = jday_name),
           source    = factor(source, levels = c("Observed", legend_names))) %>%
    filter(Days <= forcing_days)
  
  return(list(curve_df = curve_df, 
              ssd_df = ssd_df))
  
}

#plot the forcing data
helper_plot_forcing_exp_test <- function(model_res_list, obs, 
                                         jday_cut, 
                                         jday_name, 
                                         forcing_days = 50, placeholder_fail = 9999,
                                    legend_names=NULL,annotate_performance = TRUE){
  
  prep_list <- helper_prepare_forcing_data(model_res_list = model_res_list, 
                                           obs = obs,
                                           jday_cut = jday_cut,
                                           jday_name = jday_name, 
                                           forcing_days = forcing_days, 
                                           placeholder_fail = placeholder_fail,
                                           legend_names = legend_names,
                                           annotate_performance = annotate_performance)
  curve_df <- prep_list$curve_df
  ssd_sf <- prep_list$ssd_df

  
  # Assign linetypes: solid for Observed, dashed variants for model runs
  n_models   <- length(legend_names)
  lty_values <- c("solid", rep(c("dashed", "dotdash", "dotted", "longdash", "twodash"),
                               length.out = n_models))
  names(lty_values) <- c("Observed", legend_names)
  
  p_forcing <- curve_df %>%
    ggplot() +
    geom_line(aes(x = Days, y = cumsum,
                  col = source, linetype = source)) +
    scale_color_discrete(name = "Data Source") +
    scale_linetype_manual(name = "Data Source", values = lty_values) +
    ylab("Share of Buds Reaching Budbreak (%)") +
    xlab("Days After Cutting Under Forcing Conditions") +
    coord_cartesian(xlim = c(0, forcing_days + 1)) +
    scale_y_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.25)) +
    facet_wrap(~jday_fact) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  if (annotate_performance) {
    # Total SSD label: sum across all model runs, or one per run
    total_ssd_label <- ssd_df %>%
      group_by(model_run) %>%
      summarise(total = sum(ssd), .groups = "drop") %>%
      mutate(label = paste0(model_run, ": ", format(total, big.mark = ",", nsmall = 0))) %>%
      pull(label) %>%
      paste(collapse = "  |  ")
    
    p_forcing <- p_forcing +
      ggtitle(paste("Forcing Experiment. Total SSD —", total_ssd_label)) +
      ggrepel::geom_label_repel(
        data = ssd_df,
        aes(x = 1, y = 0.9,
            label = paste0(model_run, "\nSSD: ", format(ssd, nsmall = 0, big.mark = ",")),
            col = model_run),
        inherit.aes = FALSE,
        direction    = "both",
        min.segment.length = Inf,
        box.padding  = 0.8,
        point.padding = 0.8,
        segment.color = "grey40",
        force         = 5,
        show.legend   = FALSE
      )
  }
  
  return(p_forcing)
}

helper_plot_forcing_exp <- function(model_res_list, obs, jday_cut, jday_name, forcing_days = 50, placeholder_fail = 9999,
                                    legend_names=NULL,annotate_performance = TRUE){
  

  # Validate inputs
  if (!is.list(model_res_list)) model_res_list <- list(model_res_list)

  # Handle legend_names
  n_models <- length(model_res_list)
  if (!is.null(legend_names)) {
    if (length(legend_names) != n_models) {
      stop("'legend_names' must have the same length as 'model_res_list'")
    }
  } else {
    legend_names <- if (n_models == 1) "Modelled" else paste("Model", seq_len(n_models))
  }
  
  #prepare modeled forcing data
  mod_df <- model_res_list$exp %>% 
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
    summarise(ssd = sum(((cumsum - cumsum_mod)*100)^2) %>% round(digits = 0))
  
  # create curve dataset used for lines
  curve_df <- mod_ecdf %>% 
    rename(cumsum = cumsum_mod) %>% 
    mutate(source = "Modelled") %>% 
    rbind(
      cbind(obs[, c("Days", "cumsum", "jday_fact")],
            source = "Observed")
    ) %>% 
    mutate(jday_fact = factor(jday_fact, levels = jday_name)) %>% 
    filter(Days <= forcing_days)
  
  p_forcing <- curve_df %>%
    ggplot() +
    geom_line(aes(x = Days, y = cumsum,
                  col = source, linetype = source),
              size = 1.5) +
      scale_color_discrete(name = 'Data Source') +
      scale_linetype_discrete(name = 'Data Source') +
      ylab('Share of buds flowering') +
      xlab('Days after cutting under forcing conditions') +
      coord_cartesian(xlim = c(0,forcing_days+1)) +
      facet_wrap(~jday_fact) +
      theme_bw() +
      theme(legend.position = 'bottom')
  
  if(annotate_performance){
    p_forcing <- p_forcing +
      ggtitle(paste('Forcing Experiment. Total Sum of Squared Differences (SSD):', format(sum(ssd_df$ssd), big.mark = ',', nsmall = 0))) +
      ggrepel::geom_label_repel(
        data = ssd_df,
        aes(x = 1, y = 0.9,
            label = paste('SSD:',
                          format(ssd, nsmall = 0, big.mark = ','))),
        inherit.aes = FALSE,
        direction = "both",
        min.segment.length = Inf,
        box.padding = 0.8,
        point.padding = 0.8,
        segment.color = "grey40",
        force = 5
      ) 
  }
  return(p_forcing)
}
  # mod_ecdf %>% 
  #   rename(cumsum = 'cumsum_mod') %>% 
  #   mutate(source = 'Modelled') %>% 
  #   rbind(cbind(obs[,c('Days', 'cumsum', 'jday_fact')], source = 'Observed')) %>% 
  #   mutate(jday_fact = factor(jday_fact, levels = jday_name)) %>% 
  #   filter(Days <= forcing_days) %>% 
  #   ggplot() +
  #   geom_line(aes(x = Days, y = cumsum, col = source, linetype = source),
  #             size = 1.5) +
  #   # geom_ribbon(data = ribbon_df, aes(ymin = ymin, ymax = ymax, x = x),
  #   #             fill = 'grey') +
  #   #annotate(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, nsmall = 1))), vjust = 1, hjust = 0) +
  #   #geom_text(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, digits = 2))), vjust = 1, hjust = 0) +
  #   #geom_text(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, nsmall = 0, big.mark = ','))), vjust = 1, hjust = 0) +
  #   ggrepel::geom_text_repel(data = ssd_df, aes(x = 0, y = 1, label = paste('SSD:' ,format(ssd, nsmall = 0, big.mark = ','))), vjust = 1, hjust = 0) +
  #   ggtitle(paste('Forcing Experiment. Total Sum of Squared Differences (SSD):', format(sum(ssd_df$ssd), big.mark = ',', nsmall = 0))) + 
  #   scale_color_discrete(name = 'Data Source') +
  #   scale_linetype_discrete(name = 'Data Source') +
  #   ylab('Share of buds flowering') + 
  #   xlab('Days after cutting under forcing conditions') +
  #   coord_cartesian(xlim = c(0,forcing_days+1)) +
  #   facet_wrap(~jday_fact) +
  #   theme_bw() +
  #   theme(legend.position = 'bottom')
    
  
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


# #---------------------------------#
# #OBSERVED
# 
# exp_obs <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
#                               sheet = 'T_2022_term+spur') %>% 
#   group_by(Data) %>% 
#   mutate(cumsum = cumsum(Bubble_size)) %>% 
#   mutate(Date = lubridate::parse_date_time(Data, orders = 'dmy'),
#          yday = lubridate::yday(Date),
#          h_after_cut = Days*24) %>% 
#   filter(yday >= 300 | yday < 55)
# 
# 
# #-----------------------------------#
# #cutting experiments
# 
# #the function expects to get already the position at which the forcing experiment should take place
# jday_cut <- unique(exp_obs$yday)
# jday_name <- lubridate::stamp("Nov 03", orders = '%b %d')(unique(exp_obs$Date))
# #jday_name <- c('Nov 03', 'Nov 17', 'Dec 01', 'Dec 15', 'Dec 29', 'Jan 12', 'Jan 26', 'Feb 09', 'Feb23')
# i_cut <-purrr::map_int(jday_cut, function(x){
#   floor(median(which(x == s$JDay)))
# })  
# #c++ starts counting at zero, correct for that
# i_cut <- i_cut -1
# 
# #add missing days to antons table
# miss_df <- exp_obs %>% 
#   group_by(Date, yday) %>% 
#   summarise(max = max(cumsum),
#             miss = 1 - max) %>% 
#   ungroup() %>% 
#   mutate(h_after_cut = 51 *24,
#          Bubble_size = miss,
#          Days = 51,
#          cumsum = 1) %>% 
#   select(Days, Date, yday, h_after_cut, Bubble_size, cumsum)
# 
# exp_obs <- exp_obs %>% 
#   select(Days, Date, yday, h_after_cut, Bubble_size, cumsum) %>% 
#   rbind(miss_df) %>%  
#   mutate(jday_fact = factor(yday, levels = jday_cut, 
#                             labels = jday_name)) 
# 
# temp_df = s

helper_prepare_obs_data <- function(sheet = 'T_2022_term+spur',
                                    start_yday = 300,
                                    end_yday = 55){
  
  exp_obs <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                                sheet = sheet) %>% 
    group_by(Data) %>% 
    mutate(cumsum = cumsum(Bubble_size)) %>% 
    mutate(Date = lubridate::parse_date_time(Data, orders = 'dmy'),
           yday = lubridate::yday(Date),
           h_after_cut = Days*24) 
  
  #subsetting
  
  if(start_yday >end_yday){
    exp_obs <- exp_obs[exp_obs$yday >= start_yday | exp_obs$yday <= end_yday, ]
  } else {
    exp_obs <- exp_obs[exp_obs$yday >= start_yday & exp_obs$yday <= end_yday, ]
  }
  
  if(grepl(pattern = '2020', x = sheet)) temp_file <- 'experimental/hohenheim_aug19-jul20.csv'
  if(grepl(pattern = '2022', x = sheet)) temp_file <- 'experimental/hohenheim_aug21-jul22.csv'
  #read temperature data
  s <- read.csv(temp_file, sep = ';', dec = ',') %>% 
    mutate(Date = lubridate::dmy(Tag),
           Hour = lubridate::hm(Stunde) %>% hour(),
           JDay = lubridate::yday(Date),
           Year = lubridate::year(Date),
           Temp = AVG_TA200,
    ) %>% 
    select(Temp, JDay, Year)
  
  #-----------------------------------#
  #cutting experiments
  
  #the function expects to get already the position at which the forcing experiment should take place
  jday_cut <- unique(exp_obs$yday)
  jday_name <- lubridate::stamp("03 Nov", orders = '%d %b', locale = 'C')(unique(exp_obs$Date))
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
  
  return(list(
    exp_obs = exp_obs,
    jday_cut = jday_cut,
    jday_name = jday_name,
    i_cut = i_cut,
    temp_df = s
    
  ) )
}


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

helper_plot_flowering <- function(bloom_df, 
                                  obs_df, 
                                  raw_obs = TRUE,
                                  year_select = 'all', 
                                  placeholder_fail = 9999, 
                                  annotate_performance = TRUE){
  
  years <- as.numeric(year_select)
  unique_years <- unique(bloom_df$year)
  val_years <- c(2007, 2009, 2010, 2016)
  
  if(length(year_select) == 1){
    if(year_select == 'all') years <- unique_years
    if(year_select == 'calibration') years <- unique_years[(unique_years %in% val_years) == FALSE]
    if(year_select == 'validation') years <- val_years
  }
  
  strip_fills <- ifelse(years %in% val_years,
                        "#F4CCCC",   # validation (red-ish)
                        "#D9EAD3")   # calibration (green-ish)
  
  
  if(raw_obs){
    obs_df <- obs_df %>% 
      filter(year %in% years) %>% 
      select(year, firstbloom, fullbloom) %>% 
      pivot_longer(cols = -year) %>% 
      mutate(share = ifelse(name == 'firstbloom', yes = 0.1, no = 0.5))
  }

  
  bloom_df$value_mod <- ifelse(bloom_df$value == placeholder_fail, yes = NA, no = bloom_df$value)
  
  yday_x_scale <- (floor(min(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE))):(ceiling(max(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE)))
  yday_x_scale_label <- as.Date(yday_x_scale,
                                origin = '2021-12-31') %>% 
    format("%b %d")
  
  yday_x_scale_sub <- yday_x_scale[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  yday_x_scale_label_sub <- yday_x_scale_label[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  

  
  p_flower <- bloom_df %>% 
    filter(year %in% years) %>% 
    ggplot() +
    stat_ecdf(aes(x = value,
                  color = 'Modelled'),
              geom = 'step', size = 1.5) +
    ylab('Share of buds flowering (%)') +
    xlab('Date') +
    scale_x_continuous(breaks = yday_x_scale_sub, 
                       labels = yday_x_scale_label_sub,
                       limits = c(min(yday_x_scale), 
                                  max(yday_x_scale))) +
    scale_y_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.25)) +
    geom_point(data = obs_df, aes(x = value, y = share, 
                                  col = 'Observed')) +
    ggh4x::facet_wrap2(~year, strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_fills,
                                                                                 colour = 'black'))) +
    scale_color_manual(values = c('steelblue', 'firebrick', '#E69F00')) +
    scale_shape_manual(values = 2) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(annotate_performance){
    
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
    
    p_flower <- p_flower  +
      annotate(geom = 'text', x = min(yday_x_scale), y = 1, label = 'Error',
               vjust = 1, hjust = 0) +
      geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.8, 
                                           label = paste('F1:', format(round(firstbloom, digits = 1), nsmall = 1))),
                vjust = 1, hjust = 0) +
      geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.6, 
                                           label = paste('F2:', format(round(fullbloom, digits = 1), nsmall = 1))),
                vjust = 1, hjust = 0)+
      ggtitle(paste('Predicted Bloom KOB. RMSE Firstbloom (F1):', format(rmse_df$firstbloom, nsmall = 1),
                    'RMSE Fullbloom (F2):', format(rmse_df$fullbloom, nsmall = 1) ))
      
  }
  return(p_flower)  
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

helper_plot_budbreak_orchard <- function(bloom_df, 
                                         obs_df, 
                                         year_select = 'all', 
                                         placeholder_fail = 9999, 
                                         annotate_performance = TRUE){
  years <- as.numeric(year_select)
  unique_years <- unique(bloom_df$year)
  val_years <- c(2007, 2009, 2010, 2016)
  
  if(length(year_select) == 1){
    if(year_select == 'all') years <- unique_years
    if(year_select == 'calibration') years <- unique_years[(unique_years %in% val_years) == FALSE]
    if(year_select == 'validation') years <- val_years
  }
  
  #subset observations
  obs_df <- obs_df %>% 
    filter(year %in% years) %>% 
    select(year, budbreak) |> 
    mutate(year = as.character(year))
  
  
  bloom_df$value_mod <- ifelse(bloom_df$value == placeholder_fail, yes = NA, no = bloom_df$value)
  
  yday_y_scale <- (floor(min(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE))):(ceiling(max(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE)))
  yday_y_scale_label <- as.Date(yday_y_scale,
                                origin = '2021-12-31') %>% 
    format("%b %d")
  
  yday_y_scale_sub <- yday_y_scale[seq.int(from = 1, to = length(yday_y_scale), length.out = 10)]
  yday_y_scale_label_sub <- yday_y_scale_label[seq.int(from = 1, to = length(yday_y_scale), length.out = 10)]
  
  
  ggplot() +
    geom_point(data = obs_df, aes(x = year, y = budbreak, 
                                  col = 'Observed'))
  
  
  p_budbreak <- ggplot() +
    geom_boxplot(data = bloom_df, aes(x = year, y = value_mod,
                                      col = 'Predicted'),
                 fill = 'grey70') +
    ylab('Date of budbreak (%)') +
    xlab('Year') +
    # scale_y_continuous(breaks = yday_y_scale_sub, 
    #                    labels = yday_y_scale_label_sub,
    #                    limits = c(min(yday_y_scale), 
    #                               max(yday_y_scale))) +
    geom_point(data = obs_df, aes(x = year, y = budbreak, 
                                  col = 'Observed'),
               shape = 25, size = 2, fill = 'firebrick') +
    scale_color_manual(values = c('firebrick', 'black')) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(annotate_performance){
    
    #merge predicted and observed
    rmse_df <- bloom_df %>% 
      filter(year %in% years) %>% 
      group_by(year) |> 
      summarise(budbreak = median(value)) |> 
      pivot_longer(cols = -year, values_to = 'pred') %>% 
      merge(obs_df, by = c('year')) %>% 
      summarise(rmse = RMSEP(predicted = pred, observed = budbreak) %>% round(digits = 1)) %>% 
      ungroup() 
    
    p_budbreak <- p_budbreak  +
       ggtitle(paste('Predicted Budbreak KOB. RMSE (days):', format(rmse_df$rmse, nsmall = 1)))
    
  }
  return(p_budbreak)
  }  

# 
# input <- list(
#   yc = 61.1571429,
#   yc_sd = 10.8367347,
#   zc = 185.5714286,
#   zc_sd = 10,
#   s1 = 0.1662224,
#   theta_star = 280.04901-273.15,
#   theta_c = 286.05151-273.15,
#   tau = 27.88297,
#   pi_c = 30.41736,
#   Tf = 4,
#   slope = 1.6,
#   Tb = 4,
#   Tu = 25,
#   Tc = 36
# ) 
# 
# #topaz first flowering caspersen 2025
# input <- list(
#   yc = 40.0336321,
#   yc_sd = 5,
#   zc = 181.2843981,
#   zc_sd = 10,
#   s1 = 0.1473177,
#   theta_star = 279-273.15,
#   theta_c = 285.6807267-273.15,
#   tau = 34.2354385,
#   pi_c = 32.6464321,
#   Tf = 8.0650228,
#   slope = 1.9426370,
#   Tb = 5.9122156,
#   Tu = 21.0231964,
#   Tc = 36
# )
# 
# par <- c(input$yc, input$zc, input$s1, input$Tu, input$theta_star + 273.15, input$theta_c + 273.15, input$tau, input$pi_c, input$Tf, input$Tc, input$Tb, input$slope) %>% 
#   LarsChill::convert_parameters()
# 
# #run model
# pop_out <- helper_run_pop_model(par = par, 
#                                 yc_sd = input$yc_sd, 
#                                 zc_sd = input$zc_sd, 
#                                 jday_cut = jday_cut, 
#                                 temp_df = s, 
#                                 n = 100)
# 
# bloom_list <- purrr::map(kob_season, function(s1){
#   bloom <- helper_run_pop_model(par = par, yc_sd = input$yc_sd, zc_sd =  input$zc_sd, jday_cut = jday_cut, temp_df = s1, n = 100,
#                                 basic_output = TRUE) %>% 
#     purrr::pluck('bloomindex') %>% 
#     purrr::map_dbl(helper_bloomint_to_jday, x = s) %>% 
#     return()
# })
# pop_bloom <- do.call(cbind, bloom_list) %>% 
#   as.data.frame() %>% 
#   setNames(names(kob_season)) %>% 
#   pivot_longer(cols = everything(), names_to = 'year')
# 
# p_flower <- helper_plot_flowering(bloom_df = pop_bloom[pop_bloom$year == 2022, ], obs_df = kob_bloom, year_select = 2022)
# 
# 
# lubridate::yday('2022-04-25')
# 
# yday_x_scale_label <- as.Date(103:115) %>% 
#   format("%b %d")
# 
# p_flower +
#   theme_bw(base_size = 15) +
#   scale_x_continuous(breaks = 103:115, labels = yday_x_scale_label) +
#   theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# 
#     
# ggsave('experimental/example_bloom-curve_kob-2022.jpeg', height = 15, width = 20, units = 'cm', device = 'jpeg')
#   
# 
# 
# p_force <- helper_plot_forcing_exp(model_res = pop_out, obs = exp_obs, jday_cut = jday_cut, jday_name = jday_name)
# p_force +
#   theme_bw(base_size = 15) +
#   theme(legend.position = 'bottom')
# ggsave('experimental/example_forcing-2022_topaz.jpeg', height = 15, width = 20, units = 'cm', device = 'jpeg')
# 
# 
# 
# cbind(obs[,c('Days', 'cumsum', 'jday_fact')], source = 'Observed') %>% 
#   mutate(jday_fact = factor(jday_fact, levels = jday_name)) %>% 
#   filter(Days <= forcing_days) %>% 
#   ggplot() +
#   geom_line(aes(x = Days, y = cumsum, col = source, linetype = source),
#             size = 1.5) +
#   scale_color_manual(name = 'Data Source', values =  '#00BFC4') +
#   scale_linetype_discrete(name = 'Data Source') +
#   ylab('Share of buds flowering') + 
#   xlab('Days after cutting under forcing conditions') +
#   coord_cartesian(xlim = c(0,forcing_days+1)) +
#   facet_wrap(~jday_fact) +
#   theme_bw(base_size = 15) +
#   theme(legend.position = 'bottom')
# ggsave('experimental/example_forcing-2022_topaz_blank.jpeg', height = 15, width = 20, units = 'cm', device = 'jpeg')
# 
# 
# 
# 
# input <- list(
#   yc = 40.0336321,
#   yc_sd = 0,
#   zc = 181.2843981,
#   zc_sd = 0,
#   s1 = 0.1473177,
#   theta_star = 279-273.15,
#   theta_c = 285.6807267-273.15,
#   tau = 34.2354385,
#   pi_c = 32.6464321,
#   Tf = 8.0650228,
#   slope = 1.9426370,
#   Tb = 5.9122156,
#   Tu = 21.0231964,
#   Tc = 36
# )
# 
# par <- c(input$yc, input$zc, input$s1, input$Tu, input$theta_star + 273.15, input$theta_c + 273.15, input$tau, input$pi_c, input$Tf, input$Tc, input$Tb, input$slope) %>% 
#   LarsChill::convert_parameters()
# 
# #run model
# pop_out <- helper_run_pop_model(par = par, 
#                                 yc_sd = input$yc_sd, 
#                                 zc_sd = input$zc_sd, 
#                                 jday_cut = jday_cut, 
#                                 temp_df = s, 
#                                 n = 100)
# 
# bloom_list <- purrr::map(kob_season, function(s1){
#   bloom <- helper_run_pop_model(par = par, yc_sd = input$yc_sd, zc_sd =  input$zc_sd, jday_cut = jday_cut, temp_df = s1, n = 100,
#                                 basic_output = TRUE) %>% 
#     purrr::pluck('bloomindex') %>% 
#     purrr::map_dbl(helper_bloomint_to_jday, x = s) %>% 
#     return()
# })
# pop_bloom <- do.call(cbind, bloom_list) %>% 
#   as.data.frame() %>% 
#   setNames(names(kob_season)) %>% 
#   pivot_longer(cols = everything(), names_to = 'year')
# 
# p_flower <- helper_plot_flowering(bloom_df = pop_bloom[pop_bloom$year == 2022, ], obs_df = kob_bloom, year_select = 2022)
# 
# p_flower +
#   theme_bw(base_size = 15) +
#   scale_x_continuous(breaks = 103:115, labels = yday_x_scale_label) +
#   theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# ggsave('experimental/example_bloom-curve_kob-2022_no-sd.jpeg', height = 15, width = 20, units = 'cm', device = 'jpeg')
# 
# p_force <- helper_plot_forcing_exp(model_res = pop_out, obs = exp_obs, jday_cut = jday_cut, jday_name = jday_name)
# p_force +
#   theme_bw(base_size = 15) +
#   theme(legend.position = 'bottom')
# ggsave('experimental/example_forcing-2022_topaz_no-sd.jpeg', height = 15, width = 20, units = 'cm', device = 'jpeg')

