library(shiny)
library(bslib)
library(patchwork)
options(shiny.sanitize.errors = FALSE)

source('experimental/app_pop-model/helpers_phenoflex_pop.R')

default_params <- list(
  yc = 75,
  yc_sd = 5,
  zc = 181,
  zc_sd = 10,
  s1 = 0.1473177,
  theta_star = 280.04901-273.15,
  theta_c = 286.05151-273.15,
  tau = 27.88297,
  pi_c = 30.41736,
  Tf = 4,
  slope = 1.6,
  Tb = 4,
  Tu = 25,
  Tc = 36
)

kob_study_par <- list(
  yc = 40.0336321,
  yc_sd = 5,
  zc = 181.2843981,
  zc_sd = 10,
  s1 = 0.1473177,
  theta_star = 279-273.15,
  theta_c = 285.6807267-273.15,
  tau = 34.2354385,
  pi_c = 32.6464321,
  Tf = 8.0650228,
  slope = 1.9426370,
  Tb = 5.9122156,
  Tu = 21.0231964,
  Tc = 36
)

# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  # App title ----
  title = "PhenoFlex Population Model - Forcing Experiment",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    
    # actionButton(
    #   inputId = "run_model",
    #   label = "Run population model",
    #   class = "btn-primary"
    # ),
    # radioButtons(
    #   inputId = "plot_type",
    #   label = "Choose plot type",
    #   choices = c("Forcing Experiment" = "forcing",
    #               "Flowering" = "flower",
    #               'Both' = 'forcing_and_flower'),
    #   selected = "forcing"
    # ),
    # selectInput(
    #   inputId = "plot_type",
    #   label = "Choose plot type",
    #   choices = c("Forcing Experiment" = "forcing",
    #               "Flowering" = "flower",
    #               'Forcing and Flowering' = 'forcing_and_flower',
    #               'Temperature Response' = ''),
    #   selected = "forcing"
    # ),
    checkboxInput(
      inputId = "plot_forcing",
      label = "Show Forcing Plot",
      value = TRUE
    ),
    checkboxInput(
      inputId = "plot_flowering",
      label = "Show Flowering Plot",
      value = FALSE
    ),
    checkboxInput(
      inputId = "plot_tempresponse",
      label = "Show Temperature Response Plot",
      value = FALSE
    ),
    numericInput(
      inputId = "yc",
      label = "yc\n(mean chill requirement)",
      min = 10,
      max = 80,
      value = 40.0336321,
    ),
    numericInput(
      inputId = "yc_sd",
      label = "yc sd\n(standard deviation of chill requirement)",
      min = 0,
      max = 20,
      value = 5
    ),
    numericInput(
      inputId = "zc",
      label = "zc\n(mean of heat requirement)",
      min = 100,
      max = 500,
      value = 181.2843981
    ),
    numericInput(
      inputId = "zc_sd",
      label = "zc sd\n(standard deviation of heat requirement)",
      min = 0,
      max = 20,
      value = 10
    ),
    numericInput(
      inputId = "s1",
      label = "s1\n(slope of transition function between chill and heat accumulation)",
      min = 0.01,
      max = 1.5,
      value = 0.1473177
    ),
    numericInput(
      inputId = "theta_star",
      label = "theta_star\n(optimal temperature in °C for chill accumulation)",
      min = 5,
      max = 8,
      value = 279-273.15
    ),
    numericInput(
      inputId = "theta_c",
      label = "theta_c\n(critical temperature in K°C for chill accumulation)",
      min = 12,
      max = 15,
      value = 285.6807267-273.15
    ),
    numericInput(
      inputId = "tau",
      label = "tau\n(time interval for chill accumulation under optimal conditions)",
      min = 16,
      max = 48,
      value = 34.2354385
    ),
    numericInput(
      inputId = "pi_c",
      label = "pi_c\n(time interval leading to chill negation)",
      min = 24,
      max = 40,
      value = 32.6464321
    ),
    numericInput(
      inputId = "Tf",
      label = "Tf\n(temperature for transition from PDBF to DBF)",
      min = 0,
      max = 10,
      value = 8.0650228
    ),
    numericInput(
      inputId = "slope",
      label = "slope\n(slope of transition function from PDBF to DBF)",
      min = 0.1,
      max = 15,
      value = 1.9426370
    ),
    numericInput(
      inputId = "Tb",
      label = "Tb\n(base temperature heat accumulation)",
      min = 0,
      max = 10,
      value = 5.9122156
    ),
    numericInput(
      inputId = "Tu",
      label = "Tu\n(optimal temperature heat accumulation)",
      min = 20,
      max = 35,
      value = 21.0231964
    ),
    numericInput(
      inputId = "Tc",
      label = "Tc\n(critical temperature heat accumulation)",
      min = 30,
      max = 40,
      value = 36
    ),
    actionButton(
      inputId = "default_params",
      label = "Default parameters",
      class = "btn-secondary"
    ),
    actionButton(
      inputId = "kob_params",
      label = "KOB study parameters",
      class = "btn-secondary"
    ),
    selectizeInput(
      inputId = "years_bloom",
      label = "Choose years for bloom prediction",
      choices = c('all', as.character(2004:2022)),
      selected = 'all',
      multiple = TRUE
    )
  ),
  # Output: Histogram ----
  plotOutput(outputId = "pop_plot")
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  # Run PhenoFlex Population Model for 2022 season in Bavendorf
  # Model returns results of Forcing Experiment together with observation of experiment
  #message("Starting model run")
  
  # pop_out <- eventReactive(input$run_model,{
  #   
  #   
  #   par <- c(input$yc, input$zc, input$s1, input$Tu, input$theta_star + 273.15, input$theta_c + 273.15, input$tau, input$pi_c, input$Tf, input$Tc, input$Tb, input$slope) %>% 
  #   LarsChill::convert_parameters()
  #   
  #   helper_run_pop_model(par = par, yc_sd = input$yc_sd, zc_sd = input$zc_sd, jday_cut = jday_cut, temp_df = s, n = 100)
  #   
  # })
  
  observeEvent(input$default_params, {
    
    updateNumericInput(session, "yc", value = default_params$yc)
    updateNumericInput(session, "yc_sd", value = default_params$yc_sd)
    updateNumericInput(session, "zc", value = default_params$zc)
    updateNumericInput(session, "zc_sd", value = default_params$zc_sd)
    updateNumericInput(session, "s1", value = default_params$s1)
    updateNumericInput(session, "theta_star", value = default_params$theta_star)
    updateNumericInput(session, "theta_c", value = default_params$theta_c)
    updateNumericInput(session, "tau", value = default_params$tau)
    updateNumericInput(session, "pi_c", value = default_params$pi_c)
    updateNumericInput(session, "Tf", value = default_params$Tf)
    updateNumericInput(session, "slope", value = default_params$slope)
    updateNumericInput(session, "Tb", value = default_params$Tb)
    updateNumericInput(session, "Tu", value = default_params$Tu)
    updateNumericInput(session, "Tc", value = default_params$Tc)
    
  })
  
  observeEvent(input$kob_params, {
    
    updateNumericInput(session, "yc", value = kob_study_par$yc)
    updateNumericInput(session, "yc_sd", value = kob_study_par$yc_sd)
    updateNumericInput(session, "zc", value = kob_study_par$zc)
    updateNumericInput(session, "zc_sd", value = kob_study_par$zc_sd)
    updateNumericInput(session, "s1", value = kob_study_par$s1)
    updateNumericInput(session, "theta_star", value = kob_study_par$theta_star)
    updateNumericInput(session, "theta_c", value = kob_study_par$theta_c)
    updateNumericInput(session, "tau", value = kob_study_par$tau)
    updateNumericInput(session, "pi_c", value = kob_study_par$pi_c)
    updateNumericInput(session, "Tf", value = kob_study_par$Tf)
    updateNumericInput(session, "slope", value = kob_study_par$slope)
    updateNumericInput(session, "Tb", value = kob_study_par$Tb)
    updateNumericInput(session, "Tu", value = kob_study_par$Tu)
    updateNumericInput(session, "Tc", value = kob_study_par$Tc)
    
  })
  
  par <- reactive({
    c(input$yc, input$zc, input$s1, input$Tu, input$theta_star + 273.15, input$theta_c + 273.15, input$tau, input$pi_c, input$Tf, input$Tc, input$Tb, input$slope) %>% 
      LarsChill::convert_parameters()
  })
  
  pop_out <- reactive({
    
    helper_run_pop_model(par = par(), yc_sd = input$yc_sd, zc_sd = input$zc_sd, jday_cut = jday_cut, temp_df = s, n = 100)
    
  })
  
  #calculation of the flowerings in KOB. Only run when the checkbox is checked
  pop_bloom <- reactive({
    
    req(input$plot_flowering)
    
    # test <- helper_run_pop_model(par = par(), yc_sd = input$yc_sd, zc_sd = input$zc_sd, jday_cut = jday_cut, temp_df = s, n = 100,
    #                      basic_output = TRUE)
    
    #calculate population of bloom for kob each season
    bloom_list <- purrr::map(kob_season, function(s1){
      bloom <- helper_run_pop_model(par = par(), yc_sd = input$yc_sd, zc_sd =  input$zc_sd, jday_cut = jday_cut, temp_df = s1, n = 100,
                                   basic_output = TRUE) %>% 
        purrr::pluck('bloomindex') %>% 
        purrr::map_dbl(helper_bloomint_to_jday, x = s) %>% 
        return()
    })
    do.call(cbind, bloom_list) %>% 
      as.data.frame() %>% 
      setNames(names(kob_season)) %>% 
      pivot_longer(cols = everything(), names_to = 'year')
  })
  

  message("Model finished")
  
  output$pop_plot <- renderPlot({
    req(pop_out())
    
    p_force <- p_flower <- p_tempresp <- NULL
    
    if(input$plot_forcing){
      p_force <- helper_plot_forcing_exp(model_res = pop_out(), obs = exp_obs, jday_cut = jday_cut, jday_name = jday_name)
    } 
    if(input$plot_flowering){
      p_flower <- helper_plot_flowering(bloom_df = pop_bloom(), obs_df = kob_bloom, year_select = input$years_bloom)
    }
    if(input$plot_tempresponse){
      #p_tempresp <- LarsChill::get_temp_response_plot(par = par(), temp_values = seq(from = -10, to = 40, by = 0.1))
      
      p_tempresp <- helper_combined_response_plot(par = par(), temp_values = seq(from = -10, to = 40, by = 0.1))
    }
    
    # p_force <- helper_plot_forcing_exp(model_res = pop_out, obs = test, jday_cut = jday_cut)
    # p_flower <- NULL
    # p_tempresp <- NULL
    # n_plot <- 1
    # 
    # plot_present <- c(TRUE, FALSE, FALSE)
    # plot_list <- list(p_force, p_flower, p_tempresp)
    
    plot_list <- list(p_force, p_flower, p_tempresp)
    plot_present <- c(input$plot_forcing, input$plot_flowering, input$plot_tempresponse)
    n_plot <- input$plot_forcing + input$plot_flowering + input$plot_tempresponse
    
    # width_p1 <- 1
    # if(input$plot_forcing) width_p1 <- 2.5
    # plot_list[[1]] + (plot_list[[2]] / plot_list[[3]]) +
    #   plot_layout(widths = c(width_p1, 1))
    
    
    #design the combined plots. keep it flexible, so that all combinations work
    if(n_plot == 1){
      plot_list[[which(plot_present)]]
    } else if(n_plot == 2){
      p1 <- plot_list[[which(plot_present)[1]]]
      p2 <- plot_list[[which(plot_present)[2]]]
      
      #in case forcing and 

      # plot_height <- 1
      # if(input$plot_forcing) plot_height <- 2
      p1 + p2 


    } else if(n_plot == 3){
      # plot_list[[1]] + (plot_list[[2]] / plot_list[[3]]) +
      #   plot_layout(widths = c(2.5, 1))

      design <- "AAABBB
                 AAABBB
                 CCCCCC"
      
      plot_list[[1]] + plot_list[[2]] +  plot_list[[3]] +
        plot_layout(design = design)
    }
    
    
    
    # else if(input$plot_type == 'flower'){
    #   helper_plot_flowering(model_res = pop_out(), temp_df = s)
    # } else if(input$plot_type == 'forcing_and_flower'){
    #   p1 <- helper_plot_forcing_exp(model_res = pop_out(), obs = test, jday_cut = jday_cut)
    #   p2 <- helper_plot_flowering(model_res = pop_out(), temp_df = s)
    #   library(patchwork)
    #   p1+p2+plot_layout(widths = c(3,1))
    #   
    #}
  })

  
}

shinyApp(ui = ui, server = server)

