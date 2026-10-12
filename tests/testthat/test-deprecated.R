test_that("deprecated exports warn with their own name and replacement", {
  replacements <- c(
    phenologyFitter_combined_fitting = "fit_phenology",
    eval_all_daoptim = "phenology_rss",
    eval_fixed_daoptim = "phenology_rss",
    evaluation_function_gensa = "phenology_rss",
    eval_phenoflex_single = "phenology_rss",
    eval_phenoflex_single_twofixed = "phenology_rss",
    eval_phenoflex_onlyreq = "phenology_rss",
    eval_phenoflex_combined = "phenology_rss_combined",
    evaluation_function_meigo = "phenology_rss",
    evaluation_function_meigo_vns = "phenology_rss",
    evaluation_function_meigo_nonlinear = "phenology_rss",
    evaluation_function_meigo_nonlinear_combined = "phenology_rss_combined",
    evaluation_function_meigo_nonlinear_fixed = "phenology_rss_combined",
    wrapper_seq_model = "predict_phenology",
    wrapper_parallel_model = "predict_phenology",
    wrapper_po_model = "predict_phenology",
    custom_PhenoFlex_GDHwrapper = "predict_phenology",
    custom_PhenoFlex_GDHwrapper_v2 = "predict_phenology",
    helper_run_pop_model = "phenoflex_population",
    convert_parameters = "characteristic_to_kinetic",
    convert_parameters_old_to_new = "kinetic_to_characteristic"
  )
  for (name in names(replacements)) {
    expect_true(name %in% getNamespaceExports("evalpheno"))
    # Catch the warning before argument validation or model execution.
    warning <- tryCatch(withCallingHandlers(
      eval(call(name), envir = asNamespace("evalpheno")),
      deprecatedWarning = function(w) {
        stop(structure(list(message = "Captured deprecation", warning = w),
                       class = c("deprecation_probe", "error", "condition")))
      }), deprecation_probe = function(e) e$warning)
    expect_true(inherits(warning, "deprecatedWarning"), info = name)
    expect_identical(warning$old, name, info = name)
    expect_identical(warning$new, unname(replacements[[name]]), info = name)
    expect_identical(warning$package, "evalpheno", info = name)
    expect_match(conditionMessage(warning), paste0(name, "() is deprecated"),
                 fixed = TRUE, info = name)
    expect_match(conditionMessage(warning), replacements[[name]], fixed = TRUE,
                 info = name)
  }
})

test_that("deprecated adapters retain results without nested alias warnings", {
  characteristic <- c(40, 190, .5, 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  warnings <- list()
  result <- withCallingHandlers(
    eval_all_daoptim(characteristic, function(weather, par) 100,
                     bloomJDays = 102, SeasonList = list(NULL),
                     intermed_chill = TRUE),
    deprecatedWarning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(result, 4)
  expect_length(warnings, 1)
  expect_identical(warnings[[1]]$old, "eval_all_daoptim")
})

test_that("modular prediction and conversion remain warning-free", {
  characteristic <- c(40, 190, .5, 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  expect_no_warning(kinetic <- characteristic_to_kinetic(characteristic))
  expect_no_warning(kinetic_to_characteristic(kinetic))
  model <- pheno_model("sequential", chill_dynamic("characteristic"),
                      parameters = c(yc = 0.1, zc = 5))
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24))
  expect_no_warning(observed <- predict_phenology(model, weather))
  expect_no_warning(fit <- fit_phenology(
    model, phenology_rss, lower = numeric(), upper = numeric(),
    observed = observed, seasons = list(weather)
  ))
  expect_equal(fit$value, 0)
})
