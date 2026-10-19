# Shared case definitions. No expected values are calculated here.
regression_parameters <- function() {
  characteristic <- c(20, 100, .5, 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  list(characteristic = characteristic,
       # Original Dynamic Model coefficients, separate from the converted equivalent.
       original = c(20, 100, .5, 25, 4153.5, 12888.8, 139500, 2.567e18, 4, 36, 4, 1.6))
}

run_regression_cases <- function(fixtures, api = asNamespace("evalpheno")) {
  call <- function(name, ...) get(name, envir = api)(...)
  wrapper <- get("custom_PhenoFlex_GDHwrapper", envir = api)
  detailed <- get("custom_PhenoFlex_GDHwrapper_v2", envir = api)
  rows <- list()
  add <- function(id, station, season, function_name, output, value,
                  cultivar = "", stage = "", parameterization = "kinetic") {
    value <- unlist(value, use.names = FALSE)
    stopifnot(is.numeric(value), length(value) > 0L)
    rows[[length(rows) + 1L]] <<- data.frame(case_id = id, station = station,
      season = as.character(season), function_name = function_name,
      cultivar = cultivar, stage = stage, parameterization = parameterization,
      output = output, element = seq_along(value),
      status = ifelse(is.na(value), "missing", ifelse(is.infinite(value), "Inf", "finite")),
      value = value, stringsAsFactors = FALSE)
  }
  for (station in names(fixtures$seasons)) {
    seasons <- fixtures$seasons[[station]]
    p <- fixtures$parameters$kinetic
    cp <- fixtures$parameters$characteristic
    obs <- fixtures$observations[[station]]
    for (year in names(seasons)) {
      s <- seasons[[year]]
      emit <- function(id, fn, output, value, parameterization = "kinetic", stage = "")
        add(id, station, year, fn, output, value, stage = stage, parameterization = parameterization)
      emit("ordinary", "custom_PhenoFlex_GDHwrapper", "bloom_jday", wrapper(s, p))
      emit("original_dynamic", "custom_PhenoFlex_GDHwrapper", "bloom_jday",
           wrapper(s, fixtures$parameters$original), "original_kinetic")
      emit("characteristic", "custom_PhenoFlex_GDHwrapper", "bloom_jday",
           wrapper(s, call("characteristic_to_kinetic", cp)), "characteristic")
      no_bloom <- p; no_bloom[1:2] <- 1e12
      emit("no_bloom", "custom_PhenoFlex_GDHwrapper", "bloom_jday", wrapper(s, no_bloom))
      emit("no_bloom_penalty", "eval_all_daoptim", "rss",
           call("eval_all_daoptim", no_bloom, wrapper, obs$A[match(year, names(seasons))], list(s)))
      out <- detailed(s, p)
      emit("detailed", "custom_PhenoFlex_GDHwrapper_v2", "bloom_jday", out$JDay)
      # Freeze daily endpoints, not just final totals: protects accumulation trajectories.
      ends <- which(!duplicated(paste(s$Year, s$JDay), fromLast = TRUE))
      emit("detailed", "custom_PhenoFlex_GDHwrapper_v2", "daily_chill", out$chill_heat[ends, 3])
      emit("detailed", "custom_PhenoFlex_GDHwrapper_v2", "daily_heat", out$chill_heat[ends, 4])
      for (zc in c(100, 190, 250)) {
        stage_par <- p; stage_par[1:2] <- c(40, zc)
        emit("stage_predictions", "custom_PhenoFlex_GDHwrapper", "bloom_jday",
             wrapper(s, stage_par), stage = paste0("heat_", zc))
      }
      sub <- p[c(5:9, 12, 11, 4, 10)]
      variants <- list(sequential = c(20, 100, sub),
                       parallel = c(20, 100, .1, sub),
                       partial_overlap = c(20, 150, 100, .1, .5, sub))
      fn <- c(sequential = "wrapper_seq_model", parallel = "wrapper_parallel_model",
              partial_overlap = "wrapper_po_model")
      for (model in names(variants)) {
        emit(model, fn[[model]], "bloom_jday", call(fn[[model]], s, variants[[model]]))
        emit(paste0(model, "_evaluation"), "eval_all_daoptim", "rss",
             call("eval_all_daoptim", variants[[model]], get(fn[[model]], api),
                  obs$A[match(year, names(seasons))], list(s)))
      }
      # Real hourly rows exercise Dec 31, Jan 1, and Feb 29 explicitly, irrespective of bloom timing.
      for (boundary in c("dec31", "jan01", "feb29")) {
        dates <- as.Date(paste0(s$Year, "-01-01")) + s$JDay - 1L
        key <- c(dec31 = "12-31", jan01 = "01-01", feb29 = "02-29")[[boundary]]
        indices <- which(format(dates, "%m-%d") == key)
        if (length(indices)) {
          selected <- indices[unique(c(1L, 12L, length(indices)))]
          emit(paste0("calendar_", boundary), "return_JDay", "jday",
               vapply(selected, function(i) call("return_JDay", i, s$JDay, s$Year), numeric(1)))
        }
      }
    }
    emit_group <- function(id, fn, output, value, cultivar = "", parameterization = "kinetic")
      add(id, station, "2008+2009", fn, output, value, cultivar, parameterization = parameterization)
    emit_group("all_parameters", "eval_all_daoptim", "rss",
      call("eval_all_daoptim", p, wrapper, obs$A, seasons))
    emit_group("characteristic_evaluation", "eval_all_daoptim", "rss",
      call("eval_all_daoptim", cp, wrapper, obs$A, seasons, intermed_chill = TRUE),
      parameterization = "characteristic")
    fixed <- call("eval_phenoflex_onlyreq", p[1:3], wrapper, obs$A, seasons, p[4:12])
    emit_group("fixed_parameters", "eval_phenoflex_onlyreq", "F", fixed$F)
    emit_group("fixed_parameters", "eval_phenoflex_onlyreq", "g", fixed$g)
    emit_group("fixed_sequential", "eval_fixed_daoptim", "rss",
      call("eval_fixed_daoptim", c(20, 100), get("wrapper_seq_model", api), obs$A, seasons))
    combined_x <- c(20, 25, 100, 150, .5, .6, 25, cp[6:9], cp[11:12])
    combined_args <- list(x = combined_x, modelfn = wrapper, bloomJDays = list(obs$A, obs$B),
                         SeasonList = list(seasons, seasons), ncult = 2)
    combined <- do.call(get("eval_phenoflex_combined", api), combined_args)
    predictions <- do.call(get("eval_phenoflex_combined", api), c(combined_args, list(return_pred = TRUE)))
    emit_group("combined", "eval_phenoflex_combined", "F", combined$F)
    emit_group("combined", "eval_phenoflex_combined", "g", combined$g)
    for (i in 1:2) emit_group("combined", "eval_phenoflex_combined", "bloom_jday",
                             predictions[(2*i-1):(2*i)], c("A", "B")[i])
    stages <- call("eval_phenoflex_three_stages", c(40, 100, 190, 250, .5, 25, cp[6:9], cp[11:12]),
                   detailed, obs$stages, seasons)
    emit_group("three_stages", "eval_phenoflex_three_stages", "F", stages$F)
    emit_group("three_stages", "eval_phenoflex_three_stages", "g", stages$g)
    if (station != "quillota") {
      # Preserve a discovered failure separately, rather than hiding it or fixing production code.
      failure <- tryCatch(call("eval_phenoflex_three_stages",
        c(20, 50, 100, 150, .5, 25, cp[6:9], cp[11:12]), detailed, obs$stages, seasons),
        error = function(e) e)
      emit_group("known_three_stage_index_error", "eval_phenoflex_three_stages",
        "raises_purrr_indexed_error", as.numeric(inherits(failure, "purrr_error_indexed")))
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

compare_regression <- function(expected, actual, atol = 1e-8, rtol = 1e-10) {
  keys <- setdiff(names(expected), "value")
  if (!identical(expected[keys], actual[keys])) stop("Regression case keys, shapes or statuses changed.")
  finite <- is.finite(expected$value) & is.finite(actual$value)
  equal <- rep(FALSE, nrow(expected))
  equal[finite] <- abs(actual$value[finite] - expected$value[finite]) <=
    atol + rtol * abs(expected$value[finite])
  equal[!finite] <- (is.na(expected$value[!finite]) & is.na(actual$value[!finite])) |
    (!is.na(expected$value[!finite]) & !is.na(actual$value[!finite]) &
       expected$value[!finite] == actual$value[!finite])
  result <- actual
  result$expected <- expected$value
  result$difference <- actual$value - expected$value
  result$passed <- equal
  result
}

