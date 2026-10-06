test_that("three-stage chill checks use hourly indices across calendar boundaries", {
  parameters <- c(40, 100, 190, 250, .5, 25, 286.1, 47.7, 28, 4, 4, 1.6)
  for (start in as.Date(c("2007-12-30", "2008-02-28"))) {
    dates <- as.Date(start, origin = "1970-01-01") + 0:2
    year <- rep(as.integer(format(dates, "%Y")), each = 24)
    jday <- rep(as.integer(format(dates, "%j")), each = 24)
    # Budburst at hour 6, first bloom at 30, full bloom at 54.
    heat <- c(rep(0, 5), rep(100, 24), rep(190, 24), rep(250, 19))
    chill <- c(rep(0, 5), rep(40, 67))
    trajectory <- cbind(Year = year, JDay = jday,
                        chill_accumulated = chill, heat_accumulated = heat)
    model <- function(x, par) list(JDay = return_JDay(30, jday, year), chill_heat = trajectory)
    # Independent calendar calculation, with hour 12 represented as an integer day.
    expected <- as.numeric(dates - as.Date(paste0(max(year), "-01-01"))) + 1
    observed <- data.frame(budbreak = expected[1] - .25,
                           firstbloom = expected[2] - .25,
                           fullbloom = expected[3] - .25)
    result <- eval_phenoflex_three_stages(parameters, model, observed, list(NULL))
    expect_equal(result$F, 0)
    expect_length(result$g, 8)
    # The threshold is reached, but chill is insufficient at that exact hour.
    trajectory[6, "chill_accumulated"] <- 39
    expect_equal(eval_phenoflex_three_stages(parameters, model, observed, list(NULL))$F,
                 sum((365 - unlist(observed))^2))
    # Never reaching budburst must give a scalar penalty, not index with 365.
    trajectory[, "heat_accumulated"] <- 0
    expect_equal(eval_phenoflex_three_stages(parameters, model, observed, list(NULL))$F,
                 sum((365 - unlist(observed))^2))
  }
})

test_that("three-stage station objectives match independently evaluated stage dates", {
  fixtures <- readRDS(test_path("fixtures", "weather-regression", "inputs.rds"))
  for (station in names(fixtures$seasons)) {
    seasons <- fixtures$seasons[[station]]
    for (requirements in list(c(40, 100, 190, 250), c(20, 50, 100, 150))) {
      p <- fixtures$parameters$kinetic
      p[1] <- requirements[1]
      predictions <- sapply(requirements[2:4], function(zc) {
        p[2] <- zc
        vapply(seasons, custom_PhenoFlex_GDHwrapper, numeric(1), par = p)
      })
      observed <- data.frame(budbreak = predictions[, 1], firstbloom = predictions[, 2],
                             fullbloom = predictions[, 3])
      cp <- fixtures$parameters$characteristic
      x <- c(requirements, .5, 25, cp[6:9], cp[11:12])
      result <- eval_phenoflex_three_stages(x, custom_PhenoFlex_GDHwrapper_v2, observed, seasons)
      expect_equal(result$F, 0, info = station)
    }
  }
})
