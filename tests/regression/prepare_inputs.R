# Run from the package root. Uses the installed evalpheno selected by R_LIBS.
source("tests/regression/cases.R")
destination <- "tests/testthat/fixtures/weather-regression"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
path <- file.path(destination, "inputs.rds")
if (file.exists(path) && !"--overwrite" %in% commandArgs(TRUE))
  stop("Frozen inputs exist. Use --overwrite only when deliberately replacing the fixtures.")
coords <- read.csv("tests/stations.csv")
ka <- new.env(); utils::data("KA_weather", package = "chillR", envir = ka)
weather <- list(cka = read.csv("tests/cka.csv"), quillota = read.csv("tests/quillota.csv"),
                chillR_KA = subset(ka$KA_weather, Year >= 2007 & Year <= 2009))
latitude <- c(cka = coords$latitude[match("Klein-Altendorf", coords$station)],
              quillota = coords$latitude[match("Quillota", coords$station)],
              chillR_KA = coords$latitude[match("Klein-Altendorf", coords$station)])
stopifnot(all(is.finite(latitude)))
seasons <- list(); metadata <- list()
for (station in names(weather)) {
  d <- weather[[station]]
  dates <- as.Date(sprintf("%04d-%02d-%02d", d$Year, d$Month, d$Day))
  stopifnot(!anyNA(dates), !anyDuplicated(dates), all(diff(dates) == 1),
            all(is.finite(d$Tmin)), all(is.finite(d$Tmax)), all(d$Tmin <= d$Tmax))
  hourly <- chillR::stack_hourly_temps(d, latitude = latitude[[station]])$hourtemps
  # stack_hourly_temps is deterministic here: there are no missing daily inputs.
  stopifnot(nrow(hourly) == nrow(d) * 24L, all(is.finite(hourly$Temp)))
  for (year in 2008:2009) {
    if (latitude[[station]] > 0) {
      s <- chillR::genSeasonList(hourly, mrange = c(8, 6), years = year)[[1]]
      start <- as.Date(paste0(year-1, "-08-01")); end <- as.Date(paste0(year, "-06-30"))
    } else {
      # genSeasonList only accepts seasons spanning two years. Southern seasons
      # here run April-November within the bloom year and need an explicit subset.
      s <- hourly[hourly$Year == year & hourly$Month %in% 4:11, c("Temp", "JDay", "Year")]
      start <- as.Date(paste0(year, "-04-01")); end <- as.Date(paste0(year, "-11-30"))
    }
    rownames(s) <- NULL
    stopifnot(nrow(s) == (as.integer(end-start)+1L)*24L)
    seasons[[station]][[as.character(year)]] <- s
    metadata[[length(metadata)+1L]] <- data.frame(station, year, latitude = latitude[[station]],
      start = as.character(start), end = as.character(end), hours = nrow(s),
      crosses_december = length(unique(s$Year)) > 1,
      contains_feb29 = any(format(as.Date(paste0(s$Year, "-01-01"))+s$JDay-1, "%m-%d") == "02-29"))
  }
}
p <- regression_parameters()
p$kinetic <- unname(evalpheno::characteristic_to_kinetic(p$characteristic))
stopifnot(is.numeric(p$kinetic), length(p$kinetic) == 12L)
# Fixed synthetic dates, independent of the model's predictions. Not measured phenology.
observations <- lapply(names(seasons), function(station) {
  base <- if (station == "quillota") 260 else 110
  list(A = c(base, base+5), B = c(base+10, base+16),
       stages = data.frame(year = 2008:2009, budbreak = c(base-10, base-5),
                           firstbloom = c(base, base+5), fullbloom = c(base+10, base+15)))
})
names(observations) <- names(seasons)
saveRDS(list(seasons = seasons, parameters = p, observations = observations), path, version = 2)
write.csv(do.call(rbind, metadata), file.path(destination, "seasons.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(p), function(name)
  data.frame(parameterization = name, position = seq_along(p[[name]]), value = p[[name]]))),
  file.path(destination, "parameters.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(observations), function(station)
  data.frame(station, year = 2008:2009, cultivar_A = observations[[station]]$A,
             cultivar_B = observations[[station]]$B, observations[[station]]$stages[-1]))),
  file.path(destination, "synthetic_observations.csv"), row.names = FALSE)
provenance <- c("tests/cka.csv", "tests/quillota.csv", "tests/stations.csv")
write.csv(data.frame(file = provenance, md5 = unname(tools::md5sum(provenance))),
          file.path(destination, "weather_sources.csv"), row.names = FALSE)
cat("Prepared six frozen seasons in", path, "\n")
