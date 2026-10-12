# Regenerate from the unchanged daily weather used by the regression tests.
# Run this script from the package root.
cka_weather <- utils::read.csv("tests/cka.csv")
stations <- utils::read.csv("tests/stations.csv")
station <- stations[stations$station == "Klein-Altendorf", ]
stopifnot(nrow(station) == 1L)
dates <- as.Date(sprintf("%04d-%02d-%02d", cka_weather$Year,
                         cka_weather$Month, cka_weather$Day))
stopifnot(!anyNA(dates), !anyDuplicated(dates), all(diff(dates) == 1),
          all(is.finite(cka_weather$Tmin)), all(is.finite(cka_weather$Tmax)),
          all(cka_weather$Tmin <= cka_weather$Tmax))
attr(cka_weather, "latitude") <- station$latitude
attr(cka_weather, "longitude") <- station$longitude
attr(cka_weather, "source") <- "tests/cka.csv; coordinates from tests/stations.csv"
save(cka_weather, file = "data/cka_weather.rda", compress = "xz", version = 2)
