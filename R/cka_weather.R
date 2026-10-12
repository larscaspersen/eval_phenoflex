#' Daily temperatures at Klein-Altendorf (CKA)
#'
#' Daily minimum and maximum temperatures for 2007--2009, taken from the
#' package's existing test weather file. Hourly temperatures are not observed
#' in this dataset; they can be estimated with chillR::stack_hourly_temps().
#' No phenology observations are included.
#'
#' @docType data
#' @usage data("cka_weather", package = "evalpheno")
#' @format A data frame with 1,096 rows and five columns:
#' \describe{
#'   \item{Year}{Calendar year.}
#'   \item{Month}{Calendar month, 1--12.}
#'   \item{Day}{Day of the month.}
#'   \item{Tmax}{Daily maximum temperature in degrees Celsius.}
#'   \item{Tmin}{Daily minimum temperature in degrees Celsius.}
#' }
#' The \code{latitude} and \code{longitude} attributes retain the station
#' coordinates supplied in the test metadata (50.4 and 6.9945 degrees).
#' @source The repository files \code{tests/cka.csv} and
#' \code{tests/stations.csv}. The original external data source is not recorded
#' in these files.
#' @examples
#' data("cka_weather", package = "evalpheno")
#' head(cka_weather)
#' range(cka_weather$Year)
"cka_weather"
