#' Convert an hourly row index to a fractional day of the year
#'
#' Convert the row at which a phenological requirement is met to a day of the
#' year with a fractional part representing the time of day. This helper is
#' used when predicting successive phenological stages.
#'
#' @param index Integer row index at which the phenological requirement is met.
#'
#' @param Jday_vec Numeric vector of days of the year, one per weather row.
#'
#' @param year_vec Vector of calendar years, one per weather row.
#'
#' @return Numeric day of the year with a fractional part for the time of day.
#' Noon corresponds to an integer value. Times after noon add up to 0.5;
#' times before noon subtract up to 0.5. In a season spanning two calendar years,
#' dates in the first year are expressed relative to the end of that year.
#'
#' @author Lars Caspersen, \email{lars.caspersen@@uni-bonn.de}
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom assertthat are_equal
#' @importFrom nleqslv nleqslv
#' @examples
#' \dontrun{
#'
#'  i <- 56
#'  Jday_vec <- rep(1:30, each = 24)
#'  year_vec <- rep(2024, length(Jday_vec))
#'  return_JDay(i, Jday_vec, year_vec)
#' }
#' @export return_JDay
return_JDay <- function(index, Jday_vec, year_vec = NULL){
  if (!is.numeric(index) || length(index) != 1L ||
      (!is.na(index) && (!is.finite(index) || index < 0 ||
                        index != floor(index) || index > length(Jday_vec)))) {
    stop("index must be a single valid row index, 0, or NA.", call. = FALSE)
  }
  if (is.na(index) || index == 0) return(NA_real_)
  if (!is.numeric(Jday_vec) || any(!is.finite(Jday_vec))) {
    stop("Jday_vec must contain finite day-of-year values.", call. = FALSE)
  }
  if (!is.null(year_vec) &&
      (!is.numeric(year_vec) || length(year_vec) != length(Jday_vec) ||
       any(!is.finite(year_vec)) || any(year_vec != floor(year_vec)))) {
    stop("year_vec must contain one finite integer year per row.", call. = FALSE)
  }
  day <- Jday_vec[index]
  same_day <- Jday_vec == day
  if (!is.null(year_vec)) same_day <- same_day & year_vec == year_vec[index]
  rows <- which(same_day)
  if (!is.null(year_vec) && length(unique(year_vec)) == 2L &&
      year_vec[index] == min(year_vec)) {
    year <- year_vec[index]
    leap <- year %% 4 == 0 & (year %% 100 != 0 | year %% 400 == 0)
    day <- day - 365 - as.integer(leap)
  }
  if (length(rows) == 1L) return(day)
  day + (match(index, rows) - ceiling(length(rows) / 2)) / length(rows)
}
