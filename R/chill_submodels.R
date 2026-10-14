#' Specify alternative chill submodels
#'
#' These constructors select the chill algorithm; configure its parameters through
#' [pheno_model()]. [parameter_schema()] lists the parameters and defaults.
#'
#' Chill hours counts temperatures inclusively between `CH_lower` (0 Celsius)
#' and `CH_upper` (7.2 Celsius). Utah uses seven zones with upper boundaries
#' 1.4, 2.4, 9.1, 12.4, 15.9, 18 and 1000 Celsius, lower-exclusive and
#' upper-inclusive, with weights 0, 0.5, 1, 0.5, 0, -0.5 and -1. Utah's
#' parameters are `utah_w1` through `utah_w7`; positive Utah uses
#' `pos_utah_w1` through `pos_utah_w7`.
#'
#' Positive Utah sums signed hourly contributions within each calendar day, adds
#' only positive daily totals, and holds accumulated chill constant within a day.
#' It also flushes the final partial day of processed intervals. Ordinary Utah
#' retains its signed cumulative total in `y_raw`; the last column `y` clips it
#' at zero for coupling with the selected structure.
#'
#' With `chill_none()`, the chosen coupling is bypassed and only heat is
#' accumulated, starting at the one-based interval row `heat_start_index`
#' (default 1). The structure has only a `zc` requirement; `yc`, `s1` and other
#' chill-dependent coupling parameters are omitted. Fractional start indices
#' enable the first interval whose row is at least that value. Indices beyond
#' the available intervals produce no heat accumulation.
#'
#' Chill requirements use each submodel's own units: hours, Utah chill units or
#' Dynamic Model portions. The illustrative starting requirements are `yc = 1000`
#' for chill hours, `yc = 1200` for Utah and positive Utah, and `yc = 40` for the
#' Dynamic Model. These are not conversions between units. Choose starting values and search
#' bounds appropriate to your data when changing chill models.
#' @return A chill specification accepted by [pheno_model()].
#' @examples
#' hours <- pheno_model(chill = chill_hours(), parameters = c(yc = 600))
#' utah <- pheno_model(chill = chill_utah())
#' positive <- pheno_model(chill = chill_positive_utah())
#' heat_only <- pheno_model(chill = chill_none(), parameters = c(heat_start_index = 25))
#' heat_only
#' @md
#' @export
chill_hours <- function() list(name = "chill_hours")

#' @rdname chill_hours
#' @export
chill_utah <- function() list(name = "utah")

#' @rdname chill_hours
#' @export
chill_positive_utah <- function() list(name = "positive_utah")

#' @rdname chill_hours
#' @export
chill_none <- function() list(name = "none")
