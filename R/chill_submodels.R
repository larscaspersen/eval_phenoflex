#' Specify alternative chill submodels
#'
#' These constructors select the chill algorithm; configure its parameters through
#' [pheno_model()]. [parameter_schema()] lists the parameters and defaults.
#'
#' @details
#' ## Chilling Hours and Utah Models
#'
#' The Chilling Hours Model is described by Weinberger (1950), and the Utah
#' Model by Richardson et al. (1974).
#'
#' Chill hours counts temperatures inclusively between `CH_lower` (0 Celsius)
#' and `CH_upper` (7.2 Celsius). Utah uses seven zones with upper boundaries
#' 1.4, 2.4, 9.1, 12.4, 15.9, 18 and 1000 Celsius, lower-exclusive and
#' upper-inclusive, with weights 0, 0.5, 1, 0.5, 0, -0.5 and -1. Utah's
#' parameters are `utah_w1` through `utah_w7`; positive Utah uses
#' `pos_utah_w1` through `pos_utah_w7`.
#'
#' ## Positive Utah and signed chill totals
#'
#' The Positive Utah modification follows Linsley-Noakes et al. (1994).
#' Positive Utah sums signed hourly contributions within each calendar day, adds
#' only positive daily totals, and holds accumulated chill constant within a day.
#' It also flushes the final partial day of processed intervals. Ordinary Utah
#' retains its signed cumulative total in `y_raw`; the last column `y` clips it
#' at zero for coupling with the selected structure.
#'
#' ## Heat-only models
#'
#' With `chill_none()`, the chosen coupling is bypassed and only heat is
#' accumulated, starting at the one-based interval row `heat_start_index`
#' (default 1). The structure has only a `zc` requirement; `yc`, `s1` and other
#' chill-dependent coupling parameters are omitted. Fractional start indices
#' enable the first interval whose row is at least that value. Indices beyond
#' the available intervals produce no heat accumulation.
#'
#' ## Chill units and starting values
#'
#' Chill requirements use each submodel's own units: hours, Utah chill units or
#' Dynamic Model portions. The illustrative starting requirements are `yc = 1000`
#' for chill hours, `yc = 1200` for Utah and positive Utah, and `yc = 40` for the
#' Dynamic Model. These are not conversions between units. Choose starting values and search
#' bounds appropriate to your data when changing chill models.
#'
#' ## Choosing a chill model
#'
#' Based on extensive model comparisons and a review of the relevant literature,
#' I do not recommend using the Chilling Hours or Utah Models, especially in warm
#' climates. The Dynamic Model (Chill Portions), although far from perfect, seems
#' more reliable. See Luedeling and Brown (2011) and Fernandez et al. (2020) for
#' comparisons, and [pheno_model()] for the Dynamic Model and structure references.
#' Model performance should still be evaluated for the species, cultivar and
#' climate of interest.
#' @references
#' Weinberger, J. H. (1950). Chilling requirements of peach varieties.
#' Proceedings of the American Society for Horticultural Science 56, 122--128.
#'
#' Richardson, E. A., Seeley, S. D. and Walker, D. R. (1974). A model for
#' estimating the completion of rest for 'Redhaven' and 'Elberta' peach trees.
#' HortScience 9(4), 331--332. \doi{10.21273/HORTSCI.9.4.331}.
#'
#' Linsley-Noakes, G. C., Allan, P. and Matthee, G. (1994). Modification of rest
#' completion prediction models for improved accuracy in South African stone
#' fruit orchards. Journal of the Southern African Society for Horticultural
#' Sciences 4, 13--15.
#'
#' Luedeling, E. and Brown, P. H. (2011). A global analysis of the comparability
#' of winter chill models for fruit and nut trees. International Journal of
#' Biometeorology 55, 411--421. \doi{10.1007/s00484-010-0352-y}.
#'
#' Fernandez, E., Whitney, C. and Luedeling, E. (2020). The importance of chill
#' model selection -- a multi-site analysis. European Journal of Agronomy 119,
#' 126103. \doi{10.1016/j.eja.2020.126103}.
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
