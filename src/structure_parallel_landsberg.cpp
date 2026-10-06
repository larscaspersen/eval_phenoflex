#include "phenology_modules.h"
#include "model_validation.h"
#include "structure_validation.h"
#include <cmath>

//' Apply parallel coupling with Landsberg exponential effectiveness
//'
//' Uses the normalized growth coefficient from Landsberg (1974), equation (5):
//' effectiveness = 1 - exp(-y/y0). Here y is the last chill column, evaluated
//' at the updated interval endpoint. y0 must use the same units as that column.
//' There is no kmin: effectiveness is zero at zero chill and approaches one.
//'
//' This adapts the coefficient to precomputed chill and potential heat. It is
//' not the complete Landsberg logistic bud-growth model, nor his original
//' chilling index or temperature sum. For Dynamic Model chill, y0 must be
//' calibrated in chill portions; Landsberg's I0 = 6.7 is not transferable.
//'
//' All chill pools continue updating, so effectiveness continues approaching one.
//' yc is retained for calling compatibility and does not affect this structure.
//' @param chill Finite matrix with at least two rows and one column. The last
//' column is non-negative accumulated chill; other columns are internal pools.
//' @param heat Finite non-negative potential heat increments, length nrow(chill)-1.
//' @param y0 Positive finite exponential scale in the units of the last chill
//' column. At y=y0, effectiveness is approximately 0.6321.
//' @param yc Positive finite compatibility argument; unused by the exponential curve.
//' @param zc Positive finite effective heat requirement.
//' @param stopatzc Stop at the first bloom if TRUE, leaving unprocessed rows zero.
//' FALSE completes the trajectory while retaining the first bloom index.
//' @param basic_output Return only bloomindex if TRUE.
//' @return A list with bloomindex (one-based input row; zero if not reached).
//' With basic_output=FALSE, also returns chill (unrestricted pools) and z (effective
//' accumulated heat). Inputs are not modified.
//' @references Landsberg, J. J. (1974). Apple fruit bud development and growth;
//' analysis and an empirical model. Annals of Botany 38(5), 1013-1023.
//' Equation (5). doi:10.1093/oxfordjournals.aob.a084891.
//' @seealso \code{\link{apply_parallel_structure}}
//' @examples
//' chill <- cbind(x = rep(0, 5), xs = rep(0, 5), y = c(0, 0, 10, 20, 40))
//' apply_parallel_landsberg_structure(chill, rep(10, 4), y0 = 10,
//'                                    yc = 30, zc = 20, basic_output = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::List apply_parallel_landsberg_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double y0, double yc = 40, double zc = 1119,
    bool stopatzc = true, bool basic_output = true) {
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(yc, "yc");
  evalpheno::require_positive(zc, "zc");
  evalpheno::require_positive(y0, "y0");
  const int N = chill.nrow(), columns = chill.ncol();
  for (int i = 0; i < N; ++i)
    if (chill(i, columns - 1) < 0)
      Rcpp::stop("Parallel coupling requires non-negative accumulated chill in the last column.");
  Rcpp::NumericMatrix effective_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    effective_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) effective_chill(0, j) = chill(0, j);
  }
  int bloomindex = 0;
  double accumulated_heat = 0;
  for (int i = 0; i < N - 1; ++i) {
    const int row = i + 1;
    const double y = chill(row, columns - 1);
    const double effectiveness = -std::expm1(-y / y0);
    accumulated_heat += effectiveness * heat[i];
    if (!R_finite(accumulated_heat)) Rcpp::stop("Effective heat accumulation overflowed.");
    if (!basic_output) {
      for (int j = 0; j < columns; ++j) effective_chill(i + 1, j) = chill(row, j);
      z[i + 1] = accumulated_heat;
    }
    if (bloomindex == 0 && accumulated_heat >= zc) {
      bloomindex = i + 2;
      if (stopatzc) break;
    }
  }
  if (basic_output) return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex);
  return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex,
                             Rcpp::Named("chill") = effective_chill, Rcpp::Named("z") = z);
}
