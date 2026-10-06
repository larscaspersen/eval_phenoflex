#include "phenology_modules.h"
#include "model_validation.h"
#include "structure_validation.h"

//' Apply sequential coupling to precomputed chill and heat
//'
//' The last chill column controls heat accumulation, regardless of its name.
//' Heat starts when the last column first reaches yc and remains enabled.
//' Chill pools continue updating. The crossing interval contributes its full
//' potential heat increment.
//' @param chill Finite numeric matrix with at least two rows and one column.
//' The last column is the accumulated chill measure. Earlier columns are
//' model-specific diagnostics or internal pools. Rows represent interval endpoints.
//' @param heat Finite non-negative potential heat increments, one per interval;
//' length must equal nrow(chill)-1. Inputs must share the same hourly time grid.
//' @param yc,zc Positive finite chill and heat requirements, in the units of the
//' supplied submodel outputs.
//' @param stopatzc Stop processing at the first bloom event. Unprocessed output
//' rows remain zero. FALSE completes the trajectory and retains the first event.
//' @param basic_output Return only bloomindex when TRUE.
//' @return A list with bloomindex (one-based input row; zero if not reached).
//' Detailed output also contains the structure-controlled chill matrix and
//' cumulative effective heat z. The supplied inputs are not modified.
//' @export
// [[Rcpp::export]]
Rcpp::List apply_sequential_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc = 40, double zc = 6000,
    bool stopatzc = true, bool basic_output = true) {
  const int N = chill.nrow();
  const int columns = chill.ncol();
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(yc, "yc");
  evalpheno::require_positive(zc, "zc");

  Rcpp::NumericMatrix effective_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    effective_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) effective_chill(0, j) = chill(0, j);
  }
  bool heat_started = chill(0, columns - 1) >= yc;
  int bloomindex = 0;
  double accumulated_heat = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    if (chill(i + 1, columns - 1) >= yc) heat_started = true;
    const int row = i + 1;
    if (heat_started) accumulated_heat += heat[i];
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
