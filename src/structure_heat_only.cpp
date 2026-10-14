#include "model_validation.h"
#include "structure_validation.h"

//' Accumulate potential heat without chill coupling
//'
//' Heat is enabled by the start-of-interval activation state in the last chill
//' column. Use calculate_chill_none() to construct this indicator. No chill
//' requirement or transition parameters are used.
//' @inheritParams apply_sequential_structure
//' @param chill Finite matrix with at least two rows. Its last column contains
//' zero/one activation states, not accumulated chill, from calculate_chill_none().
//' @param zc Positive finite heat requirement in the units of the heat increments.
//' @return A list containing the first bloomindex (zero if not reached). Detailed
//' output also contains chill activation states and cumulative heat z. Early
//' stopping leaves unprocessed rows zero. Inputs are not modified.
//' @export
// [[Rcpp::export]]
Rcpp::List apply_heat_only_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double zc = 190.0, bool stopatzc = true, bool basic_output = true) {
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(zc, "zc");
  const int N = chill.nrow(), columns = chill.ncol();
  for (int i = 0; i < N; ++i)
    if (chill(i, columns - 1) != 0.0 && chill(i, columns - 1) != 1.0)
      Rcpp::stop("Heat-only activation states must be zero or one.");
  Rcpp::NumericMatrix output_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    output_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) output_chill(0, j) = chill(0, j);
  }
  int bloomindex = 0;
  double accumulated_heat = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    if (chill(i, columns - 1) == 1.0) accumulated_heat += heat[i];
    if (!R_finite(accumulated_heat)) Rcpp::stop("Heat accumulation overflowed.");
    if (!basic_output) {
      for (int j = 0; j < columns; ++j) output_chill(i + 1, j) = chill(i + 1, j);
      z[i + 1] = accumulated_heat;
    }
    if (bloomindex == 0 && accumulated_heat >= zc) {
      bloomindex = i + 2;
      if (stopatzc) break;
    }
  }
  if (basic_output) return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex);
  return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex,
    Rcpp::Named("chill") = output_chill, Rcpp::Named("z") = z);
}
