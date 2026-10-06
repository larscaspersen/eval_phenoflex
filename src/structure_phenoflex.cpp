#include "phenology_modules.h"
#include "model_validation.h"
#include "structure_validation.h"
#include <cmath>

//' Apply PhenoFlex coupling to precomputed chill and heat
//'
//' Uses the last chill column, regardless of its name. Heat effectiveness uses
//' start-of-interval chill y, preserving the existing PhenoFlex timing:
//' logistic(s1*yc*(y-yc)/y). At y=0 effectiveness is zero; at y=yc it is 0.5.
//' As in PhenoFlex, sigmoid arguments <= -20 give zero and >= 17 give one.
//' Chill pools continue evolving throughout all processed intervals.
//' @inheritParams apply_parallel_structure
//' @param yc Positive finite chill value at which heat effectiveness is 0.5.
//' @param s1 Positive finite slope parameter for chill-dependent heat effectiveness.
//' @param zc Positive finite heat requirement in the units of heat. The default
//' 190 assumes unscaled GDH; multiply by Tu-Tb when using scaled GDH instead.
//' @return A list containing bloomindex (one-based row; zero if not reached).
//' Detailed output also includes unrestricted chill and cumulative effective
//' heat z. With stopatzc=TRUE, unprocessed rows remain zero.
//' Input matrices and vectors are not modified.
//' @examples
//' temp <- rep(8, 2000)
//' times <- seq_along(temp)
//' chill <- calculate_chill_dynamic(temp, times)
//' heat <- calculate_heat_gdh_unscaled(temp, times)
//' apply_phenoflex_structure(chill, heat, yc=40, zc=190, s1=0.5)
//' @export
// [[Rcpp::export]]
Rcpp::List apply_phenoflex_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc = 40, double zc = 190, double s1 = 0.5,
    bool stopatzc = true, bool basic_output = true) {
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(yc, "yc");
  evalpheno::require_positive(zc, "zc");
  evalpheno::require_positive(s1, "s1");
  const int N = chill.nrow(), columns = chill.ncol();
  for (int i = 0; i < N; ++i)
    if (chill(i, columns - 1) < 0)
      Rcpp::stop("PhenoFlex coupling requires non-negative accumulated chill in the last column.");
  Rcpp::NumericMatrix output_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    output_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) output_chill(0, j) = chill(0, j);
  }
  int bloomindex = 0;
  double accumulated_heat = 0;
  for (int i = 0; i < N - 1; ++i) {
    const double y = chill(i, columns - 1);
    double effectiveness = 0;
    if (y > 0) {
      // Long double avoids intermediate overflow for finite input parameters.
      const long double a = static_cast<long double>(s1) * yc *
        (static_cast<long double>(y) - yc) / y;
      if (a >= 17) effectiveness = 1;
      else if (a > -20) {
        const double e = std::exp(static_cast<double>(a));
        effectiveness = e / (1.0 + e);
      }
    }
    accumulated_heat += effectiveness * heat[i];
    if (!R_finite(accumulated_heat)) Rcpp::stop("Effective heat accumulation overflowed.");
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
