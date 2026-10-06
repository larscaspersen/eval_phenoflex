#include "phenology_modules.h"
#include "model_validation.h"
#include <cmath>

namespace {
double conversion_fraction(double T, double Tf, double slope) {
  const double value = slope * Tf * (T - Tf) / T;
  if (value >= 17) return 1.0;
  if (value <= -20) return 0.0;
  const double rate = std::exp(value);
  return rate / (1.0 + rate);
}
}

//' Precompute unrestricted Dynamic Model chill
//'
//' Computes chill throughout the supplied series, independently of a structure's
//' requirements. A structure can freeze this trajectory once its requirement is met.
//' The final column is the chill measure consumed by structure functions.
//' @param temp At least two finite hourly temperatures, above -273 Celsius or
//' zero Kelvin. The existing Celsius conversion adds 273.
//' @param times Finite time points matching temp, one hour apart (tolerance 1e-8 hours).
//' @param A0,A1 Positive finite Dynamic Model rate coefficients.
//' @param E0,E1 Finite energies satisfying 0 < E0 < E1.
//' @param slope Positive finite slope of the conversion function.
//' @param Tf Conversion transition temperature in the same units as temp.
//' @param deg_celsius TRUE for Celsius inputs; FALSE for Kelvin inputs.
//' @return A matrix with columns x, xs, y and one row per input time point.
//' The first row is zero. Row i+1 contains x and y after processing interval i;
//' xs is the equilibrium value used for that update. The last temperature is
//' an endpoint and does not generate a further interval.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_chill_dynamic(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double A0 = 6319.5, double A1 = 5.939917e13,
    double E0 = 3372.8, double E1 = 9900.3,
    double slope = 1.6, double Tf = 4, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  evalpheno::validate_dynamic(A0, A1, E0, E1, slope, Tf, deg_celsius);
  const int N = temp.size();
  Rcpp::NumericMatrix chill(N, 3);
  Rcpp::colnames(chill) = Rcpp::CharacterVector::create("x", "xs", "y");
  const double transition = Tf + (deg_celsius ? 273.0 : 0.0);
  for (int i = 0; i < N - 1; ++i) {
    const double T = temp[i] + (deg_celsius ? 273.0 : 0.0);
    const double xs = A0 / A1 * std::exp(-(E0 - E1) / T);
    const double k1 = A1 * std::exp(-E1 / T);
    double x = xs - (xs - chill(i, 0)) * std::exp(-k1 * (times[i + 1] - times[i]));
    double y = chill(i, 2);
    if (x >= 1.0) {
      const double delta = conversion_fraction(T, transition, slope) * x;
      y += delta;
      x -= delta;
    }
    if (!R_finite(xs) || !R_finite(x) || !R_finite(y))
      Rcpp::stop("Dynamic Model produced a non-finite chill state.");
    chill(i + 1, 0) = x;
    chill(i + 1, 1) = xs;
    chill(i + 1, 2) = y;
  }
  return chill;
}
