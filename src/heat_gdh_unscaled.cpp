#include "phenology_modules.h"
#include "model_validation.h"
#include <cmath>

//' Precompute unscaled GDH increments
//'
//' Uses the GDH temperature response without the Tu-Tb multiplier, as in
//' PhenoFlex. The hourly increment ranges from zero to one, peaking at Tu.
//' Heat requirements must use these normalized units. To convert a requirement
//' from scaled GDH, divide it by Tu-Tb (both b1 and b2 for partial overlap).
//' @inheritParams calculate_heat_gdh
//' @return A numeric vector of length length(temp)-1. Element i is the potential
//' heat increment over interval i to i+1, using temp[i]. It is not cumulative.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculate_heat_gdh_unscaled(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double Tu = 25, double Tb = 4, double Tc = 36, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  evalpheno::validate_gdh(Tb, Tu, Tc, deg_celsius);
  const double pi = 3.14159265358979323846;
  const double offset = deg_celsius ? 273.0 : 0.0;
  const double lower = Tb + offset, optimum = Tu + offset, upper = Tc + offset;
  Rcpp::NumericVector heat(temp.size() - 1);
  for (R_xlen_t i = 0; i < heat.size(); ++i) {
    const double T = temp[i] + offset;
    if (T >= lower && T <= optimum)
      heat[i] = 0.5 * (1.0 + std::cos(pi + pi * (T - lower) / (optimum - lower)));
    else if (T > optimum && T <= upper)
      heat[i] = 1.0 + std::cos(pi / 2.0 + pi / 2.0 * (T - optimum) / (upper - optimum));
    if (!R_finite(heat[i])) Rcpp::stop("Unscaled GDH produced a non-finite heat increment.");
  }
  return heat;
}
