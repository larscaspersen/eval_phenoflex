#include "phenology_modules.h"
#include "model_validation.h"
#include <cmath>

namespace {
double gdh(double T, double Tu, double Tb, double Tc) {
  const double pi = 3.14159265358979323846;
  if (T >= Tb && T <= Tu)
    return (Tu - Tb) / 2.0 * (1.0 + std::cos(pi + pi * (T - Tb) / (Tu - Tb)));
  if (T > Tu && T <= Tc)
    return (Tu - Tb) * (1.0 + std::cos(pi / 2.0 + pi / 2.0 * (T - Tu) / (Tc - Tu)));
  return 0.0;
}
}

//' Precompute potential GDH increments
//'
//' Calculates potential heat independently of chill and model structure. Uses
//' the GDH scaling of the sequential, parallel and partial-overlap kernels.
//' @param temp At least two finite hourly temperatures, above -273 Celsius or zero Kelvin.
//' @param times Finite time points matching temp, one hour apart (tolerance 1e-8 hours).
//' @param Tu,Tb,Tc Finite temperature thresholds satisfying Tb < Tu < Tc,
//' in the same units as temp.
//' @param deg_celsius TRUE for Celsius inputs; FALSE for Kelvin inputs.
//' @return A numeric vector of length length(temp)-1. Element i is the potential
//' heat increment over interval i to i+1, using temp[i]. It is not cumulative.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculate_heat_gdh(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double Tu = 25, double Tb = 4, double Tc = 36, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  evalpheno::validate_gdh(Tb, Tu, Tc, deg_celsius);
  const double offset = deg_celsius ? 273.0 : 0.0;
  Rcpp::NumericVector heat(temp.size() - 1);
  for (R_xlen_t i = 0; i < heat.size(); ++i) {
    heat[i] = gdh(temp[i] + offset, Tu + offset, Tb + offset, Tc + offset);
    if (!R_finite(heat[i])) Rcpp::stop("GDH produced a non-finite heat increment.");
  }
  return heat;
}
