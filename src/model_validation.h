#ifndef EVALPHENO_MODEL_VALIDATION_H
#define EVALPHENO_MODEL_VALIDATION_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace evalpheno {

inline void require_positive(double value, const char* name) {
  if (!R_finite(value) || value <= 0.0)
    Rcpp::stop("%s must be finite and positive.", name);
}

inline void require_nonnegative(double value, const char* name) {
  if (!R_finite(value) || value < 0.0)
    Rcpp::stop("%s must be finite and non-negative.", name);
}

inline void require_temperature(double value, bool deg_celsius, const char* name) {
  // Match the existing kernels' Celsius-to-Kelvin offset of 273.
  const double lower = deg_celsius ? -273.0 : 0.0;
  if (!R_finite(value) || value <= lower)
    Rcpp::stop("%s must be finite and greater than %g in the supplied temperature units.",
               name, lower);
}

inline void validate_hourly_weather(const Rcpp::NumericVector& temp,
                                    const Rcpp::NumericVector& times,
                                    bool deg_celsius) {
  if (temp.size() < 2 || times.size() != temp.size())
    Rcpp::stop("temp and times must have the same length, with at least two values.");
  if (temp.size() > std::numeric_limits<int>::max())
    Rcpp::stop("Input length exceeds the supported integer index range.");
  for (R_xlen_t i = 0; i < temp.size(); ++i) {
    require_temperature(temp[i], deg_celsius, "temp");
    if (!R_finite(times[i]))
      Rcpp::stop("times must contain only finite values.");
    // Absolute tolerance in hours permits floating-point roundoff, not gaps.
    if (i > 0 && std::abs((times[i] - times[i - 1]) - 1.0) > 1e-8)
      Rcpp::stop("times must be consecutive hourly values (one-hour intervals).");
  }
}

inline void validate_dynamic(double A0, double A1, double E0, double E1,
                              double slope, double Tf, bool deg_celsius) {
  require_positive(A0, "A0");
  require_positive(A1, "A1");
  require_positive(E0, "E0");
  require_positive(E1, "E1");
  if (E1 <= E0) Rcpp::stop("Chill energies must satisfy 0 < E0 < E1.");
  require_positive(slope, "slope");
  require_temperature(Tf, deg_celsius, "Tf");
}

inline void validate_gdh(double Tb, double Tu, double Tc, bool deg_celsius) {
  require_temperature(Tb, deg_celsius, "Tb");
  require_temperature(Tu, deg_celsius, "Tu");
  require_temperature(Tc, deg_celsius, "Tc");
  if (!(Tb < Tu && Tu < Tc))
    Rcpp::stop("Heat thresholds must satisfy Tb < Tu < Tc.");
}

inline void validate_dynamic_gdh(double A0, double A1, double E0, double E1,
                                 double slope, double Tf, double Tb, double Tu,
                                 double Tc, double Delta, bool deg_celsius) {
  validate_dynamic(A0, A1, E0, E1, slope, Tf, deg_celsius);
  validate_gdh(Tb, Tu, Tc, deg_celsius);
  require_positive(Delta, "Delta");
}

} // namespace evalpheno

#endif
