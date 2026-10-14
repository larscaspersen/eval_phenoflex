#include "model_validation.h"
#include <algorithm>

namespace {
void validate_weights(const Rcpp::NumericVector& weights) {
  if (weights.size() != 7 || weights.hasAttribute("dim"))
    Rcpp::stop("weights must be a numeric vector of seven finite zone weights.");
  for (double value : weights)
    if (!R_finite(value)) Rcpp::stop("weights must contain only finite values.");
}

double zone_response(double temperature, const Rcpp::NumericVector& weights, double offset) {
  const double upper[] = {1.4, 2.4, 9.1, 12.4, 15.9, 18.0, 1000.0};
  double lower = -1000.0 + offset;
  for (int zone = 0; zone < 7; ++zone) {
    if (temperature > lower && temperature <= upper[zone] + offset) return weights[zone];
    lower = upper[zone] + offset;
  }
  return 0.0;
}
}

//' Precompute chill hours
//'
//' Counts each hourly interval whose start temperature is inclusively between
//' CH_lower and CH_upper. The first row is zero; the final temperature is an
//' endpoint and does not contribute a further interval.
//' @inheritParams calculate_chill_dynamic
//' @param CH_lower,CH_upper Finite lower and upper thresholds in the same units
//' as temp, satisfying -273 < CH_lower < CH_upper in Celsius (zero in Kelvin).
//' NA uses defaults of 0 and 7.2 Celsius, or 273 and 280.2 Kelvin.
//' @return A matrix with a cumulative chill-hours column y and one row per time point.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_chill_hours(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double CH_lower = NA_REAL, double CH_upper = NA_REAL, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  if (R_IsNA(CH_lower)) CH_lower = deg_celsius ? 0.0 : 273.0;
  if (R_IsNA(CH_upper)) CH_upper = deg_celsius ? 7.2 : 280.2;
  evalpheno::require_temperature(CH_lower, deg_celsius, "CH_lower");
  evalpheno::require_temperature(CH_upper, deg_celsius, "CH_upper");
  if (CH_lower >= CH_upper) Rcpp::stop("CH_lower must be less than CH_upper.");
  Rcpp::NumericMatrix chill(temp.size(), 1);
  Rcpp::colnames(chill) = Rcpp::CharacterVector::create("y");
  for (int i = 0; i < temp.size() - 1; ++i)
    chill(i + 1, 0) = chill(i, 0) + (temp[i] >= CH_lower && temp[i] <= CH_upper ? 1.0 : 0.0);
  return chill;
}

//' Precompute Utah chill units
//'
//' Uses the supplied seven zone weights at upper boundaries 1.4, 2.4, 9.1,
//' 12.4, 15.9, 18 and 1000 Celsius. Intervals are lower-exclusive and
//' upper-inclusive; temperatures outside the zones contribute zero.
//' The first row is zero and the final temperature is an endpoint.
//' @inheritParams calculate_chill_dynamic
//' @param weights Seven finite zone weights, in increasing temperature order.
//' @return A matrix with signed accumulated units y_raw and non-negative coupling
//' units y = max(0, y_raw). Negative contributions are retained in y_raw and
//' must be recovered by later positive contributions before y increases.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_chill_utah(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    Rcpp::NumericVector weights = Rcpp::NumericVector::create(0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0),
    bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  validate_weights(weights);
  Rcpp::NumericMatrix chill(temp.size(), 2);
  Rcpp::colnames(chill) = Rcpp::CharacterVector::create("y_raw", "y");
  for (int i = 0; i < temp.size() - 1; ++i) {
    const double value = chill(i, 0) + zone_response(temp[i], weights, deg_celsius ? 0.0 : 273.0);
    if (!R_finite(value)) Rcpp::stop("Utah chill accumulation overflowed.");
    chill(i + 1, 0) = value;
    chill(i + 1, 1) = std::max(0.0, value);
  }
  return chill;
}

//' Precompute positive Utah chill units
//'
//' Adds max(0, daily sum) at each day boundary, retaining the signed hourly
//' contributions within that day. Final processed intervals are flushed even
//' when the last day is partial. The final temperature is an endpoint.
//' @inheritParams calculate_chill_utah
//' @param day_group Optional numeric vector of finite, non-decreasing integer
//' day identifiers, one per temperature row. NULL groups consecutive blocks of
//' 24 rows from the first row. Supply calendar groups when starting within a day.
//' @return A matrix with accumulated positive daily chill units y. The first row
//' is zero and y is held constant between day boundaries.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_chill_positive_utah(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    Rcpp::NumericVector weights = Rcpp::NumericVector::create(0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0),
    Rcpp::Nullable<Rcpp::NumericVector> day_group = R_NilValue, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  validate_weights(weights);
  const int N = temp.size();
  Rcpp::NumericVector groups(N);
  if (day_group.isNotNull()) {
    Rcpp::NumericVector supplied(day_group.get());
    if (supplied.size() != N || supplied.hasAttribute("dim"))
      Rcpp::stop("day_group must have one integer day identifier per temperature row.");
    for (int i = 0; i < N; ++i) {
      const double value = supplied[i];
      if (!R_finite(value) || value != std::floor(value) || (i > 0 && value < supplied[i - 1]))
        Rcpp::stop("day_group must contain finite, non-decreasing integer day identifiers.");
      groups[i] = value;
    }
  } else {
    for (int i = 0; i < N; ++i) groups[i] = i / 24;
  }
  Rcpp::NumericMatrix chill(N, 1);
  Rcpp::colnames(chill) = Rcpp::CharacterVector::create("y");
  double daily_sum = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    daily_sum += zone_response(temp[i], weights, deg_celsius ? 0.0 : 273.0);
    if (!R_finite(daily_sum)) Rcpp::stop("Positive Utah daily accumulation overflowed.");
    double value = chill(i, 0);
    if (i == N - 2 || groups[i + 1] != groups[i]) {
      value += std::max(0.0, daily_sum);
      daily_sum = 0.0;
    }
    if (!R_finite(value)) Rcpp::stop("Positive Utah chill accumulation overflowed.");
    chill(i + 1, 0) = value;
  }
  return chill;
}

//' Precompute the activation state for a heat-only model
//'
//' This is an activation indicator, not accumulated chill.
//' @inheritParams calculate_chill_dynamic
//' @param heat_start_index Finite one-based interval row at which heat starts;
//' must be at least one. Fractional values enable the next interval row.
//' @return A matrix with column y equal to one at rows greater than or equal
//' to heat_start_index, and zero before it.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_chill_none(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double heat_start_index = 1.0, bool deg_celsius = true) {
  evalpheno::validate_hourly_weather(temp, times, deg_celsius);
  if (!R_finite(heat_start_index) || heat_start_index < 1.0)
    Rcpp::stop("heat_start_index must be finite and at least one.");
  Rcpp::NumericMatrix chill(temp.size(), 1);
  Rcpp::colnames(chill) = Rcpp::CharacterVector::create("y");
  for (int i = 0; i < temp.size(); ++i) chill(i, 0) = i + 1 >= heat_start_index ? 1.0 : 0.0;
  return chill;
}
