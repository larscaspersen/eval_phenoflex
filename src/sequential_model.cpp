#include <Rcpp.h>
#include "phenology_modules.h"
#include "model_validation.h"
using namespace Rcpp;

//' @title Run a sequential phenology model
//' @description Combine Dynamic Model chill accumulation and GDH heat
//' accumulation in a sequential model.
//' Heat uses the updated chill state: the interval reaching yc contributes its
//' full heat increment. States are recorded at interval endpoints.
//' Dynamic chill and potential GDH are precomputed in separate modules.
//' All chill pools continue updating after yc is reached.
//'
//' @param temp At least two finite hourly temperatures, greater than -273 Celsius
//' (or zero Kelvin), matching the legacy temperature conversion.
//' @param times Finite time points in hours, matching temp in length. Consecutive
//' values must be one hour apart (absolute tolerance 1e-8 hours).
//' @param yc numeric. Chill threshold for full heat effectiveness
//' @param zc numeric. Critical value of z determining the end of heat accumulation
//' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
//' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
//' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
//' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
//' @param slope numeric. Slope parameter for sigmoidal function
//' @param Tf numeric. Transition temperature in the same units as temp.
//' @param Tb numeric. GDH base temperature (lower threshold)
//' @param Tu numeric. GDH optimal temperature
//' @param Tc numeric. GDH upper temperature (upper threshold)
//' @param Delta Positive compatibility argument; unused by the GDH heat model.
//' @param stopatzc logical. If TRUE, stop when the bloom requirement is reached;
//' otherwise process the complete series while preserving the first bloom index.
//' @param deg_celsius If TRUE, temperatures are supplied in Celsius; otherwise Kelvin.
//' @param basic_output logical. If `TRUE`, return a list containing only bloomindex.
//' @useDynLib evalpheno
//' @author Lars Caspersen <lcaspers@uni-bonn.de>
//' @return
//' A list with `bloomindex`, the row index at which bloom occurs.
//' With `basic_output = FALSE`, also return numeric vectors `x`, `y`, `z` and `xs`,
//' each with the same length as the input temperature vector `temp`.
//' @examples
//' # Synthetic winter (90 days at 6 C), followed by spring (45 days at 18 C).
//' temperatures <- c(rep(6, 90 * 24), rep(18, 45 * 24))
//' times <- seq_along(temperatures) - 1 # Elapsed hours
//' # Illustrative requirements: 40 chill portions and 6000 scaled GDH.
//' # Heat contributes to development only after the chill requirement is met.
//' result <- seq_model(
//'   temp = temperatures, times = times, yc = 40, zc = 6000,
//'   stopatzc = FALSE, basic_output = FALSE
//' )
//' bloom_day <- times[result$bloomindex] / 24
//' bloom_day # Elapsed days until bloom, rather than a calendar day of the year
//' data.frame(chill = result$y[result$bloomindex],
//'            heat = result$z[result$bloomindex])
//' plot(times / 24, result$z, type = "l",
//'      xlab = "Days since start", ylab = "Accumulated effective heat (GDH)")
//' abline(h = 6000, lty = 2)
//' abline(v = bloom_day, lty = 3)
//' @export
// [[Rcpp::export]]
 List seq_model(NumericVector temp,
               NumericVector times,
               const double yc=40,
               const double zc=6000,
               const double A0=6319.5,
               const double A1=5.939917e13,
               const double E0=3372.8,
               const double E1=9900.3,
               const double slope=1.6,
               const double Tf=4,
               const double Tu=25,
               const double Tb=4,
               const double Tc=36,
               const double Delta=4,
               bool stopatzc = true,
               bool deg_celsius = true,
               bool basic_output = true) {
   evalpheno::validate_hourly_weather(temp, times, deg_celsius);
   evalpheno::validate_dynamic_gdh(A0, A1, E0, E1, slope, Tf, Tb, Tu, Tc, Delta, deg_celsius);
   evalpheno::require_positive(yc, "yc");
   evalpheno::require_positive(zc, "zc");

   // Precompute unrestricted submodel outputs, then apply structure rules.
   NumericMatrix chill = calculate_chill_dynamic(temp, times, A0, A1, E0, E1,
                                                  slope, Tf, deg_celsius);
   NumericVector heat = calculate_heat_gdh(temp, times, Tu, Tb, Tc, deg_celsius);
   List result = apply_sequential_structure(chill, heat, yc, zc, stopatzc, basic_output);
   if (basic_output) return result;

   NumericMatrix effective_chill = result["chill"];
   NumericVector x = effective_chill(_, 0);
   NumericVector y = effective_chill(_, effective_chill.ncol() - 1);
   NumericVector z = result["z"];
   const int bloomindex = result["bloomindex"];
   const int end = stopatzc && bloomindex > 0 ? bloomindex : temp.size();
   NumericVector xs(temp.size());
   // Preserve the old xs diagnostic convention: interval-start equilibrium,
   // the final/unprocessed values remain zero.
   for (int i = 0; i < end - 1; ++i) xs[i] = chill(i + 1, 1);
   return List::create(Named("x") = x, Named("y") = y, Named("z") = z,
                       Named("xs") = xs, Named("bloomindex") = bloomindex);
 }
