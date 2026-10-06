#include <Rcpp.h>
#include "phenology_modules.h"
#include "model_validation.h"
using namespace Rcpp;

//' @title parallel_model
//' @description Parallel model, combining dynamic model for chill accumulation and the GDH model
//' Uses precomputed chill and heat modules. All chill pools continue updating.
//'
//' @param temp At least two finite hourly temperatures, greater than -273 Celsius
//' (or zero Kelvin), matching the legacy temperature conversion.
//' @param times Finite time points in hours, matching temp in length. Consecutive
//' values must be one hour apart (absolute tolerance 1e-8 hours).
//' @param yc numeric. Chill threshold for full heat effectiveness
//' @param zc numeric. Critical value of z determining the end of heat accumulation
//' @param kmin numeric. Minimum relative heat effectiveness at zero accumulated chill.
//' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
//' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
//' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
//' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
//' @param slope numeric. Slope parameter for sigmoidal function
//' @param Tf numeric. Transition temperature in the same units as temp.
//' @param Tb numeric. GDH base temperature (lower threshold) 
//' @param Tu numeric. GDH optimal temperature 
//' @param Tc numeric. GDH upper temperature (upper threshold)
//' @param Delta numeric. Positive compatibility argument; unused by the GDH heat model.
//' @param stopatzc boolean. If TRUE, stop when the bloom requirement is reached;
//' otherwise process the complete series while preserving the first bloom index.
//' @param deg_celsius boolean. If set `TRUE` function assumes degree celsius temperature parameters, otherwise kelvin.
//' @param basic_output boolean. If `TRUE`, only the bloomindex is returned as a named element of the return list.
//' @useDynLib evalpheno
//' @author Lars Caspersen <lcaspers@uni-bonn.de>
//' @return
//' A list is returned with named element `bloomindex`, which is the index at which blooming occurs. When `basic_output=FALSE` also `x`, `y`, `z` and `xs` are
//' returned as named element of this list, which are numeric vectors of the same length as the input vector `temp` containing the hourly temperatures.
//' @examples
//' temperatures <- rep(15, 240)
//' result <- parallel_model(temp = temperatures, times = seq_along(temperatures))
//' result$bloomindex
//' @references Hänninen, H. and Kramer, K. (2007). A framework for modelling the annual cycle of trees in boreal and temperate regions. Silva Fennica 41(1), 167-205. Equation B4b. doi:10.14214/sf.313.
//' @export
// [[Rcpp::export]]
 List parallel_model(NumericVector temp,
               NumericVector times,
               const double yc=40,
               const double zc=1119,
               const double kmin=0.1,
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
   if (!R_finite(kmin) || kmin < 0.0 || kmin > 1.0)
     Rcpp::stop("kmin must be finite and between 0 and 1.");

   // Precompute unrestricted submodel outputs, then apply structure rules.
   NumericMatrix chill = calculate_chill_dynamic(temp, times, A0, A1, E0, E1,
                                                  slope, Tf, deg_celsius);
   NumericVector heat = calculate_heat_gdh(temp, times, Tu, Tb, Tc, deg_celsius);
   List result = apply_parallel_structure(chill, heat, yc, zc, kmin, stopatzc, basic_output);
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
