#include <Rcpp.h>
#include "phenology_modules.h"
#include "model_validation.h"
using namespace Rcpp;

//' @title partial_overlap_model
//' @description Partial Overlap Model, combining the dynamic model for chill accumulation and the GDH model
//' Uses precomputed chill and heat modules. Chill continues throughout. Once heat at
//' the start of an interval reaches b1*ol, additional chill no longer changes the heat requirement.
//'
//' Additional chill is measured from the first row reaching yc, including any
//' overshoot in the baseline. The Dynamic Model pools retain their total values.
//'
//' @param temp At least two finite hourly temperatures, greater than -273 Celsius
//' (or zero Kelvin), matching the legacy temperature conversion.
//' @param times Finite time points in hours, matching temp in length. Consecutive
//' values must be one hour apart (absolute tolerance 1e-8 hours).
//' @param yc numeric, Chill requirement for heat accumulation to begin
//' @param b1 Positive finite baseline in the heat requirement b1 + b2*exp(-b3*Ca), where Ca is additional chill after yc is first reached.
//' @param b2 Non-negative finite additive heat requirement at zero additional chill.
//' @param b3 Non-negative finite decay coefficient for additional chill Ca.
//' @param ol Non-negative overlap multiplier. Additional chill stops affecting the heat requirement when heat reaches b1*ol.
//' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
//' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
//' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
//' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
//' @param slope numeric. Slope parameter for sigmoidal function
//' @param Tf numeric. Transition temperature in the same units as temp.
//' @param Tb numeric. GDH base temperature (lower threshold)  
//' @param Tu numeric. GDH optimal temperature 
//' @param Tc numeric. GDH critical temperature 
//' @param Delta numeric. Positive compatibility argument; unused by the GDH heat model.
//' @param stopatzc boolean. If TRUE, stop when the bloom requirement is reached;
//' otherwise process the complete series while preserving the first bloom index.
//' @param deg_celsius If TRUE, temperatures are supplied in Celsius; otherwise Kelvin.
//' @param basic_output boolean. If `TRUE`, only the bloomindex is returned as a named element of the return list.
//' @useDynLib evalpheno
//' @author Lars Caspersen <lcaspers@uni-bonn.de>
//' @return
//' A list is returned with named element `bloomindex`, which is the index at which blooming occurs. When `basic_output=FALSE` also `x`, `y`, `z` and `xs` are
//' returned as named element of this list, which are numeric vectors of the same length as the input vector `temp` containing the hourly temperatures.
//' @examples
//' temperatures <- rep(15, 240)
//' result <- po_model(temp = temperatures, times = seq_along(temperatures))
//' result$bloomindex
//' @export
// [[Rcpp::export]]
 List po_model(NumericVector temp,
                NumericVector times,
                const double yc=40,
                const double b1=1119,
                const double b2=8677,
                const double b3=0.01119,
                const double ol=0.75,
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
   evalpheno::require_positive(b1, "b1");
   evalpheno::require_nonnegative(b2, "b2");
   evalpheno::require_nonnegative(b3, "b3");
   evalpheno::require_nonnegative(ol, "ol");

   // Precompute unrestricted submodel outputs, then apply structure rules.
   NumericMatrix chill = calculate_chill_dynamic(temp, times, A0, A1, E0, E1,
                                                  slope, Tf, deg_celsius);
   NumericVector heat = calculate_heat_gdh(temp, times, Tu, Tb, Tc, deg_celsius);
   List result = apply_partial_overlap_structure(chill, heat, yc, b1, b2, b3, ol, stopatzc, basic_output);
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
