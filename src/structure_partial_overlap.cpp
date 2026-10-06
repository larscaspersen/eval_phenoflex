#include "phenology_modules.h"
#include "model_validation.h"
#include "structure_validation.h"
#include <cmath>

//' Apply partial-overlap coupling to precomputed chill and heat
//'
//' The last chill column supplies y, regardless of its name. Heat starts in the
//' interval where updated y reaches yc. Once effective heat at the start of an
//' interval reaches b1*ol, the chill value used for coupling is held constant.
//' All returned chill pools continue along their unrestricted trajectories.
//' Bloom requires effective heat
//' to reach b1 + b2*exp(-b3*Ca), where Ca is additional chill after the first
//' yc crossing. The chill value at that crossing is the baseline, so Ca starts
//' at zero even when the hourly update overshoots yc. If chill already meets yc
//' in the initial row, that row is the baseline. Ca is bounded below by zero.
//' The Dynamic Model pools are not reset; returned chill retains total chill.
//' @param chill Finite matrix with at least two rows and one column. The last
//' column is accumulated chill; earlier columns are pools or diagnostics.
//' Rows represent interval endpoints.
//' @param heat Finite non-negative potential heat increments on the same hourly
//' grid, with length nrow(chill)-1.
//' @param yc Positive finite chill requirement for heat accumulation to begin.
//' @param b1 Positive finite baseline heat requirement, also used in b1*ol.
//' @param b2 Non-negative finite additive heat requirement at Ca=0 (no additional chill).
//' @param b3 Non-negative finite coefficient in exp(-b3*Ca).
//' @param ol Non-negative finite multiplier defining when additional chill stops
//' affecting the heat requirement. With zero overlap, Ca remains zero after
//' initial chilling completes, so the heat requirement is b1+b2.
//' @param stopatzc Stop at first bloom if TRUE, leaving unprocessed output rows
//' zero. FALSE completes the trajectory while retaining the first bloom index.
//' @param basic_output Return only bloomindex if TRUE.
//' @return A list containing bloomindex (one-based row; zero if not reached).
//' Detailed output also contains unrestricted chill and effective heat z,
//' with unprocessed rows zero when stopatzc is TRUE.
//' Input matrices and vectors are not modified.
//' @references Pope, K.S., Da Silva, D., Brown, P.H. and DeJong, T.M. (2014).
//' A biologically based approach to modeling spring phenology in temperate
//' deciduous trees. Agricultural and Forest Meteorology 198-199, 15-23.
//' Equation (3), sections 2.2.1-2.2.4. doi:10.1016/j.agrformet.2014.07.009.
//' @export
// [[Rcpp::export]]
Rcpp::List apply_partial_overlap_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc = 40, double b1 = 1119, double b2 = 8677,
    double b3 = 0.01119, double ol = 0.75,
    bool stopatzc = true, bool basic_output = true) {
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(yc, "yc");
  evalpheno::require_positive(b1, "b1");
  evalpheno::require_nonnegative(b2, "b2");
  evalpheno::require_nonnegative(b3, "b3");
  evalpheno::require_nonnegative(ol, "ol");
  const double z_chill_stop = b1 * ol;
  if (!R_finite(z_chill_stop)) Rcpp::stop("b1*ol must be finite.");
  const int N = chill.nrow(), columns = chill.ncol();
  Rcpp::NumericMatrix effective_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    effective_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) effective_chill(0, j) = chill(0, j);
  }
  int requirement_row = -1, bloomindex = 0;
  bool chill_completed = chill(0, columns - 1) >= yc;
  double chill_at_completion = chill_completed ? chill(0, columns - 1) : 0.0;
  double accumulated_heat = 0;
  for (int i = 0; i < N - 1; ++i) {
    // The interval crossing b1*ol still contributes additional chill. Afterward,
    // hold only the coupling value; the reported pools continue accumulating.
    // Even with zero overlap, initial chilling must complete first.
    if (chill_completed && requirement_row < 0 && accumulated_heat >= z_chill_stop)
      requirement_row = i;
    const int row = requirement_row >= 0 ? requirement_row : i + 1;
    const double y = chill(row, columns - 1);
    double z_crit = 0;
    if (y >= yc) {
      if (!chill_completed) {
        chill_at_completion = y;
        chill_completed = true;
      }
      // Count increments after completion, capped at the end of overlap.
      const double additional_chill = y > chill_at_completion ? y - chill_at_completion : 0.0;
      accumulated_heat += heat[i];
      z_crit = b1 + b2 * std::exp(-b3 * additional_chill);
      if (!R_finite(z_crit)) Rcpp::stop("Partial-overlap heat requirement overflowed.");
    }
    if (!R_finite(accumulated_heat)) Rcpp::stop("Effective heat accumulation overflowed.");
    if (!basic_output) {
      for (int j = 0; j < columns; ++j) effective_chill(i + 1, j) = chill(i + 1, j);
      z[i + 1] = accumulated_heat;
    }
    if (bloomindex == 0 && z_crit > 0 && accumulated_heat >= z_crit) {
      bloomindex = i + 2;
      if (stopatzc) break;
    }
  }
  if (basic_output) return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex);
  return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex,
                             Rcpp::Named("chill") = effective_chill, Rcpp::Named("z") = z);
}
