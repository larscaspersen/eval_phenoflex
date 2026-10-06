#include "phenology_modules.h"
#include "model_validation.h"
#include "structure_validation.h"

//' Apply parallel coupling to precomputed chill and heat
//'
//' Uses the last chill column, regardless of its name. Heat effectiveness is
//' kmin + (1-kmin)*y/yc below yc and one at or above yc, using updated chill.
//' All chill columns continue updating after yc; effectiveness is capped at one.
//' @param chill Finite matrix with at least two rows and one column. The last
//' column is non-negative accumulated chill; earlier columns are internal pools
//' or diagnostics. Rows represent interval endpoints.
//' @param heat Finite non-negative potential heat increments on the same hourly
//' grid, with length nrow(chill)-1.
//' @param yc,zc Positive finite chill and heat requirements in the submodel units.
//' @param kmin Finite minimum heat effectiveness between zero and one.
//' @param stopatzc Stop at first bloom if TRUE, leaving unprocessed output rows
//' zero. FALSE completes the trajectory while retaining the first bloom index.
//' @param basic_output Return only bloomindex if TRUE.
//' @return A list containing bloomindex (one-based row; zero if not reached).
//' Detailed output also includes the unrestricted chill matrix and effective heat z.
//' Input matrices and vectors are not modified.
//' @references Hänninen, H. and Kramer, K. (2007). A framework for modelling the annual cycle of trees in boreal and temperate regions. Silva Fennica 41(1), 167-205. Equation B4b. doi:10.14214/sf.313.
//' @export
// [[Rcpp::export]]
Rcpp::List apply_parallel_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc = 40, double zc = 1119, double kmin = 0.1,
    bool stopatzc = true, bool basic_output = true) {
  evalpheno::validate_structure_inputs(chill, heat);
  evalpheno::require_positive(yc, "yc");
  evalpheno::require_positive(zc, "zc");
  if (!R_finite(kmin) || kmin < 0 || kmin > 1)
    Rcpp::stop("kmin must be finite and between 0 and 1.");
  const int N = chill.nrow(), columns = chill.ncol();
  for (int i = 0; i < N; ++i)
    if (chill(i, columns - 1) < 0)
      Rcpp::stop("Parallel coupling requires non-negative accumulated chill in the last column.");
  Rcpp::NumericMatrix effective_chill(basic_output ? 0 : N, basic_output ? 0 : columns);
  Rcpp::NumericVector z(basic_output ? 0 : N);
  if (!basic_output) {
    effective_chill.attr("dimnames") = chill.attr("dimnames");
    for (int j = 0; j < columns; ++j) effective_chill(0, j) = chill(0, j);
  }
  int bloomindex = 0;
  double accumulated_heat = 0;
  for (int i = 0; i < N - 1; ++i) {
    const int row = i + 1;
    const double y = chill(row, columns - 1);
    const double effectiveness = y >= yc ? 1.0 : kmin + ((1.0 - kmin) / yc) * y;
    accumulated_heat += effectiveness * heat[i];
    if (!R_finite(accumulated_heat)) Rcpp::stop("Effective heat accumulation overflowed.");
    if (!basic_output) {
      for (int j = 0; j < columns; ++j) effective_chill(i + 1, j) = chill(row, j);
      z[i + 1] = accumulated_heat;
    }
    if (bloomindex == 0 && accumulated_heat >= zc) {
      bloomindex = i + 2;
      if (stopatzc) break;
    }
  }
  if (basic_output) return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex);
  return Rcpp::List::create(Rcpp::Named("bloomindex") = bloomindex,
                             Rcpp::Named("chill") = effective_chill, Rcpp::Named("z") = z);
}
