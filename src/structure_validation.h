#ifndef EVALPHENO_STRUCTURE_VALIDATION_H
#define EVALPHENO_STRUCTURE_VALIDATION_H
#include <Rcpp.h>

namespace evalpheno {
inline void validate_structure_inputs(const Rcpp::NumericMatrix& chill,
                                      const Rcpp::NumericVector& heat) {
  if (chill.nrow() < 2 || chill.ncol() < 1 || heat.size() != chill.nrow() - 1)
    Rcpp::stop("chill needs at least two rows and one column; heat must have nrow(chill)-1 values.");
  for (double value : chill)
    if (!R_finite(value)) Rcpp::stop("chill must contain only finite values.");
  for (double value : heat)
    if (!R_finite(value) || value < 0.0)
      Rcpp::stop("heat increments must be finite and non-negative.");
}
}
#endif
