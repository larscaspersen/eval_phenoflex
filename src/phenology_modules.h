#ifndef EVALPHENO_PHENOLOGY_MODULES_H
#define EVALPHENO_PHENOLOGY_MODULES_H

#include <Rcpp.h>

Rcpp::NumericMatrix calculate_chill_dynamic(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double A0, double A1, double E0, double E1, double slope, double Tf,
    bool deg_celsius);

Rcpp::NumericVector calculate_heat_gdh(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double Tu, double Tb, double Tc, bool deg_celsius);

Rcpp::List apply_sequential_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc, double zc, bool stopatzc, bool basic_output);

Rcpp::NumericVector calculate_heat_gdh_unscaled(
    Rcpp::NumericVector temp, Rcpp::NumericVector times,
    double Tu, double Tb, double Tc, bool deg_celsius);

Rcpp::List apply_phenoflex_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc, double zc, double s1, bool stopatzc, bool basic_output);

Rcpp::List apply_parallel_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc, double zc, double kmin, bool stopatzc, bool basic_output);

Rcpp::List apply_parallel_landsberg_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double y0, double yc, double zc, bool stopatzc, bool basic_output);

Rcpp::List apply_partial_overlap_structure(
    Rcpp::NumericMatrix chill, Rcpp::NumericVector heat,
    double yc, double b1, double b2, double b3, double ol,
    bool stopatzc, bool basic_output);

#endif
