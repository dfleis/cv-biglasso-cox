#include <Rcpp.h>

using namespace Rcpp;

double norm( double x, double y ) {
  return sqrt( x*x + y*y );
}

RcppExport SEXP norm_wrapper(SEXP x_, SEXP y_) {
  // step 0: convert input to C++ types
  double x = as<double>(x_), y = as<double>(y_);
  // step 1: call the underlying C++ function
  double res = norm(x, y);
  // step 2: return the result as a SEXP
  return wrap(res);
}