#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
double smooth_bspline(rowvec pt1, rowvec pt2, double support);

double bspline3(const double& x);

double tensorProductSpline(const rowvec& x);

double get_bspline(const rowvec& x, const rowvec& y, const double support) ;
