#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

double smooth_gaussian(double distance, double sigma);

double smooth_laplacian(double distance, double sigma);

double smooth_exponential(double distance, double sigma);

rowvec smooth_field_at_point(const mat &displacement,const rowvec &distance, double sigma, double gamma, int smoothtype = 0, const std::vector<float> &use = std::vector<float>());
rowvec smooth_field_at_point_bspline(const rowvec &pt, const mat &clostPoints, const mat &displacement,double support, double gamma, const std::vector<float>  &use = std::vector<float>());
