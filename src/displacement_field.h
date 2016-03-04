#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

double smooth_gaussian(double distance, double sigma);

double smooth_laplacian(double distance, double sigma);

double smooth_exponential(double distance, double sigma);

rowvec smooth_field_at_point(mat displacement, rowvec distance, double sigma, double gamma, int smoothtype = 0, std::vector<float> use = std::vector<float>());
rowvec smooth_field_at_point_bspline(rowvec pt, mat clostPoints, mat displacement,double support, double gamma, std::vector<float> use = std::vector<float>());
