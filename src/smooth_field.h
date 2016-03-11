#include "RcppArmadillo.h"
#include "kernels.h"
using namespace mesheR;

rowvec smooth_field_at_point_with_distance(const mat &displacement, const rowvec &distance, ScalarValuedKernel<rowvec>* scalarkernel, double gamma, const std::vector<float>  &use = std::vector<float>());

rowvec smooth_field_at_point(const rowvec &pt, const mat &clostPoints, const mat &displacement, ScalarValuedKernel<rowvec>* scalarkernel, double gamma, const std::vector<float>  &use = std::vector<float>());
