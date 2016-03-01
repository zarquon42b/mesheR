#include "RcppArmadillo.h"

double smooth_gaussian(double distance, double sigma);

double smooth_laplacian(double distance, double sigma);

double smooth_exponential(double distance, double sigma);
