#include "RcppArmadillo.h"
#include "bspline.h"
#include "smooth_field.h"
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


RcppExport SEXP smoothField(SEXP evaluatePoints_, SEXP clostPoints_, SEXP displacement_, SEXP sigma_, SEXP gamma_, SEXP closestInds_, SEXP distances_, SEXP iterations_, SEXP threads_= wrap(1), SEXP smoothtype_ = wrap(0)){
  try {
    mat evaluatePoints = as<mat>(evaluatePoints_);
    mat clostPoints = as<mat>(clostPoints_);
    mat displacement = as<mat>(displacement_);
    double sigma = as<double>(sigma_);
    double gamma = as<double>(gamma_);
    imat closestInds = as<imat>(closestInds_);
    mat distances = as<mat>(distances_);
    int iterations = as<int>(iterations_);
    int threads = as<int>(threads_);
    int smoothtype = as<int>(smoothtype_);
    mat out = evaluatePoints;
    ScalarValuedKernel<rowvec>* gk;
    if (smoothtype == 0)
      gk = new GaussianKernel(sigma);
    else if (smoothtype == 1)
      gk = new LaplacianKernel(sigma);
    else if (smoothtype == 2)
      gk = new ExponentialKernel(sigma);
    else if (smoothtype == 3)
      gk = new BsplineKernel(sigma);

    for (int iter = 0; iter < iterations; iter++) {
#pragma omp parallel for schedule(static) num_threads(threads)
      for (uint i = 0; i < evaluatePoints.n_rows; i++) {
	rowvec pt = evaluatePoints.row(i);
	uvec tmpW = conv_to<uvec>::from(closestInds.row(i));
	//rowvec distances_at_pt = distances.row(i);
	if (gk->CanUseDistance())
	  out.row(i) = smooth_field_at_point_with_distance(displacement.rows(tmpW),distances.row(i),gk, gamma,smoothtype);
	else
	  out.row(i) = smooth_field_at_point(pt, clostPoints.rows(tmpW),displacement.rows(tmpW),gk, gamma);
      }
      displacement = out;
    }
    delete gk;
    return wrap(out);

  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
    
