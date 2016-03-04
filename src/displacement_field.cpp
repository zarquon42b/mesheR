#include "RcppArmadillo.h"
#include "bspline.h"
#include "displacement_field.h"
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// distance is the distance between two points and sigma the bandwidth of the kernel
double smooth_gaussian(double distance, double sigma) {
  sigma = 2*sigma*sigma;
  distance *= distance;
  double weight = exp(-distance/sigma);
  //else
  //weight = 1e12;
  return weight;
}

//exponential kernel
double smooth_exponential(double distance, double sigma) {
  sigma = 2*sigma*sigma;
  double weight = exp(-distance/sigma);
  return weight;
}

//laplacian kernel
double smooth_laplacian(double distance, double sigma) {
  double weight = exp(-distance/sigma);
  return weight;
}



rowvec smooth_field_at_point(mat displacement, rowvec distance, double sigma, double gamma, int smoothtype,  std::vector<float>  use) {
  rowvec outpt(3); outpt.zeros();
  float g1sum = 0;
  for (uint i = 0; i < distance.size(); i++) {
    float g1tmp = 0;
    if (smoothtype == 1)
      g1tmp = smooth_laplacian(distance[i],sigma);
    else if (smoothtype == 2) 
      g1tmp = smooth_exponential(distance[i],sigma);
    else
      g1tmp = smooth_gaussian(distance[i],sigma);
    if (use.size() != distance.size()) {
      outpt += g1tmp*displacement.row(i);
      g1sum += g1tmp;
    } else {
      outpt += g1tmp*displacement.row(i)*use[i];
      g1sum += g1tmp*use[i];
    }
  }
  if (g1sum > 0)
    outpt /= g1sum;

  outpt /= gamma;
  return outpt;
}

rowvec smooth_field_at_point_bspline(rowvec pt, mat clostPoints, mat displacement,double support, double gamma,  std::vector<float>  use) {
  rowvec outpt(3); outpt.zeros();
  float g1sum = 0;
  for (uint i = 0; i < clostPoints.n_rows; i++) {
    float g1tmp = 0;
    g1tmp = smooth_bspline(pt,clostPoints.row(i),support);
    if (use.size() != clostPoints.n_rows) {
      outpt += g1tmp*displacement.row(i);
      g1sum += g1tmp;
    } else {
      outpt += g1tmp*displacement.row(i)*use[i];
      g1sum += g1tmp*use[i];
    }
  }
  if (g1sum > 0)
    outpt /= g1sum;
  
  outpt /= gamma;
  return outpt;
}

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

    for (int iter = 0; iter < iterations; iter++) {
#pragma omp parallel for schedule(static) num_threads(threads)
      for (uint i = 0; i < evaluatePoints.n_rows; i++) {
	rowvec pt = evaluatePoints.row(i);
	uvec tmpW = conv_to<uvec>::from(closestInds.row(i));
	//rowvec distances_at_pt = distances.row(i);
	if (smoothtype != 3)
	  out.row(i) = smooth_field_at_point(displacement.rows(tmpW),distances.row(i),sigma, gamma,smoothtype);
	else
	  out.row(i) = smooth_field_at_point_bspline(pt, clostPoints.rows(tmpW),displacement.rows(tmpW),sigma, gamma);
      }
      displacement = out;
    }

    return wrap(out);

  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
    
