#include "RcppArmadillo.h"

#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double smooth_weight(double distance, double sigma) {
  double weight = 0;
  //if (distance > 0)
    weight = exp(-distance/sigma);
    //else
    //weight = 1e12;
  return weight;
}

rowvec smooth_field_at_point(mat displacement, rowvec distance, double sigma, double gamma) {
  rowvec outpt(3); outpt.zeros();
  float g1sum = 0;
  for (uint i = 0; i < distance.size(); i++) {
    float g1tmp = smooth_weight(distance[i],sigma);
    outpt += g1tmp*displacement.row(i);
    g1sum += g1tmp;
  }
  if (g1sum > 0)
    outpt /= g1sum;

  outpt /= gamma;
  return outpt;
}

RcppExport SEXP smoothField(SEXP evaluatePoints_, SEXP displacement_, SEXP sigma_, SEXP gamma_, SEXP closestInds_, SEXP distances_, SEXP threads_= wrap(1)){
  mat evaluatePoints = as<mat>(evaluatePoints_);
  mat displacement = as<mat>(displacement_);
  double sigma = as<double>(sigma_);
  double gamma = as<double>(gamma_);
  imat closestInds = as<imat>(closestInds_);
  mat distances = as<mat>(distances_);
  int threads = as<int>(threads_);
   mat out = evaluatePoints;
#ifdef SUPPORT_OPENMP
  omp_set_num_threads(threads);
#endif
#pragma omp parallel for schedule(static)
  for (uint i = 0; i < evaluatePoints.n_rows; i++) {
    //rowvec pt = iomat.row(i);
    uvec tmpW = conv_to<uvec>::from(closestInds.row(i));
    //rowvec distances_at_pt = distances.row(i);
    out.row(i) = smooth_field_at_point(displacement.rows(tmpW),distances.row(i),sigma, gamma);
    
  }
  
  return wrap(out);
}
