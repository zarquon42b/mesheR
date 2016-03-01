#include "RcppArmadillo.h"
#include "displacement_field.h"
#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

typedef unsigned int uint;



  
rowvec relax_pt(mat D1, mat D2,rowvec distance1, rowvec distance2, double sigma, double gamma, bool oneway, vec dif1, vec dif2, int smoothtype) {
  
  rowvec outpt(3); outpt.zeros();
  rowvec outpt2(3); outpt2.zeros();
  float g1sum = 0;

  float g2sum = 0;
  
  for (uint i = 0; i < distance1.size(); i++) {
    float g1tmp = 0;
    if (smoothtype == 1)
      g1tmp = smooth_laplacian(distance1[i],sigma);
    else if (smoothtype == 2) 
      g1tmp = smooth_exponential(distance1[i],sigma);
    else
      g1tmp = smooth_gaussian(distance1[i],sigma);
    outpt += g1tmp*D1.row(i)*dif1(i);
    g1sum += g1tmp*dif1(i);

  }
  if (g1sum > 0)
    outpt /= g1sum;
  if (!oneway) {
    for (uint i = 0; i < distance2.size(); i++) {
      float g2tmp = 0;
    if (smoothtype == 1)
      g2tmp = smooth_laplacian(distance2[i],sigma);
    else if (smoothtype == 2) 
      g2tmp = smooth_exponential(distance2[i],sigma);
    else
      g2tmp = smooth_gaussian(distance2[i],sigma);
						   
    outpt2 += g2tmp*D2.row(i)*dif2(i);
    g2sum += g2tmp*dif2(i);
    }
     if (g2sum > 0)
       outpt2 /= g2sum;
     
     outpt = (outpt-outpt2)/gamma;
  }
  else 
    outpt /= gamma;
    
  return outpt;
    
}


RcppExport SEXP displaceGauss(SEXP iomat_, SEXP D1_, SEXP D2_, SEXP sigma_, SEXP gamma_, SEXP clIW_, SEXP clIWdistances_, SEXP clIP_, SEXP clIPdistances_, SEXP tol_, SEXP rt0_, SEXP rt1_, SEXP rc_, SEXP oneway_,SEXP type_, SEXP threads_= wrap(1)){
  try {
  mat out = as<mat>(iomat_);
  
    //int cores = as<int>(cores_);

  mat armaD1 = as<mat>(D1_);
  mat armaD2 = as<mat>(D2_);
  
  double sigma = as<double>(sigma_);
  double gamma = as<double>(gamma_);
  imat armaclIW = as<imat>(clIW_);
  mat clIWdistances = as<mat>(clIWdistances_);
  imat armaclIP = as<imat>(clIP_);
  mat clIPdistances = as<mat>(clIPdistances_);
  NumericVector rt0(rt0_);
  NumericVector rt1(rt1_);
  double rc = as<double>(rc_);
  double tol = as<double>(tol_);
  bool oneway = as<bool>(oneway_);
  int type = as<int>(type_);
  int threads = as<int>(threads_);
  //omp_set_num_threads(8);
  vec diff1(armaD1.n_rows);
  diff1.fill(1);
  for (uint i=0; i < armaD1.n_rows;i++) {
    //diff1[i] = sum(square(armaD1.row(i)));
    if (tol > 0 && diff1[i] > tol)
      diff1[i] = 0;
    if (rc > 0 &&  rt0[i] > rc)
      diff1[i] = 0;
  }
  
  vec diff2(armaD2.n_rows);
  diff2.fill(1);
  for (uint i=0; i < armaD2.n_rows;i++) {
    if (tol > 0 && diff2[i] > tol)
      diff2[i] = 0;
    if (rc > 0 &&  rt1[i] > rc)
      diff2[i] = 0;
  }
  
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(threads);
#endif
#pragma omp parallel for schedule(static)
  for (uint i = 0; i < out.n_rows; i++) {
    //rowvec pt = iomat.row(i);
    uvec tmpW = conv_to<uvec>::from(armaclIW.row(i));
    uvec tmpP = conv_to<uvec>::from(armaclIP.row(i));

    out.row(i) = relax_pt(armaD1.rows(tmpW),armaD2.rows(tmpP),clIWdistances.row(i), clIPdistances.row(i),sigma, gamma, oneway,diff1(tmpW),diff2(tmpP),type);
    
  }
  
 
  return wrap(out);
 }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
    
