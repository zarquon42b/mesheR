#include "RcppArmadillo.h"

#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

typedef unsigned int uint;

double gauss_smooth(rowvec pt, rowvec tarpt, double sigma){
  rowvec diffvec = pt - tarpt;
  double k = arma::dot(diffvec, diffvec);
  double out = exp(-k/sigma);
  return out;
}

  
rowvec relax_pt(rowvec pt, mat Wvb, mat  Pvb, mat D1, mat D2, double sigma, double gamma, bool oneway, vec dif1, vec dif2) {
  
  rowvec outpt(3); outpt.zeros();
  rowvec outpt2(3); outpt2.zeros();
  float g1sum = 0;

  float g2sum = 0;
  
  for (uint i = 0; i < Wvb.n_rows; i++) {
    float g1tmp = gauss_smooth(pt, Wvb.row(i),sigma);
    outpt += g1tmp*D1.row(i)*dif1(i);
    g1sum += g1tmp*dif1(i);

  }
  if (g1sum > 0)
    outpt /= g1sum;
  if (!oneway) {
    for (uint i = 0; i < Pvb.n_rows; i++) {
      float g2tmp = gauss_smooth(pt, Pvb.row(i),sigma);
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


RcppExport SEXP displaceGauss(SEXP iomat_, SEXP Wvb_, SEXP Pvb_, SEXP D1_, SEXP D2_, SEXP sigma_, SEXP gamma_, SEXP clIW_, SEXP clIP_, SEXP tol_, SEXP rt0_, SEXP rt1_, SEXP rc_, SEXP oneway_,SEXP threads_= wrap(1)){
  try {
  NumericMatrix iomat(iomat_);
  NumericMatrix Wvb(Wvb_);
  mat armaWvb(Wvb.begin(), Wvb.nrow(), Wvb.ncol());
  NumericMatrix Pvb(Pvb_);
  mat armaPvb(Pvb.begin(), Pvb.nrow(), Pvb.ncol());
  //int cores = as<int>(cores_);

  NumericMatrix D1(D1_);
  mat armaD1(D1.begin(), D1.nrow(), D1.ncol());
  
  NumericMatrix D2(D2_);
  mat armaD2(D2.begin(), D2.nrow(), D2.ncol());
  
  double sigma = as<double>(sigma_);
  double gamma = as<double>(gamma_);
  IntegerMatrix clIW(clIW_);
  imat armaclIW(clIW.begin(),clIW.nrow(),clIW.ncol());
 
  IntegerMatrix clIP(clIP_);
  imat armaclIP(clIP.begin(),clIP.nrow(),clIP.ncol());
  NumericVector rt0(rt0_);
  NumericVector rt1(rt1_);
  double rc = as<double>(rc_);
  double tol = as<double>(tol_);
  bool oneway = as<bool>(oneway_);
  int threads = as<int>(threads_);
  //omp_set_num_threads(8);
  vec diff1(D1.nrow());
  for (uint i=0; i < D1.nrow();i++) {
    diff1[i] = sum(D1(i,_)*D1(i,_));
    if (tol > 0 && diff1[i] > tol)
      diff1[i] = 0;
    if (rc > 0 &&  rt0[i] > rc)
      diff1[i] = 0;
  }
  
  vec diff2(D2.nrow());
  for (uint i=0; i < D2.nrow();i++) {

    diff2[i] = sum(D2(i,_)*D2(i,_));
    if (tol > 0 && diff2[i] > tol)
      diff2[i] = 0;
    if (rc > 0 &&  rt1[i] > rc)
      diff2[i] = 0;
  }
  mat out(iomat.begin(), iomat.nrow(), iomat.ncol());
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(threads);
#endif
#pragma omp parallel for schedule(static)
  for (uint i = 0; i < iomat.nrow(); i++) {
    rowvec pt = iomat(i,_);
    uvec tmpW = conv_to<uvec>::from(armaclIW.row(i));
    uvec tmpP = conv_to<uvec>::from(armaclIP.row(i));

    out.row(i) = relax_pt(pt, armaWvb.rows(tmpW),armaPvb.rows(tmpP),armaD1.rows(tmpW),armaD2.rows(tmpP),sigma, gamma, oneway,diff1(tmpW),diff2(tmpP));
    
  }
  
 
  return wrap(out);
 }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
    
