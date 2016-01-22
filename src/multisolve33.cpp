#ifndef ARMA_DONT_PRINT_ERRORS 
#define ARMA_DONT_PRINT_ERRORS 
#endif
#include "RcppArmadillo.h"
#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
RcppExport SEXP multisolve3(SEXP A_, SEXP trans_) {
  try {
    bool trans = as<bool>(trans_);
    mat armaA = as<mat>(A_);
    int k = armaA.n_rows;
    int iterlen = k/3;
    uvec v02; v02 << 0 << 1 << 2 << endr;
    mat outmat = armaA;
    // #pragma omp parallel for schedule(static)
   
    for (unsigned int i=0; i < iterlen; i++) {
      mat Asolve(3,3);
      Asolve.zeros();
      mat Aslice = armaA.rows(v02);
      
      if (trans) 
	Aslice = Aslice.t();
      bool check = inv(Asolve, Aslice);
      if (!check) {
	Aslice = armaA.rows(v02);
	Asolve.resize(3,3); Asolve.zeros();
	check = pinv(Asolve, Aslice);
	outmat.rows(v02) = Asolve;
      } else 
	outmat.rows(v02) = Asolve;
      
      v02 += 3;
    }
    return wrap(outmat);
  
  }  catch (std::exception& e) {
    //::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
    return wrap(1);
  }
}
    
