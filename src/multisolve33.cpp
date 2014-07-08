#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
RcppExport SEXP multisolve3(SEXP A_, SEXP trans_) {
  try {
    NumericMatrix A(A_);
    bool trans = as<bool>(trans_);
    mat armaA(A.begin(), A.nrow(), A.ncol());
    int k = A.nrow();
    int iterlen = k/3;
    uvec v02; v02 << 0 << 1 << 2 << endr;
    mat outmat = armaA;
    // #pragma omp parallel for schedule(static)
  
    for (unsigned int i=0; i < iterlen; i++) {
      mat Aslice = armaA.rows(v02);
      if (trans) 
	Aslice = Aslice.t();
      try {
	Aslice = inv(Aslice);
      } catch (std::runtime_error& e) {
	::Rf_warning( e.what());
	Aslice = Aslice*0;
      }
      outmat.rows(v02) = Aslice;
      v02 += 3;
    }
    return wrap(outmat);
  
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
    return wrap(1);
  }
}
    
