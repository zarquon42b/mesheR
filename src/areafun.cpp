#include "RcppArmadillo.h"
#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

//#include "angcheck.h"
using namespace Rcpp;
using namespace arma;
vec armacross(vec& x, vec& y) {
  vec z(3); z.fill(1);
  
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
  return z;
}



RcppExport SEXP areafun(SEXP A_) {
  try {
    mat armaA = as<mat>(A_);
    uvec v02; v02 << 0 << 1 << 2 << endr;
    uvec v35; v35 << 3 << 4 << 5 << endr;
    uvec v68; v68 << 6 << 7 << 8 << endr;
    uvec v911; v911 << 9 << 10 << 11 << endr;
    vec out(armaA.n_rows);
#pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < armaA.n_rows;i++) {
    uvec ui(1); ui[0] = i;
    mat ac = armaA.submat(ui,v02) - armaA.submat(ui,v35) ;
    mat bd = armaA.submat(ui,v68) - armaA.submat(ui,v911);
        
    vec acvec = conv_to<vec>::from(ac.row(0));
    vec bdvec = conv_to<vec>::from(bd.row(0));
    vec z = armacross(acvec,bdvec);
    out[i] = 0.5*arma::norm(z,2);
    // ac = ac.elem(v02);
  }
  return wrap(out);
}  catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
    return wrap(1);
  }
}
