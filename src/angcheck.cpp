#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

// double dot(vec c,vec d)
//    {     int nc= c.size();
//          double e=0;
//          for (int i = 0; i < nc; i++)
//          {  e=e+c(i)*d(i); }
//          return e;
//    }

// double norm(Rcpp::NumericVector c)
//    {     
//      double e = sqrt(dot(c,c));
     
//          return e;
//    }


double anglecalc(vec a, vec b) {
  a = normalise(a);
  b = normalise(b);
  vec diffvec = a -b;
   double angle = acos((dot(diffvec,diffvec)-2)/-2);
   return angle;
}


RcppExport SEXP angcheck(SEXP mat1_, SEXP mat2_, SEXP threads_) {
  try {
    mat mat1 = as<mat>(mat1_);
    mat mat2 = as<mat>(mat2_);
    int threads = as<int>(threads_);
    arma::vec angles(mat1.n_cols);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int i = 0; i < static_cast<int>(mat1.n_cols);i++) {
      angles(i) = anglecalc(mat1.col(i), mat2.col(i));
    }
    
    return wrap(angles);
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

