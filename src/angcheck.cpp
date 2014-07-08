#include "Rcpp.h"
using namespace Rcpp;

double dot(Rcpp::NumericVector c,Rcpp::NumericVector d)
   {     int nc= c.size();
         double e=0;
         for (int i = 0; i < nc; i++)
         {  e=e+c(i)*d(i); }
         return e;
   }

double norm(Rcpp::NumericVector c)
   {     
     double e = sqrt(dot(c,c));
     
         return e;
   }


double anglecalc(NumericVector a, NumericVector b) {
   double alen = norm(a);
   double blen = norm(b);
   if (alen > 0)
     a = a/alen;
   if (blen > 0)
     b = b/blen;
   NumericVector diffvec = a -b;
   double angle = acos((dot(diffvec,diffvec)-2)/-2);
   return angle;
}


RcppExport SEXP angcheck(SEXP mat1_, SEXP mat2_) {
  try {
    NumericMatrix mat1(mat1_);
    NumericMatrix mat2(mat2_);
    NumericVector angles(mat1.ncol());
    for (int i = 0; i < mat1.ncol();i++) {
      angles[i] = anglecalc(mat1(_,i), mat2(_,i));
    }
    
    return angles;
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
    return wrap(1);
  }
}

