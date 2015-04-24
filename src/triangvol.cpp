#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

typedef unsigned int uint;

RcppExport SEXP trianvol(SEXP vb1_, SEXP vb2_, SEXP it_) {
  double V = 0.0;
  
  try {
  NumericMatrix vb1(vb1_);
  NumericMatrix vb2(vb2_);
  mat armavb1(vb1.begin(),vb1.nrow(),vb1.ncol());
  mat armavb2(vb2.begin(),vb2.nrow(),vb2.ncol());
  IntegerMatrix it(it_);
  imat armacIT(it.begin(),it.nrow(),it.ncol());
  umat uIT = conv_to<umat>::from(armacIT);
  mat allvol(4,6);
  umat subvols(3,4);
  uvec tmp; tmp  << 0 << 1 << 2 << 3 << endr;
  subvols.row(0) = tmp.t();
  tmp  << 0 << 1 << 3 << 4 << endr;
  subvols.row(1) = tmp.t();
  tmp  <<  1 << 3 << 4 << 5 << endr;
  subvols.row(2) = tmp.t();
  uvec tmp02; tmp02 << 0 << 1 << 2 << endr;
  uvec tmp453; tmp453 << 4 << 5 << 3 << endr;
  int dimit = it.ncol();
  for (unsigned int i = 0; i < dimit; i++) {
    double Vtmp = 0.0;
    uvec slice0 = uIT.col(i);
    allvol.cols(tmp02) = armavb1.cols(uIT.col(i));
    allvol.cols(tmp453) = armavb2.cols(uIT.col(i));
    for (unsigned int j = 0; j < 3; j++) {
      mat allvolsub = allvol.cols(subvols.row(j));
      allvolsub = allvolsub.t();
      double tmpdet = det(allvolsub);
      Vtmp += std::abs(tmpdet);
    }
    V += Vtmp/6;
  }
  } catch (...)
  {
  V=NA_REAL;
  }
  
  return wrap(V);
}
