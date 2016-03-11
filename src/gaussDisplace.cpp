#include "RcppArmadillo.h"
#include "smooth_field.h"

#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

typedef unsigned int uint;


RcppExport SEXP displaceGauss(SEXP points_,SEXP clostPointsM_, SEXP D1_, SEXP D2_, SEXP sigma_, SEXP gamma_, SEXP clIW_, SEXP clIWdistances_, SEXP clIP_, SEXP clIPdistances_, SEXP tol_, SEXP rt0_, SEXP rt1_, SEXP rc_, SEXP oneway_,SEXP smoothtype_, SEXP threads_= wrap(1)){
  try {
    mat points = as<mat>(points_);
    mat out = points*0;
    mat toward = out*0;
    mat back = out*0;
    mat clostPointsW = points;
    mat clostPointsM = as<mat>(clostPointsM_);
   
    // displacement field D1 = towards, D2 = back
    mat D1 = as<mat>(D1_);
    mat D2 = as<mat>(D2_);
    double sigma = as<double>(sigma_);
    double gamma = as<double>(gamma_);
    imat armaclIW = as<imat>(clIW_);
    mat clIWdistances = as<mat>(clIWdistances_);
    imat armaclIP = as<imat>(clIP_);
    mat clIPdistances = as<mat>(clIPdistances_);
    // results from angle checking
    NumericVector rt0(rt0_);
    NumericVector rt1(rt1_);
    double rc = as<double>(rc_);
    double tol = as<double>(tol_);
    bool oneway = as<bool>(oneway_);
    int smoothtype = as<int>(smoothtype_);
    int threads = as<int>(threads_);
    vec use1(D1.n_rows);
    use1.fill(1);
    for (uint i=0; i < D1.n_rows;i++) {
      if (tol > 0 && sum(square(D1.row(i))) > tol)
	use1[i] = 0;
      if (rc > 0 &&  rt0[i] > rc)
	use1[i] = 0;
    }
    vec use2(D2.n_rows);
    if (!oneway) {
      use2.fill(1);
      for (uint i=0; i < D2.n_rows;i++) {
	if (tol > 0 && sum(square(D2.row(i))) > tol)
	  use2[i] = 0;
	if (rc > 0 &&  rt1[i] > rc)
	  use2[i] = 0;
      }
    }
    ScalarValuedKernel<rowvec>* gk;
    if (smoothtype == 0)
      gk = new mesheRrow::GaussianKernel(sigma);
    else if (smoothtype == 1)
      gk = new mesheRrow::LaplacianKernel(sigma);
    else if (smoothtype == 2)
      gk = new mesheRrow::ExponentialKernel(sigma);
    else if (smoothtype == 3)
      gk = new mesheRrow::BsplineKernel(sigma);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (uint i = 0; i < out.n_rows; i++) {
      rowvec direction_towards(3), direction_backwards(3);
      direction_towards.zeros();direction_backwards.zeros();
      rowvec pt = points.row(i);
      uvec tmpW = conv_to<uvec>::from(armaclIW.row(i));
      
      std::vector<float> use1tmp = conv_to< std::vector<float> >::from(use1(tmpW));
      if (gk->CanUseDistance()) {
	direction_towards = smooth_field_at_point_with_distance(D1.rows(tmpW),clIWdistances.row(i),gk, gamma,use1tmp);
	if (!oneway) {
	  uvec tmpP = conv_to<uvec>::from(armaclIP.row(i));
	  std::vector<float> use2tmp = conv_to< std::vector<float> >::from(use2(tmpP));
	  direction_backwards = smooth_field_at_point_with_distance(D2.rows(tmpP),clIPdistances.row(i),gk, gamma,use2tmp);
	}
      } else {
	direction_towards = smooth_field_at_point(pt,clostPointsW.rows(tmpW),D1.rows(tmpW),gk, gamma,use1tmp);
	if (!oneway) {
	  uvec tmpP = conv_to<uvec>::from(armaclIP.row(i));
	  std::vector<float> use2tmp = conv_to< std::vector<float> >::from(use2(tmpP));
	  direction_backwards = smooth_field_at_point(pt,clostPointsM.rows(tmpP),D2.rows(tmpP),gk, gamma,use2tmp);
	}
      }
      out.row(i) = direction_towards-direction_backwards;
    
    }
    delete gk;
    return wrap(out);
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
    
