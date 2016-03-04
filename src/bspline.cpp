#include "bspline.h"

double smooth_bspline(rowvec pt1, rowvec pt2, double support) {
  double out = get_bspline(pt1,pt2,support);
  return out;
}

double bspline3(const double& x) {
    const double absX = std::fabs(x);
    const double absXSquared = absX * absX;
    const double absXCube = absXSquared * absX;
    const double twoMinAbsX = 2.0 - absX;
    const double twoMinAbsXCube = twoMinAbsX * twoMinAbsX * twoMinAbsX;
    const double twoByThree = 2.0 / 3.0;

    double splineValue = 0;
    if (absX >= 0 && absX < 1) {
      splineValue = twoByThree - absXSquared + 0.5 * absXCube;
    }
    else if (absX >= 1 && absX < 2) {
      splineValue = twoMinAbsXCube / 6.0;
    }
    else {
      splineValue = 0;
    }
    return splineValue;
  }
double tensorProductSpline(const rowvec& x) {
    double prod = 1;
    for (unsigned d = 0; d < 3; ++d) {
      prod *= bspline3(x[d]);
    }
    return prod;
  }
double get_bspline(const rowvec& x, const rowvec& y, const double support) {

    const unsigned dim =3;
    const double supportBasisFunction = 4.0;
    const double scale = -1.0 * std::log(support / supportBasisFunction) / std::log(2.0);
      
    rowvec xScaled = x * std::pow(2.0, scale);
    rowvec yScaled = y * std::pow(2.0, scale);
    
    std::vector<int> kLower(dim);
    std::vector<int> kUpper(dim);
    for (unsigned d = 0; d < dim; ++d) {
      kLower[d] = static_cast<int>(std::ceil(std::max(xScaled[d], yScaled[d]) - 0.5 * supportBasisFunction));
      kUpper[d] = static_cast<int>(std::floor(std::min(xScaled[d], yScaled[d]) + 0.5 * supportBasisFunction));
    }
      

    // We need to generate the cartesian product k_1 x ... x k_d, where k_i goes through all the integers
    // within the given bounds. A non-recursive solution requires d loops. Here we just write down the cases
    // for 1 2 and 3D
    double sum = 0.0;
    double kx = kLower[0];
    while (kx <= kUpper[0]) {
      double ky = kLower[1];
      while (ky <= kUpper[1]) {
	double kz = kLower[2];
	while (kz <= kUpper[2]) {
	  rowvec k(3);
	  k[0] = kx; k[1] = ky; k[2] = kz;
	  sum += (tensorProductSpline(xScaled - k) * tensorProductSpline(yScaled- k));
	  kz += 1;
	}
	ky += 1;
      }
      kx += 1;
    }
    
    return sum;
   
}
