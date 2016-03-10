#ifndef __KERNELS_H
#define __KERNELS_H
#include "RcppArmadillo.h"
#include "bspline.h"

namespace mesheR {
  template<class TPoint>
  class ScalarValuedKernel {
  public:

    /**
     * Create a new scalar valued kernel.
     */
    ScalarValuedKernel() {	}

    virtual ~ScalarValuedKernel() {
    }

    /**
     * Evaluate the kernel function at the points x and y
     */
    virtual double getWeightFromVectors(const TPoint& x, const TPoint& y) const = 0;

    /* Evaluate the kernel function based on a precomputed distance
     */
    virtual double getWeightFromDistance(const double& distance) const = 0;
    
    /**
     * Return whether the kernel supports distance based functions
     */
    virtual bool CanUseDistance() const = 0;

  };
}
#endif // __KERNELS_H


using namespace Rcpp;
using namespace mesheR;
using namespace arma;
typedef rowvec VectorType;
class GaussianKernel: public ScalarValuedKernel<rowvec> {
public:
  GaussianKernel(double sigma) : m_sigma(sigma), m_sigma2(2*sigma * sigma){
  }

  inline double getWeightFromVectors(const rowvec& x, const rowvec& y) const {
    VectorType r = x-y;
    return exp(-dot(r,r) / m_sigma2);
  }
  inline double getWeightFromDistance(const double& distance) const {
    double weight = exp(-(distance*distance)/m_sigma2);
    return weight;
  }
  inline bool CanUseDistance() const {
    return m_CanUseDistance;
  }
  
private:

  double m_sigma;
  double m_sigma2;
  const static bool m_CanUseDistance = true;
};

class LaplacianKernel: public ScalarValuedKernel<rowvec> {
public:
  LaplacianKernel(double sigma) : m_sigma(sigma), m_sigma2(sigma * sigma){
  }

  inline double getWeightFromVectors(const rowvec& x, const rowvec& y) const {
    VectorType r = x-y;
    return exp(-sqrt(dot(r,r)) / m_sigma);
  }
  inline double getWeightFromDistance(const double& distance) const {
    double weight = exp(-distance/m_sigma);
    return weight;
  }
  inline bool CanUseDistance() const {
    return m_CanUseDistance;
  }
  
  
private:

  double m_sigma;
  double m_sigma2;
  const static bool m_CanUseDistance=true;
};

class ExponentialKernel: public ScalarValuedKernel<rowvec> {
public:
  ExponentialKernel(double sigma) : m_sigma(sigma), m_sigma2(2*sigma * sigma){
  }

  inline double getWeightFromVectors(const rowvec& x, const rowvec& y) const {
    VectorType r = x-y;
    return exp(-sqrt(dot(r,r)) / m_sigma2);
  }
  inline double getWeightFromDistance(const double& distance) const {
    double weight = exp(-distance/m_sigma2);
    return weight;
  }
  
  inline bool CanUseDistance() const {
    return m_CanUseDistance;
  }
  
  
private:

  double m_sigma;
  double m_sigma2;
  const static bool m_CanUseDistance=true;

};

class BsplineKernel: public ScalarValuedKernel<rowvec> {
public:
  BsplineKernel(double sigma) : m_sigma(sigma) {
  }

  inline double getWeightFromVectors(const rowvec& x, const rowvec& y) const {
    double weight = get_bspline(x,y,m_sigma);
    return weight;
  }
  inline double getWeightFromDistance(const double& distance) const {
    return 0;
  }
  inline bool CanUseDistance() const {
    return m_CanUseDistance;
  }
  
  
private:

  double m_sigma;
  const static bool m_CanUseDistance=false;
};
