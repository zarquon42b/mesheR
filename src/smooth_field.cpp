#include "smooth_field.h"

rowvec smooth_field_at_point_with_distance(const mat &displacement, const rowvec &distance, ScalarValuedKernel<rowvec>* scalarkernel, double gamma, int smoothtype, const std::vector<float>  &use) {
  rowvec outpt(3); outpt.zeros();
  float g1sum = 0;
  for (uint i = 0; i < distance.size(); i++) {
    float g1tmp = 0;
    
    g1tmp = scalarkernel->getWeightFromDistance(distance[i]);
    if (use.size() != distance.size()) {
      outpt += g1tmp*displacement.row(i);
      g1sum += g1tmp;
    } else {
      outpt += g1tmp*displacement.row(i)*use[i];
      g1sum += g1tmp*use[i];
    }
  }
  if (g1sum > 0)
    outpt /= g1sum;

  outpt /= gamma;
  return outpt;
}
rowvec smooth_field_at_point (const rowvec &pt, const mat &clostPoints, const mat &displacement,ScalarValuedKernel<rowvec>* scalarkernel, double gamma, const std::vector<float>  &use) {
  rowvec outpt(3); outpt.zeros();
  float g1sum = 0;
  for (uint i = 0; i < clostPoints.n_rows; i++) {
    float g1tmp = 0;
    g1tmp = scalarkernel->getWeightFromVectors(pt,clostPoints.row(i));
    if (use.size() != clostPoints.n_rows) {
      outpt += g1tmp*displacement.row(i);
      g1sum += g1tmp;
    } else {
      outpt += g1tmp*displacement.row(i)*use[i];
      g1sum += g1tmp*use[i];
    }
  }
  if (g1sum > 0)
    outpt /= g1sum;
  
  outpt /= gamma;
  return outpt;
}
