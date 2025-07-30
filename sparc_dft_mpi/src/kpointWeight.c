/**
 * @brief   Calculate the weight of a given k-point.
 */
double kpointWeight(double kx, double ky, double kz) {
  /*
   * Appropriate weights for k-point sampling
   * we find the weight of the x-direction k point. If the kx has 0, then the
   * weight is 1.0 else the weight is 2.0
   */
  double w;
  if (fabs(kx) > 0) {
    w = 2.0;
  } else {
    w = 1.0;
  }
  return w;
}
