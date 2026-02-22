#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "initialization.h"
#include "tools.h"

/**
 * @brief   Find equivalent mesh size to a given Ecut.
 *
 * @param Ecut  Energy cutoff used in plane-wave codes, in Hartree.
 * @param FDn   Finite difference order divided by 2.
 */
double Ecut2h(double Ecut, int FDn) {
  double epsilon = 0.1;
  double *w2;
  w2 = (double *)malloc((FDn + 1) * sizeof(double));

  // 2nd derivative weights
  w2[0] = 0;
  for (int p = 1; p < FDn + 1; p++) {
    w2[0] -= (2.0 / (p * p));
    w2[p] = (2 * (p % 2) - 1) * 2 * fract(FDn, p) / (p * p);
  }

  // k grid within interval (0,pi]
  int N = 1000;
  double dk = M_PI / (double)N;
  double *kk = (double *)malloc(N * sizeof(double));
  double *y_fd = (double *)malloc(N * sizeof(double));
  double *k2 = (double *)malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    kk[i] = (i + 1) * dk;
    k2[i] = kk[i] * kk[i];
  }

  // y_fd(k) = -w[0] + sum_{p=1}^{FDn} -2*w[p]*cos(k*p) (assuming h = 1)
  for (int i = 0; i < N; i++) {
    y_fd[i] = -w2[0];
    for (int p = 1; p < FDn + 1; p++) {
      y_fd[i] -= 2 * w2[p] * cos(kk[i] * p);
    }
  }

  // find out at which point |k^2 - y_fd| > epsilon
  int ind = 0;
  for (int i = 0; i < N; i++) {
    if (fabs(y_fd[i] - k2[i]) > epsilon) {
      ind = i;
      break;
    }
  }

  // k_cutoff
  double k_cutoff = kk[ind];

  double h = k_cutoff / sqrt(2.0 * Ecut);

  free(kk);
  free(y_fd);
  free(k2);
  free(w2);

  return h;
}

