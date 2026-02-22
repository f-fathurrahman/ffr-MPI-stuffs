#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"


/**
 * @brief: function to convert gradients along lattice directions to cartesian
 * gradients
 */
void nonCart2Cart_grad(SPARC_OBJ *pSPARC, double *x, double *y, double *z) {
  double x1, x2, x3;
  x1 = pSPARC->gradT[0] * (*x) + pSPARC->gradT[3] * (*y) + pSPARC->gradT[6] * (*z);
  x2 = pSPARC->gradT[1] * (*x) + pSPARC->gradT[4] * (*y) + pSPARC->gradT[7] * (*z);
  x3 = pSPARC->gradT[2] * (*x) + pSPARC->gradT[5] * (*y) + pSPARC->gradT[8] * (*z);
  *x = x1;
  *y = x2;
  *z = x3;
}

