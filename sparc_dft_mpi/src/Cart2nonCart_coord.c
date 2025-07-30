#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"

/**
 * @brief: function to convert cartesian to non cartesian coordinates
 */
void Cart2nonCart_coord(const SPARC_OBJ *pSPARC, double *x, double *y, double *z) {
  double x1, x2, x3;
  x1 = pSPARC->gradT[0] * (*x) + pSPARC->gradT[1] * (*y) + pSPARC->gradT[2] * (*z);
  x2 = pSPARC->gradT[3] * (*x) + pSPARC->gradT[4] * (*y) + pSPARC->gradT[5] * (*z);
  x3 = pSPARC->gradT[6] * (*x) + pSPARC->gradT[7] * (*y) + pSPARC->gradT[8] * (*z);
  *x = x1;
  *y = x2;
  *z = x3;
}

