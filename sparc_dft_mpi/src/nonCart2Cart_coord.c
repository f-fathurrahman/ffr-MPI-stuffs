#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"

/**
 * @ brief: function to convert non cartesian to cartesian coordinates
 */
void nonCart2Cart_coord(const SPARC_OBJ *pSPARC, double *x, double *y, double *z) {
  double x1, x2, x3;
  x1 = pSPARC->LatUVec[0] * (*x) + pSPARC->LatUVec[3] * (*y) + pSPARC->LatUVec[6] * (*z);
  x2 = pSPARC->LatUVec[1] * (*x) + pSPARC->LatUVec[4] * (*y) + pSPARC->LatUVec[7] * (*z);
  x3 = pSPARC->LatUVec[2] * (*x) + pSPARC->LatUVec[5] * (*y) + pSPARC->LatUVec[8] * (*z);
  *x = x1;
  *y = x2;
  *z = x3;
}
