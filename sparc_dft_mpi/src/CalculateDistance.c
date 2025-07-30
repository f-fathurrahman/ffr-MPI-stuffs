#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"

/**
 * @brief: function to calculate the distance btween two points
 */
void CalculateDistance(SPARC_OBJ *pSPARC, double x, double y, double z, double xref, double yref, double zref,
                       double *d) {
  if (pSPARC->cell_typ == 0) {
    *d = sqrt(pow((x - xref), 2.0) + pow((y - yref), 2.0) + pow((z - zref), 2.0));
  } else if (pSPARC->cell_typ > 10 && pSPARC->cell_typ < 20) {
    double xx = x - xref;
    double yy = y - yref;
    double zz = z - zref;
    *d = sqrt(pSPARC->metricT[0] * (xx * xx) + pSPARC->metricT[1] * (xx * yy) + pSPARC->metricT[2] * (xx * zz) +
              pSPARC->metricT[4] * (yy * yy) + pSPARC->metricT[5] * (yy * zz) + pSPARC->metricT[8] * (zz * zz));
  }
}

