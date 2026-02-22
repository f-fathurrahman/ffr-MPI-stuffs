#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "initialization.h"

/**
 * @brief   Call Spline to calculate derivatives of the tabulated functions and
 *          store them for later use (during interpolation).
 */
void Calculate_SplineDerivRadFun(SPARC_OBJ *pSPARC) {
  int ityp, l, lcount, lcount2, np, ppl_sum, psd_len;
  for (ityp = 0; ityp < pSPARC->Ntypes; ityp++) {
    int lloc = pSPARC->localPsd[ityp];
    psd_len = pSPARC->psd[ityp].size;
    pSPARC->psd[ityp].SplinerVlocD = (double *)malloc(sizeof(double) * psd_len);
    pSPARC->psd[ityp].SplineFitIsoAtomDen = (double *)malloc(sizeof(double) * psd_len);
    pSPARC->psd[ityp].SplineRhocD = (double *)malloc(sizeof(double) * psd_len);
    assert(pSPARC->psd[ityp].SplinerVlocD != NULL);
    assert(pSPARC->psd[ityp].SplineFitIsoAtomDen != NULL);
    assert(pSPARC->psd[ityp].SplineRhocD != NULL);
    getYD_gen(pSPARC->psd[ityp].RadialGrid, pSPARC->psd[ityp].rVloc, pSPARC->psd[ityp].SplinerVlocD, psd_len);
    getYD_gen(pSPARC->psd[ityp].RadialGrid, pSPARC->psd[ityp].rhoIsoAtom, pSPARC->psd[ityp].SplineFitIsoAtomDen,
              psd_len);
    getYD_gen(pSPARC->psd[ityp].RadialGrid, pSPARC->psd[ityp].rho_c_table, pSPARC->psd[ityp].SplineRhocD, psd_len);
    // note we neglect lloc
    ppl_sum = 0;
    for (l = 0; l <= pSPARC->psd[ityp].lmax; l++) {
      // if (l == pSPARC->localPsd[ityp]) continue; // this fails under -O3, -O2
      // optimization
      if (l == lloc)
        continue;
      ppl_sum += pSPARC->psd[ityp].ppl[l];
    }
    pSPARC->psd[ityp].SplineFitUdV = (double *)malloc(sizeof(double) * psd_len * ppl_sum);
    if (pSPARC->psd[ityp].SplineFitUdV == NULL) {
      printf("Memory allocation failed!\n");
      exit(EXIT_FAILURE);
    }
    for (l = lcount = lcount2 = 0; l <= pSPARC->psd[ityp].lmax; l++) {
      if (l == lloc) {
        lcount2 += pSPARC->psd[ityp].ppl[l];
        continue;
      }
      for (np = 0; np < pSPARC->psd[ityp].ppl[l]; np++) {
        // note that UdV is of size (psd_len, lmax+1), while SplineFitUdV has
        // size (psd_len, lmax)
        getYD_gen(pSPARC->psd[ityp].RadialGrid, pSPARC->psd[ityp].UdV + lcount2 * psd_len,
                  pSPARC->psd[ityp].SplineFitUdV + lcount * psd_len, psd_len);
        lcount++;
        lcount2++;
      }
    }
    if (pSPARC->psd[ityp].pspsoc) {
      ppl_sum = 0;
      for (l = 1; l <= pSPARC->psd[ityp].lmax; l++) {
        // if (l == pSPARC->localPsd[ityp]) continue; // this fails under -O3,
        // -O2 optimization
        if (l == lloc)
          continue;
        ppl_sum += pSPARC->psd[ityp].ppl_soc[l - 1];
      }
      pSPARC->psd[ityp].SplineFitUdV_soc = (double *)malloc(sizeof(double) * psd_len * ppl_sum);
      assert(pSPARC->psd[ityp].SplineFitUdV_soc != NULL);
      lcount = lcount2 = 0;
      for (l = 1; l <= pSPARC->psd[ityp].lmax; l++) {
        if (l == lloc) {
          lcount2 += pSPARC->psd[ityp].ppl_soc[l - 1];
          continue;
        }
        for (np = 0; np < pSPARC->psd[ityp].ppl_soc[l - 1]; np++) {
          // note that UdV is of size (psd_len, lmax+1), while SplineFitUdV has
          // size (psd_len, lmax)
          getYD_gen(pSPARC->psd[ityp].RadialGrid, pSPARC->psd[ityp].UdV_soc + lcount2 * psd_len,
                    pSPARC->psd[ityp].SplineFitUdV_soc + lcount * psd_len, psd_len);
          lcount++;
          lcount2++;
        }
      }
    }
  }
}

