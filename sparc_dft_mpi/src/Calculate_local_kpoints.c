#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"

void Calculate_local_kpoints(SPARC_OBJ *pSPARC) {
  // Initialize the weights of k points
  // by going over all the k points from the Monkhorst Pack grid
  int k, kstart, kend;
  kstart = pSPARC->kpt_start_indx;
  kend = pSPARC->kpt_end_indx;
  int nk;
  k = 0;
  for (nk = kstart; nk <= kend; nk++) {
    pSPARC->k1_loc[k] = pSPARC->k1[nk];
    pSPARC->k2_loc[k] = pSPARC->k2[nk];
    pSPARC->k3_loc[k] = pSPARC->k3[nk];
    pSPARC->kptWts_loc[k] = pSPARC->kptWts[nk];
    k++;
  }
}
