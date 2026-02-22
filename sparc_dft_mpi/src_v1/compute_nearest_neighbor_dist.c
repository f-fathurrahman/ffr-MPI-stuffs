#include <stdlib.h>
#include <stdio.h>

#include "initialization.h"

/**
 * @brief   Computing nearest neighbohr distance
 */
double compute_nearest_neighbor_dist(SPARC_OBJ *pSPARC, char CorN) {
#ifdef DEBUG
  double t1, t2;
  t1 = MPI_Wtime();
#endif
  int atm1, atm2;
  double nn, dist = 0.0;
  nn = 100000000000;
  if (CorN == 'N') { // Non-Cartesian coordinates
    for (atm1 = 0; atm1 < pSPARC->n_atom - 1; atm1++) {
      for (atm2 = atm1 + 1; atm2 < pSPARC->n_atom; atm2++) {
        CalculateDistance(pSPARC, pSPARC->atom_pos[3 * atm1], pSPARC->atom_pos[3 * atm1 + 1],
                          pSPARC->atom_pos[3 * atm1 + 2], pSPARC->atom_pos[3 * atm2], pSPARC->atom_pos[3 * atm2 + 1],
                          pSPARC->atom_pos[3 * atm2 + 2], &dist);
        if (dist < nn)
          nn = dist;
      }
    }
  } else if (CorN == 'C') { // Cartesian coordinates
    for (atm1 = 0; atm1 < pSPARC->n_atom - 1; atm1++) {
      for (atm2 = atm1 + 1; atm2 < pSPARC->n_atom; atm2++) {
        dist = fabs(sqrt(pow(pSPARC->atom_pos[3 * atm1] - pSPARC->atom_pos[3 * atm2], 2.0) +
                         pow(pSPARC->atom_pos[3 * atm1 + 1] - pSPARC->atom_pos[3 * atm2 + 1], 2.0) +
                         pow(pSPARC->atom_pos[3 * atm1 + 2] - pSPARC->atom_pos[3 * atm2 + 2], 2.0)));
        if (dist < nn)
          nn = dist;
      }
    }
  } else {
    printf("ERROR: please use 'N' for non-cartesian coordinates and 'C' for "
           "cartesian coordinates in compute_nearest_neighbor_dist function.");
    exit(-1);
  }

#ifdef DEBUG
  t2 = MPI_Wtime();
  printf("\nComputing nearest neighbor distance (%.3f Bohr) takes %.3f ms\n", nn, (t2 - t1) * 1000);
#endif
  return nn;
}
