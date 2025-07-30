#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "initialization.h"

#define TEMP_TOL 1e-12

/**
 * @ brief: function to calculate cell type, the det(Jacobian) for integration
 *          and transformation matrices for distance, gradient, and laplacian
 *          in a non cartesian coordinate system.
 **/
void Cart2nonCart_transformMat(SPARC_OBJ *pSPARC) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int i, j, k;

  // Construct LatUVec;
  double mag;
  for (i = 0; i < 3; i++) {
    mag = sqrt(pow(pSPARC->LatVec[3 * i], 2.0) + pow(pSPARC->LatVec[3 * i + 1], 2.0) +
               pow(pSPARC->LatVec[3 * i + 2], 2.0));
    pSPARC->LatUVec[3 * i] = pSPARC->LatVec[3 * i] / mag;
    pSPARC->LatUVec[3 * i + 1] = pSPARC->LatVec[3 * i + 1] / mag;
    pSPARC->LatUVec[3 * i + 2] = pSPARC->LatVec[3 * i + 2] / mag;
  }

  // determinant of 3x3 Jacobian
  pSPARC->Jacbdet = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if (i != j && j != k && k != i)
          pSPARC->Jacbdet += ((i - j) * (j - k) * (k - i) / 2) * pSPARC->LatUVec[3 * i] * pSPARC->LatUVec[3 * j + 1] *
                             pSPARC->LatUVec[3 * k + 2];
      }
    }
  }

  if (pSPARC->Jacbdet <= 0) {
    if (rank == 0)
      printf("ERROR: Volume(det(jacobian)) %lf is <= 0\n", pSPARC->Jacbdet);
    exit(EXIT_FAILURE);
  }

  // transformation matrix for distance
  for (i = 0; i < 9; i++)
    pSPARC->metricT[i] = 0.0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        pSPARC->metricT[3 * i + j] += pSPARC->LatUVec[3 * i + k] * pSPARC->LatUVec[3 * j + k];
      }
    }
  }

  pSPARC->metricT[1] = 2 * pSPARC->metricT[1];
  pSPARC->metricT[2] = 2 * pSPARC->metricT[2];
  pSPARC->metricT[5] = 2 * pSPARC->metricT[5];

  // transformation matrix for gradient
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      pSPARC->gradT[3 * j + i] =
          (pSPARC->LatUVec[3 * ((j + 1) % 3) + (i + 1) % 3] * pSPARC->LatUVec[3 * ((j + 2) % 3) + (i + 2) % 3] -
           pSPARC->LatUVec[3 * ((j + 1) % 3) + (i + 2) % 3] * pSPARC->LatUVec[3 * ((j + 2) % 3) + (i + 1) % 3]) /
          pSPARC->Jacbdet;
    }
  }

  // transformation matrix for laplacian
  for (i = 0; i < 9; i++)
    pSPARC->lapcT[i] = 0.0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        pSPARC->lapcT[3 * i + j] += pSPARC->gradT[3 * i + k] * pSPARC->gradT[3 * j + k];
      }
    }
  }

  /* Different cell types for laplacian */
  if (fabs(pSPARC->lapcT[1]) > TEMP_TOL && fabs(pSPARC->lapcT[2]) < TEMP_TOL && fabs(pSPARC->lapcT[5]) < TEMP_TOL)
    pSPARC->cell_typ = 11;
  else if (fabs(pSPARC->lapcT[1]) < TEMP_TOL && fabs(pSPARC->lapcT[2]) > TEMP_TOL && fabs(pSPARC->lapcT[5]) < TEMP_TOL)
    pSPARC->cell_typ = 12;
  else if (fabs(pSPARC->lapcT[1]) < TEMP_TOL && fabs(pSPARC->lapcT[2]) < TEMP_TOL && fabs(pSPARC->lapcT[5]) > TEMP_TOL)
    pSPARC->cell_typ = 13;
  else if (fabs(pSPARC->lapcT[1]) > TEMP_TOL && fabs(pSPARC->lapcT[2]) > TEMP_TOL && fabs(pSPARC->lapcT[5]) < TEMP_TOL)
    pSPARC->cell_typ = 14;
  else if (fabs(pSPARC->lapcT[1]) < TEMP_TOL && fabs(pSPARC->lapcT[2]) > TEMP_TOL && fabs(pSPARC->lapcT[5]) > TEMP_TOL)
    pSPARC->cell_typ = 15;
  else if (fabs(pSPARC->lapcT[1]) > TEMP_TOL && fabs(pSPARC->lapcT[2]) < TEMP_TOL && fabs(pSPARC->lapcT[5]) > TEMP_TOL)
    pSPARC->cell_typ = 16;
  else if (fabs(pSPARC->lapcT[1]) > TEMP_TOL && fabs(pSPARC->lapcT[2]) > TEMP_TOL && fabs(pSPARC->lapcT[5]) > TEMP_TOL)
    pSPARC->cell_typ = 17;
#ifdef DEBUG
  if (!rank)
    printf("\n\nCELL_TYP: %d\n\n", pSPARC->cell_typ);
#endif
  /* transform the coefficiens of lapacian*/
  // int p, FDn = pSPARC->order/2;
  // double dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv;
  // dx_inv = 1.0 / (pSPARC->delta_x);
  // dy_inv = 1.0 / (pSPARC->delta_y);
  // dz_inv = 1.0 / (pSPARC->delta_z);
  // dx2_inv = 1.0 / (pSPARC->delta_x * pSPARC->delta_x);
  // dy2_inv = 1.0 / (pSPARC->delta_y * pSPARC->delta_y);
  // dz2_inv = 1.0 / (pSPARC->delta_z * pSPARC->delta_z);
  // for (p = 0; p < FDn + 1; p++) {
  //     pSPARC->D2_stencil_coeffs_x[p] = pSPARC->lapcT[0] *
  //     pSPARC->FDweights_D2[p] * dx2_inv; pSPARC->D2_stencil_coeffs_y[p] =
  //     pSPARC->lapcT[4] * pSPARC->FDweights_D2[p] * dy2_inv;
  //     pSPARC->D2_stencil_coeffs_z[p] = pSPARC->lapcT[8] *
  //     pSPARC->FDweights_D2[p] * dz2_inv; pSPARC->D2_stencil_coeffs_xy[p] = 2
  //     * pSPARC->lapcT[1] * pSPARC->FDweights_D1[p] * dx_inv; // 2*T_12
  //     d/dx(df/dy) pSPARC->D2_stencil_coeffs_xz[p] = 2 * pSPARC->lapcT[2] *
  //     pSPARC->FDweights_D1[p] * dx_inv; // 2*T_13 d/dx(df/dz)
  //     pSPARC->D2_stencil_coeffs_yz[p] = 2 * pSPARC->lapcT[5] *
  //     pSPARC->FDweights_D1[p] * dy_inv; // 2*T_23 d/dy(df/dz)
  //     pSPARC->D1_stencil_coeffs_xy[p] = 2 * pSPARC->lapcT[1] *
  //     pSPARC->FDweights_D1[p] * dy_inv; // d/dx(2*T_12 df/dy) used in
  //     d/dx(2*T_12 df/dy + 2*T_13 df/dz) pSPARC->D1_stencil_coeffs_yx[p] = 2 *
  //     pSPARC->lapcT[1] * pSPARC->FDweights_D1[p] * dx_inv; // d/dy(2*T_12
  //     df/dx) used in d/dy(2*T_12 df/dx + 2*T_23 df/dz)
  //     pSPARC->D1_stencil_coeffs_xz[p] = 2 * pSPARC->lapcT[2] *
  //     pSPARC->FDweights_D1[p] * dz_inv; // d/dx(2*T_13 df/dz) used in
  //     d/dx(2*T_12 df/dy + 2*T_13 df/dz) pSPARC->D1_stencil_coeffs_zx[p] = 2 *
  //     pSPARC->lapcT[2] * pSPARC->FDweights_D1[p] * dx_inv; // d/dz(2*T_13
  //     df/dx) used in d/dz(2*T_13 df/dz + 2*T_23 df/dy)
  //     pSPARC->D1_stencil_coeffs_yz[p] = 2 * pSPARC->lapcT[5] *
  //     pSPARC->FDweights_D1[p] * dz_inv; // d/dy(2*T_23 df/dz) used in
  //     d/dy(2*T_12 df/dx + 2*T_23 df/dz) pSPARC->D1_stencil_coeffs_zy[p] = 2 *
  //     pSPARC->lapcT[5] * pSPARC->FDweights_D1[p] * dy_inv; // d/dz(2*T_23
  //     df/dy) used in d/dz(2*T_12 df/dx + 2*T_23 df/dy)
  // }
  // TODO: Find maximum eigenvalue of Hamiltionian (= max. eigvalue of -0.5 lap)
  // for non orthogonal periodic systems
}
