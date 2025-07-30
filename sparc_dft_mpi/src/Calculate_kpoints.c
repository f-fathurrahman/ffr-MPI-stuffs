#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "initialization.h"

#define TEMP_TOL 1e-12

/**
 * @brief   Calculate k-points and the associated weights.
 */
void Calculate_kpoints(SPARC_OBJ *pSPARC) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Initialize the weights of k points
  // by going over all the k points from the Monkhorst Pack grid
  int k, nk1, nk2, nk3, k_hf, k_hf_rd;
  double k1, k2, k3;
  double Lx = pSPARC->range_x;
  double Ly = pSPARC->range_y;
  double Lz = pSPARC->range_z;
  k = k_hf = k_hf_rd = 0;

  if (pSPARC->usefock == 1) {
    // constrains on EXX_DOWNSAMPLING
    if (pSPARC->EXXDownsampling[0] < 0 || pSPARC->EXXDownsampling[1] < 0 || pSPARC->EXXDownsampling[2] < 0) {
      if (rank == 0)
        printf(RED "ERROR: EXX_DOWNSAMPLING must be non-negative.\n" RESET);
      exit(EXIT_FAILURE);
    }

    if ((pSPARC->EXXDownsampling[0] > 0 && pSPARC->Kx % pSPARC->EXXDownsampling[0]) ||
        (pSPARC->EXXDownsampling[1] > 0 && pSPARC->Ky % pSPARC->EXXDownsampling[1]) ||
        (pSPARC->EXXDownsampling[2] > 0 && pSPARC->Kz % pSPARC->EXXDownsampling[2])) {
      if (rank == 0)
        printf(RED "ERROR: Number of kpoints must be divisible by EXX_DOWNSAMPLING "
                   "in all directions if EXX_DOWNSAMPLING is positive.\n" RESET);
      exit(EXIT_FAILURE);
    }

    if (pSPARC->EXXDownsampling[0] == 0) {
      pSPARC->Nkpts_hf = 1;
      pSPARC->Kx_hf = 1;
    } else {
      pSPARC->Kx_hf = pSPARC->Kx / pSPARC->EXXDownsampling[0];
      pSPARC->Nkpts_hf = pSPARC->Kx / pSPARC->EXXDownsampling[0];
    }

    if (pSPARC->EXXDownsampling[1] == 0) {
      pSPARC->Nkpts_hf *= 1;
      pSPARC->Ky_hf = 1;
    } else {
      pSPARC->Ky_hf = pSPARC->Ky / pSPARC->EXXDownsampling[1];
      pSPARC->Nkpts_hf *= (pSPARC->Ky / pSPARC->EXXDownsampling[1]);
    }

    if (pSPARC->EXXDownsampling[2] == 0) {
      pSPARC->Nkpts_hf *= 1;
      pSPARC->Kz_hf = 1;
    } else {
      pSPARC->Kz_hf = pSPARC->Kz / pSPARC->EXXDownsampling[2];
      pSPARC->Nkpts_hf *= (pSPARC->Kz / pSPARC->EXXDownsampling[2]);
    }
    pSPARC->k1_hf = (double *)calloc(sizeof(double), pSPARC->Nkpts_hf);
    pSPARC->k2_hf = (double *)calloc(sizeof(double), pSPARC->Nkpts_hf);
    pSPARC->k3_hf = (double *)calloc(sizeof(double), pSPARC->Nkpts_hf);
    pSPARC->kpthf_ind = (int *)calloc(sizeof(int), pSPARC->Nkpts_hf);
    pSPARC->kpthf_ind_red = (int *)calloc(sizeof(int), pSPARC->Nkpts_hf);
    pSPARC->kpthf_pn = (int *)calloc(sizeof(int), pSPARC->Nkpts_hf);
    pSPARC->kptWts_hf = 1.0 / pSPARC->Nkpts_hf;
  }

  double sumx = 2.0 * M_PI / Lx;
  double sumy = 2.0 * M_PI / Ly;
  double sumz = 2.0 * M_PI / Lz;

  int nk, flag;
  int flag_0x, flag_0y, flag_0z; // flag for finding 0 k-point in 3 directions
                                 // when exx_downsampling is 0
  int flag_cx, flag_cy, flag_cz; // flag for correct slot in x,y,z direction
  flag_0x = flag_0y = flag_0z = 0;

  // calculate M-P grid similar to that in ABINIT
  int nk1_s = -floor((pSPARC->Kx - 1) / 2);
  int nk1_e = nk1_s + pSPARC->Kx;
  int nk2_s = -floor((pSPARC->Ky - 1) / 2);
  int nk2_e = nk2_s + pSPARC->Ky;
  int nk3_s = -floor((pSPARC->Kz - 1) / 2);
  int nk3_e = nk3_s + pSPARC->Kz;

  for (nk1 = nk1_s; nk1 < nk1_e; nk1++) {
    for (nk2 = nk2_s; nk2 < nk2_e; nk2++) {
      for (nk3 = nk3_s; nk3 < nk3_e; nk3++) {
        double k1_red, k2_red, k3_red;
        // calculate Monkhorst-Pack k points (reduced) using Monkhorst pack grid
        k1_red = nk1 * 1.0 / pSPARC->Kx;
        k2_red = nk2 * 1.0 / pSPARC->Ky;
        k3_red = nk3 * 1.0 / pSPARC->Kz;
        k1_red = fmod(k1_red + pSPARC->kptshift[0] / pSPARC->Kx + 0.5 - TEMP_TOL, 1.0) - 0.5 + TEMP_TOL;
        k2_red = fmod(k2_red + pSPARC->kptshift[1] / pSPARC->Ky + 0.5 - TEMP_TOL, 1.0) - 0.5 + TEMP_TOL;
        k3_red = fmod(k3_red + pSPARC->kptshift[2] / pSPARC->Kz + 0.5 - TEMP_TOL, 1.0) - 0.5 + TEMP_TOL;
#ifdef DEBUG
        if (!rank)
          printf("[k1_red,k2_red,k3_red] = %8.4f %8.4f %8.4f\n", k1_red, k2_red, k3_red);
#endif
        k1 = k1_red * 2.0 * M_PI / Lx;
        k2 = k2_red * 2.0 * M_PI / Ly;
        k3 = k3_red * 2.0 * M_PI / Lz;

        flag = 1;
        for (nk = 0; nk < k; nk++) {
          if ((fabs(k1 + pSPARC->k1[nk]) < TEMP_TOL || fabs(k1 + pSPARC->k1[nk] - sumx) < TEMP_TOL) &&
              (fabs(k2 + pSPARC->k2[nk]) < TEMP_TOL || fabs(k2 + pSPARC->k2[nk] - sumy) < TEMP_TOL) &&
              (fabs(k3 + pSPARC->k3[nk]) < TEMP_TOL || fabs(k3 + pSPARC->k3[nk] - sumz) < TEMP_TOL)) {
            flag = 0;
            break;
          }
        }

        if (flag) {
          pSPARC->k1[k] = k1;
          pSPARC->k2[k] = k2;
          pSPARC->k3[k] = k3;
          pSPARC->kptWts[k] = 1.0;
          k++;
        } else {
          pSPARC->kptWts[nk] = 2.0;
        }

        if (pSPARC->usefock == 1) {
          if (pSPARC->EXXDownsampling[0] == 0) {
            flag_cx = (fabs(k1) < TEMP_TOL);
            if (flag_cx)
              flag_0x = 1;
          } else {
            flag_cx = !((nk1 - nk1_s + 1) % pSPARC->EXXDownsampling[0]);
          }

          if (pSPARC->EXXDownsampling[1] == 0) {
            flag_cy = (fabs(k2) < TEMP_TOL);
            if (flag_cy)
              flag_0y = 1;
          } else {
            flag_cy = !((nk2 - nk2_s + 1) % pSPARC->EXXDownsampling[1]);
          }

          if (pSPARC->EXXDownsampling[2] == 0) {
            flag_cz = (fabs(k3) < TEMP_TOL);
            if (flag_cz)
              flag_0z = 1;
          } else {
            flag_cz = !((nk3 - nk3_s + 1) % pSPARC->EXXDownsampling[2]);
          }

          if (flag_cx && flag_cy && flag_cz) {
            pSPARC->k1_hf[k_hf] = k1;
            pSPARC->k2_hf[k_hf] = k2;
            pSPARC->k3_hf[k_hf] = k3;
            k_hf++;
          }
        }
      }
    }
  }

  pSPARC->Nkpts_sym = k; // update number of k points after symmetry reduction

  if (pSPARC->usefock == 1) {
    if (!pSPARC->EXXDownsampling[0] && !flag_0x) {
      if (rank == 0)
        printf(RED "ERROR: Gamma point is not one of the k-vectors. Please use "
                   "positive EXX_DOWNSAMPLING or change k-point grid in the "
                   "first direction.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (!pSPARC->EXXDownsampling[1] && !flag_0y) {
      if (rank == 0)
        printf(RED "ERROR: Gamma point is not one of the k-vectors. Please use "
                   "positive EXX_DOWNSAMPLING or change k-point grid in the "
                   "second direction.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (!pSPARC->EXXDownsampling[2] && !flag_0z) {
      if (rank == 0)
        printf(RED "ERROR: Gamma point is not one of the k-vectors. Please use "
                   "positive EXX_DOWNSAMPLING or change k-point grid in the "
                   "third direction.\n" RESET);
      exit(EXIT_FAILURE);
    }

    for (k_hf = 0; k_hf < pSPARC->Nkpts_hf; k_hf++) {
      for (k = 0; k < pSPARC->Nkpts_sym; k++) {
        if ((fabs(pSPARC->k1[k] - pSPARC->k1_hf[k_hf]) < TEMP_TOL) &&
            (fabs(pSPARC->k2[k] - pSPARC->k2_hf[k_hf]) < TEMP_TOL) &&
            (fabs(pSPARC->k3[k] - pSPARC->k3_hf[k_hf]) < TEMP_TOL)) {
          pSPARC->kpthf_ind[k_hf] = k; // index w.r.t. Nkpts_sym
          pSPARC->kpthf_pn[k_hf] = 1;  // 1 -> k, 0 -> -k
          break;
        }
      }
      if (pSPARC->kpthf_pn[k_hf] == 1)
        continue;
      for (k = 0; k < pSPARC->Nkpts_sym; k++) {
        if (((fabs(pSPARC->k1[k] + pSPARC->k1_hf[k_hf]) < TEMP_TOL) ||
             (fabs(pSPARC->k1[k] + pSPARC->k1_hf[k_hf] - sumx) < TEMP_TOL)) &&
            ((fabs(pSPARC->k2[k] + pSPARC->k2_hf[k_hf]) < TEMP_TOL) ||
             (fabs(pSPARC->k2[k] + pSPARC->k2_hf[k_hf] - sumy) < TEMP_TOL)) &&
            ((fabs(pSPARC->k3[k] + pSPARC->k3_hf[k_hf]) < TEMP_TOL) ||
             (fabs(pSPARC->k3[k] + pSPARC->k3_hf[k_hf] - sumz) < TEMP_TOL))) {
          pSPARC->kpthf_ind[k_hf] = k; // index w.r.t. Nkpts_sym
          pSPARC->kpthf_pn[k_hf] = 0;  // 1 -> k, 0 -> -k
        }
      }
    }

    // find list of k-points after reduce for HF
    pSPARC->Nkpts_hf_red = 1;
    pSPARC->kpts_hf_red_list = (int *)calloc((unsigned int)pSPARC->Nkpts_sym, sizeof(int));
    pSPARC->kpts_hf_red_list[0] = pSPARC->kpthf_ind[0];
    for (k_hf = 1; k_hf < pSPARC->Nkpts_hf; k_hf++) {
      flag = 0;
      for (k = 0; k < k_hf; k++) {
        if (pSPARC->kpthf_ind[k] == pSPARC->kpthf_ind[k_hf]) {
          flag = 1;
          break;
        }
      }
      if (flag)
        continue;
      pSPARC->kpts_hf_red_list[pSPARC->Nkpts_hf_red++] = pSPARC->kpthf_ind[k_hf];
    }

    pSPARC->kpthfred2kpthf = (int(*)[3])calloc(sizeof(int[3]), pSPARC->Nkpts_hf_red);
    for (k = 0; k < pSPARC->Nkpts_hf_red; k++)
      pSPARC->kpthfred2kpthf[k][0] = 0;
    // update the reduced index to be w.r.t. Nkpts_hf_red
    // find inverse mapping from Nkpts_hf_red to Nkpts_hf
    for (k_hf = 0; k_hf < pSPARC->Nkpts_hf; k_hf++) {
      for (k = 0; k < pSPARC->Nkpts_hf_red; k++) {
        if (pSPARC->kpthf_ind[k_hf] == pSPARC->kpts_hf_red_list[k]) {
          pSPARC->kpthf_ind_red[k_hf] = k; // index w.r.t. Nkpts_hf_red
          pSPARC->kpthfred2kpthf[k][0]++;
          int indx = pSPARC->kpthfred2kpthf[k][0];
          pSPARC->kpthfred2kpthf[k][indx] = k_hf;
          break;
        }
      }
    }
  }

#ifdef DEBUG
  if (!rank)
    printf("After symmetry reduction, Nkpts_sym = %d\n", pSPARC->Nkpts_sym);
  for (nk = 0; nk < pSPARC->Nkpts_sym; nk++) {
    double tpiblx = 2 * M_PI / Lx;
    double tpibly = 2 * M_PI / Ly;
    double tpiblz = 2 * M_PI / Lz;
    if (!rank)
      printf("k1[%2d]: %8.4f, k2[%2d]: %8.4f, k3[%2d]: %8.4f, kptwt[%2d]: %.3f \n", nk, pSPARC->k1[nk] / tpiblx, nk,
             pSPARC->k2[nk] / tpibly, nk, pSPARC->k3[nk] / tpiblz, nk, pSPARC->kptWts[nk]);
  }

  if (pSPARC->usefock == 1) {
    if (!rank)
      printf("K-points for Hartree-Fock operator after downsampling, Nkpts_hf "
             "%d\n",
             pSPARC->Nkpts_hf);
    for (nk = 0; nk < pSPARC->Nkpts_hf; nk++) {
      double tpiblx = 2 * M_PI / Lx;
      double tpibly = 2 * M_PI / Ly;
      double tpiblz = 2 * M_PI / Lz;
      if (!rank)
        printf("k1_hf[%2d]: %8.4f, k2_hf[%2d]: %8.4f, k3_hf[%2d]: %8.4f, "
               "kptwt[%2d]: %.3f, kpthf_ind[%2d]: %2d, kpthf_ind_red[%2d]: "
               "%2d, kpthf_pn[%2d]: %2d \n",
               nk, pSPARC->k1_hf[nk] / tpiblx, nk, pSPARC->k2_hf[nk] / tpibly, nk, pSPARC->k3_hf[nk] / tpiblz, nk,
               pSPARC->kptWts_hf, nk, pSPARC->kpthf_ind[nk], nk, pSPARC->kpthf_ind_red[nk], nk, pSPARC->kpthf_pn[nk]);
    }
    if (!rank)
      printf("K-points for Hartree-Fock operator after downsampling mapping "
             "into kpts_system, Nkpts_hf_red %d\n",
             pSPARC->Nkpts_hf_red);
    for (nk = 0; nk < pSPARC->Nkpts_hf_red; nk++) {
      if (!rank)
        printf("kpts_hf_red_list[%d]: %d\n", nk, pSPARC->kpts_hf_red_list[nk]);
    }
  }
#endif
}

