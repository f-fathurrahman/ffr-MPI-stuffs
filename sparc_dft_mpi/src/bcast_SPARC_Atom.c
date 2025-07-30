
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "initialization.h"

/**
 * @brief Broadcast Atom info. in SPARC using MPI_Pack & MPI_Unpack.
 */
void bcast_SPARC_Atom(SPARC_OBJ *pSPARC) {
  int i, l, rank, position, l_buff, Ntypes, n_atom, nproj, lmax_sum, size_sum, nprojsize_sum, nproj_sum;
  int *tempbuff, *is_r_uniformv, *pspxcv, *pspsocv, *lmaxv, *sizev, *pplv, *pplv_soc, *ppl_sdispl, *ppl_soc_sdispl;
  ppl_soc_sdispl = NULL, pplv_soc = NULL;
  char *buff;

#ifdef DEBUG
  double t1, t2;
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Ntypes = pSPARC->Ntypes;
  /* preparation for broadcasting the structure */
  lmaxv = (int *)malloc(Ntypes * sizeof(int));
  sizev = (int *)malloc(Ntypes * sizeof(int));
  is_r_uniformv = (int *)malloc(Ntypes * sizeof(int));
  pspxcv = (int *)malloc(Ntypes * sizeof(int));
  pspsocv = (int *)malloc(Ntypes * sizeof(int));
  ppl_sdispl = (int *)malloc((Ntypes + 1) * sizeof(int));
  tempbuff = (int *)malloc((5 * Ntypes + 3) * sizeof(int));
  assert(lmaxv != NULL && sizev != NULL && is_r_uniformv != NULL && pspxcv != NULL && pspsocv != NULL &&
         ppl_sdispl != NULL && tempbuff != NULL);

  // send n_atom, lmax, size
  if (rank == 0) {
    // pack the info. into temp. buffer
    tempbuff[0] = pSPARC->n_atom;
    tempbuff[1] = pSPARC->Nspinor;
    tempbuff[2] = pSPARC->SOC_Flag;
    for (i = 0; i < Ntypes; i++) {
      lmaxv[i] = pSPARC->psd[i].lmax;
      sizev[i] = pSPARC->psd[i].size;
      pspxcv[i] = pSPARC->psd[i].pspxc;
      is_r_uniformv[i] = pSPARC->psd[i].is_r_uniform;
      pspsocv[i] = pSPARC->psd[i].pspsoc;
      // pplv[i] = pSPARC->psd[i].ppl;
      tempbuff[i + 3] = lmaxv[i];
      tempbuff[i + Ntypes + 3] = sizev[i];
      tempbuff[i + 2 * Ntypes + 3] = pspxcv[i];
      tempbuff[i + 3 * Ntypes + 3] = is_r_uniformv[i];
      tempbuff[i + 4 * Ntypes + 3] = pspsocv[i];

      // tempbuff[i+2*Ntypes+1] = pplv[i];
    }
    MPI_Bcast(tempbuff, 5 * Ntypes + 3, MPI_INT, 0, MPI_COMM_WORLD);

    // pack psd[i].ppl[l] and bcast
    ppl_sdispl[0] = 0;
    for (i = 0; i < Ntypes; i++) {
      ppl_sdispl[i + 1] = ppl_sdispl[i] + lmaxv[i] + 1;
    }
    pplv = (int *)malloc(ppl_sdispl[Ntypes] * sizeof(int));
    for (i = 0; i < Ntypes; i++) {
      for (l = 0; l <= lmaxv[i]; l++) {
        pplv[ppl_sdispl[i] + l] = pSPARC->psd[i].ppl[l];
      }
    }
    MPI_Bcast(pplv, ppl_sdispl[Ntypes], MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    // allocate psd array for receiver process
    pSPARC->psd = (PSD_OBJ *)malloc(pSPARC->Ntypes * sizeof(PSD_OBJ));
    assert(pSPARC->psd != NULL);
#ifdef DEBUG
    t1 = MPI_Wtime();
#endif
    MPI_Bcast(tempbuff, 5 * Ntypes + 3, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    t2 = MPI_Wtime();
    if (rank == 0)
      printf("Bcast pre-info. took %.3f ms\n", (t2 - t1) * 1000);
#endif
    // unpack info.
    pSPARC->n_atom = tempbuff[0];
    pSPARC->Nspinor = tempbuff[1];
    pSPARC->SOC_Flag = tempbuff[2];
    for (i = 0; i < Ntypes; i++) {
      lmaxv[i] = tempbuff[i + 3];
      sizev[i] = tempbuff[i + Ntypes + 3];
      pspxcv[i] = tempbuff[i + 2 * Ntypes + 3];
      is_r_uniformv[i] = tempbuff[i + 3 * Ntypes + 3];
      pspsocv[i] = tempbuff[i + 4 * Ntypes + 3];
      // pplv[i] = tempbuff[i+2*Ntypes+1];
      pSPARC->psd[i].lmax = lmaxv[i];
      pSPARC->psd[i].size = sizev[i];
      pSPARC->psd[i].pspxc = pspxcv[i];
      pSPARC->psd[i].is_r_uniform = is_r_uniformv[i];
      pSPARC->psd[i].pspsoc = pspsocv[i];
      // pSPARC->psd[i].ppl = pplv[i];
    }

    // bcast psd[i].ppl[l] and unpack
    ppl_sdispl[0] = 0;
    for (i = 0; i < Ntypes; i++) {
      ppl_sdispl[i + 1] = ppl_sdispl[i] + lmaxv[i] + 1;
    }

    pplv = (int *)malloc(ppl_sdispl[Ntypes] * sizeof(int));
    MPI_Bcast(pplv, ppl_sdispl[Ntypes], MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < Ntypes; i++) {
      pSPARC->psd[i].ppl = (int *)malloc((lmaxv[i] + 1) * sizeof(int));
      for (l = 0; l <= lmaxv[i]; l++) {
        pSPARC->psd[i].ppl[l] = pplv[ppl_sdispl[i] + l];
      }
    }
  }
  n_atom = pSPARC->n_atom;

  if (pSPARC->SOC_Flag == 1) {
    ppl_soc_sdispl = (int *)malloc((Ntypes + 1) * sizeof(int));
    // pack psd[i].ppl[l] and bcast
    ppl_soc_sdispl[0] = 0;
    for (i = 0; i < Ntypes; i++) {
      ppl_soc_sdispl[i + 1] = ppl_soc_sdispl[i] + pspsocv[i] * lmaxv[i]; // l from 1 to lmax
    }
    pplv_soc = (int *)malloc(ppl_soc_sdispl[Ntypes] * sizeof(int));

    if (rank == 0) {
      for (i = 0; i < Ntypes; i++) {
        if (!pspsocv[i])
          continue;
        for (l = 1; l <= lmaxv[i]; l++) {
          pplv_soc[ppl_soc_sdispl[i] + l - 1] = pSPARC->psd[i].ppl_soc[l - 1];
        }
      }
      MPI_Bcast(pplv_soc, ppl_soc_sdispl[Ntypes], MPI_INT, 0, MPI_COMM_WORLD);
    } else {
      MPI_Bcast(pplv_soc, ppl_soc_sdispl[Ntypes], MPI_INT, 0, MPI_COMM_WORLD);
      for (i = 0; i < Ntypes; i++) {
        pSPARC->psd[i].ppl_soc = (int *)malloc(lmaxv[i] * sizeof(int));
        for (l = 1; l <= lmaxv[i]; l++) {
          pSPARC->psd[i].ppl_soc[l - 1] = pplv_soc[ppl_soc_sdispl[i] + l - 1];
        }
      }
    }
  }

  // allocate memory for buff, the extra 16*(Ntypes+3*n_atom) byte is spare
  // memory
  lmax_sum = 0;
  size_sum = 0;
  nproj_sum = 0;
  nprojsize_sum = 0;
  for (i = 0; i < Ntypes; i++) {
    lmax_sum += lmaxv[i] + 1;
    size_sum += sizev[i];
    // sizelmax_sum += ((lmaxv[i]+1) * sizev[i]);
    // ppllmax_sum += pplv[i] * (lmaxv[i]+1);
    nproj = 0;
    for (l = 0; l <= lmaxv[i]; l++) {
      nproj += pSPARC->psd[i].ppl[l];
    }
    nproj_sum += nproj;
    nprojsize_sum += nproj * sizev[i];
  }

  int nprojso, nprojso_sum, nprojsosize_sum;
  nprojso_sum = nprojsosize_sum = 0;
  if (pSPARC->SOC_Flag == 1) {
    for (i = 0; i < Ntypes; i++) {
      if (!pSPARC->psd[i].pspsoc)
        continue;
      nprojso = 0;
      for (l = 1; l <= lmaxv[i]; l++) {
        nprojso += pSPARC->psd[i].ppl_soc[l - 1];
      }
      nprojso_sum += nprojso;
      nprojsosize_sum += nprojso * sizev[i];
    }
  }

  l_buff = (3 * Ntypes + 3 * n_atom + 4 * size_sum + lmax_sum + nproj_sum + nprojsize_sum + nprojsosize_sum +
            nprojsosize_sum + n_atom) *
               sizeof(double) +
           (6 * Ntypes + 3 * n_atom) * sizeof(int) + Ntypes * (L_PSD + L_ATMTYPE) * sizeof(char) +
           0 * (Ntypes + 3 * n_atom) * 16; // last term is spare memory in case
  buff = (char *)malloc(l_buff * sizeof(char));
  assert(buff != NULL);

  if (rank == 0) {
    // pack the variables
    position = 0;
    MPI_Pack(pSPARC->localPsd, Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->Zatom, Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->Znucl, Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->nAtomv, Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->Mass, Ntypes, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->atomType, Ntypes * L_ATMTYPE, MPI_CHAR, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->psdName, Ntypes * L_PSD, MPI_CHAR, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->atom_pos, 3 * n_atom, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->IsFrac, pSPARC->Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->mvAtmConstraint, 3 * n_atom, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->IsSpin, pSPARC->Ntypes, MPI_INT, buff, l_buff, &position, MPI_COMM_WORLD);
    MPI_Pack(pSPARC->atom_spin, n_atom, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
    for (i = 0; i < Ntypes; i++) {
      nproj = 0;
      for (l = 0; l <= lmaxv[i]; l++) {
        nproj += pSPARC->psd[i].ppl[l];
      }
      MPI_Pack(pSPARC->psd[i].RadialGrid, sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].UdV, nproj * sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(&pSPARC->psd[i].Vloc_0, 1, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(&pSPARC->psd[i].fchrg, 1, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].rVloc, sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].rhoIsoAtom, sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].rc, lmaxv[i] + 1, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].Gamma, nproj, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      MPI_Pack(pSPARC->psd[i].rho_c_table, sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);

      if (pSPARC->psd[i].pspsoc) {
        nprojso = 0;
        for (l = 1; l <= lmaxv[i]; l++) {
          nprojso += pSPARC->psd[i].ppl_soc[l - 1];
        }
        MPI_Pack(pSPARC->psd[i].Gamma_soc, nprojso, MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
        MPI_Pack(pSPARC->psd[i].UdV_soc, nprojso * sizev[i], MPI_DOUBLE, buff, l_buff, &position, MPI_COMM_WORLD);
      }
    }
    // broadcast the packed buffer
    MPI_Bcast(buff, l_buff, MPI_PACKED, 0, MPI_COMM_WORLD);
  } else {
    /* allocate memory for receiver processes */
    pSPARC->localPsd = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->Zatom = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->Znucl = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->nAtomv = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->IsFrac = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->IsSpin = (int *)malloc(Ntypes * sizeof(int));
    pSPARC->Mass = (double *)malloc(Ntypes * sizeof(double));
    pSPARC->atomType = (char *)malloc(Ntypes * L_ATMTYPE * sizeof(char));
    pSPARC->psdName = (char *)malloc(Ntypes * L_PSD * sizeof(char));

    if (pSPARC->localPsd == NULL || pSPARC->Zatom == NULL || pSPARC->Znucl == NULL || pSPARC->nAtomv == NULL ||
        pSPARC->Mass == NULL || pSPARC->atomType == NULL || pSPARC->psdName == NULL) {
      printf("\nmemory cannot be allocated3\n");
      exit(EXIT_FAILURE);
    }

    // values to receive
    pSPARC->atom_pos = (double *)malloc(3 * n_atom * sizeof(double));
    pSPARC->mvAtmConstraint = (int *)malloc(3 * n_atom * sizeof(int));
    pSPARC->atom_spin = (double *)malloc(n_atom * sizeof(double));
    if (pSPARC->atom_pos == NULL || pSPARC->mvAtmConstraint == NULL || pSPARC->atom_spin == NULL) {
      printf("\nmemory cannot be allocated4\n");
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < Ntypes; i++) {
      nproj = 0;
      for (l = 0; l <= lmaxv[i]; l++) {
        nproj += pSPARC->psd[i].ppl[l];
      }
      pSPARC->psd[i].RadialGrid = (double *)malloc(sizev[i] * sizeof(double));
      pSPARC->psd[i].UdV = (double *)malloc(nproj * sizev[i] * sizeof(double));
      pSPARC->psd[i].rVloc = (double *)malloc(sizev[i] * sizeof(double));
      pSPARC->psd[i].rhoIsoAtom = (double *)malloc(sizev[i] * sizeof(double));
      pSPARC->psd[i].rc = (double *)malloc((lmaxv[i] + 1) * sizeof(double));
      pSPARC->psd[i].Gamma = (double *)malloc(nproj * sizeof(double));
      pSPARC->psd[i].rho_c_table = (double *)malloc(sizev[i] * sizeof(double));
      // check if memory is allocated successfully!
      if (pSPARC->psd[i].RadialGrid == NULL || pSPARC->psd[i].UdV == NULL || pSPARC->psd[i].rVloc == NULL ||
          pSPARC->psd[i].rhoIsoAtom == NULL || pSPARC->psd[i].rc == NULL || pSPARC->psd[i].Gamma == NULL ||
          pSPARC->psd[i].rho_c_table == NULL) {
        printf("\nmemory cannot be allocated5\n");
        exit(EXIT_FAILURE);
      }
      if (pSPARC->psd[i].pspsoc) {
        nprojso = 0;
        for (l = 1; l <= lmaxv[i]; l++) {
          nprojso += pSPARC->psd[i].ppl_soc[l - 1];
        }
        pSPARC->psd[i].Gamma_soc = (double *)malloc(nprojso * sizeof(double));
        pSPARC->psd[i].UdV_soc = (double *)malloc(nprojso * sizev[i] * sizeof(double));
      }
    }
#ifdef DEBUG
    t1 = MPI_Wtime();
#endif
    // broadcast the packed buffer
    MPI_Bcast(buff, l_buff, MPI_PACKED, 0, MPI_COMM_WORLD);
#ifdef DEBUG
    t2 = MPI_Wtime();
    if (rank == 0)
      printf("MPI_Bcast packed buff of length %d took %.3f ms\n", l_buff, (t2 - t1) * 1000);
#endif
    // unpack the variables
    position = 0;
    MPI_Unpack(buff, l_buff, &position, pSPARC->localPsd, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->Zatom, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->Znucl, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->nAtomv, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->Mass, Ntypes, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->atomType, Ntypes * L_ATMTYPE, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->psdName, Ntypes * L_PSD, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->atom_pos, 3 * n_atom, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->IsFrac, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->mvAtmConstraint, 3 * n_atom, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->IsSpin, Ntypes, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, l_buff, &position, pSPARC->atom_spin, n_atom, MPI_DOUBLE, MPI_COMM_WORLD);
    for (i = 0; i < Ntypes; i++) {
      nproj = 0;
      for (l = 0; l <= lmaxv[i]; l++) {
        nproj += pSPARC->psd[i].ppl[l];
      }
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].RadialGrid, sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].UdV, nproj * sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, &pSPARC->psd[i].Vloc_0, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, &pSPARC->psd[i].fchrg, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].rVloc, sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].rhoIsoAtom, sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].rc, lmaxv[i] + 1, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].Gamma, nproj, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].rho_c_table, sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      if (pSPARC->psd[i].pspsoc) {
        nprojso = 0;
        for (l = 1; l <= lmaxv[i]; l++) {
          nprojso += pSPARC->psd[i].ppl_soc[l - 1];
        }
        MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].Gamma_soc, nprojso, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(buff, l_buff, &position, pSPARC->psd[i].UdV_soc, nprojso * sizev[i], MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  // deallocate memory
  free(tempbuff);
  free(ppl_sdispl);
  free(is_r_uniformv);
  free(pspxcv);
  free(pspsocv);
  free(pplv);
  free(lmaxv);
  free(sizev);
  free(buff);
  if (pSPARC->SOC_Flag == 1) {
    free(ppl_soc_sdispl);
    free(pplv_soc);
  }
}

