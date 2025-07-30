
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "initialization.h"

/**
 * @brief   Estimate the memory required for the simulation.
 */
double estimate_memory(const SPARC_OBJ *pSPARC) {
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (pSPARC->SQFlag == 1) {
    SQ_OBJ *pSQ = pSPARC->pSQ;
    // TODO: add accurate memory estimation
    double mem_PR, mem_Rcut, mem_phi, mem_chi;
    mem_PR = (double)sizeof(double) * pSQ->DMnd_PR * nproc;
    mem_Rcut = 0;
    mem_phi = (double)sizeof(double) * pSPARC->Nd_d * (6 + pSPARC->SQ_npl_g * 2 + 1) * nproc;
    mem_chi = (double)sizeof(double) * pSPARC->Nd_d * pSQ->Nd_loc * 0.3 * nproc;
    if (pSPARC->SQ_gauss_mem == 1) { // save vectors for all FD nodes
      mem_Rcut += (double)sizeof(double) * pSPARC->Nd_d * pSQ->Nd_loc * pSPARC->SQ_npl_g * nproc;
      mem_phi += (double)sizeof(double) * pSPARC->Nd_d * pSPARC->SQ_npl_g * pSPARC->SQ_npl_g * nproc;
    } else {
      mem_Rcut += (double)sizeof(double) * pSQ->Nd_loc * pSPARC->SQ_npl_g * nproc;
    }

    double memory_usage = 0.0;
    memory_usage = mem_PR + mem_Rcut + mem_phi + mem_chi;

#ifdef DEBUG
    if (rank == 0) {
      char mem_str[32];
      printf("----------------------\n");
      printf("Estimated memory usage\n");
      formatBytes(memory_usage, 32, mem_str);
      printf("Total: %s\n", mem_str);
      formatBytes(mem_PR, 32, mem_str);
      printf("Vectors in P.R. domain : %s\n", mem_str);
      formatBytes(mem_Rcut, 32, mem_str);
      printf("Vectors in Rcut domain : %s\n", mem_str);
      formatBytes(mem_phi, 32, mem_str);
      printf("Vectors in phi domain : %s\n", mem_str);
      formatBytes(mem_chi, 32, mem_str);
      printf("All saved nonlocal projectors: %s\n", mem_str);
      printf("----------------------------------------------\n");
      formatBytes(memory_usage / nproc, 32, mem_str);
      printf("Estimated memory usage per processor: %s\n", mem_str);
    }
#endif
    return memory_usage;
  }
  int Nd = pSPARC->Nd;
  int Ns = pSPARC->Nstates;
  int Nspin = pSPARC->Nspin;
  int Nkpts_sym = pSPARC->Nkpts_sym;
  int m = pSPARC->MixingHistory;
  int npspin = pSPARC->npspin;
  int npkpt = pSPARC->npkpt;
  int npNd = pSPARC->npNd;

  int type_size;
  if (pSPARC->isGammaPoint) {
    type_size = sizeof(double);
  } else {
    type_size = sizeof(double _Complex);
  }

  // orbitals (dominant)
  double ncpy_orbitals; // extra copies required during CheFSI
  if (pSPARC->npband > 1) {
    // MKL pdgemr2d internally creates ~2.5 copies during pdgemr2d in projection
    // + Yorb, Yorb_BLCYC, HY_BLCYC
    ncpy_orbitals = 5.5;
  } else {
    // when there's no band parallelization (domain only), then pdgemr2d is not
    // needed for projection moreover, the block cyclic formats Yorb_BLCYC,
    // HY_BLCYC (also YQ_BLCYC during rotation) are not needed so the necessary
    // copies are: Yorb, Ynew (during Chebyshev filtering) sometimes dgemm (MKL)
    // also might internally create ~0.5 copy of the orbital, we add 0.5 for
    // safety
    ncpy_orbitals = 2.5;
  }
#ifdef USE_DP_SUBEIG
  ncpy_orbitals = 6; // DP requires 4 extra copies, Yorb, and a temp copy
                     // (Chebyshev filtering and Projection)
#endif
  double memory_orbitals = (double)Nd * Ns * (ncpy_orbitals * npspin * npkpt + Nspin * Nkpts_sym) * type_size;

  // vectors: rho, phi, Veff, mixing history vectors, etc.
  int ncpy_vectors = 6 + 4 * Nspin + 2 * m * Nspin + 3 * (2 * Nspin - 1) + 1;
  double memory_vectors = (double)ncpy_vectors * Nd * sizeof(double);

  // subspace matrix: Hs, Ms, Q
  int ncpy_matrices = 3 * npNd;
#ifdef USE_DP_SUBEIG
  ncpy_matrices = 3 * nproc; // DP stores Hp_local,Mp_local,eigvecs in (almost) every process
#endif
  double memory_matrices = (double)ncpy_matrices * Ns * Ns * sizeof(double);

  // total memory
  double memory_usage = memory_orbitals + memory_vectors + memory_matrices;

  // add some buffer for other memory
  const double buf_rat = 0.10; // add 10% more memory
  memory_usage *= (1.0 + buf_rat);

  // memory for Exact Exchange part
  double memory_exx = 0.0;
  if (pSPARC->usefock > 0) {
    memory_exx = estimate_memory_exx(pSPARC);
    memory_usage += memory_exx;
  }

#ifdef DEBUG
  if (rank == 0) {
    char mem_str[32];
    printf("----------------------\n");
    printf("Estimated memory usage\n");
    formatBytes(memory_usage, 32, mem_str);
    printf("Total: %s\n", mem_str);
    formatBytes(memory_orbitals, 32, mem_str);
    printf("orbitals             : %s\n", mem_str);
    formatBytes(memory_vectors, 32, mem_str);
    printf("global sized vectors : %s\n", mem_str);
    formatBytes(memory_matrices, 32, mem_str);
    printf("subspace matrices : %s\n", mem_str);
    if (pSPARC->usefock > 0) {
      formatBytes(memory_exx, 32, mem_str);
      printf("global exact exchange memory : %s\n", mem_str);
    }
    formatBytes(memory_usage * buf_rat / (1.0 + buf_rat), 32, mem_str);
    printf("others : %s\n", mem_str);
    printf("----------------------------------------------\n");
    formatBytes(memory_usage / nproc, 32, mem_str);
    printf("Estimated memory usage per processor: %s\n", mem_str);
  }
#endif
  return memory_usage;
}
