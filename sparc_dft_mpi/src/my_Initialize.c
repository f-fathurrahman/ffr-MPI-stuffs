#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "isddft.h"

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

void my_Initialize(SPARC_OBJ *pSPARC, int argc, char *argv[]) {

#ifdef DEBUG
  double t1, t2;
#endif
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Request req;

  // these two structs are for reading info. and broadcasting
  SPARC_INPUT_OBJ SPARC_Input;

#ifdef DEBUG
  t1 = MPI_Wtime();
#endif

  /* Create new MPI struct datatype SPARC_INPUT_MPI (for broadcasting) */
  MPI_Datatype SPARC_INPUT_MPI;
  SPARC_Input_MPI_create(&SPARC_INPUT_MPI);

#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("\nCreating SPARC_INPUT_MPI datatype took %.3f ms\n", (t2 - t1) * 1000);
#endif

  if (rank == 0) {
#ifdef DEBUG
    printf("Initializing ...\n");
    t1 = MPI_Wtime();
#endif
    // check input arguments and read filename
    check_inputs(&SPARC_Input, argc, argv);

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nChecking inputs parsed by commandline took %.3f ms\n", (t2 - t1) * 1000);
    t1 = MPI_Wtime();
#endif
    set_defaults(&SPARC_Input, pSPARC); // set default values

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nSet default values took %.3f ms\n", (t2 - t1) * 1000);
    t1 = MPI_Wtime();
#endif
    read_input(&SPARC_Input, pSPARC); // read input file

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nReading input file took %.3f ms\n", (t2 - t1) * 1000);
#endif

    // broadcast the parameters read from the input files
    MPI_Bcast(&SPARC_Input, 1, SPARC_INPUT_MPI, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    t1 = MPI_Wtime();
#endif
    read_ion(&SPARC_Input, pSPARC); // read ion file

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nReading ion file took %.3f ms\n", (t2 - t1) * 1000);
    t1 = MPI_Wtime();
#endif

    // broadcast Ntypes read from ion file
    MPI_Ibcast(&pSPARC->Ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD, &req);

    // disable default pseudopotential path and name
    if (pSPARC->is_default_psd) {
      printf("\n"
             "Default path and names of pseudopotentials are currently disabled,\n"
             "Please specify pseudopotential filename!\n");
      exit(EXIT_FAILURE);
    }

    // read_pseudopotential_TM(&SPARC_Input, pSPARC); // read TM format
    // pseudopotential file
    read_pseudopotential_PSP(&SPARC_Input,
                             pSPARC); // read psp format pseudopotential file

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nReading pseudopotential file took %.3f ms\n", (t2 - t1) * 1000);
#endif

  } else {

#ifdef DEBUG
    t1 = MPI_Wtime();
#endif
    MPI_Bcast(&SPARC_Input, 1, SPARC_INPUT_MPI, 0, MPI_COMM_WORLD);
#ifdef DEBUG
    t2 = MPI_Wtime();
    if (rank == 0)
      printf("Broadcasting the input parameters took %.3f ms\n", (t2 - t1) * 1000);
#endif
    // broadcast Ntypes read from ion file
    MPI_Ibcast(&pSPARC->Ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD, &req);
  }
#ifdef DEBUG
  t1 = MPI_Wtime();
#endif

  MPI_Type_free(&SPARC_INPUT_MPI); // free the new MPI datatype

#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("\nFreeing SPARC_INPUT_MPI datatype took %.3f ms\n", (t2 - t1) * 1000);
  t1 = MPI_Wtime();
#endif

  // Ntypes is no longer read from .inpt file
  // pSPARC->Ntypes = SPARC_Input.Ntypes;
  // make sure Ntypes is broadcasted
  MPI_Wait(&req, MPI_STATUS_IGNORE);

  // broadcast SPARC members regarding Atom info. using MPI_Pack & MPI_Unpack
  bcast_SPARC_Atom(pSPARC);

#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("Broadcasting Atom info. using MPI_Pack & MPI_Unpack in SPARC took "
           "%.3f ms\n",
           (t2 - t1) * 1000);

  t1 = MPI_Wtime();
#endif

  // copy the data read from input files into struct SPARC
  SPARC_copy_input(pSPARC, &SPARC_Input);

#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("\nrank = %d, Copying data from SPARC_Input into SPARC & set up "
           "subcomm took %.3f ms\n",
           rank, (t2 - t1) * 1000);
#endif

  // set up sub-communicators
  if (pSPARC->SQFlag == 1) {
    Setup_Comms_SQ(pSPARC);
  } else {
    Setup_Comms(pSPARC);

#ifdef USE_DP_SUBEIG
#if !defined(USE_MKL) && !defined(USE_SCALAPACK)
    if (pSPARC->useLAPACK == 0) {
#ifdef DEBUG
      if (rank == 0)
        printf("[WARNING] ScaLAPACK is not compiled and Nstates > MAX_NS, "
               "subspace eigen-problem will be solved in sequential.\n");
#endif
      pSPARC->useLAPACK = 1;
    }
#endif

// SPARCX_ACCEL_NOTE Need to add this. Make sure it is always within "#ifdef
// USE_DP_SUBEIG" branch
// --- BEGIN. Alfredo Metere
#ifdef ACCEL // Activating flag for using hardware acceleration at compile time.
    pSPARC->useACCEL = 1;
    //	#else
    //	pSPARC->useACCEL = 0;

    if (rank == 0) {
      char *hwaccel[2] = {"DISABLED", "ENABLED"};
      printf("[INFO] Hardware acceleration is %s\n", hwaccel[pSPARC->useACCEL]);
    }
#endif // ACCEL
    // --- END. Alfredo Metere

    pSPARC->DP_CheFSI = NULL;
    pSPARC->DP_CheFSI_kpt = NULL;
    if (pSPARC->isGammaPoint)
      init_DP_CheFSI(pSPARC);
    else
      init_DP_CheFSI_kpt(pSPARC);
#endif

    // calculate maximum number of processors for eigenvalue solver
    if (pSPARC->useLAPACK == 0) {
      if (pSPARC->eig_paral_maxnp < 0) {
        char RorC, SorG;
        RorC = (pSPARC->isGammaPoint) ? 'R' : 'C';
        SorG = (pSPARC->StandardEigenFlag) ? 'S' : 'G';
        pSPARC->eig_paral_maxnp = parallel_eigensolver_max_processor(pSPARC->Nstates, RorC, SorG);
      }

      int gridsizes[2] = {pSPARC->Nstates, pSPARC->Nstates}, ierr = 1, size_blacscomm = 0;
      if (pSPARC->blacscomm != MPI_COMM_NULL)
        MPI_Comm_size(pSPARC->blacscomm, &size_blacscomm);
      SPARC_Dims_create(min(size_blacscomm, pSPARC->eig_paral_maxnp), 2, gridsizes, 1, pSPARC->eig_paral_subdims,
                        &ierr);
      if (ierr)
        pSPARC->eig_paral_subdims[0] = pSPARC->eig_paral_subdims[1] = 1;

#ifdef DEBUG
      if (rank == 0)
        printf("\nMaximun number of processors for eigenvalue solver is %d\n", pSPARC->eig_paral_maxnp);
      if (rank == 0)
        printf("The dimension of subgrid for eigen sovler is (%d x %d).\n", pSPARC->eig_paral_subdims[0],
               pSPARC->eig_paral_subdims[1]);
#endif
    }
  }

  // Allocate memory space for Exx methods
  if (pSPARC->usefock == 1) {
    init_exx(pSPARC);
  }

#ifdef DEBUG
  t1 = MPI_Wtime();
#endif

  // calculate spline derivatives for interpolation
  Calculate_SplineDerivRadFun(pSPARC);
#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("\nCalculate_SplineDerivRadFun took %.3f ms\n", (t2 - t1) * 1000);
#endif

  if (pSPARC->SQFlag == 1) {
    init_SQ(pSPARC);
  } else {
    // calculate indices for storing nonlocal inner product
    CalculateNonlocalInnerProductIndex(pSPARC);
    if (pSPARC->SOC_Flag == 1)
      CalculateNonlocalInnerProductIndexSOC(pSPARC);
  }

#ifdef DEBUG
  t1 = MPI_Wtime();
#endif
  // calculate pseudocharge density cutoff ("rb")
  Calculate_PseudochargeCutoff(pSPARC);
#ifdef DEBUG
  t2 = MPI_Wtime();
  if (rank == 0)
    printf("\nCalculating rb for all atom types took %.3f ms\n", (t2 - t1) * 1000);
#endif

  // initialize DFT-D3
  if (pSPARC->d3Flag == 1) {
    if ((strcmpi(pSPARC->XC, "GGA_PBE") != 0) && (strcmpi(pSPARC->XC, "GGA_PBEsol") != 0) &&
        (strcmpi(pSPARC->XC, "GGA_RPBE") != 0) && (strcmpi(pSPARC->XC, "PBE0") != 0) &&
        (strcmpi(pSPARC->XC, "HSE") != 0)) {
      if (rank == 0)
        printf(RED "ERROR: Cannot find D3 coefficients for this functional. "
                   "DFT-D3 correction calculation canceled!\n" RESET);
      exit(EXIT_FAILURE);
    } else {
      set_D3_coefficients(pSPARC); // this function is moved from electronicGroundState.c
    }
  }

  // initialize vdW-DF
  if (pSPARC->vdWDFFlag != 0) {
    vdWDF_initial_read_kernel(pSPARC); // read kernel function and 2nd derivative of spline functions
    // printf("rank %d, d2 of kernel function vdWDFd2Phidk2[2][4]=%.9e\n", rank,
    // pSPARC->vdWDFd2Phidk2[2][4]); // to verify it
  }

  // initialize metaGGA
  if (pSPARC->mGGAflag == 1) {
    initialize_MGGA(pSPARC);
  }

  // estimate memory usage
  pSPARC->memory_usage = estimate_memory(pSPARC);

  // write initialized parameters into output file
  if (rank == 0) {
    write_output_init(pSPARC);
  }
}
