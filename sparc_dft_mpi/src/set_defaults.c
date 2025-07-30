#include <stdio.h>

#include "initialization.h"

/**
 * @brief   Set default values for SPARC members.
 */
void set_defaults(SPARC_INPUT_OBJ *pSPARC_Input, SPARC_OBJ *pSPARC) {
  /* init strings to null */
  memset(pSPARC_Input->filename_out, '\0', sizeof(pSPARC_Input->filename_out));
  memset(pSPARC_Input->SPARCROOT, '\0', sizeof(pSPARC_Input->SPARCROOT));
  memset(pSPARC_Input->MDMeth, '\0', sizeof(pSPARC_Input->MDMeth));
  memset(pSPARC_Input->RelaxMeth, '\0', sizeof(pSPARC_Input->RelaxMeth));
  memset(pSPARC_Input->XC, '\0', sizeof(pSPARC_Input->XC));

  /* default file names (unless otherwise supplied) */
  snprintf(pSPARC_Input->filename_out, L_STRING, "%s", pSPARC_Input->filename);
  snprintf(pSPARC_Input->SPARCROOT, L_STRING, "%s", "UNDEFINED");

  /* Parallelizing parameters */
  pSPARC_Input->npspin = 0;              // number of processes for paral. over k-points
  pSPARC_Input->npkpt = 0;               // number of processes for paral. over k-points
  pSPARC_Input->npband = 0;              // number of processes for paral. over bands
  pSPARC_Input->npNdx = 0;               // number of processes for paral. over domain in x-dir
  pSPARC_Input->npNdy = 0;               // number of processes for paral. over domain in y-dir
  pSPARC_Input->npNdz = 0;               // number of processes for paral. over domain in z-dir
  pSPARC_Input->npNdx_phi = 0;           // number of processes for calculating phi in
                                         // paral. over domain in x-dir
  pSPARC_Input->npNdy_phi = 0;           // number of processes for calculating phi in
                                         // paral. over domain in y-dir
  pSPARC_Input->npNdz_phi = 0;           // number of processes for calculating phi in
                                         // paral. over domain in z-dir
  pSPARC_Input->eig_serial_maxns = 1500; // maximum Nstates for solving the subspace eigenproblem in serial
                                         // by default, for Nstates greater than this value, a parallel
                                         // methods will be used instead, unless ScaLAPACK is not compiled or
                                         // useLAPACK is turned off.
  pSPARC_Input->eig_paral_blksz = 128;   // block size for distributing the subspace eigenproblem
  pSPARC_Input->eig_paral_orfac = 0.0;   // no reorthogonalization when using p?syevx or p?sygvx
  pSPARC_Input->eig_paral_maxnp = -1;    // using default value from linear fitting model

  /* default spin_typ */
  pSPARC_Input->spin_typ = 0; // Default is spin unpolarized calculation

  /* default MD and Relaxation options */
  pSPARC_Input->MDFlag = 0; // default: no MD
  strncpy(pSPARC_Input->MDMeth, "NVT_NH",
          sizeof(pSPARC_Input->MDMeth)); // default MD method: NVT
  pSPARC_Input->RelaxFlag = 0;           // default: no relaxation
  strncpy(pSPARC_Input->RelaxMeth, "LBFGS",
          sizeof(pSPARC_Input->RelaxMeth)); // default relax method: LBFGS

  pSPARC_Input->RestartFlag = 0; // default: no retart
  /* default finite difference scheme info. */
  pSPARC_Input->order = 12; // default FD order: 12th

  /* default poisson solver */
  pSPARC_Input->Poisson_solver = 0; // default AAR solver
  /* Iterations: tolerances and max_iters */
  pSPARC_Input->FixRandSeed = 0;                 // default flag for fixing random numbers for MPI paralllelization
  pSPARC_Input->accuracy_level = -1;             // default accuracy level (2 - 'medium', 1e-3 in energy and forces)
  pSPARC_Input->scf_err_type = 0;                // default scf error definition: relative error in rho or veff
  pSPARC_Input->MAXIT_SCF = 100;                 // default maximum number of SCF iterations
  pSPARC_Input->MINIT_SCF = 2;                   // default minimum number of SCF iterations
  pSPARC_Input->MAXIT_POISSON = 3000;            // default maximum number of iterations for Poisson solver
  pSPARC_Input->Relax_Niter = 300;               // default maximum number of relaxation
                                                 // iterations, for RelaxFlag = 1 only
  pSPARC_Input->target_energy_accuracy = -1.0;   // default target energy accuracy
  pSPARC_Input->target_force_accuracy = -1.0;    // default target force accuracy
  pSPARC_Input->TOL_SCF = -1.0;                  // default SCF tolerance
  pSPARC_Input->TOL_RELAX = 5e-4;                // default Relaxation tolerance
  pSPARC_Input->TOL_POISSON = -1.0;              // default Poisson solve tolerance (will be set up later)
  pSPARC_Input->TOL_LANCZOS = 1e-2;              // default Lanczos tolerance
  pSPARC_Input->TOL_PSEUDOCHARGE = -1.0;         // default tolerance for calculating pseudocharge density radius
                                                 // (will be set up later)
  pSPARC_Input->TOL_PRECOND = -1.0;              // default Kerker tolerance will be set up
                                                 // later, depending on the mesh size
  pSPARC_Input->precond_kerker_kTF = 1.0;        // Thomas-Fermi screening length in the Kerker preconditioner
  pSPARC_Input->precond_kerker_thresh = 0.1;     // Threshold for the truncated Kerker preconditioner
  pSPARC_Input->precond_kerker_kTF_mag = 1.0;    // Thomas-Fermi screening length in the Kerker preconditioner
  pSPARC_Input->precond_kerker_thresh_mag = 0.1; // Threshold for the truncated Kerker preconditioner
  pSPARC_Input->precond_resta_q0 = 1.36;
  pSPARC_Input->precond_resta_Rs = 2.76;

  pSPARC_Input->REFERENCE_CUTOFF = 0.5; // default reference cutoff radius for nonlocal pseudopotential

  /* default mixing */
  pSPARC_Input->MixingVariable = -1;             // default mixing variabl (will be set to 'density' mixing)
  pSPARC_Input->MixingPrecond = -1;              // default mixing preconditioner (will be set later)
  pSPARC_Input->MixingPrecondMag = -1;           // default mixing preconditioner for magnetization density/potential
                                                 // (will be set later)
  pSPARC_Input->MixingParameter = 0.3;           // default mixing parameter
  pSPARC_Input->MixingParameterSimple = -1.0;    // default mixing parameter for simple mixing step (will be set up
                                                 // later)
  pSPARC_Input->MixingParameterMag = -1.0;       // default mixing parameter for magnetization density/potential
  pSPARC_Input->MixingParameterSimpleMag = -1.0; // default mixing parameter for magnetization density/potential in
                                                 // simple mixing step (will be set up later)
  pSPARC_Input->MixingHistory = 7;               // default mixing history
  pSPARC_Input->PulayFrequency = 1;              // default Pulay frequency
  pSPARC_Input->PulayRestartFlag = 0;            // default Pulay restart flag
  pSPARC_Input->precond_fitpow = 2;              // default fit power for the real-space preconditioner

  /* default k-points info. */
  pSPARC_Input->Kx = 1;
  pSPARC_Input->Ky = 1;
  pSPARC_Input->Kz = 1;
  pSPARC_Input->Nkpts = 1;
  pSPARC_Input->NkptsGroup = 1; // unused
  pSPARC_Input->kctr = 1;       // unused
  pSPARC_Input->kptshift[0] = 0.0;
  pSPARC_Input->kptshift[1] = 0.0;
  pSPARC_Input->kptshift[2] = 0.0;

  // default lattice vector array
  for (int i = 0; i < 9; i++) {
    if (i % 4 == 0)
      pSPARC_Input->LatVec[i] = 1.0;
    else
      pSPARC_Input->LatVec[i] = 0.0;
  }

  /* discretization */
  pSPARC_Input->mesh_spacing = -1.0;
  pSPARC_Input->ecut = -1.0;
  pSPARC_Input->numIntervals_x = -1;
  pSPARC_Input->numIntervals_y = -1;
  pSPARC_Input->numIntervals_z = -1;

  /* default system info. */
  pSPARC_Input->range_x = -1.0;
  pSPARC_Input->range_y = -1.0;
  pSPARC_Input->range_z = -1.0;
  pSPARC_Input->Flag_latvec_scale = 0;
  pSPARC_Input->latvec_scale_x = -1.0;
  pSPARC_Input->latvec_scale_y = -1.0;
  pSPARC_Input->latvec_scale_z = -1.0;
  pSPARC_Input->BC = -1; // default BC will be set up after reading input
  pSPARC_Input->BCx = -1;
  pSPARC_Input->BCy = -1;
  pSPARC_Input->BCz = -1;
  // pSPARC_Input->Beta = 1000;                // electronic smearing
  // (1/(k_B*T)) [1/Ha] pSPARC_Input->elec_T = 315.7751307269723; // default
  // electronic temperature in Kelvin
  pSPARC_Input->Beta = -1.0;   // electronic smearing (1/(k_B*T)) [1/Ha], will be specified later
  pSPARC_Input->elec_T = -1.0; // default electronic temperature in Kelvin, will be specified later
  pSPARC_Input->Ntypes = -1;   // default number of atom types
  pSPARC_Input->Nstates = -1;  // default number of states
  pSPARC_Input->NetCharge = 0; // default net charge: 0
  // strncpy(pSPARC_Input->XC, "LDA",sizeof(pSPARC_Input->XC));          //
  // default exchange-correlation approx: LDA
  strncpy(pSPARC_Input->XC, "UNDEFINED",
          sizeof(pSPARC_Input->XC)); // default: UNDEFINED

  /* default Chebyshev filter */
  pSPARC_Input->ChebDegree = -1;      // default chebyshev polynomial degree (will be
                                      // automatically found based on spectral width)
  pSPARC_Input->CheFSI_Optmz = 0;     // default is off
  pSPARC_Input->chefsibound_flag = 0; // default is to find bound using Lanczos on H in the first SCF of each
                                      // MD/Relax only
  pSPARC_Input->rhoTrigger = 4;       // default step to start updating electron
                                      // density, later will be subtracted by 1

  /* default smearing */
  pSPARC_Input->elec_T_type = 1; // default smearing type: 1 - gaussian smearing
                                 // (the other option is 0 - fermi-dirac)

  /* Default MD parameters */
  pSPARC_Input->MD_dt = 1.0;                 // default MD time step: 1.0 femtosecond
  pSPARC_Input->MD_Nstep = 10000000;         // default MD maximum steps
  pSPARC_Input->ion_T = -1.0;                // default ionic temperature in Kelvin
  pSPARC_Input->thermos_Tf = -1.0;           // Final temperature of the thermostat
  pSPARC_Input->ion_elec_eqT = 0;            // default ionic and electronic temp will be different throughout MD
  pSPARC_Input->ion_vel_dstr = 2;            // default initial velocity distribution is Maxwell-Boltzmann
  pSPARC_Input->ion_vel_dstr_rand = 0;       // default initial velocity are fixed
                                             // (different runs give the same answer)
  pSPARC_Input->qmass = 40.0 * CONST_FS2ATU; // default values of thermostat parameter mass (a.u.)
  pSPARC_Input->TWtime = 1000000000;         // default value of walltime in min
  pSPARC_Input->NPTscaleVecs[0] = 1;
  pSPARC_Input->NPTscaleVecs[1] = 1;
  pSPARC_Input->NPTscaleVecs[2] = 1; // default lattice vectors to be rescaled in NPT
  pSPARC_Input->NPT_NHnnos = 0;      // default amount of thermo variable for NPT_NH. If MDMeth is this but
                                     // nnos is 0, program will stop
  for (int subscript_NPTNH_qmass = 0; subscript_NPTNH_qmass < L_QMASS; subscript_NPTNH_qmass++) {
    pSPARC_Input->NPT_NHqmass[subscript_NPTNH_qmass] = 0.0;
  } // default mass of thermo variables for NPT_NH. If MDMeth is this but one of
    // qmass is 0, program will stop
  pSPARC_Input->NPT_NHbmass = 0.0; // default mass of baro variable for NPT_NH. If MDMeth is this but
                                   // bmass is 0, program will stop
  pSPARC_Input->prtarget = 0.0;    // default target pressure for NPT_NH.

  pSPARC_Input->NPT_NP_qmass = 0.0; // default mass of thermo variables for NPT_NP. If MDMeth is this but
                                    // qmass is 0, program will stop
  pSPARC_Input->NPT_NP_bmass = 0.0; // default mass of thermo variables for NPT_NP. If MDMeth is this but
                                    // bmass is 0, program will stop

  /* Default Relax parameters */
  pSPARC_Input->NLCG_sigma = 0.5;
  pSPARC_Input->L_history = 20;
  pSPARC_Input->L_finit_stp = 5e-3;
  pSPARC_Input->L_maxmov = 0.2;
  pSPARC_Input->L_autoscale = 1;
  pSPARC_Input->L_lineopt = 1;
  pSPARC_Input->L_icurv = 1.0;
  pSPARC_Input->FIRE_dt = 1.0;
  pSPARC_Input->FIRE_mass = 1.0;
  pSPARC_Input->FIRE_maxmov = 0.2;

  /* Default cell relaxation parameters*/
  pSPARC_Input->max_dilatation = 1.06; // maximum lattice dilatation
  pSPARC_Input->TOL_RELAX_CELL = 1e-2; // in GPa (all periodic)

  /* Default DFT-D3 correction */
  pSPARC_Input->d3Flag = 0;
  pSPARC_Input->d3Rthr = 1600.0;
  pSPARC_Input->d3Cn_thr = 625.0;

  /* Default stress flags*/
  pSPARC_Input->Calc_stress = 0;
  pSPARC_Input->Calc_pres = 0;

  /* print options */
  pSPARC_Input->Verbosity = 1;         // Flag for specifying the amount of output in .out file
  pSPARC_Input->PrintForceFlag = 1;    // flag for printing forces
  pSPARC_Input->PrintAtomPosFlag = 1;  // flag for printing atomic positions
  pSPARC_Input->PrintAtomVelFlag = 1;  // flag for printing atom velocities in case of MD/relax
  pSPARC_Input->PrintElecDensFlag = 0; // flag for printing final electron density
  pSPARC_Input->PrintEigenFlag = 0;    // Flag for printing final eigenvalues and occupations
  pSPARC_Input->PrintMDout = 1;        // Flag for printing MD output in a .aimd file
  pSPARC_Input->PrintRelaxout = 1;     // Flag for printing relax output in a .relax file
  pSPARC_Input->Printrestart = 1;      // Flag for printing output needed for restarting a simulation
  pSPARC_Input->Printrestart_fq = 1;   // Steps after which the output is written in the restart file
  pSPARC_Input->PrintPsiFlag[0] = 0;   // Flag for printing Kohn-Sham orbitals
  for (int i = 1; i < 7; i++)
    pSPARC_Input->PrintPsiFlag[i] = -1;  // defualt spin, kpt, band start and end index for printing psi
  pSPARC_Input->PrintEnergyDensFlag = 0; // flag for printing kinetic energy density

  /* Default pSPARC members */
  pSPARC->is_default_psd = 0; // default pseudopotential path is disabled

  /* Eigenvalue Problem */
  pSPARC_Input->StandardEigenFlag = 0; // default using standard eigenvalue problem is disabled

  /* Default Exact exchange potential */
  pSPARC_Input->TOL_FOCK = -1.0;         // default tolerance for Fock operator
  pSPARC_Input->TOL_SCF_INIT = -1.0;     // default tolerance for first PBE SCF
  pSPARC_Input->MAXIT_FOCK = 20;         // default maximum number of iterations for Hartree-Fock outer loop
  pSPARC_Input->MINIT_FOCK = 2;          // default minimum number of iterations for Hartree-Fock outer loop
  pSPARC_Input->EXXMeth_Flag = 0;        // default method to solve Poisson's equation
                                         // of Exact Exchange in Fourier space
  pSPARC_Input->ACEFlag = 1;             // default setting for not using ACE operator
  pSPARC_Input->EXXMem_batch = 20;       // default setting for using high memory option
  pSPARC_Input->EXXACEVal_state = 3;     // default setting for using high memory option
  pSPARC_Input->EXXDownsampling[0] = 1;  // default setting for downsampling, using full k-points
  pSPARC_Input->EXXDownsampling[1] = 1;  // default setting for downsampling, using full k-points
  pSPARC_Input->EXXDownsampling[2] = 1;  // default setting for downsampling, using full k-points
  pSPARC_Input->EXXDiv_Flag = -1;        // default setting for singularity in exact
                                         // exchange, default spherical trucation
  pSPARC_Input->hyb_range_fock = 0.1587; // default using VASP's HSE03 value
  pSPARC_Input->hyb_range_pbe = 0.1587;  // default using VASP's HSE03 value
  pSPARC_Input->exx_frac = -1;           // default exx_frac

  /* Default parameter for spin-orbit coupling */
  pSPARC->Nspinor = 1;
  pSPARC->SOC_Flag = 0;

  /* Default SQ method option */
  pSPARC_Input->SQFlag = 0;
  pSPARC_Input->SQ_gauss_mem = 0; // default not saving Lanczos vectors and eigenvectors
  pSPARC_Input->SQ_npl_g = -1;
  pSPARC_Input->SQ_rcut = -1;
  pSPARC_Input->SQ_tol_occ = 1e-6;
  pSPARC_Input->npNdx_SQ = 0;
  pSPARC_Input->npNdy_SQ = 0;
  pSPARC_Input->npNdz_SQ = 0;
}
