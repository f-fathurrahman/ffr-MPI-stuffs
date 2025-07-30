#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// For F_OK
#include <unistd.h> 

#include "initialization.h"
#include "eigenSolver.h"
#include "tools.h"

#define TEMP_TOL 1e-12

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))


/**
 * @brief   Copy the data read from input files into struct SPARC.
 */
void SPARC_copy_input(SPARC_OBJ *pSPARC, SPARC_INPUT_OBJ *pSPARC_Input) {
  int rank, nproc, i, p, FDn, Ntypes;
#ifdef DEBUG
  double t1, t2;
#endif
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  /* copy input values from input struct */
  // int type values
  pSPARC->num_node = pSPARC_Input->num_node;
  pSPARC->num_cpu_per_node = pSPARC_Input->num_cpu_per_node;
  pSPARC->num_acc_per_node = pSPARC_Input->num_acc_per_node;
  pSPARC->npspin = pSPARC_Input->npspin;
  pSPARC->npkpt = pSPARC_Input->npkpt;
  pSPARC->npband = pSPARC_Input->npband;
  pSPARC->npNdx = pSPARC_Input->npNdx;
  pSPARC->npNdy = pSPARC_Input->npNdy;
  pSPARC->npNdz = pSPARC_Input->npNdz;
  pSPARC->npNdx_phi = pSPARC_Input->npNdx_phi;
  pSPARC->npNdy_phi = pSPARC_Input->npNdy_phi;
  pSPARC->npNdz_phi = pSPARC_Input->npNdz_phi;
  pSPARC->eig_serial_maxns = pSPARC_Input->eig_serial_maxns;
  pSPARC->eig_paral_blksz = pSPARC_Input->eig_paral_blksz;
  pSPARC->spin_typ = pSPARC_Input->spin_typ;
  pSPARC->MDFlag = pSPARC_Input->MDFlag;
  pSPARC->RelaxFlag = pSPARC_Input->RelaxFlag;
  pSPARC->RestartFlag = pSPARC_Input->RestartFlag;
  pSPARC->Flag_latvec_scale = pSPARC_Input->Flag_latvec_scale;
  pSPARC->numIntervals_x = pSPARC_Input->numIntervals_x;
  pSPARC->numIntervals_y = pSPARC_Input->numIntervals_y;
  pSPARC->numIntervals_z = pSPARC_Input->numIntervals_z;
  pSPARC->BC = pSPARC_Input->BC;
  pSPARC->BCx = pSPARC_Input->BCx;
  pSPARC->BCy = pSPARC_Input->BCy;
  pSPARC->BCz = pSPARC_Input->BCz;
  pSPARC->Nstates = pSPARC_Input->Nstates;
  // pSPARC->Ntypes = pSPARC_Input->Ntypes;
  pSPARC->NetCharge = pSPARC_Input->NetCharge;
  pSPARC->order = pSPARC_Input->order;
  pSPARC->ChebDegree = pSPARC_Input->ChebDegree;
  pSPARC->CheFSI_Optmz = pSPARC_Input->CheFSI_Optmz;
  pSPARC->chefsibound_flag = pSPARC_Input->chefsibound_flag;
  pSPARC->rhoTrigger = pSPARC_Input->rhoTrigger; // pSPARC->rhoTrigger--;
  pSPARC->FixRandSeed = pSPARC_Input->FixRandSeed;
  pSPARC->accuracy_level = pSPARC_Input->accuracy_level;
  pSPARC->scf_err_type = pSPARC_Input->scf_err_type;
  pSPARC->MAXIT_SCF = pSPARC_Input->MAXIT_SCF;
  pSPARC->MINIT_SCF = pSPARC_Input->MINIT_SCF;
  pSPARC->MAXIT_POISSON = pSPARC_Input->MAXIT_POISSON;
  pSPARC->Relax_Niter = pSPARC_Input->Relax_Niter;
  pSPARC->MixingVariable = pSPARC_Input->MixingVariable;
  pSPARC->MixingPrecond = pSPARC_Input->MixingPrecond;
  pSPARC->MixingPrecondMag = pSPARC_Input->MixingPrecondMag;
  pSPARC->MixingHistory = pSPARC_Input->MixingHistory;
  pSPARC->PulayFrequency = pSPARC_Input->PulayFrequency;
  pSPARC->PulayRestartFlag = pSPARC_Input->PulayRestartFlag;
  pSPARC->precond_fitpow = pSPARC_Input->precond_fitpow;
  pSPARC->Nkpts = pSPARC_Input->Nkpts;
  pSPARC->Kx = pSPARC_Input->Kx;
  pSPARC->Ky = pSPARC_Input->Ky;
  pSPARC->Kz = pSPARC_Input->Kz;
  pSPARC->NkptsGroup = pSPARC_Input->NkptsGroup;
  pSPARC->kctr = pSPARC_Input->kctr;
  pSPARC->Verbosity = pSPARC_Input->Verbosity;
  pSPARC->PrintForceFlag = pSPARC_Input->PrintForceFlag;
  pSPARC->PrintAtomPosFlag = pSPARC_Input->PrintAtomPosFlag;
  pSPARC->PrintAtomVelFlag = pSPARC_Input->PrintAtomVelFlag;
  pSPARC->PrintEigenFlag = pSPARC_Input->PrintEigenFlag;
  pSPARC->PrintElecDensFlag = pSPARC_Input->PrintElecDensFlag;
  pSPARC->PrintMDout = pSPARC_Input->PrintMDout;
  pSPARC->PrintRelaxout = pSPARC_Input->PrintRelaxout;
  pSPARC->Printrestart = pSPARC_Input->Printrestart;
  pSPARC->Printrestart_fq = pSPARC_Input->Printrestart_fq;
  pSPARC->elec_T_type = pSPARC_Input->elec_T_type;
  pSPARC->MD_Nstep = pSPARC_Input->MD_Nstep;
  pSPARC->NPTscaleVecs[0] = pSPARC_Input->NPTscaleVecs[0];
  pSPARC->NPTscaleVecs[1] = pSPARC_Input->NPTscaleVecs[1];
  pSPARC->NPTscaleVecs[2] = pSPARC_Input->NPTscaleVecs[2];
  pSPARC->NPT_NHnnos = pSPARC_Input->NPT_NHnnos;
  pSPARC->ion_elec_eqT = pSPARC_Input->ion_elec_eqT;
  pSPARC->ion_vel_dstr = pSPARC_Input->ion_vel_dstr;
  pSPARC->ion_vel_dstr_rand = pSPARC_Input->ion_vel_dstr_rand;
  pSPARC->L_history = pSPARC_Input->L_history;
  pSPARC->L_autoscale = pSPARC_Input->L_autoscale;
  pSPARC->L_lineopt = pSPARC_Input->L_lineopt;
  pSPARC->Calc_stress = pSPARC_Input->Calc_stress;
  pSPARC->Calc_pres = pSPARC_Input->Calc_pres;
  pSPARC->d3Flag = pSPARC_Input->d3Flag;
  pSPARC->MAXIT_FOCK = pSPARC_Input->MAXIT_FOCK;
  pSPARC->EXXMeth_Flag = pSPARC_Input->EXXMeth_Flag;
  pSPARC->ACEFlag = pSPARC_Input->ACEFlag;
  pSPARC->EXXMem_batch = pSPARC_Input->EXXMem_batch;
  pSPARC->EXXACEVal_state = pSPARC_Input->EXXACEVal_state;
  pSPARC->EXXDiv_Flag = pSPARC_Input->EXXDiv_Flag;
  pSPARC->EXXDownsampling[0] = pSPARC_Input->EXXDownsampling[0];
  pSPARC->EXXDownsampling[1] = pSPARC_Input->EXXDownsampling[1];
  pSPARC->EXXDownsampling[2] = pSPARC_Input->EXXDownsampling[2];
  pSPARC->MINIT_FOCK = pSPARC_Input->MINIT_FOCK;
  pSPARC->SQFlag = pSPARC_Input->SQFlag;
  pSPARC->SQ_gauss_mem = pSPARC_Input->SQ_gauss_mem;
  pSPARC->SQ_npl_g = pSPARC_Input->SQ_npl_g;
  pSPARC->npNdx_SQ = pSPARC_Input->npNdx_SQ;
  pSPARC->npNdy_SQ = pSPARC_Input->npNdy_SQ;
  pSPARC->npNdz_SQ = pSPARC_Input->npNdz_SQ;
  for (i = 0; i < 7; i++)
    pSPARC->PrintPsiFlag[i] = pSPARC_Input->PrintPsiFlag[i];
  pSPARC->PrintEnergyDensFlag = pSPARC_Input->PrintEnergyDensFlag;
  pSPARC->StandardEigenFlag = pSPARC_Input->StandardEigenFlag;

  // double type values
  pSPARC->range_x = pSPARC_Input->range_x;
  pSPARC->range_y = pSPARC_Input->range_y;
  pSPARC->range_z = pSPARC_Input->range_z;
  if (pSPARC->range_x <= 0.0 || pSPARC->range_y <= 0.0 || pSPARC->range_z <= 0.0) {
    if (!rank)
      printf("\nERROR: Please specify valid CELL dimensions!\n");
    exit(EXIT_FAILURE);
  }
  pSPARC->latvec_scale_x = pSPARC_Input->latvec_scale_x;
  pSPARC->latvec_scale_y = pSPARC_Input->latvec_scale_y;
  pSPARC->latvec_scale_z = pSPARC_Input->latvec_scale_z;
  // allocate the memory for lattice vector array
  for (i = 0; i < 9; i++)
    pSPARC->LatVec[i] = pSPARC_Input->LatVec[i];
  pSPARC->mesh_spacing = pSPARC_Input->mesh_spacing;
  pSPARC->ecut = pSPARC_Input->ecut;
  pSPARC->kptshift[0] = pSPARC_Input->kptshift[0];
  pSPARC->kptshift[1] = pSPARC_Input->kptshift[1];
  pSPARC->kptshift[2] = pSPARC_Input->kptshift[2];
  pSPARC->target_energy_accuracy = pSPARC_Input->target_energy_accuracy;
  pSPARC->target_force_accuracy = pSPARC_Input->target_force_accuracy;
  pSPARC->TOL_SCF = pSPARC_Input->TOL_SCF;
  pSPARC->TOL_RELAX = pSPARC_Input->TOL_RELAX;
  pSPARC->TOL_POISSON = pSPARC_Input->TOL_POISSON;
  pSPARC->TOL_LANCZOS = pSPARC_Input->TOL_LANCZOS;
  pSPARC->TOL_PSEUDOCHARGE = pSPARC_Input->TOL_PSEUDOCHARGE;
  pSPARC->TOL_PRECOND = pSPARC_Input->TOL_PRECOND;
  pSPARC->POISSON_SOLVER = pSPARC_Input->Poisson_solver;
  pSPARC->precond_kerker_kTF = pSPARC_Input->precond_kerker_kTF;
  pSPARC->precond_kerker_thresh = pSPARC_Input->precond_kerker_thresh;
  pSPARC->precond_kerker_kTF_mag = pSPARC_Input->precond_kerker_kTF_mag;
  pSPARC->precond_kerker_thresh_mag = pSPARC_Input->precond_kerker_thresh_mag;
  pSPARC->precond_resta_q0 = pSPARC_Input->precond_resta_q0;
  pSPARC->precond_resta_Rs = pSPARC_Input->precond_resta_Rs;
  pSPARC->REFERENCE_CUTOFF = pSPARC_Input->REFERENCE_CUTOFF;
  pSPARC->Beta = pSPARC_Input->Beta;
  pSPARC->elec_T = pSPARC_Input->elec_T;
  pSPARC->MixingParameter = pSPARC_Input->MixingParameter;
  pSPARC->MixingParameterSimple = pSPARC_Input->MixingParameterSimple;
  pSPARC->MixingParameterMag = pSPARC_Input->MixingParameterMag;
  pSPARC->MixingParameterSimpleMag = pSPARC_Input->MixingParameterSimpleMag;
  pSPARC->MD_dt = pSPARC_Input->MD_dt;
  pSPARC->ion_T = pSPARC_Input->ion_T;
  pSPARC->thermos_Tf = pSPARC_Input->thermos_Tf;
  pSPARC->qmass = pSPARC_Input->qmass;
  pSPARC->TWtime = pSPARC_Input->TWtime; // buffer time
  for (i = 0; i < pSPARC->NPT_NHnnos; i++) {
    pSPARC->NPT_NHqmass[i] = pSPARC_Input->NPT_NHqmass[i];
  }
  pSPARC->NPT_NHbmass = pSPARC_Input->NPT_NHbmass;
  pSPARC->prtarget = pSPARC_Input->prtarget;
  pSPARC->NPT_NP_bmass = pSPARC_Input->NPT_NP_bmass;
  pSPARC->NPT_NP_qmass = pSPARC_Input->NPT_NP_qmass;
  pSPARC->NLCG_sigma = pSPARC_Input->NLCG_sigma;
  pSPARC->L_finit_stp = pSPARC_Input->L_finit_stp;
  pSPARC->L_maxmov = pSPARC_Input->L_maxmov;
  pSPARC->L_icurv = pSPARC_Input->L_icurv;
  pSPARC->FIRE_dt = pSPARC_Input->FIRE_dt;
  pSPARC->FIRE_mass = pSPARC_Input->FIRE_mass;
  pSPARC->FIRE_maxmov = pSPARC_Input->FIRE_maxmov;
  pSPARC->max_dilatation = pSPARC_Input->max_dilatation;
  pSPARC->TOL_RELAX_CELL = pSPARC_Input->TOL_RELAX_CELL;
  pSPARC->eig_paral_orfac = pSPARC_Input->eig_paral_orfac;
  pSPARC->eig_paral_maxnp = pSPARC_Input->eig_paral_maxnp;
  pSPARC->d3Rthr = pSPARC_Input->d3Rthr;
  pSPARC->d3Cn_thr = pSPARC_Input->d3Cn_thr;
  pSPARC->TOL_FOCK = pSPARC_Input->TOL_FOCK;
  pSPARC->TOL_SCF_INIT = pSPARC_Input->TOL_SCF_INIT;
  pSPARC->hyb_range_fock = pSPARC_Input->hyb_range_fock;
  pSPARC->hyb_range_pbe = pSPARC_Input->hyb_range_pbe;
  pSPARC->exx_frac = pSPARC_Input->exx_frac;
  pSPARC->SQ_rcut = pSPARC_Input->SQ_rcut;
  pSPARC->SQ_tol_occ = pSPARC_Input->SQ_tol_occ;

  // char type values
  strncpy(pSPARC->MDMeth, pSPARC_Input->MDMeth, sizeof(pSPARC->MDMeth));
  strncpy(pSPARC->RelaxMeth, pSPARC_Input->RelaxMeth, sizeof(pSPARC->RelaxMeth));
  strncpy(pSPARC->XC, pSPARC_Input->XC, sizeof(pSPARC->XC));
  if (strcmp(pSPARC->XC, "UNDEFINED") == 0) {
    if (!rank)
      printf("\nERROR: Please specify XC type!\n");
    exit(EXIT_FAILURE);
  }
  strncpy(pSPARC->filename, pSPARC_Input->filename, sizeof(pSPARC->filename));
  strncpy(pSPARC->filename_out, pSPARC_Input->filename_out, sizeof(pSPARC->filename_out));
  strncpy(pSPARC->SPARCROOT, pSPARC_Input->SPARCROOT, sizeof(pSPARC->SPARCROOT));

  // check XC compatibility with pseudopotential
  pSPARC->usefock = 0; // default: no fock operator
  int ixc = 0;         // input XC
  if (strcmpi(pSPARC->XC, "LDA_PZ") == 0) {
    ixc = 2;
  } else if (strcmpi(pSPARC->XC, "LDA_PW") == 0) {
    ixc = 7;
  } else if (strcmpi(pSPARC->XC, "GGA_PBE") == 0) {
    ixc = 11;
  } else if (strcmpi(pSPARC->XC, "GGA_RPBE") == 0) {
    ixc = 15;
  } else if (strcmpi(pSPARC->XC, "GGA_PBEsol") == 0) {
    ixc = 116;
  } else if (strcmpi(pSPARC->XC, "HF") == 0) {
    ixc = 40;
    pSPARC->usefock = 1;
    if (pSPARC->exx_frac < 0)
      pSPARC->exx_frac = 1;
    if (fabs(1 - pSPARC->exx_frac) > TEMP_TOL) {
      if (!rank)
        printf("ERROR: HF functional could be only defined with 1.0 EXX_FRAC.\n");
    }
  } else if (strcmpi(pSPARC->XC, "PBE0") == 0) {
    ixc = 41;
    pSPARC->usefock = 1;
    if (pSPARC->exx_frac < 0)
      pSPARC->exx_frac = 0.25;
    if (fabs(0.25 - pSPARC->exx_frac) > TEMP_TOL) {
      if (!rank)
        printf("Note: You are using PBE0 with %.5g exact exchange.\n", pSPARC->exx_frac);
    }
  } else if (strcmpi(pSPARC->XC, "HSE") == 0) {
    ixc = 427;
    pSPARC->usefock = 1;
    if (pSPARC->exx_frac < 0)
      pSPARC->exx_frac = 0.25;
    if (fabs(0.25 - pSPARC->exx_frac) > TEMP_TOL) {
      if (!rank)
        printf("Note: You are using HSE with %.5g exact exchange.\n", pSPARC->exx_frac);
    }
  } else if (strcmpi(pSPARC->XC, "SCAN") == 0) {
    ixc = -263267;
  } else if (strcmpi(pSPARC->XC, "vdWDF1") == 0) {
    ixc = -102; // this is the index of Zhang-Yang revPBE exchange in Libxc
  } else if (strcmpi(pSPARC->XC, "vdWDF2") == 0) {
    ixc = -108; // this is the index of PW86 exchange in Libxc
  }
  for (int ityp = 0; ityp < pSPARC->Ntypes; ityp++) {
    if (!pSPARC->usefock && pSPARC->psd[ityp].pspxc != ixc) {
      if (!rank)
        printf(YEL "\nWARNING: Pseudopotential file for atom type %s has pspxc = %d,\n"
                   "not equal to input ixc = %d (%s). Be careful with the "
                   "result.\n" RESET,
               &pSPARC->atomType[ityp * L_ATMTYPE], pSPARC->psd[ityp].pspxc, ixc, pSPARC->XC);
    }
    if (pSPARC->usefock && pSPARC->psd[ityp].pspxc != 11) {
      if (!rank)
        printf(YEL "\nWARNING: Pseudopotential file for atom type %s has pspxc = %d,\n"
                   "while hybrid calculation needs a PBE pseudopotential. Be careful "
                   "with the result.\n" RESET,
               &pSPARC->atomType[ityp * L_ATMTYPE], pSPARC->psd[ityp].pspxc);
    }
  }

  if (strcmpi(pSPARC->XC, "HSE") == 0) {
    // pSPARC->hyb_range_fock = 0.106;     // QE's value
    // pSPARC->hyb_range_pbe = 0.106;      // QE's value
    // pSPARC->hyb_range_fock = 0.106066017177982;     // ABINIT's value
    // pSPARC->hyb_range_pbe = 0.188988157484231;      // ABINIT's value
    if (!rank) {
      printf("Careful: You are using HSE with range-separation parameter "
             "omega_HF = %.6f (1/Bohr) and omega_PBE = %.6f (1/Bohr)\n",
             pSPARC->hyb_range_fock, pSPARC->hyb_range_pbe);
      printf("If you want to change it, please use EXX_RANGE_FOCK and "
             "EXX_RANGE_PBE input options.\n");
    }
  } else {
    pSPARC->hyb_range_fock = -1;
    pSPARC->hyb_range_pbe = -1;
  }

  // check MDMeth availability
  if ((strcmpi(pSPARC->MDMeth, "NVT_NH") && strcmpi(pSPARC->MDMeth, "NVE") && strcmpi(pSPARC->MDMeth, "NVK_G") &&
       strcmpi(pSPARC->MDMeth, "NPT_NH") && strcmpi(pSPARC->MDMeth, "NPT_NP")) != 0) {
    if (!rank) {
      printf("\nCannot recognize MDMeth = \"%s\"\n", pSPARC->MDMeth);
      printf("MDMeth (MD Method) must be one of the following:\n\tNVT_NH\t "
             "NVE\t NVK_G\t NPT_NH\t NPT_NP\n");
    }
    exit(EXIT_FAILURE);
  }

  /* process the data read from input files */
  Ntypes = pSPARC->Ntypes;

  // calculate total number of electrons
  pSPARC->Nelectron = 0;
  for (i = 0; i < Ntypes; i++) {
    pSPARC->Nelectron += pSPARC->Znucl[i] * pSPARC->nAtomv[i];
  }
  pSPARC->Nelectron -= pSPARC->NetCharge;

  // check if NLCC is present
  int NLCC_flag = 0;
  for (int ityp = 0; ityp < Ntypes; ityp++) {
    if (pSPARC->psd[ityp].fchrg > TEMP_TOL) {
      NLCC_flag = 1;
      break;
    }
  }
  pSPARC->NLCC_flag = NLCC_flag;

  // check if exchange-correlation functional is metaGGA
  pSPARC->mGGAflag = 0;
  if (strcmpi(pSPARC->XC, "SCAN") == 0) { // it can be expand, such as adding r2SCAN
    if (pSPARC->NLCC_flag) {
      if (!rank)
        printf("\nERROR: currently SCAN functional does not support applying "
               "NLCC pseudopotential!\n");
      exit(EXIT_FAILURE);
    }
    pSPARC->mGGAflag = 1;
  }
  // check if exchange-correlation functional is vdW-DF1 or vdW-DF2
  pSPARC->vdWDFFlag = 0;
  if (strcmpi(pSPARC->XC, "vdWDF1") == 0) {
    pSPARC->vdWDFFlag = 1;
  }
  if (strcmpi(pSPARC->XC, "vdWDF2") == 0) {
    pSPARC->vdWDFFlag = 2;
  }

  // initialize energy values
  pSPARC->Esc = 0.0;
  pSPARC->Efermi = 0.0;
  pSPARC->Exc = 0.0;
  pSPARC->Eband = 0.0;
  pSPARC->Entropy = 0.0;
  pSPARC->Etot = 0.0;

  // estimate Nstates if not provided
  if (pSPARC->Nstates == -1) {
    // estimate Nstates using the linear function y = 1.2 * x + 5
    pSPARC->Nstates = (int)((pSPARC->Nelectron / 2) * 1.2 + 5) * pSPARC->Nspinor;
  }

  if (pSPARC->Nspinor == 1) {
    if (pSPARC->Nstates < (pSPARC->Nelectron / 2) && !pSPARC->SQFlag) {
      if (!rank)
        printf("\nERROR: number of states is less than Nelectron/2!\n");
      exit(EXIT_FAILURE);
    }
  } else if (pSPARC->Nspinor == 2) {
    if (pSPARC->Nstates < pSPARC->Nelectron && !pSPARC->SQFlag) {
      if (!rank)
        printf("\nERROR: number of states is less than Nelectron!\n");
      exit(EXIT_FAILURE);
    }
  }

  // filenames
  if (rank == 0) {
    snprintf(pSPARC->OutFilename, L_STRING, "%s.out", pSPARC->filename_out);
    snprintf(pSPARC->StaticFilename, L_STRING, "%s.static", pSPARC->filename_out);
    snprintf(pSPARC->AtomFilename, L_STRING, "%s.atom", pSPARC->filename_out);
    snprintf(pSPARC->EigenFilename, L_STRING, "%s.eigen", pSPARC->filename_out);
    snprintf(pSPARC->MDFilename, L_STRING, "%s.aimd", pSPARC->filename_out);
    snprintf(pSPARC->RelaxFilename, L_STRING, "%s.geopt", pSPARC->filename_out);
    snprintf(pSPARC->restart_Filename, L_STRING, "%s.restart", pSPARC->filename_out);
    snprintf(pSPARC->restartC_Filename, L_STRING, "%s.restart-0", pSPARC->filename_out);
    snprintf(pSPARC->restartP_Filename, L_STRING, "%s.restart-1", pSPARC->filename_out);
    snprintf(pSPARC->DensTCubFilename, L_STRING, "%s.dens", pSPARC->filename_out);
    snprintf(pSPARC->DensUCubFilename, L_STRING, "%s.densUp", pSPARC->filename_out);
    snprintf(pSPARC->DensDCubFilename, L_STRING, "%s.densDwn", pSPARC->filename_out);
    snprintf(pSPARC->OrbitalsFilename, L_STRING, "%s.psi", pSPARC->filename_out);
    snprintf(pSPARC->KinEnDensTCubFilename, L_STRING, "%s.kedens", pSPARC->filename_out);
    snprintf(pSPARC->KinEnDensUCubFilename, L_STRING, "%s.kedensUp", pSPARC->filename_out);
    snprintf(pSPARC->KinEnDensDCubFilename, L_STRING, "%s.kedensDwn", pSPARC->filename_out);
    snprintf(pSPARC->XcEnDensCubFilename, L_STRING, "%s.xcedens", pSPARC->filename_out);
    snprintf(pSPARC->ExxEnDensTCubFilename, L_STRING, "%s.exxedens", pSPARC->filename_out);
    snprintf(pSPARC->ExxEnDensUCubFilename, L_STRING, "%s.exxedensUp", pSPARC->filename_out);
    snprintf(pSPARC->ExxEnDensDCubFilename, L_STRING, "%s.exxedensDwn", pSPARC->filename_out);

    // check if the name for out file exits
    char temp_outfname[L_STRING];
    snprintf(temp_outfname, L_STRING, "%s", pSPARC->OutFilename);
#ifdef DEBUG
    t1 = MPI_Wtime();
#endif
    int MAX_OUTPUT = 100; // max number of output files allowed before we
                          // overwrite existing output files
    i = 0;
    while ((access(temp_outfname, F_OK) != -1) && i <= MAX_OUTPUT) {
      i++;
      snprintf(temp_outfname, L_STRING, "%s_%02d", pSPARC->OutFilename, i);
    }
    pSPARC->suffixNum = i; // note that this is only known to rank 0!

#ifdef DEBUG
    t2 = MPI_Wtime();
    printf("\nChecking existence of (%d) out file(s) took %.3f ms\n", i, (t2 - t1) * 1000);
#endif
    if (i >= (int)(MAX_OUTPUT * 0.75) && i <= MAX_OUTPUT) {
      printf("\n#WARNING: There's a limit on total number of output files "
             "allowed!\n"
             "         After the limit is reached, SPARC will start using  the\n"
             "         name provided exactly (without attached index) and old\n"
             "         results will be overwritten! Please move some existing\n"
             "         results to some other directories if you want to keep \n"
             "         storing results in different files!\n"
             "         Current usage: %.1f %%\n\n",
             i / (double)MAX_OUTPUT * 100);
    }

    if (i > MAX_OUTPUT) {
      printf("\n#WARNING: The limit of total number of output files is reached! \n"
             "         Current output name (without suffix): %s\n\n",
             pSPARC->filename_out);
    } else if (i > 0) {
      char tempchar[L_STRING];
      snprintf(tempchar, L_STRING, "%s", pSPARC->OutFilename);
      snprintf(pSPARC->OutFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->StaticFilename);
      snprintf(pSPARC->StaticFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->AtomFilename);
      snprintf(pSPARC->AtomFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->EigenFilename);
      snprintf(pSPARC->EigenFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->MDFilename);
      snprintf(pSPARC->MDFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->RelaxFilename);
      snprintf(pSPARC->RelaxFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->DensTCubFilename);
      snprintf(pSPARC->DensTCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->DensUCubFilename);
      snprintf(pSPARC->DensUCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->DensDCubFilename);
      snprintf(pSPARC->DensDCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->OrbitalsFilename);
      snprintf(pSPARC->OrbitalsFilename, L_STRING, "%s_%02d", tempchar, i);
      // energy density files
      snprintf(tempchar, L_STRING, "%s", pSPARC->KinEnDensTCubFilename);
      snprintf(pSPARC->KinEnDensTCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->KinEnDensUCubFilename);
      snprintf(pSPARC->KinEnDensUCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->KinEnDensDCubFilename);
      snprintf(pSPARC->KinEnDensDCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->XcEnDensCubFilename);
      snprintf(pSPARC->XcEnDensCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->ExxEnDensTCubFilename);
      snprintf(pSPARC->ExxEnDensTCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->ExxEnDensUCubFilename);
      snprintf(pSPARC->ExxEnDensUCubFilename, L_STRING, "%s_%02d", tempchar, i);
      snprintf(tempchar, L_STRING, "%s", pSPARC->ExxEnDensDCubFilename);
      snprintf(pSPARC->ExxEnDensDCubFilename, L_STRING, "%s_%02d", tempchar, i);
    }
  }
  // Not only rank 0 printing orbitals
  MPI_Bcast(pSPARC->OrbitalsFilename, L_STRING, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Initialize MD/relax variables
  pSPARC->RelaxCount = 0;   // initialize current relaxation step
  pSPARC->StressCount = 0;  // initialize current stress relaxation step (used in full relaxation)
  pSPARC->elecgs_Count = 0; // count the number of times forces are
                            // calculated(useful in relaxation)
  pSPARC->MDCount = 0;      // initialize current MD step
  pSPARC->StopCount = 0;
  pSPARC->restartCount = 0;
  pSPARC->amu2au = CONST_AMU2AU; // 1 au = 9.10938356e-31 Kg; 1 amu =  1.660539040e-27 Kg;
  pSPARC->fs2atu = CONST_FS2ATU; // 1atu = 2.418884326509e-17 s;
  pSPARC->kB = CONST_KB;         // Boltzmann constant in Ha/K

  if (pSPARC->RelaxFlag == 2 || pSPARC->RelaxFlag == 3) {
    pSPARC->Printrestart = 0;
  }

  // Value of xc_rhotol used in xc functional
  pSPARC->xc_rhotol = 1e-14;

  if (pSPARC->Beta < 0) {
    if (pSPARC->elec_T_type == 1) { // gaussian
      // The electronic temperature corresponding to 0.2 eV is 2320.904422 K
      pSPARC->Beta = CONST_EH / 0.2; // smearing = 0.2 eV = 0.00734986450 Ha, Beta := 1 / smearing
    } else {                         // fermi-dirac
      // The electronic temperature corresponding to 0.1 eV is 1160.452211 K
      pSPARC->Beta = CONST_EH / 0.1; // smearing = 0.1 eV = 0.00367493225 Ha, Beta := 1 / smearing
    }
    pSPARC->elec_T = 1. / (pSPARC->kB * pSPARC->Beta);
  }

  // Check the cell typ
  double mult;
  int j;
  pSPARC->cell_typ = 0; // orthogonal cell by default
  for (i = 0; i < 2; i++) {
    for (j = i + 1; j < 3; j++) {
      mult =
          fabs(pSPARC->LatVec[3 * i] * pSPARC->LatVec[3 * j] + pSPARC->LatVec[3 * i + 1] * pSPARC->LatVec[3 * j + 1] +
               pSPARC->LatVec[3 * i + 2] * pSPARC->LatVec[3 * j + 2]);
      if (mult > TEMP_TOL) {
        pSPARC->cell_typ = 1;
        i = j = 3;
      }
    }
  }

  // determine boundary conditions in each direction
  // BCx = 0 -> periodic, BCx = 1 -> dirichlet
  if (pSPARC->BC > 0) {
    if (pSPARC->BC == 1) {
      // dirichlet boundary
      pSPARC->BCx = 1;
      pSPARC->BCy = 1;
      pSPARC->BCz = 1;
    } else if (pSPARC->BC == 2 || pSPARC->BC == 0) {
      // periodic in all three directions
      pSPARC->BCx = 0;
      pSPARC->BCy = 0;
      pSPARC->BCz = 0;
    } else if (pSPARC->BC == 3) {
      // periodic in x and y directions
      pSPARC->BCx = 0;
      pSPARC->BCy = 0;
      pSPARC->BCz = 1;
    } else if (pSPARC->BC == 4) {
      // periodic in x direction
      pSPARC->BCx = 1;
      pSPARC->BCy = 1;
      pSPARC->BCz = 0;
    } else if (pSPARC->BC > 7) {
      exit(EXIT_FAILURE);
    }
  } else if (pSPARC->BCx >= 0 && pSPARC->BCy >= 0 && pSPARC->BCz >= 0) {
    int n_Dirichlet = pSPARC->BCx + pSPARC->BCy + pSPARC->BCz;
    switch (n_Dirichlet) {
    case 0:
      pSPARC->BC = 2;
      break;
    case 1:
      pSPARC->BC = 3;
      break;
    case 2:
      pSPARC->BC = 4;
      break;
    case 3:
      pSPARC->BC = 1;
      break;
    default:
      printf("Error in BC values\n");
      break;
    }
  } else {
    // if user does not provide any BC, set default to periodic in all
    // directions
    pSPARC->BC = 2;
    pSPARC->BCx = pSPARC->BCy = pSPARC->BCz = 0;
  }

  FDn = pSPARC->order / 2; // half the FD order
  // calculate number of finite-difference intervals in case it's provided
  // indirectly
  if (pSPARC->ecut > 0) {
    double h_ecut = Ecut2h(pSPARC->ecut, FDn);
    pSPARC->numIntervals_x = max(ceil(pSPARC->range_x / h_ecut), FDn);
    pSPARC->numIntervals_y = max(ceil(pSPARC->range_y / h_ecut), FDn);
    pSPARC->numIntervals_z = max(ceil(pSPARC->range_z / h_ecut), FDn);
  } else if (pSPARC->mesh_spacing > 0) {
    pSPARC->numIntervals_x = max(ceil(pSPARC->range_x / pSPARC->mesh_spacing), FDn);
    pSPARC->numIntervals_y = max(ceil(pSPARC->range_y / pSPARC->mesh_spacing), FDn);
    pSPARC->numIntervals_z = max(ceil(pSPARC->range_z / pSPARC->mesh_spacing), FDn);
  }

  // calculate number of nodes in each direction
  pSPARC->Nx = pSPARC->numIntervals_x + pSPARC->BCx;
  pSPARC->Ny = pSPARC->numIntervals_y + pSPARC->BCy;
  pSPARC->Nz = pSPARC->numIntervals_z + pSPARC->BCz;
  pSPARC->Nd = pSPARC->Nx * pSPARC->Ny * pSPARC->Nz;

  // mesh size
  pSPARC->delta_x = pSPARC->range_x / (pSPARC->numIntervals_x);
  pSPARC->delta_y = pSPARC->range_y / (pSPARC->numIntervals_y);
  pSPARC->delta_z = pSPARC->range_z / (pSPARC->numIntervals_z);
  pSPARC->dV = pSPARC->delta_x * pSPARC->delta_y * pSPARC->delta_z; // will be multiplied by the Jacobian later

  pSPARC->Jacbdet = 1.0;
  // Compute transformation matrices needed in non cartesian coordinate system
  // and perform atomic coordinate transformation from cartesian to cell
  // coordinates
  if (pSPARC->cell_typ < 20) {
    Cart2nonCart_transformMat(pSPARC);
  }

  // Convert cartesian coordinates of atom positions into cell coordinates
  int atm, count = 0;
  if (pSPARC->cell_typ != 0) {
    // Cart2nonCart_transformMat(pSPARC);
    for (i = 0; i < pSPARC->Ntypes; i++) {
      if (pSPARC->IsFrac[i] == 0) {
        for (atm = 0; atm < pSPARC->nAtomv[i]; atm++) {
          Cart2nonCart_coord(pSPARC, &pSPARC->atom_pos[3 * count], &pSPARC->atom_pos[3 * count + 1],
                             &pSPARC->atom_pos[3 * count + 2]);
          count++;
        }
      } else {
        count += pSPARC->nAtomv[i];
      }
    }
  }

  pSPARC->xin = 0.0; // starting coordinate (global) of the cell in the x-direction

  // Provide number of spin
  if (pSPARC->spin_typ == 0)
    pSPARC->Nspin = 1;
  else
    pSPARC->Nspin = 2;

  // Provide default spin if not already provided
  if (pSPARC->spin_typ != 0) { // spin polarized calculation
    srand(1);                  // TODO: provide this as a user input
    count = 0;
    for (i = 0; i < pSPARC->Ntypes; i++) {
      if (pSPARC->IsSpin[i] == 0) {
        for (atm = 0; atm < pSPARC->nAtomv[i]; atm++) {
          pSPARC->atom_spin[count] = -pSPARC->Znucl[i] + 2 * pSPARC->Znucl[i] * ((double)rand() / RAND_MAX);
          count++;
        }
      } else {
        count += pSPARC->nAtomv[i];
      }
    }
  }

  // Update the volume of the volume element bounded by finite difference grid
  pSPARC->dV *= pSPARC->Jacbdet;

  // map positions back into the domain if necessary
  if (pSPARC->BCx == 0) {
    double x;
    for (i = 0; i < pSPARC->n_atom; i++) {
      x = *(pSPARC->atom_pos + 3 * i);
      if (x < 0 || x > pSPARC->range_x) {
        x = fmod(x, pSPARC->range_x);
        *(pSPARC->atom_pos + 3 * i) = x + (x < 0.0) * pSPARC->range_x;
      }
    }
  } else if (rank == nproc - 1) {
    // for Dirichlet BC, exit if atoms goes out of the domain
    double x;
    for (i = 0; i < pSPARC->n_atom; i++) {
      x = *(pSPARC->atom_pos + 3 * i);
      if (x < pSPARC->xin || x > pSPARC->range_x + pSPARC->xin) {
        printf("\nERROR: position of atom # %d is out of the domain!\n", i + 1);
        exit(EXIT_FAILURE);
      }
    }
  }

  if (pSPARC->BCy == 0) {
    double y;
    for (i = 0; i < pSPARC->n_atom; i++) {
      y = *(pSPARC->atom_pos + 1 + 3 * i);
      if (y < 0 || y > pSPARC->range_y) {
        y = fmod(y, pSPARC->range_y);
        *(pSPARC->atom_pos + 1 + 3 * i) = y + (y < 0.0) * pSPARC->range_y;
      }
    }
  } else if (rank == nproc - 1) {
    // for Dirichlet BC, exit if atoms goes out of the domain
    double y;
    for (i = 0; i < pSPARC->n_atom; i++) {
      y = *(pSPARC->atom_pos + 1 + 3 * i);
      if (y < 0 || y > pSPARC->range_y) {
        printf("\nERROR: position of atom # %d is out of the domain!\n", i + 1);
        exit(EXIT_FAILURE);
      }
    }
  }

  if (pSPARC->BCz == 0) {
    double z;
    for (i = 0; i < pSPARC->n_atom; i++) {
      z = *(pSPARC->atom_pos + 2 + 3 * i);
      if (z < 0 || z > pSPARC->range_z) {
        z = fmod(z, pSPARC->range_z);
        *(pSPARC->atom_pos + 2 + 3 * i) = z + (z < 0.0) * pSPARC->range_z;
      }
    }
  } else if (rank == nproc - 1) {
    // for Dirichlet BC, exit if atoms goes out of the domain
    double z;
    for (i = 0; i < pSPARC->n_atom; i++) {
      z = *(pSPARC->atom_pos + 2 + 3 * i);
      if (z < 0 || z > pSPARC->range_z) {
        printf("\nERROR: position of atom # %d is out of the domain!\n", i + 1);
        exit(EXIT_FAILURE);
      }
    }
  }

#ifdef DEBUG
  if (!rank) {
    printf("\nRange:\n%12.6f\t%12.6f\t%12.6f\n", pSPARC->range_x, pSPARC->range_y, pSPARC->range_z);
    printf("\nCOORD AFTER MAPPING:\n");
    for (i = 0; i < 3 * pSPARC->n_atom; i++) {
      printf("%12.6f\t", pSPARC->atom_pos[i]);
      if (i % 3 == 2 && i > 0)
        printf("\n");
    }
    printf("\n");
  }
#endif
  /* find finite difference weights for first & second derivatives */
  pSPARC->FDweights_D1 = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->FDweights_D2 = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D1_stencil_coeffs_x = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D1_stencil_coeffs_y = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D1_stencil_coeffs_z = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D2_stencil_coeffs_x = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D2_stencil_coeffs_y = (double *)malloc((FDn + 1) * sizeof(double));
  pSPARC->D2_stencil_coeffs_z = (double *)malloc((FDn + 1) * sizeof(double));
  if (pSPARC->FDweights_D1 == NULL || pSPARC->FDweights_D1 == NULL || pSPARC->D1_stencil_coeffs_x == NULL ||
      pSPARC->D1_stencil_coeffs_y == NULL || pSPARC->D1_stencil_coeffs_z == NULL ||
      pSPARC->D2_stencil_coeffs_x == NULL || pSPARC->D2_stencil_coeffs_y == NULL ||
      pSPARC->D2_stencil_coeffs_z == NULL) {
    printf("\nmemory cannot be allocated6\n");
    exit(EXIT_FAILURE);
  }

  // 1st derivative weights
  pSPARC->FDweights_D1[0] = 0;
  for (p = 1; p < FDn + 1; p++) {
    pSPARC->FDweights_D1[p] = (2 * (p % 2) - 1) * fract(FDn, p) / p;
  }
  // 2nd derivative weights
  pSPARC->FDweights_D2[0] = 0;
  for (p = 1; p < FDn + 1; p++) {
    pSPARC->FDweights_D2[0] -= (2.0 / (p * p));
    pSPARC->FDweights_D2[p] = (2 * (p % 2) - 1) * 2 * fract(FDn, p) / (p * p);
  }

  // 1st derivative weights including mesh
  double dx_inv, dy_inv, dz_inv;
  dx_inv = 1.0 / pSPARC->delta_x;
  dy_inv = 1.0 / pSPARC->delta_y;
  dz_inv = 1.0 / pSPARC->delta_z;
  for (p = 1; p < FDn + 1; p++) {
    pSPARC->D1_stencil_coeffs_x[p] = pSPARC->FDweights_D1[p] * dx_inv;
    pSPARC->D1_stencil_coeffs_y[p] = pSPARC->FDweights_D1[p] * dy_inv;
    pSPARC->D1_stencil_coeffs_z[p] = pSPARC->FDweights_D1[p] * dz_inv;
  }

  // 2nd derivative weights including mesh
  double dx2_inv, dy2_inv, dz2_inv;
  dx2_inv = 1.0 / (pSPARC->delta_x * pSPARC->delta_x);
  dy2_inv = 1.0 / (pSPARC->delta_y * pSPARC->delta_y);
  dz2_inv = 1.0 / (pSPARC->delta_z * pSPARC->delta_z);

  // Stencil coefficients for mixed derivatives
  if (pSPARC->cell_typ == 0) {
    for (p = 0; p < FDn + 1; p++) {
      pSPARC->D2_stencil_coeffs_x[p] = pSPARC->FDweights_D2[p] * dx2_inv;
      pSPARC->D2_stencil_coeffs_y[p] = pSPARC->FDweights_D2[p] * dy2_inv;
      pSPARC->D2_stencil_coeffs_z[p] = pSPARC->FDweights_D2[p] * dz2_inv;
    }
  } else if (pSPARC->cell_typ > 10 && pSPARC->cell_typ < 20) {
    pSPARC->D2_stencil_coeffs_xy = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D2_stencil_coeffs_yz = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D2_stencil_coeffs_xz = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_xy = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_yx = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_xz = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_zx = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_yz = (double *)malloc((FDn + 1) * sizeof(double));
    pSPARC->D1_stencil_coeffs_zy = (double *)malloc((FDn + 1) * sizeof(double));
    if (pSPARC->D2_stencil_coeffs_xy == NULL || pSPARC->D2_stencil_coeffs_yz == NULL ||
        pSPARC->D2_stencil_coeffs_xz == NULL || pSPARC->D1_stencil_coeffs_xy == NULL ||
        pSPARC->D1_stencil_coeffs_yx == NULL || pSPARC->D1_stencil_coeffs_xz == NULL ||
        pSPARC->D1_stencil_coeffs_zx == NULL || pSPARC->D1_stencil_coeffs_yz == NULL ||
        pSPARC->D1_stencil_coeffs_zy == NULL) {
      printf("\nmemory cannot be allocated\n");
      exit(EXIT_FAILURE);
    }
    for (p = 0; p < FDn + 1; p++) {
      pSPARC->D2_stencil_coeffs_x[p] = pSPARC->lapcT[0] * pSPARC->FDweights_D2[p] * dx2_inv;
      pSPARC->D2_stencil_coeffs_y[p] = pSPARC->lapcT[4] * pSPARC->FDweights_D2[p] * dy2_inv;
      pSPARC->D2_stencil_coeffs_z[p] = pSPARC->lapcT[8] * pSPARC->FDweights_D2[p] * dz2_inv;
      pSPARC->D2_stencil_coeffs_xy[p] = 2 * pSPARC->lapcT[1] * pSPARC->FDweights_D1[p] * dx_inv; // 2*T_12 d/dx(df/dy)
      pSPARC->D2_stencil_coeffs_xz[p] = 2 * pSPARC->lapcT[2] * pSPARC->FDweights_D1[p] * dx_inv; // 2*T_13 d/dx(df/dz)
      pSPARC->D2_stencil_coeffs_yz[p] = 2 * pSPARC->lapcT[5] * pSPARC->FDweights_D1[p] * dy_inv; // 2*T_23 d/dy(df/dz)
      pSPARC->D1_stencil_coeffs_xy[p] =
          2 * pSPARC->lapcT[1] * pSPARC->FDweights_D1[p] * dy_inv; // d/dx(2*T_12 df/dy) used in d/dx(2*T_12 df/dy +
                                                                   // 2*T_13 df/dz)
      pSPARC->D1_stencil_coeffs_yx[p] =
          2 * pSPARC->lapcT[1] * pSPARC->FDweights_D1[p] * dx_inv; // d/dy(2*T_12 df/dx) used in d/dy(2*T_12 df/dx +
                                                                   // 2*T_23 df/dz)
      pSPARC->D1_stencil_coeffs_xz[p] =
          2 * pSPARC->lapcT[2] * pSPARC->FDweights_D1[p] * dz_inv; // d/dx(2*T_13 df/dz) used in d/dx(2*T_12 df/dy +
                                                                   // 2*T_13 df/dz)
      pSPARC->D1_stencil_coeffs_zx[p] =
          2 * pSPARC->lapcT[2] * pSPARC->FDweights_D1[p] * dx_inv; // d/dz(2*T_13 df/dx) used in d/dz(2*T_13 df/dz +
                                                                   // 2*T_23 df/dy)
      pSPARC->D1_stencil_coeffs_yz[p] =
          2 * pSPARC->lapcT[5] * pSPARC->FDweights_D1[p] * dz_inv; // d/dy(2*T_23 df/dz) used in d/dy(2*T_12 df/dx +
                                                                   // 2*T_23 df/dz)
      pSPARC->D1_stencil_coeffs_zy[p] =
          2 * pSPARC->lapcT[5] * pSPARC->FDweights_D1[p] * dy_inv; // d/dz(2*T_23 df/dy) used in d/dz(2*T_12 df/dx +
                                                                   // 2*T_23 df/dy)
    }
  }

#ifdef DEBUG
  t1 = MPI_Wtime();
#endif
  pSPARC->MaxEigVal_mhalfLap = 0.0; // initialize to 0 to avoid accessing uninitialized variable
  if (pSPARC->cell_typ == 0) {
    // maximum eigenvalue of -0.5 * Lap (only accurate with periodic boundary
    // conditions)
    pSPARC->MaxEigVal_mhalfLap =
        pSPARC->D2_stencil_coeffs_x[0] + pSPARC->D2_stencil_coeffs_y[0] + pSPARC->D2_stencil_coeffs_z[0];
    double scal_x, scal_y, scal_z;
    scal_x = (pSPARC->Nx - pSPARC->Nx % 2) / (double)pSPARC->Nx;
    scal_y = (pSPARC->Ny - pSPARC->Ny % 2) / (double)pSPARC->Ny;
    scal_z = (pSPARC->Nz - pSPARC->Nz % 2) / (double)pSPARC->Nz;
    for (p = 1; p < FDn + 1; p++) {
      pSPARC->MaxEigVal_mhalfLap += 2.0 * (pSPARC->D2_stencil_coeffs_x[p] * cos(M_PI * p * scal_x) +
                                           pSPARC->D2_stencil_coeffs_y[p] * cos(M_PI * p * scal_y) +
                                           pSPARC->D2_stencil_coeffs_z[p] * cos(M_PI * p * scal_z));
    }
    pSPARC->MaxEigVal_mhalfLap *= -0.5;
  } else if (pSPARC->cell_typ > 10 && pSPARC->cell_typ < 20) {
    // WARNING: for non-orthogonal cells, this is only a rough estimate, which
    // gives a lowerbound of the accurate answer maximum eigenvalue of -0.5 *
    // Lap (non-orthogonal lattice)
    pSPARC->MaxEigVal_mhalfLap = pSPARC->D2_stencil_coeffs_x[0] // note lapcT[0] is already multiplied in
                                                                // D2_stencil_coeffs_x above
                                 + pSPARC->D2_stencil_coeffs_y[0] + pSPARC->D2_stencil_coeffs_z[0];
    double scal_x, scal_y, scal_z;
    scal_x = (pSPARC->Nx - pSPARC->Nx % 2) / (double)pSPARC->Nx;
    scal_y = (pSPARC->Ny - pSPARC->Ny % 2) / (double)pSPARC->Ny;
    scal_z = (pSPARC->Nz - pSPARC->Nz % 2) / (double)pSPARC->Nz;
    for (p = 1; p < FDn + 1; p++) {
      pSPARC->MaxEigVal_mhalfLap += 2.0 * pSPARC->D2_stencil_coeffs_x[p] * cos(M_PI * p * scal_x) +
                                    2.0 * pSPARC->D2_stencil_coeffs_y[p] * cos(M_PI * p * scal_y) +
                                    2.0 * pSPARC->D2_stencil_coeffs_z[p] * cos(M_PI * p * scal_z);
    }
    // for mixed terms (actually it's a better approx. if we neglect these
    // terms!)
    double sx, sy, sz;
    sx = sy = sz = 0.0;
    for (p = 1; p < FDn + 1; p++) {
      sx += 2.0 * pSPARC->D1_stencil_coeffs_x[p] * sin(M_PI * p * scal_x); // very close to 0 (exactly 0 for even Nx)
      sy += 2.0 * pSPARC->D1_stencil_coeffs_y[p] * sin(M_PI * p * scal_y); // very close to 0 (exactly 0 for even Ny)
      sz += 2.0 * pSPARC->D1_stencil_coeffs_z[p] * sin(M_PI * p * scal_z); // very close to 0 (exactly 0 for even Nz)
    }
    sx = sy = sz = 0.0;                                             // forcing the mixed terms to be zero here!
    pSPARC->MaxEigVal_mhalfLap += 2.0 * pSPARC->lapcT[1] * sx * sy; // x,y
    pSPARC->MaxEigVal_mhalfLap += 2.0 * pSPARC->lapcT[2] * sx * sz; // x,z
    pSPARC->MaxEigVal_mhalfLap += 2.0 * pSPARC->lapcT[5] * sy * sz; // y,z
    pSPARC->MaxEigVal_mhalfLap *= -0.5;
  }
#ifdef DEBUG
  t2 = MPI_Wtime();
  if (!rank)
    printf("Max eigenvalue of -0.5*Lap is %.13f, time taken: %.3f ms\n", pSPARC->MaxEigVal_mhalfLap, (t2 - t1) * 1e3);
#endif

  // find Chebyshev polynomial degree based on max eigenvalue (spectral width)
  if (pSPARC->ChebDegree < 0) {
    double h_eff = 0.0;
    if (pSPARC->cell_typ == 0) {
      if (fabs(pSPARC->delta_x - pSPARC->delta_y) < 1E-12 && fabs(pSPARC->delta_y - pSPARC->delta_z) < 1E-12) {
        h_eff = pSPARC->delta_x;
      } else {
        // find effective mesh s.t. it has same spectral width
        h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
      }
    } else if (pSPARC->cell_typ > 10 && pSPARC->cell_typ < 20) {
      // use max eigenvalue of -1/2*Lap to estimate the effective mesh size for
      // orthogonal Laplacian
      const double lambda_ref = 6.8761754299116333; // max eigval for 1D orthogonal -1.0*Lap for h_eff
                                                    // = 1.0
      h_eff = sqrt(3.0 * lambda_ref / (2.0 * pSPARC->MaxEigVal_mhalfLap));
    }
    pSPARC->ChebDegree = Mesh2ChebDegree(h_eff);
#ifdef DEBUG
    if (!rank && h_eff < 0.1) {
      printf("#WARNING: for mesh less than 0.1, the default Chebyshev "
             "polynomial degree might not be enought!\n");
    }
    if (!rank)
      printf("h_eff = %.2f, npl = %d\n", h_eff, pSPARC->ChebDegree);
#endif
  }
#ifdef DEBUG
  else {
    if (!rank)
      printf("Chebyshev polynomial degree (provided by user): npl = %d\n", pSPARC->ChebDegree);
  }
#endif

  // set default simple (linear) mixing parameter to be the same as for pulay
  // mixing
  if (pSPARC->MixingParameterSimple < 0.0) {
    pSPARC->MixingParameterSimple = pSPARC->MixingParameter;
  }

  // set default mixing parameter for magnetization density to the same as
  // mixing parameter for total density/potential
  if (pSPARC->MixingParameterMag < 0.0) {
    pSPARC->MixingParameterMag = pSPARC->MixingParameter;
  }

  // set default simple (linear) mixing parameter for magnetization density to
  // be the same as for pulay mixing
  if (pSPARC->MixingParameterSimpleMag < 0.0) {
    pSPARC->MixingParameterSimpleMag = pSPARC->MixingParameterMag;
  }

  if (pSPARC->MixingVariable < 0) { // mixing variable not provided
    pSPARC->MixingVariable = 0;     // set default to 'density'
  }

  if (pSPARC->MixingPrecond < 0) { // mixing preconditioner not provided
    pSPARC->MixingPrecond = 1;     // set default to 'Kerker' preconditioner
  }

  if (pSPARC->MixingPrecondMag < 0) { // mixing preconditioner for magnetization not provided
    pSPARC->MixingPrecondMag = 0;     // set default to 'none'
  }

  // set up real space preconditioner coefficients
  if (pSPARC->MixingPrecond == 2) { // coeff for Resta preconditioner
    pSPARC->precondcoeff_n = 1;
    pSPARC->precondcoeff_a = (double _Complex *)malloc(pSPARC->precondcoeff_n * sizeof(double _Complex));
    pSPARC->precondcoeff_lambda_sqr = (double _Complex *)malloc(pSPARC->precondcoeff_n * sizeof(double _Complex));
    pSPARC->precondcoeff_a[0] = 0.820368124615603 - 0.330521220406859 * I;
    pSPARC->precondcoeff_lambda_sqr[0] = 1.539184309857566 + 1.454446707757472 * I;
    pSPARC->precondcoeff_k = 0.179473117041025;

    // check env for what system we're running to set the coeffs,
    // if not provided, the default will be used
    // TODO: remove the following after check!
    char system[L_STRING] = "default";
    char *precond_system = getenv("PRECOND_SYSTEM");
    if (precond_system != NULL)
      snprintf(system, L_STRING, "%s", precond_system);

    if (strcmpi(system, "MoS2") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for MoS2 system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif

        pSPARC->precondcoeff_n = 1;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 0.897326075519806 - 0.837703986538753 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.328766315380339 + 0.183508748834006 * I;
        pSPARC->precondcoeff_k = 0.102576229227011;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for MoS2 system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 3;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 0.797410061005122 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = -0.000029265620523 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[2] = 0.103239343777798 - 0.003381206211038 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.601186842883198 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = -0.256060441488722 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[2] = -0.104178068950438 + 0.493948292725977 * I;
        pSPARC->precondcoeff_k = 0.099398800263940;
      }
    } else if (strcmpi(system, "Si") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for Si system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 1;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 0.914678024418436 - 1.055347015597097 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.238671535971552 + 0.106323808659314 * I;
        pSPARC->precondcoeff_k = 0.085289070702772;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for Si system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 3;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = -0.000124974499632 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = 0.822613437367865 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[2] = 0.094666235811611 - 0.004627781592542 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = -1.072175758908308 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.420975552998538 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[2] = -0.054999300909744 + 0.349588273989346 * I;
        pSPARC->precondcoeff_k = 0.082856817316465;
      }
    } else if (strcmpi(system, "C") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for C system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 1;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 0.8206 - 0.3427 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.4284 + 0.4019 * I;
        pSPARC->precondcoeff_k = 0.1793;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "WARNING: coeffs for C system with fitpow %d are not "
                     "set!\n" RESET,
                 pSPARC->precond_fitpow);
#endif
        // pSPARC->precondcoeff_n = 3;
        // // reallocate memory
        // pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a,
        // pSPARC->precondcoeff_n * sizeof(double _Complex));
        // pSPARC->precondcoeff_lambda_sqr =
        // realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n *
        // sizeof(double _Complex)); pSPARC->precondcoeff_a[0]          =
        // -0.000124974499632 + 0.000000000000000 * I; pSPARC->precondcoeff_a[1]
        // = 0.822613437367865 + 0.000000000000000 * I;
        // pSPARC->precondcoeff_a[2]          = 0.094666235811611 -
        // 0.004627781592542 * I; pSPARC->precondcoeff_lambda_sqr[0] =
        // -1.072175758908308 + 0.000000000000000 * I;
        // pSPARC->precondcoeff_lambda_sqr[1] = 0.420975552998538 +
        // 0.000000000000000 * I; pSPARC->precondcoeff_lambda_sqr[2] =
        // -0.054999300909744 + 0.349588273989346 * I; pSPARC->precondcoeff_k =
        // 0.082856817316465;
      }
    }
  } else if (pSPARC->MixingPrecond == 3) { // coeff for truncated Kerker preconditioner
    pSPARC->precondcoeff_n = 1;
    pSPARC->precondcoeff_a = (double _Complex *)malloc(pSPARC->precondcoeff_n * sizeof(double _Complex));
    pSPARC->precondcoeff_lambda_sqr = (double _Complex *)malloc(pSPARC->precondcoeff_n * sizeof(double _Complex));
    pSPARC->precondcoeff_a[0] = 0.740197283447608 - 2.187940485005530 * I;
    pSPARC->precondcoeff_lambda_sqr[0] = 0.515764278984552 + 0.261718938132583 * I;
    pSPARC->precondcoeff_k = 0.259680843800232;

    // check env for what system we're running to set the coeffs,
    // it not provided, the default will be used
    // TODO: remove the following after check!
    char system[L_STRING] = "default";
    char *precond_system = getenv("PRECOND_SYSTEM");
    if (precond_system != NULL)
      snprintf(system, L_STRING, "%s", precond_system);

    if (strcmpi(system, "MoS2") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for MoS2 system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 2;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 1.069131757115932 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = -0.171827850593795 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.261519729188790 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.024058288033320 + 0.000000000000000 * I;
        pSPARC->precondcoeff_k = 0.102669136088733;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for MoS2 system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 3;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 0.000011385765477 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = 0.994255001880647 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[2] = -0.093994967542657 - 0.006240439304379 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = -0.580143676837624 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.281390031341584 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[2] = -0.005192385910338 + 0.009670637051448 * I;
        pSPARC->precondcoeff_k = 0.099729735832187;
      }
    } else if (strcmpi(system, "Si") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for Si system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 2;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 1.045423322787217 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = -0.130145326907590 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.267115428215830 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.019530203373891 + 0.000000000000000 * I;
        pSPARC->precondcoeff_k = 0.084702403406033;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for Si system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 3;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = -0.000450002447564 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[1] = 0.991616958994114 + 0.000000000000000 * I;
        pSPARC->precondcoeff_a[2] = -0.074468796694241 - 0.014060128507695 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 3.578501584073372 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.283063390321347 + 0.000000000000000 * I;
        pSPARC->precondcoeff_lambda_sqr[2] = -0.004905277505535 + 0.011599970024290 * I;
        pSPARC->precondcoeff_k = 0.083301273707655;
      }
    } else if (strcmpi(system, "C") == 0) {
      if (pSPARC->precond_fitpow == 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "Using coeffs for C system with fitpow %d\n" RESET, pSPARC->precond_fitpow);
#endif
        pSPARC->precondcoeff_n = 2;
        // reallocate memory
        pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_lambda_sqr =
            realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n * sizeof(double _Complex));
        pSPARC->precondcoeff_a[0] = 1.2926 + 0.0000 * I;
        pSPARC->precondcoeff_a[1] = -0.4780 + 0.0000 * I;
        pSPARC->precondcoeff_lambda_sqr[0] = 0.2310 + 0.0000 * I;
        pSPARC->precondcoeff_lambda_sqr[1] = 0.0552 + 0.0000 * I;
        pSPARC->precondcoeff_k = 0.1854;
      } else if (pSPARC->precond_fitpow == 4) {
#ifdef DEBUG
        if (!rank)
          printf(RED "WARNING: coeffs for C system with fitpow %d are not "
                     "set!\n" RESET,
                 pSPARC->precond_fitpow);
#endif
        // pSPARC->precondcoeff_n = 3;
        // // reallocate memory
        // pSPARC->precondcoeff_a = realloc(pSPARC->precondcoeff_a,
        // pSPARC->precondcoeff_n * sizeof(double _Complex));
        // pSPARC->precondcoeff_lambda_sqr =
        // realloc(pSPARC->precondcoeff_lambda_sqr, pSPARC->precondcoeff_n *
        // sizeof(double _Complex)); pSPARC->precondcoeff_a[0]          =
        // -0.000124974499632 + 0.000000000000000 * I; pSPARC->precondcoeff_a[1]
        // = 0.822613437367865 + 0.000000000000000 * I;
        // pSPARC->precondcoeff_a[2]          = 0.094666235811611 -
        // 0.004627781592542 * I; pSPARC->precondcoeff_lambda_sqr[0] =
        // -1.072175758908308 + 0.000000000000000 * I;
        // pSPARC->precondcoeff_lambda_sqr[1] = 0.420975552998538 +
        // 0.000000000000000 * I; pSPARC->precondcoeff_lambda_sqr[2] =
        // -0.054999300909744 + 0.349588273989346 * I; pSPARC->precondcoeff_k =
        // 0.082856817316465;
      }
    }
  }

  // scf error type, 0 - default, 1 - QE (conv_thr)
  if (pSPARC->scf_err_type != 0 && pSPARC->scf_err_type != 1) {
    if (!rank)
      printf("Cannot recognize SCF error type!\n");
    exit(EXIT_FAILURE);
  }

  // for evaluating QE scf error, we need to perform some extra calculations
  // e.g., an extra Poisson solve, this timer keeps track of the extra time
  // spent
  if (pSPARC->scf_err_type == 1) {
    pSPARC->t_qe_extra = 0.0;
  }

  // default SCF tolerance based on accuracy_level
  // we use a model curve to correlate scf tolerance and energy and force
  // accuracy
  //     log(y) = a * log(x) + b
  // y could be either energy accuracy (Ha/atom) or force accuracy  (Ha/Bohr)
  // if scf tol is not set, we'll use accuracy_level to find scf tol
  // All models satisfy 90% dataset
  if (pSPARC->TOL_SCF < 0.0) {
    if (pSPARC->MDFlag != 0) {
      // in case of MD, using 1E-3 Ha/Bohr force accuracy as target
      double target_force_accuracy = 1E-3;
      const double a = 1.025;
      const double b = 1.368;
      double log_target = log(target_force_accuracy);
      pSPARC->TOL_SCF = exp((log_target - b) / a);
    } else if (pSPARC->RelaxFlag != 0) {
      // in case of relaxation, using TOL_RELAX/5 force accuracy as target
      double target_force_accuracy = pSPARC->TOL_RELAX / 5;
      const double a = 1.025;
      const double b = 1.468;
      double log_target = log(target_force_accuracy);
      pSPARC->TOL_SCF = exp((log_target - b) / a);
    } else {
      // in case of single point calculation, using 1E-5 Ha/atom energy accuracy
      // as target
      double target_force_accuracy = -1.0;
      double target_energy_accuracy = -1.0;

      // accuracy_levels      : 0 - minimal | 1 - low  | 2 - medium | 3 - high |
      // 4 - extreme target force accuracy: 0 - 1e-1    | 1 - 1e-2 | 2 - 1e-3 |
      // 3 - 1e-4 | 4 - 1e-5
      if (pSPARC->accuracy_level >= 0) {
        target_force_accuracy = pow(10.0, -(pSPARC->accuracy_level + 1.0));
      } else if (pSPARC->target_force_accuracy > 0.0) {
        target_force_accuracy = pSPARC->target_force_accuracy;
      } else if (pSPARC->target_energy_accuracy > 0.0) {
        target_energy_accuracy = pSPARC->target_energy_accuracy;
      }

      // if none of the accuracy levels are specified, set energy_accuracy to
      // 1e-5 Ha/atom
      if (target_force_accuracy < 0 && target_energy_accuracy < 0) {
        target_energy_accuracy = 1e-5;
      }

      // calculate SCF TOL based on specified target accuracy
      if (target_energy_accuracy > 0.0) { // find scf tol based on target energy accuracy
        const double a = 1.502;
        const double b = 1.165;
        double log_target = log(target_energy_accuracy);
        pSPARC->TOL_SCF = exp((log_target - b) / a);
      } else if (target_force_accuracy > 0.0) { // find scf tol based on target force accuracy
        const double a = 1.025;
        const double b = 1.368;
        double log_target = log(target_force_accuracy);
        pSPARC->TOL_SCF = exp((log_target - b) / a);
      }
    }
  }

  // default Kerker tolerance
  if (pSPARC->TOL_PRECOND < 0.0) { // kerker tol not provided by user
    double h_eff = 0.0;
    if (fabs(pSPARC->delta_x - pSPARC->delta_y) < 1E-12 && fabs(pSPARC->delta_y - pSPARC->delta_z) < 1E-12) {
      h_eff = pSPARC->delta_x;
    } else {
      // find effective mesh s.t. it has same spectral width
      h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
    }
    pSPARC->TOL_PRECOND = (h_eff * h_eff) * 1e-3;
  }

  // default poisson tolerance
  // In most cases this is enough. The error in energy/atom due to poisson
  // error is usually 1~2 orders of magnitude smaller than the poisson tol.
  // And the error in forces due to poisson error is of 1 order of magnitude
  // larger than poisson tol.
  // Another issue about using moderate poisson tolerance is SCF convergence.
  // we found that usually keeping poisson tol as 1 order of magnitude lower
  // than scf_tol, the SCF convergence won't be affected. However, if SCF
  // stagers near convergence, consider reducing poisson tolerance further
  // (e.g. scf_tol*0.01).
  if (pSPARC->TOL_POISSON < 0.0) { // if poisson tol not provided by user
    pSPARC->TOL_POISSON = pSPARC->TOL_SCF * 0.01;
  }

  // default rb tolerance
  // Note that rb tolerance is the absolute tolerance.
  // The error in energy/atom due to rb error is usually 1~2 orders of
  // magnitude smaller than rb tol. While the error in force due to rb
  // error is roughly the same order as rb tol.
  if (pSPARC->TOL_PSEUDOCHARGE < 0.0) { // if rb tol not provided by user
    pSPARC->TOL_PSEUDOCHARGE = pSPARC->TOL_SCF * 0.001;
  }

  // The following will override the user-provided FixRandSeed
  // check env for FixRandSeed variable: 0 - off, 1 - on
  int FixRandSeed = pSPARC->FixRandSeed;
  char *env_var = getenv("FIX_RAND");
  if (env_var != NULL)
    FixRandSeed = atoi(env_var);
  pSPARC->FixRandSeed = FixRandSeed;

  // allocate memory for pseudocharge cutoff radius
  pSPARC->CUTOFF_x = (double *)malloc(Ntypes * sizeof(double));
  pSPARC->CUTOFF_y = (double *)malloc(Ntypes * sizeof(double));
  pSPARC->CUTOFF_z = (double *)malloc(Ntypes * sizeof(double));

  if (pSPARC->CUTOFF_x == NULL || pSPARC->CUTOFF_y == NULL || pSPARC->CUTOFF_z == NULL) {
    printf("\nmemory cannot be allocated7\n");
    exit(EXIT_FAILURE);
  }

  // number of k-points after symmetry reduction (currently only
  // takes into account the time reversal symmetry)
  // WARNING: Time-reversal symmetry only if there is no magnetic field applied
  // this won't work with kpt shift
  // pSPARC->Nkpts_sym = ceil(pSPARC->Kx*pSPARC->Ky*pSPARC->Kz/2.0);
  pSPARC->Nkpts_sym = pSPARC->Nkpts; // will be updated after actual symmetry reduction

  // at this point symmetry reduction is not done yet
  pSPARC->kptWts = (double *)malloc(pSPARC->Nkpts_sym * sizeof(double));
  pSPARC->k1 = (double *)malloc(pSPARC->Nkpts_sym * sizeof(double));
  pSPARC->k2 = (double *)malloc(pSPARC->Nkpts_sym * sizeof(double));
  pSPARC->k3 = (double *)malloc(pSPARC->Nkpts_sym * sizeof(double));
  // calculate the k point and weights (shift in kpt may apply)
  Calculate_kpoints(pSPARC);

  // flag to indicate if it is a gamma-point calculation
  pSPARC->isGammaPoint =
      (int)(pSPARC->Nkpts_sym == 1 && fabs(pSPARC->k1[0]) < TEMP_TOL && fabs(pSPARC->k2[0]) < TEMP_TOL &&
            fabs(pSPARC->k3[0]) < TEMP_TOL && pSPARC->SOC_Flag == 0);

  if (pSPARC->vdWDFFlag != 0) {
    if ((pSPARC->BCx) || (pSPARC->BCy) || (pSPARC->BCz)) {
      if (rank == 0)
        printf(RED "ERROR: vdW-DF does not support Dirichlet boundary "
                   "condition!\n" RESET);
      exit(EXIT_FAILURE);
    }
  }

#if !defined(USE_MKL) && !defined(USE_FFTW)
  if (pSPARC->vdWDFFlag != 0) {
    if (rank == 0)
      printf(RED "ERROR: To use vdW-DF, please turn on MKL or FFTW in makefile!\n"
                 "Or you can stop using vdW-DF by setting other "
                 "exchange-correlation functionals.\n" RESET);
    exit(EXIT_FAILURE);
  }
  if (pSPARC->usefock == 1) {
    if (rank == 0)
      printf(RED "ERROR: To use hybrid functionals like PBE0 or HF, please "
                 "turn on MKL or FFTW in makefile!\n" RESET);
    exit(EXIT_FAILURE);
  }
#endif // #if !defined(USE_MKL) && !defined(USE_FFTW)

  if (pSPARC->mGGAflag == 1) {
    // if (pSPARC->spin_typ != 0) {
    //     if (rank == 0)
    //         printf(RED "ERROR: currently SCAN does not support spin
    //         polarization!\n" RESET);
    //     exit(EXIT_FAILURE);
    // }
    if (pSPARC->SOC_Flag || pSPARC->usefock || pSPARC->SQFlag) {
      if (!rank)
        printf(RED "ERROR: Spin-orbit coupling, hybrid and SQ are not "
                   "supported in this version of SCAN implementation.\n" RESET);
      exit(EXIT_FAILURE);
    }
  }

#if !defined(USE_MKL) && !defined(USE_SCALAPACK)
  if (pSPARC->usefock == 1 && pSPARC->ACEFlag == 1) {
    if (rank == 0)
      printf(RED "ERROR: To use hybrid functional with ACE method, please turn "
                 "on MKL or SCALAPACK in makefile!\n" RESET);
    exit(EXIT_FAILURE);
  }
#endif // #if defined(USE_MKL) || defined(USE_SCALAPACK)

  if (pSPARC->usefock == 1) {
    if (pSPARC->SOC_Flag || pSPARC->SQFlag) {
      if (!rank)
        printf(RED "ERROR: Spin-orbit coupling and SQ are not supported in "
                   "this version of hybrid implementation.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC->EXXDiv_Flag < 0) {
      if (pSPARC->BC > 2) {
#ifdef DEBUG
        if (!rank)
          printf(RED "For Wire and Slab with hybrid funcitonal, deafults to "
                     "use auxiliary functioin method.\n" RESET);
#endif
        pSPARC->EXXDiv_Flag = 1;
      } else {
        if (strcmpi(pSPARC->XC, "HSE") == 0) {
#ifdef DEBUG
          if (!rank)
            printf(RED "For Bulk and Cluster with HSE hybrid funcitonal, "
                       "deafults to use ERFC method.\n" RESET);
#endif
          pSPARC->EXXDiv_Flag = 2;
        } else {
#ifdef DEBUG
          if (!rank)
            printf(RED "For Bulk and Cluster with hybrid funcitonal, deafults "
                       "to use spherical truncation method.\n" RESET);
#endif
          pSPARC->EXXDiv_Flag = 0;
        }
      }
    } else {
      if (strcmpi(pSPARC->XC, "HSE") != 0 && pSPARC->EXXDiv_Flag == 2) {
        printf(RED "ERROR: ERFC method could only be used with HSE "
                   "functional.\n" RESET);
        exit(EXIT_FAILURE);
      }
    }

    if (pSPARC->TOL_FOCK < 0.0) {
      // TODO: This model is based on the results of 5 tests. Do more tests to
      // improve the robustness. default FOCK outer loop tolerance based on
      // accuracy_level we use a model curve to correlate fock tolerance and
      // energy and force accuracy
      //     log10(y) = a * log10(x) + b, a = 1, b = 0.3
      // when related to TOL_SCF, this could be simplified as 0.2 * TOL_SCF for
      // the same force accuracy.
      pSPARC->TOL_FOCK = 0.2 * pSPARC->TOL_SCF;
    }

    // If initial PBE SCF tolerance is not defined, use default
    // 10*pSPARC->TOL_FOCK
    if (pSPARC->TOL_SCF_INIT < 0.0) {
      pSPARC->TOL_SCF_INIT = max(10 * pSPARC->TOL_FOCK, 1e-3);
    }
    pSPARC->MAXIT_FOCK = max(1, pSPARC->MAXIT_FOCK);
    pSPARC->MINIT_FOCK = max(1, pSPARC->MINIT_FOCK);

    if (pSPARC->EXXMem_batch < 0) {
      // use default EXXMem_batch if it's negative
      pSPARC->EXXMem_batch = 0;
    }

    // If using ACE operator, only do domain parallelization
    if (pSPARC->ACEFlag == 1) {
      if (pSPARC->EXXACEVal_state < 0) {
        // Use default EXXACEVal_state if it's negative
        pSPARC->EXXACEVal_state = 3;
      }
    } else {
      pSPARC->EXXACEVal_state = 0;
    }
  } else {
    pSPARC->ACEFlag = 0;
    pSPARC->EXXMem_batch = 0;
    pSPARC->EXXACEVal_state = 0;
  }

  // constraints on SOC
  if (pSPARC->SOC_Flag == 1) {
    if (pSPARC->usefock || pSPARC->mGGAflag || pSPARC->SQFlag) {
      if (rank == 0)
        printf(RED "ERROR: Hybrid functional, SCAN and SQ are not supported in "
                   "this version of spin-orbit coupling implementation.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC->spin_typ == 1) {
      if (rank == 0)
        printf(RED "ERROR: Spin-polarized calculation is not supported in this "
                   "version of spin-orbit coupling implementation.\n" RESET);
      exit(EXIT_FAILURE);
    }
  }

  // constraints on SQ
  if (pSPARC->SQFlag == 1) {
    if (pSPARC->BCx || pSPARC->BCy || pSPARC->BCz) {
      if (!rank)
        printf(RED "ERROR: SQ method only supports periodic boundary "
                   "conditions.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC->cell_typ != 0) {
      if (!rank)
        printf(RED "ERROR: SQ method only supports orthogonal systems in this "
                   "version.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC->isGammaPoint != 1 || pSPARC->spin_typ == 1) {
      if (rank == 0)
        printf(RED "ERROR: Polarized calculation and Kpoint options are not "
                   "supported in this version of SQ implementation.\n" RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC_Input->Nstates != -1) {
      if (rank == 0)
        printf(RED "ERROR: NSTATES is not vaild in SQ method.\n" RESET);
      exit(EXIT_FAILURE);
    }

    if (pSPARC->SOC_Flag || pSPARC->usefock || pSPARC->mGGAflag) {
      if (!rank)
        printf(RED "ERROR: Hybrid functional, spin-orbit coupling, and SCAN are "
                   "not supported in this version of SQ implementation." RESET);
      exit(EXIT_FAILURE);
    }
    if (pSPARC->SQ_rcut < 0) {
      if (!rank)
        printf(RED "ERROR: SQ_RCUT must be provided when SQ method is turned "
                   "on.\n" RESET);
      exit(EXIT_FAILURE);
    }

    if (pSPARC->SQ_npl_g <= 0) {
      if (!rank)
        printf(RED "ERROR: SQ_NPL_G must be provided a positive integer when Gauss "
                   "Quadrature method is turned on in SQ method.\n" RESET);
      exit(EXIT_FAILURE);
    }

    if (pSPARC->PrintEigenFlag > 0) {
      if (!rank)
        printf(RED "ERROR: PRINT_EIGEN is not valid in SQ method.\n" RESET);
      exit(EXIT_FAILURE);
    }
    pSPARC->SQ_correction = 0; // The correction term in energy and forces
                               // hasn't been implemented in this version.
  }

  if (pSPARC->PrintPsiFlag[0] == 1 && pSPARC->PrintPsiFlag[1] < 0) {
    pSPARC->PrintPsiFlag[1] = 0;
    pSPARC->PrintPsiFlag[2] = pSPARC->Nspin - 1; // spin start/end index
    pSPARC->PrintPsiFlag[3] = 0;
    pSPARC->PrintPsiFlag[4] = pSPARC->Nkpts - 1; // k-point start/end index
    pSPARC->PrintPsiFlag[5] = 0;
    pSPARC->PrintPsiFlag[6] = pSPARC->Nstates - 1; // band start/end index
  }
}
