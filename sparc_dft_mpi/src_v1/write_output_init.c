#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "initialization.h"


/**
 * @brief   Write the initialized parameters into the output file.
 */
void write_output_init(SPARC_OBJ *pSPARC) {
  int i, j, nproc, count;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // time_t current_time = time(NULL);
  // char *c_time_str = ctime(&current_time);
  time_t current_time;
  time(&current_time);
  char *c_time_str = ctime(&current_time);
  // ctime includes a newline char '\n', remove manually
  if (c_time_str[strlen(c_time_str) - 1] == '\n')
    c_time_str[strlen(c_time_str) - 1] = '\0';

  FILE *output_fp = fopen(pSPARC->OutFilename, "w");
  if (output_fp == NULL) {
    printf("\nCannot open file \"%s\"\n", pSPARC->OutFilename);
    exit(EXIT_FAILURE);
  }

  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "*                       SPARC (version May 11, 2023)     "
                     "                 *\n");
  fprintf(output_fp, "*   Copyright (c) 2020 Material Physics & Mechanics "
                     "Group, Georgia Tech   *\n");
  fprintf(output_fp, "*           Distributed under GNU General Public License "
                     "3 (GPL)          *\n");
  fprintf(output_fp, "*                   Start time: %s                  *\n", c_time_str);
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "                           Input parameters              "
                     "                  \n");
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  if (pSPARC->Flag_latvec_scale == 0) {
    fprintf(output_fp, "CELL: %.15g %.15g %.15g \n", pSPARC->range_x, pSPARC->range_y, pSPARC->range_z);
  } else {
    fprintf(output_fp, "LATVEC_SCALE: %.15g %.15g %.15g \n", pSPARC->latvec_scale_x, pSPARC->latvec_scale_y,
            pSPARC->latvec_scale_z);
  }
  fprintf(output_fp, "LATVEC:\n");
  fprintf(output_fp, "%.15g %.15g %.15g \n", pSPARC->LatVec[0], pSPARC->LatVec[1], pSPARC->LatVec[2]);
  fprintf(output_fp, "%.15g %.15g %.15g \n", pSPARC->LatVec[3], pSPARC->LatVec[4], pSPARC->LatVec[5]);
  fprintf(output_fp, "%.15g %.15g %.15g \n", pSPARC->LatVec[6], pSPARC->LatVec[7], pSPARC->LatVec[8]);
  fprintf(output_fp, "FD_GRID: %d %d %d\n", pSPARC->numIntervals_x, pSPARC->numIntervals_y, pSPARC->numIntervals_z);
  fprintf(output_fp, "FD_ORDER: %d\n", pSPARC->order);
  fprintf(output_fp, "BC:");
  fprintf(output_fp, " %s", pSPARC->BCx == 0 ? "P" : "D");
  fprintf(output_fp, " %s", pSPARC->BCy == 0 ? "P" : "D");
  fprintf(output_fp, " %s", pSPARC->BCz == 0 ? "P" : "D");
  fprintf(output_fp, "\n");
  if (pSPARC->BC > 1 && !pSPARC->SQFlag) {
    fprintf(output_fp, "KPOINT_GRID: %d %d %d\n", pSPARC->Kx, pSPARC->Ky, pSPARC->Kz);
    fprintf(output_fp, "KPOINT_SHIFT: %.15g %.15g %.15g\n", pSPARC->kptshift[0], pSPARC->kptshift[1],
            pSPARC->kptshift[2]);
  }
  fprintf(output_fp, "SPIN_TYP: %d\n", pSPARC->spin_typ);
  if (pSPARC->elec_T_type == 0)
    fprintf(output_fp, "ELEC_TEMP_TYPE: Fermi-Dirac\n");
  else if (pSPARC->elec_T_type == 1)
    fprintf(output_fp, "ELEC_TEMP_TYPE: Gaussian\n");

  if (pSPARC->MDFlag == 1) {
    fprintf(output_fp, "ELEC_TEMP: %.15g\n", pSPARC->elec_T);
  } else {
    fprintf(output_fp, "SMEARING: %.10g\n", 1. / pSPARC->Beta);
  }
  fprintf(output_fp, "EXCHANGE_CORRELATION: %s\n", pSPARC->XC);
  if (strcmpi(pSPARC->XC, "HSE") == 0) {
    fprintf(output_fp, "EXX_RANGE_FOCK: %.6f\n", pSPARC->hyb_range_fock);
    fprintf(output_fp, "EXX_RANGE_PBE: %.6f\n", pSPARC->hyb_range_pbe);
  }
  if (pSPARC->SQFlag == 1) {
    fprintf(output_fp, "SQ_FLAG: %d\n", pSPARC->SQFlag);
    fprintf(output_fp, "SQ_RCUT: %.10g\n", pSPARC->SQ_rcut);
    fprintf(output_fp, "SQ_NPL_G: %d\n", pSPARC->SQ_npl_g);
    if (pSPARC->SQ_gauss_mem == 1) {
      fprintf(output_fp, "SQ_GAUSS_MEM: HIGH\n");
    } else {
      fprintf(output_fp, "SQ_GAUSS_MEM: LOW\n");
    }
    fprintf(output_fp, "SQ_TOL_OCC: %.2E\n", pSPARC->SQ_tol_occ);
  } else {
    fprintf(output_fp, "NSTATES: %d\n", pSPARC->Nstates);
    // this should depend on temperature and preconditoner used
    if (pSPARC->Nstates <
        (int)(1.2 * (pSPARC->Nelectron / 2) + 5) * pSPARC->Nspinor) { // with kerker a factor of 1.1 might be needed
      printf("#WARNING: Number of bands may be insufficient for efficient SCF "
             "convergence.\n");
    }
    fprintf(output_fp, "CHEB_DEGREE: %d\n", pSPARC->ChebDegree);
    if (pSPARC->CheFSI_Optmz == 1) {
      fprintf(output_fp, "CHEFSI_OPTMZ: %d\n", pSPARC->CheFSI_Optmz);
    }
    fprintf(output_fp, "CHEFSI_BOUND_FLAG: %d\n", pSPARC->chefsibound_flag);
  }

  if (pSPARC->RelaxFlag >= 1) {
    fprintf(output_fp, "RELAX_FLAG: %d\n", pSPARC->RelaxFlag);
  }

  if (pSPARC->RelaxFlag == 1 || pSPARC->RelaxFlag == 3) {
    fprintf(output_fp, "RELAX_METHOD: %s\n", pSPARC->RelaxMeth);
    fprintf(output_fp, "RELAX_NITER: %d\n", pSPARC->Relax_Niter);
    if (strcmpi(pSPARC->RelaxMeth, "LBFGS") == 0) {
      fprintf(output_fp, "L_HISTORY: %d\n", pSPARC->L_history);
      fprintf(output_fp, "L_FINIT_STP: %.15g\n", pSPARC->L_finit_stp);
      fprintf(output_fp, "L_MAXMOV: %.15g\n", pSPARC->L_maxmov);
      fprintf(output_fp, "L_AUTOSCALE: %d\n", pSPARC->L_autoscale);
      fprintf(output_fp, "L_LINEOPT: %d\n", pSPARC->L_lineopt);
      fprintf(output_fp, "L_ICURV: %.15g\n", pSPARC->L_icurv);
    } else if (strcmpi(pSPARC->RelaxMeth, "NLCG") == 0) {
      fprintf(output_fp, "NLCG_sigma: %.15g\n", pSPARC->NLCG_sigma);
    } else if (strcmpi(pSPARC->RelaxMeth, "FIRE") == 0) {
      fprintf(output_fp, "FIRE_dt: %.15g\n", pSPARC->FIRE_dt);
      fprintf(output_fp, "FIRE_mass: %.15g\n", pSPARC->FIRE_mass);
      fprintf(output_fp, "FIRE_maxmov: %.15g\n", pSPARC->FIRE_maxmov);
    }
  }

  fprintf(output_fp, "CALC_STRESS: %d\n", pSPARC->Calc_stress);
  if (pSPARC->Calc_stress == 0)
    fprintf(output_fp, "CALC_PRES: %d\n", pSPARC->Calc_pres);
  if (pSPARC->MDFlag == 1 || pSPARC->RelaxFlag == 1)
    fprintf(output_fp, "TWTIME: %G\n", pSPARC->TWtime);
  if (pSPARC->MDFlag == 1) {
    fprintf(output_fp, "MD_FLAG: %d\n", pSPARC->MDFlag);
    fprintf(output_fp, "MD_METHOD: %s\n", pSPARC->MDMeth);
    fprintf(output_fp, "MD_TIMESTEP: %.15g\n", pSPARC->MD_dt);
    fprintf(output_fp, "MD_NSTEP: %d\n", pSPARC->MD_Nstep);
    // fprintf(output_fp,"ION_ELEC_EQT: %d\n",pSPARC->ion_elec_eqT);
    fprintf(output_fp, "ION_VEL_DSTR: %d\n", pSPARC->ion_vel_dstr);
    fprintf(output_fp, "ION_VEL_DSTR_RAND: %d\n", pSPARC->ion_vel_dstr_rand);
    fprintf(output_fp, "ION_TEMP: %.15g\n", pSPARC->ion_T);
    if (strcmpi(pSPARC->MDMeth, "NVT_NH") == 0) {
      fprintf(output_fp, "ION_TEMP_END: %.15g\n", pSPARC->thermos_Tf);
      fprintf(output_fp, "QMASS: %.15g\n", pSPARC->qmass);
    }
    if (strcmpi(pSPARC->MDMeth, "NPT_NH") == 0) {
      // fprintf(output_fp,"AMOUNT_THERMO_VARIABLE: %d\n",pSPARC->NPT_NHnnos);
      fprintf(output_fp, "NPT_SCALE_VECS:");
      if (pSPARC->NPTscaleVecs[0] == 1)
        fprintf(output_fp, " 1");
      if (pSPARC->NPTscaleVecs[1] == 1)
        fprintf(output_fp, " 2");
      if (pSPARC->NPTscaleVecs[2] == 1)
        fprintf(output_fp, " 3");
      fprintf(output_fp, "\n");
      fprintf(output_fp, "NPT_NH_QMASS:");
      fprintf(output_fp, " %d", pSPARC->NPT_NHnnos);
      for (i = 0; i < pSPARC->NPT_NHnnos; i++) {
        if (i % 5 == 0) {
          fprintf(output_fp, "\n");
        }
        fprintf(output_fp, " %.15g", pSPARC->NPT_NHqmass[i]);
      }
      fprintf(output_fp, "\n");
      fprintf(output_fp, "NPT_NH_BMASS: %.15g\n", pSPARC->NPT_NHbmass);
      fprintf(output_fp, "TARGET_PRESSURE: %.15g GPa\n", pSPARC->prtarget);
    }
    if (strcmpi(pSPARC->MDMeth, "NPT_NP") == 0) {
      // fprintf(output_fp,"AMOUNT_THERMO_VARIABLE: %d\n",pSPARC->NPT_NHnnos);
      fprintf(output_fp, "NPT_NP_QMASS: %.15g\n", pSPARC->NPT_NP_qmass);
      fprintf(output_fp, "NPT_NP_BMASS: %.15g\n", pSPARC->NPT_NP_bmass);
      fprintf(output_fp, "TARGET_PRESSURE: %.15g GPa\n", pSPARC->prtarget);
    }
  }

  if (pSPARC->RestartFlag == 1) {
    fprintf(output_fp, "RESTART_FLAG: %d\n", pSPARC->RestartFlag);
  }
  if (pSPARC->NetCharge != 0) {
    fprintf(output_fp, "NET_CHARGE: %d\n", pSPARC->NetCharge);
  }
  fprintf(output_fp, "MAXIT_SCF: %d\n", pSPARC->MAXIT_SCF);
  fprintf(output_fp, "MINIT_SCF: %d\n", pSPARC->MINIT_SCF);
  fprintf(output_fp, "MAXIT_POISSON: %d\n", pSPARC->MAXIT_POISSON);
  if (pSPARC->scf_err_type == 0) {
    fprintf(output_fp, "TOL_SCF: %.2E\n", pSPARC->TOL_SCF);
  } else if (pSPARC->scf_err_type == 1) {
    fprintf(output_fp, "TOL_SCF_QE: %.2E\n", pSPARC->TOL_SCF);
    if (pSPARC->Nspin > 1) {
      fprintf(output_fp, "#WARNING: TOL_SCF_QE is not appropriatly set up for "
                         "spin-polarized systems\n");
    }
    if (pSPARC->MixingVariable == 1) {
      fprintf(output_fp, "#WARNING: TOL_SCF_QE is not appropriatly set up for "
                         "potential mixing\n");
    }
  }
  if (pSPARC->POISSON_SOLVER == 0) {
    fprintf(output_fp, "POISSON_SOLVER: AAR\n");
  } else {
    fprintf(output_fp, "POISSON_SOLVER: CG\n");
  }
  fprintf(output_fp, "TOL_POISSON: %.2E\n", pSPARC->TOL_POISSON);
  fprintf(output_fp, "TOL_LANCZOS: %.2E\n", pSPARC->TOL_LANCZOS);
  fprintf(output_fp, "TOL_PSEUDOCHARGE: %.2E\n", pSPARC->TOL_PSEUDOCHARGE);
  if (pSPARC->MixingVariable == 0) {
    fprintf(output_fp, "MIXING_VARIABLE: density\n");
  } else if (pSPARC->MixingVariable == 1) {
    fprintf(output_fp, "MIXING_VARIABLE: potential\n");
  }

  if (pSPARC->MixingPrecond == 0) {
    fprintf(output_fp, "MIXING_PRECOND: none\n");
  } else if (pSPARC->MixingPrecond == 1) {
    fprintf(output_fp, "MIXING_PRECOND: kerker\n");
  } else if (pSPARC->MixingPrecond == 2) {
    fprintf(output_fp, "MIXING_PRECOND: resta\n");
  } else if (pSPARC->MixingPrecond == 3) {
    fprintf(output_fp, "MIXING_PRECOND: truncated_kerker\n");
  }

  if (pSPARC->Nspin > 1) {
    if (pSPARC->MixingPrecondMag == 0) {
      fprintf(output_fp, "MIXING_PRECOND_MAG: none\n");
    } else if (pSPARC->MixingPrecondMag == 1) {
      fprintf(output_fp, "MIXING_PRECOND_MAG: kerker\n");
    } else if (pSPARC->MixingPrecondMag == 2) {
      fprintf(output_fp, "MIXING_PRECOND_MAG: resta\n");
    } else if (pSPARC->MixingPrecondMag == 3) {
      fprintf(output_fp, "MIXING_PRECOND_MAG: truncated_kerker\n");
    }
  }

  // for large periodic systems, give warning if preconditioner is not chosen
  if (pSPARC->BC == 2) {
    double Lx, Ly, Lz, L_diag;
    Lx = pSPARC->range_x;
    Ly = pSPARC->range_y;
    Lz = pSPARC->range_z;
    L_diag = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
    if (L_diag > 20.0 && pSPARC->MixingPrecond == 0) {
      fprintf(output_fp, "#WARNING: the preconditioner for SCF has been turned off, this \n"
                         "#might lead to slow SCF convergence. To specify SCF preconditioner, "
                         "\n"
                         "#use 'MIXING_PRECOND' in the .inpt file\n");
    }
  }

  if (pSPARC->MixingPrecond != 0) {
    fprintf(output_fp, "TOL_PRECOND: %.2E\n", pSPARC->TOL_PRECOND);
  }

  if (pSPARC->MixingPrecond == 1) { // kerker
    fprintf(output_fp, "PRECOND_KERKER_KTF: %.10G\n", pSPARC->precond_kerker_kTF);
    fprintf(output_fp, "PRECOND_KERKER_THRESH: %.10G\n", pSPARC->precond_kerker_thresh);
  } else if (pSPARC->MixingPrecond == 2) { // resta
    fprintf(output_fp, "PRECOND_RESTA_Q0: %.3f\n", pSPARC->precond_resta_q0);
    fprintf(output_fp, "PRECOND_RESTA_RS: %.3f\n", pSPARC->precond_resta_Rs);
    fprintf(output_fp, "PRECOND_FITPOW: %d\n", pSPARC->precond_fitpow);
  } else if (pSPARC->MixingPrecond == 3) { // truncated kerker
    fprintf(output_fp, "PRECOND_KERKER_KTF: %.10G\n", pSPARC->precond_kerker_kTF);
    fprintf(output_fp, "PRECOND_KERKER_THRESH: %.10G\n", pSPARC->precond_kerker_thresh);
    fprintf(output_fp, "PRECOND_FITPOW: %d\n", pSPARC->precond_fitpow);
  }

  if (pSPARC->Nspin > 1) {
    if (pSPARC->MixingPrecondMag == 1) {
      fprintf(output_fp, "PRECOND_KERKER_KTF_MAG: %.10G\n", pSPARC->precond_kerker_kTF_mag);
      fprintf(output_fp, "PRECOND_KERKER_THRESH_MAG: %.10G\n", pSPARC->precond_kerker_thresh_mag);
    }
  }

  fprintf(output_fp, "MIXING_PARAMETER: %.10G\n", pSPARC->MixingParameter);
  if (pSPARC->PulayFrequency > 1) {
    fprintf(output_fp, "MIXING_PARAMETER_SIMPLE: %.10G\n", pSPARC->MixingParameterSimple);
  }
  if (pSPARC->Nspin > 1) {
    fprintf(output_fp, "MIXING_PARAMETER_MAG: %.10G\n", pSPARC->MixingParameterMag);
    if (pSPARC->PulayFrequency > 1) {
      fprintf(output_fp, "MIXING_PARAMETER_SIMPLE_MAG: %.10G\n", pSPARC->MixingParameterSimpleMag);
    }
  }
  fprintf(output_fp, "MIXING_HISTORY: %d\n", pSPARC->MixingHistory);
  fprintf(output_fp, "PULAY_FREQUENCY: %d\n", pSPARC->PulayFrequency);
  fprintf(output_fp, "PULAY_RESTART: %d\n", pSPARC->PulayRestartFlag);
  fprintf(output_fp, "REFERENCE_CUTOFF: %.10g\n", pSPARC->REFERENCE_CUTOFF);
  fprintf(output_fp, "RHO_TRIGGER: %d\n", pSPARC->rhoTrigger);
  fprintf(output_fp, "FIX_RAND: %d\n", pSPARC->FixRandSeed);
  if (pSPARC->StandardEigenFlag == 1)
    fprintf(output_fp, "STANDARD_EIGEN: %d\n", pSPARC->StandardEigenFlag);
  fprintf(output_fp, "VERBOSITY: %d\n", pSPARC->Verbosity);
  fprintf(output_fp, "PRINT_FORCES: %d\n", pSPARC->PrintForceFlag);
  fprintf(output_fp, "PRINT_ATOMS: %d\n", pSPARC->PrintAtomPosFlag);
  fprintf(output_fp, "PRINT_EIGEN: %d\n", pSPARC->PrintEigenFlag);
  fprintf(output_fp, "PRINT_DENSITY: %d\n", pSPARC->PrintElecDensFlag);
  if (pSPARC->MDFlag == 1)
    fprintf(output_fp, "PRINT_MDOUT: %d\n", pSPARC->PrintMDout);
  if (pSPARC->MDFlag == 1 || pSPARC->RelaxFlag >= 1) {
    fprintf(output_fp, "PRINT_VELS: %d\n", pSPARC->PrintAtomVelFlag);
    fprintf(output_fp, "PRINT_RESTART: %d\n", pSPARC->Printrestart);
    if (pSPARC->Printrestart == 1)
      fprintf(output_fp, "PRINT_RESTART_FQ: %d\n", pSPARC->Printrestart_fq);
  }
  if (pSPARC->PrintPsiFlag[0] == 1) {
    fprintf(output_fp, "PRINT_ORBITAL: %d %d %d %d %d %d %d\n", pSPARC->PrintPsiFlag[0], pSPARC->PrintPsiFlag[1],
            pSPARC->PrintPsiFlag[2], pSPARC->PrintPsiFlag[3], pSPARC->PrintPsiFlag[4], pSPARC->PrintPsiFlag[5],
            pSPARC->PrintPsiFlag[6]);
  }
  fprintf(output_fp, "PRINT_ENERGY_DENSITY: %d\n", pSPARC->PrintEnergyDensFlag);

  if (pSPARC->RelaxFlag == 1) {
    fprintf(output_fp, "TOL_RELAX: %.2E\n", pSPARC->TOL_RELAX);
    fprintf(output_fp, "PRINT_RELAXOUT: %d\n", pSPARC->PrintRelaxout);
  } else if (pSPARC->RelaxFlag == 2) {
    fprintf(output_fp, "TOL_RELAX_CELL: %.2E\n", pSPARC->TOL_RELAX_CELL);
    fprintf(output_fp, "RELAX_MAXDILAT: %.2E\n", pSPARC->max_dilatation);
    fprintf(output_fp, "PRINT_RELAXOUT: %d\n", pSPARC->PrintRelaxout);
  } else if (pSPARC->RelaxFlag == 3) {
    fprintf(output_fp, "TOL_RELAX: %.2E\n", pSPARC->TOL_RELAX);
    fprintf(output_fp, "TOL_RELAX_CELL: %.2E\n", pSPARC->TOL_RELAX_CELL);
    fprintf(output_fp, "RELAX_MAXDILAT: %.2E\n", pSPARC->max_dilatation);
    fprintf(output_fp, "PRINT_RELAXOUT: %d\n", pSPARC->PrintRelaxout);
  }
  if (pSPARC->usefock == 1) {
    fprintf(output_fp, "EXX_FRAC: %.5g\n", pSPARC->exx_frac);
    fprintf(output_fp, "TOL_FOCK: %.2E\n", pSPARC->TOL_FOCK);
    fprintf(output_fp, "TOL_SCF_INIT: %.2E\n", pSPARC->TOL_SCF_INIT);
    fprintf(output_fp, "MAXIT_FOCK: %d\n", pSPARC->MAXIT_FOCK);
    fprintf(output_fp, "MINIT_FOCK: %d\n", pSPARC->MINIT_FOCK);
    if (pSPARC->EXXMeth_Flag == 0)
      fprintf(output_fp, "EXX_METHOD: FOURIER_SPACE\n");
    else
      fprintf(output_fp, "EXX_METHOD: REAL_SPACE\n");
    if (pSPARC->EXXDiv_Flag == 0)
      fprintf(output_fp, "EXX_DIVERGENCE: SPHERICAL\n");
    else if (pSPARC->EXXDiv_Flag == 1)
      fprintf(output_fp, "EXX_DIVERGENCE: AUXILIARY\n");
    else if (pSPARC->EXXDiv_Flag == 2)
      fprintf(output_fp, "EXX_DIVERGENCE: ERFC\n");
    fprintf(output_fp, "EXX_MEM: %d\n", pSPARC->EXXMem_batch);
    fprintf(output_fp, "ACE_FLAG: %d\n", pSPARC->ACEFlag);
    if (pSPARC->ACEFlag == 1) {
      fprintf(output_fp, "EXX_ACE_VALENCE_STATES: %d\n", pSPARC->EXXACEVal_state);
    }
    fprintf(output_fp, "EXX_DOWNSAMPLING: %d %d %d\n", pSPARC->EXXDownsampling[0], pSPARC->EXXDownsampling[1],
            pSPARC->EXXDownsampling[2]);
  }
  if (pSPARC->d3Flag == 1) {
    fprintf(output_fp, "D3_FLAG: %d\n", pSPARC->d3Flag);
    fprintf(output_fp, "D3_RTHR: %.15G\n", pSPARC->d3Rthr);
    fprintf(output_fp, "D3_CN_THR: %.15G\n", pSPARC->d3Cn_thr);
  }

  fprintf(output_fp, "OUTPUT_FILE: %s\n", pSPARC->filename_out);
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "                                Cell                     "
                     "                  \n");
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "Lattice vectors (Bohr):\n");
  fprintf(output_fp, "%.15f %.15f %.15f \n", pSPARC->LatUVec[0] * pSPARC->range_x, pSPARC->LatUVec[1] * pSPARC->range_x,
          pSPARC->LatUVec[2] * pSPARC->range_x);
  fprintf(output_fp, "%.15f %.15f %.15f \n", pSPARC->LatUVec[3] * pSPARC->range_y, pSPARC->LatUVec[4] * pSPARC->range_y,
          pSPARC->LatUVec[5] * pSPARC->range_y);
  fprintf(output_fp, "%.15f %.15f %.15f \n", pSPARC->LatUVec[6] * pSPARC->range_z, pSPARC->LatUVec[7] * pSPARC->range_z,
          pSPARC->LatUVec[8] * pSPARC->range_z);
  fprintf(output_fp, "Volume: %-.10E (Bohr^3)\n",
          pSPARC->range_x * pSPARC->range_y * pSPARC->range_z * pSPARC->Jacbdet);
  fprintf(output_fp, "Density: %-.10E (amu/Bohr^3), %-.10E (g/cc)\n",
          pSPARC->TotalMass / (pSPARC->range_x * pSPARC->range_y * pSPARC->range_z * pSPARC->Jacbdet),
          pSPARC->TotalMass / (pSPARC->range_x * pSPARC->range_y * pSPARC->range_z * pSPARC->Jacbdet) *
              CONST_AMU_BOHR3_GCC);
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "                           Parallelization               "
                     "                  \n");
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  if (pSPARC->SQFlag == 1) {
    fprintf(output_fp, "NP_DOMAIN_SQ_PARAL: %d %d %d\n", pSPARC->npNdx_SQ, pSPARC->npNdy_SQ, pSPARC->npNdz_SQ);
    fprintf(output_fp, "NP_DOMAIN_PHI_PARAL: %d %d %d\n", pSPARC->npNdx_phi, pSPARC->npNdy_phi, pSPARC->npNdz_phi);
  } else {
    if (pSPARC->num_node || pSPARC->num_cpu_per_node || pSPARC->num_acc_per_node) {
      fprintf(output_fp, "# Command line arguments: ");
      if (pSPARC->num_node)
        fprintf(output_fp, "-n %d ", pSPARC->num_node);
      if (pSPARC->num_cpu_per_node)
        fprintf(output_fp, "-c %d ", pSPARC->num_cpu_per_node);
      if (pSPARC->num_acc_per_node)
        fprintf(output_fp, "-a %d ", pSPARC->num_acc_per_node);
      fprintf(output_fp, "\n");
    }
    fprintf(output_fp, "NP_SPIN_PARAL: %d\n", pSPARC->npspin);
    fprintf(output_fp, "NP_KPOINT_PARAL: %d\n", pSPARC->npkpt);
    fprintf(output_fp, "NP_BAND_PARAL: %d\n", pSPARC->npband);
    fprintf(output_fp, "NP_DOMAIN_PARAL: %d %d %d\n", pSPARC->npNdx, pSPARC->npNdy, pSPARC->npNdz);
    fprintf(output_fp, "NP_DOMAIN_PHI_PARAL: %d %d %d\n", pSPARC->npNdx_phi, pSPARC->npNdy_phi, pSPARC->npNdz_phi);
    fprintf(output_fp, "EIG_SERIAL_MAXNS: %d\n", pSPARC->eig_serial_maxns);
    if (pSPARC->useLAPACK == 0) {
      fprintf(output_fp, "EIG_PARAL_BLKSZ: %d\n", pSPARC->eig_paral_blksz);
      fprintf(output_fp, "EIG_PARAL_ORFAC: %.1e\n", pSPARC->eig_paral_orfac);
      fprintf(output_fp, "EIG_PARAL_MAXNP: %d\n", pSPARC->eig_paral_maxnp);
    }
  }
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "                             Initialization              "
                     "                  \n");
  fprintf(output_fp, "*********************************************************"
                     "******************\n");
  fprintf(output_fp, "Number of processors               :  %d\n", nproc);

  if ((fabs(pSPARC->delta_x - pSPARC->delta_y) <= 1e-12) && (fabs(pSPARC->delta_x - pSPARC->delta_z) <= 1e-12) &&
      (fabs(pSPARC->delta_y - pSPARC->delta_z) <= 1e-12)) {
    fprintf(output_fp, "Mesh spacing                       :  %.6g (Bohr)\n", pSPARC->delta_x);
  } else {
    fprintf(output_fp, "Mesh spacing in x-direction        :  %.6g (Bohr)\n", pSPARC->delta_x);
    fprintf(output_fp, "Mesh spacing in y-direction        :  %.6g (Bohr)\n", pSPARC->delta_y);
    fprintf(output_fp, "Mesh spacing in z-direction        :  %.6g (Bohr)\n", pSPARC->delta_z);
  }

  if (pSPARC->BC == 2 || pSPARC->BC == 3 || pSPARC->BC == 4) {
    fprintf(output_fp, "Number of symmetry adapted k-points:  %d\n", pSPARC->Nkpts_sym);
  }

  fprintf(output_fp, "Output printed to                  :  %s\n", pSPARC->OutFilename);

  // if (pSPARC->PrintAtomPosFlag==1) {
  //    fprintf(output_fp,"Atom positions printed to          :
  //    %s\n",pSPARC->AtomFilename);
  //}

  // if (pSPARC->PrintForceFlag==1) {
  //    fprintf(output_fp,"Forces printed to                  :
  //    %s\n",pSPARC->StaticFilename);
  //}

  if (pSPARC->PrintEigenFlag == 1) {
    fprintf(output_fp, "Final eigenvalues printed to       :  %s\n", pSPARC->EigenFilename);
  }

  if (pSPARC->MDFlag == 1 && pSPARC->PrintMDout == 1) {
    fprintf(output_fp, "MD output printed to               :  %s\n", pSPARC->MDFilename);
  }

  if (pSPARC->RelaxFlag == 1 && pSPARC->PrintRelaxout == 1) {
    fprintf(output_fp, "Relax output printed to            :  %s\n", pSPARC->RelaxFilename);
  }

  fprintf(output_fp, "Total number of atom types         :  %d\n", pSPARC->Ntypes);
  fprintf(output_fp, "Total number of atoms              :  %d\n", pSPARC->n_atom);
  fprintf(output_fp, "Total number of electrons          :  %d\n", pSPARC->Nelectron);

  // count = 0;
  for (i = 0; i < pSPARC->Ntypes; i++) {
    fprintf(output_fp, "Atom type %-2d (valence electrons)   :  %s %d\n", i + 1, &pSPARC->atomType[L_ATMTYPE * i],
            pSPARC->Znucl[i]);
    fprintf(output_fp, "Pseudopotential                    :  %s\n", pSPARC->psdName + i * L_PSD);
    // fprintf(output_fp,"lloc                               :
    // %d\n",pSPARC->localPsd[i]);
    fprintf(output_fp, "Atomic mass                        :  %.15g\n", pSPARC->Mass[i]);
    fprintf(output_fp,
            "Pseudocharge radii of atom type %-2d :  %.2f %.2f %.2f (x, y, z "
            "dir)\n",
            i + 1, pSPARC->CUTOFF_x[i], pSPARC->CUTOFF_y[i], pSPARC->CUTOFF_z[i]);
    fprintf(output_fp, "Number of atoms of type %-2d         :  %d\n", i + 1, pSPARC->nAtomv[i]);
    // if (pSPARC->PrintAtomPosFlag == 1 && pSPARC->MDFlag == 0 &&
    // pSPARC->RelaxFlag == 0) {
    //     fprintf(output_fp,"Fractional coordinates of atoms of type %-2d
    //     :\n",i+1); for (j = 0; j < pSPARC->nAtomv[i]; j++) {
    //         fprintf(output_fp,"%18.10f %18.10f
    //         %18.10f\n",pSPARC->atom_pos[3*count]/pSPARC->range_x,pSPARC->atom_pos[3*count+1]/pSPARC->range_y,pSPARC->atom_pos[3*count+2]/pSPARC->range_z);
    //         count++;
    //     }
    // }
  }

  char mem_str[32];
  formatBytes(pSPARC->memory_usage, 32, mem_str);
  fprintf(output_fp, "Estimated total memory usage       :  %s\n", mem_str);
  formatBytes(pSPARC->memory_usage / nproc, 32, mem_str);
  fprintf(output_fp, "Estimated memory per processor     :  %s\n", mem_str);

  fclose(output_fp);

  // write .static file
  if ((pSPARC->PrintAtomPosFlag == 1 || pSPARC->PrintForceFlag == 1) && pSPARC->MDFlag == 0 && pSPARC->RelaxFlag == 0) {
    FILE *static_fp = fopen(pSPARC->StaticFilename, "w");
    if (static_fp == NULL) {
      printf("\nCannot open file \"%s\"\n", pSPARC->StaticFilename);
      exit(EXIT_FAILURE);
    }

    // print atoms
    if (pSPARC->PrintAtomPosFlag == 1) {
      fprintf(static_fp, "*****************************************************"
                         "**********************\n");
      fprintf(static_fp, "                            Atom positions           "
                         "                      \n");
      fprintf(static_fp, "*****************************************************"
                         "**********************\n");
      count = 0;
      for (i = 0; i < pSPARC->Ntypes; i++) {
        fprintf(static_fp, "Fractional coordinates of %s:\n", &pSPARC->atomType[L_ATMTYPE * i]);
        for (j = 0; j < pSPARC->nAtomv[i]; j++) {
          fprintf(static_fp, "%18.10f %18.10f %18.10f\n", pSPARC->atom_pos[3 * count] / pSPARC->range_x,
                  pSPARC->atom_pos[3 * count + 1] / pSPARC->range_y, pSPARC->atom_pos[3 * count + 2] / pSPARC->range_z);
          count++;
        }
      }
    }
    fclose(static_fp);
  }
}

