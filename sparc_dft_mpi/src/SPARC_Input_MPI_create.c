#include <stdio.h>
#include "initialization.h"

#define N_MEMBR 162

/**
 * @brief Create MPI struct type SPARC_INPUT_MPI for broadcasting.
 */
void SPARC_Input_MPI_create(MPI_Datatype *pSPARC_INPUT_MPI) {
  SPARC_INPUT_OBJ sparc_input_tmp;

  MPI_Datatype SPARC_types[N_MEMBR] = {
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,
      MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR,   MPI_CHAR,   MPI_CHAR,   MPI_CHAR,   MPI_CHAR,   MPI_CHAR};
  int blens[N_MEMBR] = {3,       3,  7, /* int array */
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* int */
                        9,       3,  L_QMASS,                                                         /* double array */
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1,       1,  1,       1,        1,        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* double */
                        32,      32, 32,      L_STRING, L_STRING,                                           /* char */
                        L_STRING};

  // calculating offsets in an architecture independent manner
  MPI_Aint addr[N_MEMBR], disps[N_MEMBR], base;
  int i = 0;
  MPI_Get_address(&sparc_input_tmp, &base);
  // int array type
  MPI_Get_address(&sparc_input_tmp.NPTscaleVecs, addr + i++);
  MPI_Get_address(&sparc_input_tmp.EXXDownsampling, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintPsiFlag, addr + i++);
  // int type
  MPI_Get_address(&sparc_input_tmp.num_node, addr + i++);
  MPI_Get_address(&sparc_input_tmp.num_cpu_per_node, addr + i++);
  MPI_Get_address(&sparc_input_tmp.num_acc_per_node, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npspin, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npkpt, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npband, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdx, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdy, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdz, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdx_phi, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdy_phi, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdz_phi, addr + i++);
  MPI_Get_address(&sparc_input_tmp.eig_serial_maxns, addr + i++);
  MPI_Get_address(&sparc_input_tmp.eig_paral_blksz, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MDFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.spin_typ, addr + i++);
  MPI_Get_address(&sparc_input_tmp.RelaxFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.RestartFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Flag_latvec_scale, addr + i++);
  MPI_Get_address(&sparc_input_tmp.numIntervals_x, addr + i++);
  MPI_Get_address(&sparc_input_tmp.numIntervals_y, addr + i++);
  MPI_Get_address(&sparc_input_tmp.numIntervals_z, addr + i++);
  MPI_Get_address(&sparc_input_tmp.BC, addr + i++);
  MPI_Get_address(&sparc_input_tmp.BCx, addr + i++);
  MPI_Get_address(&sparc_input_tmp.BCy, addr + i++);
  MPI_Get_address(&sparc_input_tmp.BCz, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Nstates, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Ntypes, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NetCharge, addr + i++);
  MPI_Get_address(&sparc_input_tmp.order, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ChebDegree, addr + i++);
  MPI_Get_address(&sparc_input_tmp.CheFSI_Optmz, addr + i++);
  MPI_Get_address(&sparc_input_tmp.chefsibound_flag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.rhoTrigger, addr + i++);
  MPI_Get_address(&sparc_input_tmp.FixRandSeed, addr + i++);
  MPI_Get_address(&sparc_input_tmp.accuracy_level, addr + i++);
  MPI_Get_address(&sparc_input_tmp.scf_err_type, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MAXIT_SCF, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MINIT_SCF, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MAXIT_POISSON, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Relax_Niter, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingVariable, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingPrecond, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingPrecondMag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingHistory, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PulayFrequency, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PulayRestartFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_fitpow, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Nkpts, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Kx, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Ky, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Kz, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NkptsGroup, addr + i++);
  MPI_Get_address(&sparc_input_tmp.kctr, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Verbosity, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintForceFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintAtomPosFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintAtomVelFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintEigenFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintElecDensFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintMDout, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintRelaxout, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Printrestart, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Printrestart_fq, addr + i++);
  MPI_Get_address(&sparc_input_tmp.elec_T_type, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MD_Nstep, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ion_elec_eqT, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ion_vel_dstr, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ion_vel_dstr_rand, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_history, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_autoscale, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_lineopt, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Calc_stress, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Calc_pres, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Poisson_solver, addr + i++);
  MPI_Get_address(&sparc_input_tmp.d3Flag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NPT_NHnnos, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MAXIT_FOCK, addr + i++);
  MPI_Get_address(&sparc_input_tmp.EXXMeth_Flag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ACEFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.EXXMem_batch, addr + i++);
  MPI_Get_address(&sparc_input_tmp.EXXACEVal_state, addr + i++);
  MPI_Get_address(&sparc_input_tmp.EXXDiv_Flag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MINIT_FOCK, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SQFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SQ_gauss_mem, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SQ_npl_g, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdx_SQ, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdy_SQ, addr + i++);
  MPI_Get_address(&sparc_input_tmp.npNdz_SQ, addr + i++);
  MPI_Get_address(&sparc_input_tmp.PrintEnergyDensFlag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.eig_paral_maxnp, addr + i++);
  MPI_Get_address(&sparc_input_tmp.StandardEigenFlag, addr + i++);
  // double array type
  MPI_Get_address(&sparc_input_tmp.LatVec, addr + i++);
  MPI_Get_address(&sparc_input_tmp.kptshift, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NPT_NHqmass, addr + i++);
  // double type
  MPI_Get_address(&sparc_input_tmp.range_x, addr + i++);
  MPI_Get_address(&sparc_input_tmp.range_y, addr + i++);
  MPI_Get_address(&sparc_input_tmp.range_z, addr + i++);
  MPI_Get_address(&sparc_input_tmp.latvec_scale_x, addr + i++);
  MPI_Get_address(&sparc_input_tmp.latvec_scale_y, addr + i++);
  MPI_Get_address(&sparc_input_tmp.latvec_scale_z, addr + i++);
  MPI_Get_address(&sparc_input_tmp.mesh_spacing, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ecut, addr + i++);
  MPI_Get_address(&sparc_input_tmp.target_energy_accuracy, addr + i++);
  MPI_Get_address(&sparc_input_tmp.target_force_accuracy, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_SCF, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_RELAX, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_POISSON, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_LANCZOS, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_PSEUDOCHARGE, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_PRECOND, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_kerker_kTF, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_kerker_thresh, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_kerker_kTF_mag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_kerker_thresh_mag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_resta_q0, addr + i++);
  MPI_Get_address(&sparc_input_tmp.precond_resta_Rs, addr + i++);
  MPI_Get_address(&sparc_input_tmp.REFERENCE_CUTOFF, addr + i++);
  MPI_Get_address(&sparc_input_tmp.Beta, addr + i++);
  MPI_Get_address(&sparc_input_tmp.elec_T, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingParameter, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingParameterSimple, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingParameterMag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MixingParameterSimpleMag, addr + i++);
  MPI_Get_address(&sparc_input_tmp.MD_dt, addr + i++);
  MPI_Get_address(&sparc_input_tmp.ion_T, addr + i++);
  MPI_Get_address(&sparc_input_tmp.thermos_Tf, addr + i++);
  MPI_Get_address(&sparc_input_tmp.qmass, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TWtime, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NLCG_sigma, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_finit_stp, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_maxmov, addr + i++);
  MPI_Get_address(&sparc_input_tmp.L_icurv, addr + i++);
  MPI_Get_address(&sparc_input_tmp.FIRE_dt, addr + i++);
  MPI_Get_address(&sparc_input_tmp.FIRE_mass, addr + i++);
  MPI_Get_address(&sparc_input_tmp.FIRE_maxmov, addr + i++);
  MPI_Get_address(&sparc_input_tmp.max_dilatation, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_RELAX_CELL, addr + i++);
  MPI_Get_address(&sparc_input_tmp.eig_paral_orfac, addr + i++);
  MPI_Get_address(&sparc_input_tmp.d3Rthr, addr + i++);
  MPI_Get_address(&sparc_input_tmp.d3Cn_thr, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NPT_NHbmass, addr + i++);
  MPI_Get_address(&sparc_input_tmp.prtarget, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NPT_NP_qmass, addr + i++);
  MPI_Get_address(&sparc_input_tmp.NPT_NP_bmass, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_FOCK, addr + i++);
  MPI_Get_address(&sparc_input_tmp.TOL_SCF_INIT, addr + i++);
  MPI_Get_address(&sparc_input_tmp.hyb_range_fock, addr + i++);
  MPI_Get_address(&sparc_input_tmp.hyb_range_pbe, addr + i++);
  MPI_Get_address(&sparc_input_tmp.exx_frac, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SQ_rcut, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SQ_tol_occ, addr + i++);
  // char type
  MPI_Get_address(&sparc_input_tmp.MDMeth, addr + i++);
  MPI_Get_address(&sparc_input_tmp.RelaxMeth, addr + i++);
  MPI_Get_address(&sparc_input_tmp.XC, addr + i++);
  MPI_Get_address(&sparc_input_tmp.filename, addr + i++);
  MPI_Get_address(&sparc_input_tmp.filename_out, addr + i++);
  MPI_Get_address(&sparc_input_tmp.SPARCROOT, addr + i++);

  for (i = 0; i < N_MEMBR; i++) {
    disps[i] = addr[i] - base;
  }

  MPI_Type_create_struct(N_MEMBR, blens, disps, SPARC_types, pSPARC_INPUT_MPI);
  MPI_Type_commit(pSPARC_INPUT_MPI);

  // MPI_Aint extend = sizeof(sparc_input_tmp);
  // MPI_Datatype SPARC_INPUT_MPI_tmp;
  // MPI_Type_create_struct(N_MEMBR, blens, disps, SPARC_types,
  // &SPARC_INPUT_MPI_tmp); MPI_Type_create_resized(SPARC_INPUT_MPI_tmp, 0,
  // extend, pSPARC_INPUT_MPI); MPI_Type_commit(pSPARC_INPUT_MPI);
}
