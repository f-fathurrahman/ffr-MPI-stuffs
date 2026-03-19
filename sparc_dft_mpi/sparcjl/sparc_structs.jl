
struct _D2D_OBJ
    n_target::Cint
    target_ranks::Ptr{Cint}
end

const D2D_OBJ = _D2D_OBJ

struct _PSD_OBJ
    rVloc::Ptr{Cdouble}
    UdV::Ptr{Cdouble}
    rhoIsoAtom::Ptr{Cdouble}
    RadialGrid::Ptr{Cdouble}
    SplinerVlocD::Ptr{Cdouble}
    SplineFitUdV::Ptr{Cdouble}
    SplineFitIsoAtomDen::Ptr{Cdouble}
    SplineRhocD::Ptr{Cdouble}
    rc::Ptr{Cdouble}
    Gamma::Ptr{Cdouble}
    rho_c_table::Ptr{Cdouble}
    fchrg::Cdouble
    Vloc_0::Cdouble
    is_r_uniform::Cint
    pspxc::Cint
    lmax::Cint
    size::Cint
    ppl::Ptr{Cint}
    pspsoc::Cint
    ppl_soc::Ptr{Cint}
    Gamma_soc::Ptr{Cdouble}
    UdV_soc::Ptr{Cdouble}
    SplineFitUdV_soc::Ptr{Cdouble}
end

const PSD_OBJ = _PSD_OBJ

struct _ATOM_INFLUENCE_OBJ
    n_atom::Cint
    coords::Ptr{Cdouble}
    atom_spin::Ptr{Cdouble}
    atom_index::Ptr{Cint}
    xs::Ptr{Cint}
    xe::Ptr{Cint}
    ys::Ptr{Cint}
    ye::Ptr{Cint}
    zs::Ptr{Cint}
    ze::Ptr{Cint}
end

const ATOM_INFLUENCE_OBJ = _ATOM_INFLUENCE_OBJ

struct _ATOM_NLOC_INFLUENCE_OBJ
    n_atom::Cint
    coords::Ptr{Cdouble}
    atom_index::Ptr{Cint}
    xs::Ptr{Cint}
    xe::Ptr{Cint}
    ys::Ptr{Cint}
    ye::Ptr{Cint}
    zs::Ptr{Cint}
    ze::Ptr{Cint}
    ndc::Ptr{Cint}
    grid_pos::Ptr{Ptr{Cint}}
end

const ATOM_NLOC_INFLUENCE_OBJ = _ATOM_NLOC_INFLUENCE_OBJ

struct _NLOC_PROJ_OBJ
    nproj::Cint
    Chi::Ptr{Ptr{Cdouble}}
    Chi_c::Ptr{Ptr{ComplexF32}}
    nprojso::Cint
    Chiso::Ptr{Ptr{ComplexF32}}
    nprojso_ext::Cint
    Chisowt0::Ptr{Ptr{ComplexF32}}
    Chisowtl::Ptr{Ptr{ComplexF32}}
    Chisowtnl::Ptr{Ptr{ComplexF32}}
    Chi_cyclix::Ptr{Ptr{Cdouble}}
    Chi_c_cyclix::Ptr{Ptr{ComplexF32}}
    Chiso_cyclix::Ptr{Ptr{ComplexF32}}
    Chisowt0_cyclix::Ptr{Ptr{ComplexF32}}
    Chisowtl_cyclix::Ptr{Ptr{ComplexF32}}
    Chisowtnl_cyclix::Ptr{Ptr{ComplexF32}}
end

const NLOC_PROJ_OBJ = _NLOC_PROJ_OBJ

struct _SPARC_OBJ
    SPARCROOT::NTuple{512, Cchar}
    filename::NTuple{512, Cchar}
    filename_out::NTuple{512, Cchar}
    OutFilename::NTuple{512, Cchar}
    StaticFilename::NTuple{512, Cchar}
    AtomFilename::NTuple{512, Cchar}
    EigenFilename::NTuple{512, Cchar}
    MDFilename::NTuple{512, Cchar}
    RelaxFilename::NTuple{512, Cchar}
    restart_Filename::NTuple{512, Cchar}
    restartC_Filename::NTuple{512, Cchar}
    restartP_Filename::NTuple{512, Cchar}
    DensTCubFilename::NTuple{512, Cchar}
    DensDCubFilename::NTuple{512, Cchar}
    DensUCubFilename::NTuple{512, Cchar}
    OrbitalsFilename::NTuple{512, Cchar}
    KinEnDensTCubFilename::NTuple{512, Cchar}
    KinEnDensUCubFilename::NTuple{512, Cchar}
    KinEnDensDCubFilename::NTuple{512, Cchar}
    XcEnDensCubFilename::NTuple{512, Cchar}
    ExxEnDensTCubFilename::NTuple{512, Cchar}
    ExxEnDensUCubFilename::NTuple{512, Cchar}
    ExxEnDensDCubFilename::NTuple{512, Cchar}
    num_node::Cint
    num_cpu_per_node::Cint
    num_acc_per_node::Cint
    npspin::Cint
    npkpt::Cint
    npband::Cint
    npNdx::Cint
    npNdy::Cint
    npNdz::Cint
    npNd::Cint
    npNdx_phi::Cint
    npNdy_phi::Cint
    npNdz_phi::Cint
    npNdx_kptcomm::Cint
    npNdy_kptcomm::Cint
    npNdz_kptcomm::Cint
    useDefaultParalFlag::Cint
    FixRandSeed::Cint
    spincomm::Cint
    spin_bridge_comm::Cint
    kptcomm::Cint
    kptcomm_topo::Cint
    kptcomm_topo_excl::Cint
    kptcomm_inter::Cint
    kpt_bridge_comm::Cint
    bandcomm::Cint
    dmcomm::Cint
    blacscomm::Cint
    ictxt_blacs::Cint
    ictxt_blacs_topo::Cint
    nprow_ictxt_blacs_topo::Cint
    npcol_ictxt_blacs_topo::Cint
    desc_orbitals::NTuple{9, Cint}
    desc_orb_BLCYC::NTuple{9, Cint}
    desc_Hp_BLCYC::NTuple{9, Cint}
    desc_Mp_BLCYC::NTuple{9, Cint}
    desc_Q_BLCYC::NTuple{9, Cint}
    dmcomm_phi::Cint
    comm_dist_graph_phi::Cint
    comm_dist_graph_psi::Cint
    kptcomm_topo_dist_graph::Cint
    spincomm_index::Cint
    Nspin_spincomm::Cint
    Nspinor_spincomm::Cint
    kptcomm_index::Cint
    bandcomm_index::Cint
    Nkpts_kptcomm::Cint
    Nband_bandcomm::Cint
    band_start_indx::Cint
    band_end_indx::Cint
    spin_typ::Cint
    Nspin::Cint
    mag::Ptr{Cdouble}
    mag_at::Ptr{Cdouble}
    netM::NTuple{4, Cdouble}
    spin_start_indx::Cint
    spinor_start_indx::Cint
    spin_end_indx::Cint
    spinor_end_indx::Cint
    Nspinor::Cint
    SOC_Flag::Cint
    Nspinor_eig::Cint
    Nspdentd::Cint
    Nspdend::Cint
    Nspden::Cint
    Nmag::Cint
    MDFlag::Cint
    RelaxFlag::Cint
    MDMeth::NTuple{32, Cchar}
    RelaxMeth::NTuple{32, Cchar}
    RestartFlag::Cint
    TWtime::Cdouble
    cell_typ::Cint
    Flag_latvec_scale::Cint
    numIntervals_x::Cint
    numIntervals_y::Cint
    numIntervals_z::Cint
    Nx::Cint
    Ny::Cint
    Nz::Cint
    Nd::Cint
    range_x::Cdouble
    range_y::Cdouble
    range_z::Cdouble
    latvec_scale_x::Cdouble
    latvec_scale_y::Cdouble
    latvec_scale_z::Cdouble
    LatVec::NTuple{9, Cdouble}
    delta_x::Cdouble
    delta_y::Cdouble
    delta_z::Cdouble
    dV::Cdouble
    mesh_spacing::Cdouble
    ecut::Cdouble
    Jacbdet::Cdouble
    LatUVec::NTuple{9, Cdouble}
    metricT::NTuple{9, Cdouble}
    gradT::NTuple{9, Cdouble}
    lapcT::NTuple{9, Cdouble}
    DMVertices::NTuple{6, Cint}
    Nx_d::Cint
    Ny_d::Cint
    Nz_d::Cint
    Nd_d::Cint
    DMVertices_kptcomm::NTuple{6, Cint}
    Nx_d_kptcomm::Cint
    Ny_d_kptcomm::Cint
    Nz_d_kptcomm::Cint
    Nd_d_kptcomm::Cint
    DMVertices_dmcomm::NTuple{6, Cint}
    Nx_d_dmcomm::Cint
    Ny_d_dmcomm::Cint
    Nz_d_dmcomm::Cint
    Nd_d_dmcomm::Cint
    Atom_Influence_local::Ptr{ATOM_INFLUENCE_OBJ}
    isRbOut::NTuple{3, Cint}
    Atom_Influence_nloc::Ptr{ATOM_NLOC_INFLUENCE_OBJ}
    Atom_Influence_nloc_kptcomm::Ptr{ATOM_NLOC_INFLUENCE_OBJ}
    nlocProj::Ptr{NLOC_PROJ_OBJ}
    nlocProj_kptcomm::Ptr{NLOC_PROJ_OBJ}
    nlocProjso::Ptr{NLOC_PROJ_OBJ}
    nlocProjso_kptcomm::Ptr{NLOC_PROJ_OBJ}
    IP_displ::Ptr{Cint}
    IP_displ_SOC::Ptr{Cint}
    order::Cint
    FDweights_D1::Ptr{Cdouble}
    FDweights_D2::Ptr{Cdouble}
    D2_stencil_coeffs_x::Ptr{Cdouble}
    D2_stencil_coeffs_y::Ptr{Cdouble}
    D2_stencil_coeffs_z::Ptr{Cdouble}
    D2_stencil_coeffs_xy::Ptr{Cdouble}
    D2_stencil_coeffs_yz::Ptr{Cdouble}
    D2_stencil_coeffs_xz::Ptr{Cdouble}
    D1_stencil_coeffs_x::Ptr{Cdouble}
    D1_stencil_coeffs_y::Ptr{Cdouble}
    D1_stencil_coeffs_z::Ptr{Cdouble}
    D1_stencil_coeffs_xy::Ptr{Cdouble}
    D1_stencil_coeffs_yx::Ptr{Cdouble}
    D1_stencil_coeffs_xz::Ptr{Cdouble}
    D1_stencil_coeffs_zx::Ptr{Cdouble}
    D1_stencil_coeffs_yz::Ptr{Cdouble}
    D1_stencil_coeffs_zy::Ptr{Cdouble}
    MaxEigVal_mhalfLap::Cdouble
    coeffs_x::Cdouble
    coeffs_y::Cdouble
    coeffs_z::Cdouble
    stencil_coeffs_x::Cdouble
    stencil_coeffs_y::Cdouble
    stencil_coeffs_z::Cdouble
    coeffs_grad_x::Cdouble
    coeffs_grad_y::Cdouble
    coeffs_grad_z::Cdouble
    POISSON_SOLVER::Cint
    MAXIT_SCF::Cint
    MINIT_SCF::Cint
    MAXIT_POISSON::Cint
    Relax_Niter::Cint
    accuracy_level::Cint
    target_force_accuracy::Cdouble
    target_energy_accuracy::Cdouble
    TOL_SCF::Cdouble
    TOL_RELAX::Cdouble
    TOL_POISSON::Cdouble
    TOL_LANCZOS::Cdouble
    TOL_PSEUDOCHARGE::Cdouble
    TOL_PRECOND::Cdouble
    precond_fitpow::Cint
    precondcoeff_n::Cint
    precond_kerker_kTF::Cdouble
    precond_kerker_thresh::Cdouble
    precond_kerker_kTF_mag::Cdouble
    precond_kerker_thresh_mag::Cdouble
    precond_resta_q0::Cdouble
    precond_resta_Rs::Cdouble
    precondcoeff_k::Cdouble
    precondcoeff_a::Ptr{ComplexF32}
    precondcoeff_lambda_sqr::Ptr{ComplexF32}
    RelaxCount::Cint
    StressCount::Cint
    REFERENCE_CUTOFF::Cdouble
    CUTOFF_x::Ptr{Cdouble}
    CUTOFF_y::Ptr{Cdouble}
    CUTOFF_z::Ptr{Cdouble}
    Lanczos_x0::Ptr{Cdouble}
    Lanczos_x0_complex::Ptr{ComplexF32}
    Veff_loc_dmcomm::Ptr{Cdouble}
    Veff_loc_dmcomm_phi::Ptr{Cdouble}
    Veff_dia_loc_dmcomm_phi::Ptr{Cdouble}
    Veff_loc_dmcomm_phi_in::Ptr{Cdouble}
    Veff_loc_kptcomm_topo::Ptr{Cdouble}
    mixing_hist_xk::Ptr{Cdouble}
    mixing_hist_xkm1::Ptr{Cdouble}
    mixing_hist_fk::Ptr{Cdouble}
    mixing_hist_fkm1::Ptr{Cdouble}
    mixing_hist_Xk::Ptr{Cdouble}
    mixing_hist_Fk::Ptr{Cdouble}
    mixing_hist_Pfk::Ptr{Cdouble}
    psdChrgDens::Ptr{Cdouble}
    psdChrgDens_ref::Ptr{Cdouble}
    Vc::Ptr{Cdouble}
    electronDens::Ptr{Cdouble}
    electronDens_at::Ptr{Cdouble}
    electronDens_core::Ptr{Cdouble}
    electronDens_in::Ptr{Cdouble}
    elecstPotential::Ptr{Cdouble}
    XCPotential::Ptr{Cdouble}
    XCPotential_nc::Ptr{Cdouble}
    e_xc::Ptr{Cdouble}
    Dxcdgrho::Ptr{Cdouble}
    xc_rhotol::Cdouble
    xc_magtol::Cdouble
    xc_sigmatol::Cdouble
    ixc::NTuple{4, Cint}
    xcoption::NTuple{2, Cint}
    isgradient::Cint
    occ::Ptr{Cdouble}
    occ_sorted::Ptr{Cdouble}
    occfac::Cdouble
    lambda::Ptr{Cdouble}
    lambda_sorted::Ptr{Cdouble}
    Xorb::Ptr{Cdouble}
    Yorb::Ptr{Cdouble}
    Xorb_BLCYC::Ptr{Cdouble}
    Yorb_BLCYC::Ptr{Cdouble}
    Xorb_kpt::Ptr{ComplexF32}
    Yorb_kpt::Ptr{ComplexF32}
    Xorb_BLCYC_kpt::Ptr{ComplexF32}
    Yorb_BLCYC_kpt::Ptr{ComplexF32}
    nr_orb_BLCYC::Cint
    nc_orb_BLCYC::Cint
    nr_Hp_BLCYC::Cint
    nc_Hp_BLCYC::Cint
    nr_Mp_BLCYC::Cint
    nc_Mp_BLCYC::Cint
    nr_Q_BLCYC::Cint
    nc_Q_BLCYC::Cint
    Hp::Ptr{Cdouble}
    Mp::Ptr{Cdouble}
    Q::Ptr{Cdouble}
    Hp_kpt::Ptr{ComplexF32}
    Mp_kpt::Ptr{ComplexF32}
    Q_kpt::Ptr{ComplexF32}
    Hp_s::Ptr{Cdouble}
    Mp_s::Ptr{Cdouble}
    useLAPACK::Cint
    eig_serial_maxns::Cint
    eig_paral_blksz::Cint
    eig_paral_orfac::Cdouble
    eig_paral_maxnp::Cint
    eig_paral_subdims::NTuple{2, Cint}
    req_veff_loc::Cint
    d2d_dmcomm_phi::D2D_OBJ
    d2d_dmcomm::D2D_OBJ
    d2d_dmcomm_lanczos::D2D_OBJ
    d2d_kptcomm_topo::D2D_OBJ
    is_phi_eq_kpt_topo::Cint
    MixingVariable::Cint
    MixingPrecond::Cint
    MixingPrecondMag::Cint
    MixingHistory::Cint
    PulayFrequency::Cint
    PulayRestartFlag::Cint
    MixingParameter::Cdouble
    MixingParameterSimple::Cdouble
    MixingParameterMag::Cdouble
    MixingParameterSimpleMag::Cdouble
    Nkpts::Cint
    Kx::Cint
    Ky::Cint
    Kz::Cint
    Nkpts_sym::Cint
    NkptsGroup::Cint
    kptParalFlag::Cint
    kctr::Cint
    kpt_start_indx::Cint
    kpt_end_indx::Cint
    isGammaPoint::Cint
    kptshift::NTuple{3, Cdouble}
    lambdakpt::Ptr{Cdouble}
    kptWts::Ptr{Cdouble}
    k1::Ptr{Cdouble}
    k2::Ptr{Cdouble}
    k3::Ptr{Cdouble}
    kptWts_loc::Ptr{Cdouble}
    k1_loc::Ptr{Cdouble}
    k2_loc::Ptr{Cdouble}
    k3_loc::Ptr{Cdouble}
    BC::Cint
    BCx::Cint
    BCy::Cint
    BCz::Cint
    Nstates::Cint
    Ntypes::Cint
    Nelectron::Cint
    NetCharge::Cint
    PosCharge::Cdouble
    NegCharge::Cdouble
    elec_T_type::Cint
    Beta::Cdouble
    elec_T::Cdouble
    rho::Cdouble
    XC::NTuple{32, Cchar}
    is_default_psd::Cint
    NLCC_flag::Cint
    n_atom::Cint
    localPsd::Ptr{Cint}
    Mass::Ptr{Cdouble}
    TotalMass::Cdouble
    atomType::Ptr{Cchar}
    Zatom::Ptr{Cint}
    Znucl::Ptr{Cint}
    nAtomv::Ptr{Cint}
    psdName::Ptr{Cchar}
    atom_pos::Ptr{Cdouble}
    IsFrac::Ptr{Cint}
    IsSpin::Ptr{Cint}
    mvAtmConstraint::Ptr{Cint}
    atom_spin::Ptr{Cdouble}
    psd::Ptr{PSD_OBJ}
    ChebDegree::Cint
    CheFSI_Optmz::Cint
    rhoTrigger::Cint
    chefsibound_flag::Cint
    eigmin::Ptr{Cdouble}
    eigmax::Ptr{Cdouble}
    npl_min::Cint
    npl_max::Cint
    Nchefsi::Cint
    Esc::Cdouble
    Efermi::Cdouble
    Exc::Cdouble
    Exc_corr::Cdouble
    Eband::Cdouble
    Entropy::Cdouble
    Etot::Cdouble
    Escc::Cdouble
    forces::Ptr{Cdouble}
    AtomMag::Ptr{Cdouble}
    Calc_stress::Cint
    Calc_pres::Cint
    stress::NTuple{6, Cdouble}
    stress_k::NTuple{6, Cdouble}
    stress_xc::NTuple{6, Cdouble}
    stress_el::NTuple{6, Cdouble}
    stress_nl::NTuple{6, Cdouble}
    stress_i::NTuple{6, Cdouble}
    stress_exx::NTuple{6, Cdouble}
    pres::Cdouble
    pres_xc::Cdouble
    pres_el::Cdouble
    pres_nl::Cdouble
    pres_i::Cdouble
    pres_exx::Cdouble
    ion_vel::Ptr{Cdouble}
    ion_accel::Ptr{Cdouble}
    ion_T::Cdouble
    PE::Cdouble
    KE::Cdouble
    TE::Cdouble
    TE_ext::Cdouble
    kB::Cdouble
    MD_dt::Cdouble
    mean_elec_T::Cdouble
    mean_ion_T::Cdouble
    mean_TE::Cdouble
    mean_KE::Cdouble
    mean_PE::Cdouble
    mean_U::Cdouble
    mean_Entropy::Cdouble
    mean_TE_ext::Cdouble
    std_elec_T::Cdouble
    std_ion_T::Cdouble
    std_TE::Cdouble
    std_KE::Cdouble
    std_PE::Cdouble
    std_U::Cdouble
    std_Entropy::Cdouble
    std_TE_ext::Cdouble
    MD_maxStep::Cint
    MD_Nstep::Cint
    dof::Cint
    MDCount::Cint
    restartCount::Cint
    StopCount::Cint
    ion_vel_dstr::Cint
    ion_vel_dstr_rand::Cint
    ion_elec_eqT::Cint
    thermos_T::Cdouble
    thermos_Ti::Cdouble
    thermos_Tf::Cdouble
    v2nose::Cdouble
    snose::Cdouble
    xi_nose::Cdouble
    qmass::Cdouble
    amu2au::Cdouble
    fs2atu::Cdouble
    NPTscaleVecs::NTuple{3, Cint}
    NPTconstraintFlag::Cint
    NPTisotropicFlag::Cint
    prtarget::Cdouble
    scale::Cdouble
    volumeCell::Cdouble
    initialLatVecLength::NTuple{3, Cdouble}
    NPT_NHnnos::Cint
    NPT_NHqmass::NTuple{60, Cdouble}
    NPT_NHbmass::Cdouble
    glogs::NTuple{60, Cdouble}
    vlogs::NTuple{60, Cdouble}
    vlogv::Cdouble
    xlogs::NTuple{60, Cdouble}
    Hamiltonian_NPT_NH::Cdouble
    maxTimeIter::Cint
    NPT_NP_qmass::Cdouble
    NPT_NP_bmass::Cdouble
    range_x_velo::Cdouble
    range_y_velo::Cdouble
    range_z_velo::Cdouble
    G_NPT_NP::NTuple{3, Cdouble}
    Pm_NPT_NP::NTuple{3, Cdouble}
    Kbaro::Cdouble
    Ubaro::Cdouble
    S_NPT_NP::Cdouble
    Sv_NPT_NP::Cdouble
    Kther::Cdouble
    Uther::Cdouble
    Hamiltonian_NPT_NP::Cdouble
    init_Hamil_NPT_NP::Cdouble
    Relax_fac::Cdouble
    elecgs_Count::Cint
    d::Ptr{Cdouble}
    NLCG_sigma::Cdouble
    L_history::Cint
    L_finit_stp::Cdouble
    L_maxmov::Cdouble
    L_autoscale::Cint
    L_lineopt::Cint
    L_icurv::Cdouble
    deltaX::Ptr{Cdouble}
    deltaG::Ptr{Cdouble}
    iys::Ptr{Cdouble}
    fold::Ptr{Cdouble}
    atom_disp::Ptr{Cdouble}
    isFD::Cint
    isReset::Cint
    step::Cint
    FIRE_dt::Cdouble
    FIRE_mass::Cdouble
    FIRE_maxmov::Cdouble
    FIRE_alpha::Cdouble
    FIRE_vel::Ptr{Cdouble}
    FIRE_dtNow::Cdouble
    FIRE_resetIter::Cint
    cellrelax_ndim::Cint
    cellrelax_dims::NTuple{3, Cint}
    max_dilatation::Cdouble
    TOL_RELAX_CELL::Cdouble
    d3Energy::NTuple{4, Cdouble}
    d3Grads::Ptr{Cdouble}
    d3Stress::NTuple{9, Cdouble}
    d3Flag::Cint
    d3Rthr::Cdouble
    d3Cn_thr::Cdouble
    atomicNumbers::Ptr{Cint}
    atomScaledR2R4::Ptr{Cdouble}
    atomScaledRcov::Ptr{Cdouble}
    atomCN::Ptr{Cdouble}
    atomMaxci::Ptr{Cint}
    nImageCN::NTuple{3, Cint}
    nImageEG::NTuple{3, Cint}
    c6ab::Ptr{Ptr{Ptr{Ptr{Ptr{Cdouble}}}}}
    lattice::NTuple{9, Cdouble}
    BCtype::NTuple{3, Cint}
    periodicBCFlag::Cint
    d3Rs6::Cdouble
    d3S18::Cdouble
    d3Output::Ptr{Cint}
    d3Sigma::NTuple{9, Cdouble}
    vdWDFnrpoints::Cint
    vdWDFnqs::Cint
    vdWDFecLinear::Ptr{Cdouble}
    vdWDFVcLinear::Ptr{Cdouble}
    vdWDFdr::Cdouble
    vdWDFdk::Cdouble
    vdWDFZab::Cdouble
    detLattice::Cdouble
    vdWDFqmesh::Ptr{Cdouble}
    vdWDFkernelPhi::Ptr{Ptr{Cdouble}}
    vdWDFd2Phidk2::Ptr{Ptr{Cdouble}}
    vdWDFd2Splineydx2::Ptr{Ptr{Cdouble}}
    Drho::Ptr{Ptr{Cdouble}}
    gradRhoLen::Ptr{Cdouble}
    vdWDFq0::Ptr{Cdouble}
    vdWDFdq0drho::Ptr{Cdouble}
    vdWDFdq0dgradrho::Ptr{Cdouble}
    vdWDFps::Ptr{Ptr{Cdouble}}
    vdWDFdpdq0s::Ptr{Ptr{Cdouble}}
    zAxisComm::Cint
    newzAxisDims::NTuple{3, Cint}
    zAxisVertices::NTuple{6, Cint}
    gatherThetasSender::D2D_OBJ
    gatherThetasRecvr::D2D_OBJ
    scatterThetasSender::D2D_OBJ
    scatterThetasRecvr::D2D_OBJ
    vdWDFthetaFTs::Ptr{Ptr{ComplexF32}}
    reciLattice::NTuple{9, Cdouble}
    timeReciLattice::Ptr{Ptr{Cint}}
    vdWDFreciLatticeGrid::Ptr{Ptr{Cdouble}}
    vdWDFreciLength::Ptr{Cdouble}
    vdWDFkernelReciPoints::Ptr{Ptr{Cdouble}}
    vdWDFuFTs::Ptr{Ptr{ComplexF32}}
    vdWDFu::Ptr{Ptr{Cdouble}}
    vdWDFenergy::Cdouble
    vdWDFpotential::Ptr{Cdouble}
    countPotentialCalculate::Cint
    KineticTauPhiDomain::Ptr{Cdouble}
    vxcMGGA1::Ptr{Cdouble}
    vxcMGGA2::Ptr{Cdouble}
    vxcMGGA3::Ptr{Cdouble}
    vxcMGGA3_loc_dmcomm::Ptr{Cdouble}
    vxcMGGA3_loc_kptcomm::Ptr{Cdouble}
    usefock::Cint
    TOL_FOCK::Cdouble
    TOL_SCF_INIT::Cdouble
    MAXIT_FOCK::Cint
    MINIT_FOCK::Cint
    exx_frac::Cdouble
    hyb_range_fock::Cdouble
    hyb_range_pbe::Cdouble
    EXXMeth_Flag::Cint
    Eexx::Cdouble
    psi_outer::Ptr{Cdouble}
    occ_outer::Ptr{Cdouble}
    psi_outer_kptcomm_topo::Ptr{Cdouble}
    psi_outer_kpt::Ptr{ComplexF32}
    psi_outer_kptcomm_topo_kpt::Ptr{ComplexF32}
    Xi::Ptr{Cdouble}
    Xi_kptcomm_topo::Ptr{Cdouble}
    Xi_kpt::Ptr{ComplexF32}
    Xi_kptcomm_topo_kpt::Ptr{ComplexF32}
    Nkpts_hf::Cint
    Kx_hf::Cint
    Ky_hf::Cint
    Kz_hf::Cint
    Nkpts_hf_red::Cint
    k1_hf::Ptr{Cdouble}
    k2_hf::Ptr{Cdouble}
    k3_hf::Ptr{Cdouble}
    kptWts_hf::Cdouble
    kpthf_ind::Ptr{Cint}
    kpthf_ind_red::Ptr{Cint}
    kpthf_pn::Ptr{Cint}
    k1_shift::Ptr{Cdouble}
    k2_shift::Ptr{Cdouble}
    k3_shift::Ptr{Cdouble}
    Nkpts_shift::Cint
    Kptshift_map::Ptr{Cint}
    Nkpts_hf_kptcomm::Cint
    kpthf_flag_kptcomm::Ptr{Cint}
    Nkpts_hf_list::Ptr{Cint}
    kpthf_start_indx::Cint
    kpthf_start_indx_list::Ptr{Cint}
    kpts_hf_red_list::Ptr{Cint}
    neg_phase::Ptr{ComplexF32}
    pos_phase::Ptr{ComplexF32}
    kpthfred2kpthf::Ptr{NTuple{3, Cint}}
    ACEtime::Cdouble
    Exxtime::Cdouble
    pois_FFT_const::Ptr{Cdouble}
    pois_FFT_const_stress::Ptr{Cdouble}
    pois_FFT_const_stress2::Ptr{Cdouble}
    pois_FFT_const_press::Ptr{Cdouble}
    ACEFlag::Cint
    Nstates_occ::Cint
    Nstates_occ_list::NTuple{2, Cint}
    EXXMem_batch::Cint
    EXXACEVal_state::Cint
    EXXDownsampling::NTuple{3, Cint}
    const_aux::Cdouble
    EXXDiv_Flag::Cint
    flag_kpttopo_dm::Cint
    flag_kpttopo_dm_type::Cint
    kpttopo_dmcomm_inter::Cint
    desc_M::NTuple{9, Cint}
    desc_Xi::NTuple{9, Cint}
    nrows_M::Cint
    ncols_M::Cint
    Nband_bandcomm_M::Cint
    SQFlag::Cint
    SQ_gauss_mem::Cint
    SQ_npl_g::Cint
    SQ_correction::Cint
    SQ_rcut::Cdouble
    SQ_tol_occ::Cdouble
    npNdx_SQ::Cint
    npNdy_SQ::Cint
    npNdz_SQ::Cint
    d2d_s2p_sq::D2D_OBJ
    d2d_s2p_phi::D2D_OBJ
    Atom_Influence_nloc_SQ::Ptr{Ptr{ATOM_NLOC_INFLUENCE_OBJ}}
    nlocProj_SQ::Ptr{Ptr{NLOC_PROJ_OBJ}}
    pSQ::Ptr{Cint}
    CyclixFlag::Cint
    twist::Cdouble
    twistpercell::Cdouble
    xin::Cdouble
    RotM_cyclix::NTuple{9, Cdouble}
    xmin_at::Cdouble
    xmax_at::Cdouble
    xvac::Cdouble
    xout::Cdouble
    Intgwt_kpttopo::Ptr{Cdouble}
    Intgwt_psi::Ptr{Cdouble}
    Intgwt_phi::Ptr{Cdouble}
    lambda_temp1::Ptr{Cdouble}
    lambda_temp2::Ptr{Cdouble}
    lambda_temp3::Ptr{Cdouble}
    vl::Ptr{Cdouble}
    vr::Ptr{Cdouble}
    lambda_temp1_kpt::Ptr{ComplexF32}
    lambda_temp2_kpt::Ptr{ComplexF32}
    vl_kpt::Ptr{ComplexF32}
    vr_kpt::Ptr{ComplexF32}
    delectronDens::Ptr{Cdouble}
    delectronDens_0dt::Ptr{Cdouble}
    delectronDens_1dt::Ptr{Cdouble}
    delectronDens_2dt::Ptr{Cdouble}
    atom_pos_nm::Ptr{Cdouble}
    atom_pos_0dt::Ptr{Cdouble}
    atom_pos_1dt::Ptr{Cdouble}
    atom_pos_2dt::Ptr{Cdouble}
    Verbosity::Cint
    PrintForceFlag::Cint
    PrintAtomPosFlag::Cint
    PrintAtomVelFlag::Cint
    PrintEigenFlag::Cint
    PrintElecDensFlag::Cint
    PrintMDout::Cint
    PrintRelaxout::Cint
    Printrestart::Cint
    Printrestart_fq::Cint
    suffixNum::Cint
    PrintPsiFlag::NTuple{7, Cint}
    PrintEnergyDensFlag::Cint
    KineticRho::Ptr{Cdouble}
    ExxRho::Ptr{Cdouble}
    ExcRho::Ptr{Cdouble}
    time_start::Cdouble
    memory_usage::Cdouble
    DP_CheFSI::Ptr{Cvoid}
    DP_CheFSI_kpt::Ptr{Cvoid}
    StandardEigenFlag::Cint
    n_kpt_line::Cint
    kredx::NTuple{512, Cdouble}
    kredy::NTuple{512, Cdouble}
    kredz::NTuple{512, Cdouble}
    k1_inpt_kpt::NTuple{512, Cdouble}
    k2_inpt_kpt::NTuple{512, Cdouble}
    k3_inpt_kpt::NTuple{512, Cdouble}
    BandStructFlag::Cint
    kpt_per_line::Cint
    InDensTCubFilename::NTuple{512, Cchar}
    InDensUCubFilename::NTuple{512, Cchar}
    InDensDCubFilename::NTuple{512, Cchar}
    densfilecount::Cint
end

const SPARC_OBJ = _SPARC_OBJ

struct _SPARC_INPUT_OBJ
    num_node::Cint
    num_cpu_per_node::Cint
    num_acc_per_node::Cint
    npspin::Cint
    npkpt::Cint
    npband::Cint
    npNdx::Cint
    npNdy::Cint
    npNdz::Cint
    npNdx_phi::Cint
    npNdy_phi::Cint
    npNdz_phi::Cint
    eig_serial_maxns::Cint
    eig_paral_blksz::Cint
    MDFlag::Cint
    RelaxFlag::Cint
    RestartFlag::Cint
    Flag_latvec_scale::Cint
    numIntervals_x::Cint
    numIntervals_y::Cint
    numIntervals_z::Cint
    spin_typ::Cint
    BC::Cint
    BCx::Cint
    BCy::Cint
    BCz::Cint
    Nstates::Cint
    Ntypes::Cint
    NetCharge::Cint
    order::Cint
    ChebDegree::Cint
    CheFSI_Optmz::Cint
    chefsibound_flag::Cint
    rhoTrigger::Cint
    Nchefsi::Cint
    FixRandSeed::Cint
    accuracy_level::Cint
    MAXIT_SCF::Cint
    MINIT_SCF::Cint
    MAXIT_POISSON::Cint
    Relax_Niter::Cint
    MixingVariable::Cint
    MixingPrecond::Cint
    MixingPrecondMag::Cint
    MixingHistory::Cint
    PulayFrequency::Cint
    PulayRestartFlag::Cint
    precond_fitpow::Cint
    Nkpts::Cint
    Kx::Cint
    Ky::Cint
    Kz::Cint
    NkptsGroup::Cint
    kctr::Cint
    Verbosity::Cint
    PrintForceFlag::Cint
    PrintAtomPosFlag::Cint
    PrintAtomVelFlag::Cint
    PrintEigenFlag::Cint
    PrintElecDensFlag::Cint
    PrintMDout::Cint
    PrintRelaxout::Cint
    Printrestart::Cint
    Printrestart_fq::Cint
    PrintPsiFlag::NTuple{7, Cint}
    PrintEnergyDensFlag::Cint
    elec_T_type::Cint
    MD_Nstep::Cint
    ion_elec_eqT::Cint
    ion_vel_dstr::Cint
    ion_vel_dstr_rand::Cint
    L_history::Cint
    L_autoscale::Cint
    L_lineopt::Cint
    Calc_stress::Cint
    Calc_pres::Cint
    Poisson_solver::Cint
    range_x::Cdouble
    range_y::Cdouble
    range_z::Cdouble
    latvec_scale_x::Cdouble
    latvec_scale_y::Cdouble
    latvec_scale_z::Cdouble
    LatVec::NTuple{9, Cdouble}
    mesh_spacing::Cdouble
    ecut::Cdouble
    kptshift::NTuple{3, Cdouble}
    target_energy_accuracy::Cdouble
    target_force_accuracy::Cdouble
    TOL_SCF::Cdouble
    TOL_RELAX::Cdouble
    TOL_POISSON::Cdouble
    TOL_LANCZOS::Cdouble
    TOL_PSEUDOCHARGE::Cdouble
    TOL_PRECOND::Cdouble
    precond_kerker_kTF::Cdouble
    precond_kerker_thresh::Cdouble
    precond_kerker_kTF_mag::Cdouble
    precond_kerker_thresh_mag::Cdouble
    precond_resta_q0::Cdouble
    precond_resta_Rs::Cdouble
    REFERENCE_CUTOFF::Cdouble
    Beta::Cdouble
    elec_T::Cdouble
    MixingParameter::Cdouble
    MixingParameterSimple::Cdouble
    MixingParameterMag::Cdouble
    MixingParameterSimpleMag::Cdouble
    MD_dt::Cdouble
    ion_T::Cdouble
    thermos_Tf::Cdouble
    qmass::Cdouble
    NPT_NHnnos::Cint
    NPTscaleVecs::NTuple{3, Cint}
    NPTconstraintFlag::Cint
    NPT_NHqmass::NTuple{60, Cdouble}
    NPT_NHbmass::Cdouble
    prtarget::Cdouble
    NPT_NP_qmass::Cdouble
    NPT_NP_bmass::Cdouble
    TWtime::Cdouble
    NLCG_sigma::Cdouble
    L_finit_stp::Cdouble
    L_maxmov::Cdouble
    L_icurv::Cdouble
    FIRE_dt::Cdouble
    FIRE_mass::Cdouble
    FIRE_maxmov::Cdouble
    max_dilatation::Cdouble
    TOL_RELAX_CELL::Cdouble
    eig_paral_orfac::Cdouble
    eig_paral_maxnp::Cint
    MDMeth::NTuple{32, Cchar}
    RelaxMeth::NTuple{32, Cchar}
    XC::NTuple{32, Cchar}
    d3Flag::Cint
    d3Rthr::Cdouble
    d3Cn_thr::Cdouble
    TOL_FOCK::Cdouble
    TOL_SCF_INIT::Cdouble
    MAXIT_FOCK::Cint
    MINIT_FOCK::Cint
    EXXMeth_Flag::Cint
    ACEFlag::Cint
    EXXMem_batch::Cint
    EXXACEVal_state::Cint
    EXXDownsampling::NTuple{3, Cint}
    EXXDiv_Flag::Cint
    hyb_range_fock::Cdouble
    hyb_range_pbe::Cdouble
    exx_frac::Cdouble
    SQFlag::Cint
    SQ_gauss_mem::Cint
    SQ_npl_g::Cint
    SQ_rcut::Cdouble
    SQ_tol_occ::Cdouble
    npNdx_SQ::Cint
    npNdy_SQ::Cint
    npNdz_SQ::Cint
    twist::Cdouble
    StandardEigenFlag::Cint
    filename::NTuple{512, Cchar}
    filename_out::NTuple{512, Cchar}
    SPARCROOT::NTuple{512, Cchar}
    n_kpt_line::Cint
    kredx::NTuple{512, Cdouble}
    kredy::NTuple{512, Cdouble}
    kredz::NTuple{512, Cdouble}
    BandStructFlag::Cint
    kpt_per_line::Cint
    InDensTCubFilename::NTuple{512, Cchar}
    InDensUCubFilename::NTuple{512, Cchar}
    InDensDCubFilename::NTuple{512, Cchar}
    densfilecount::Cint
end

const SPARC_INPUT_OBJ = _SPARC_INPUT_OBJ