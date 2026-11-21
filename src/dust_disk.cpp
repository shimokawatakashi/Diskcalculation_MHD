#include "dust_disk.hpp"
#include "alpha.hpp"
#include "energy.hpp"

#include <algorithm>
#include <limits>
#include <iostream>
// #include <omp.h>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <stdio.h>
#include <sstream>

// (shared constants sqrt2/sqrt1_2/.../coe_vbr moved to dust_disk.hpp so
// dust_growth.cpp can use them too)

// DustDisk

DustDisk::DustDisk()
{
    // default parameter
    NG = 2;
    nr_disk_ = 100;
    nz_disk_ = 100;
    rin_ = 1.0 * ASTRONOMICAL_UNIT;
    rout_ = 1.0e2 * ASTRONOMICAL_UNIT;

    amin_ = 1.0e-5; // cm = 0.1 μm = 1.0e-7 m
    amax_ = 1.0e0;  // cm = 1.0e6 μm = 1.0 m
    rhodi_ = 1.4;
    mmin_ = (4.0 * M_PI / 3.0) * rhodi_ * CUB(amin_);
    mmax_ = (4.0 * M_PI / 3.0) * rhodi_ * CUB(amax_);
    nbin_ = 180;
    nbin_lin_ = 10;

    fdg_ = 1.0e-2;

    star_mass_ = 1.0 * SOLAR_MASS;

    on_fragmentation_ = false;
    frag_slope_ = 1.83;
    frag_frac_ = 1.0e-1;

    v_coag_ = 50.0 * 1.0e2;
    v_frag_ = 50.0 * 1.0e2;

    restart_ = false;
}

DustDisk::~DustDisk()
{
    
}

void DustDisk::InputParameter(Utils::InputConfigure &input)
{
    std::cout << "DustDisk::InputParameter" << std::endl;

    outdir_ = input.Get("outdir");
    outdir_mag = outdir_ + "_B";

    // disk radius paramteer
    rin_ = std::stod(input.Get("diskin"));
    rin_ *= ASTRONOMICAL_UNIT;
    rout_ = std::stod(input.Get("diskout"));
    rout_ *= ASTRONOMICAL_UNIT;
    r_cut = std :: stod(input.Get("cut_off"));
    r_cut *= ASTRONOMICAL_UNIT;
    z_max = std :: stod(input.Get("z_max"));

    nr_disk_ = std::stoi(input.Get("nrdisk"));
    nz_disk_ = std :: stoi(input.Get("nzdisk"));

    sigmag0_ = std::stod(input.Get("sigmag0_"));

    // dust size and mass paramter
    amin_ = std::stod(input.Get("amin"));
    amax_ = std::stod(input.Get("amax"));
    rhodi_ = std::stod(input.Get("rhodi"));
    mmin_ = (4.0 * M_PI / 3.0) * rhodi_ * CUB(amin_);
    mmax_ = (4.0 * M_PI / 3.0) * rhodi_ * CUB(amax_);
    nbin_ = std::stoi(input.Get("ndust"));
    nbin_lin_ = std::stoi(input.Get("nd_lin"));
    fdg_ = std::stod(input.Get("fdg"));
    no_dust_ = (fdg_ <= 0.0); // true : skip all dust motion/growth (CalcStokesNumberAndDustScaleHeight,
                              // CalcDustVelocity_TM2018, CalcDustFlux, CalculateDust_to_Gas_massratio,
                              // DetermineVaporDust, ApplyDustDiffusionImplicit, coagulation) -- fdg=0 means
                              // sigmad_ stays exactly 0 forever (SetInitialDustDistribution), so all of that
                              // is otherwise wasted (and O(nbin^3) per cell per substep for the full-implicit
                              // coagulation solve) work computing on an all-zero field.
    B_out = std :: stod(input.Get("B_out"));
    beta_ini = std :: stod(input.Get("beta_ini"));
    eta_cap = std :: stod(input.Get("eta_max"));
    Cw_      = std :: stod(input.Get("Cw"));
    inv_Cw_ = (Cw_!=0) ? 1.0/Cw_: 1.0e50;
    alphaDW_ini = std :: stod(input.Get("alpha_DW"));
    // wind alpha model : 0 = Suzuki2016 (default) , 1 = Enomoto et al. (2026) self-consistent alpha_DW
    wind_alpha_model_ = input.Contains("wind_alpha_model") ? std :: stoi(input.Get("wind_alpha_model")) : 0;
    // magnetic lever arm lambda of the disk wind (Tabone et al. 2022):
    //   dSigma/dt|_w = alpha_phiz*rho_mid*cs^2/v_K / (lambda - 1) = 3 alpha_DW Sigma cs^2 / (4 (lambda-1) r^2 Omega)
    // default 2.0 (lambda - 1 = 1) reproduces the previous hard-coded expression.
    lambda_w_ = input.Contains("lambda_w") ? std :: stod(input.Get("lambda_w")) : 2.0;
    if (!(lambda_w_ > 1.0))
    {
        std :: cerr << "Error: lambda_w must be > 1 (got " << lambda_w_ << ")" << std :: endl;
        std :: exit(1);
    }
    inv_lambda_w_m1_ = 1.0/(lambda_w_ - 1.0);
    if (lambda_w_ < 1.5)
        std :: cerr << "Warning: lambda_w < 1.5 makes the wind heating factor (2 lambda - 3)/(2 (lambda - 1)) negative" << std :: endl;

    // Super-time-stepping for Ohmic/ambipolar diffusion of B: SUPERSEDED, see
    // note/STS_magnetic_diffusion_note.md -- diverged even with the
    // Runge-Kutta-Gegenbauer generalization (arXiv:2606.13398) because eta's
    // ~18-orders-of-magnitude spatial range exceeds any STS variant's stable
    // stage-count range. Replaced by SetImplicitMagneticDiffusion below.
    // bool use_sts = input.Contains("use_sts") ? Utils::ConvertStrToBool(input.Get("use_sts")) : false;
    // double sts_alpha = input.Contains("sts_alpha") ? std :: stod(input.Get("sts_alpha")) : 10.0;
    // SetSuperTimeStepping(use_sts, sts_alpha);

    // Implicit (Thomas-algorithm tridiagonal) treatment of Ohmic/ambipolar
    // diffusion of B -- unconditionally stable regardless of how large eta or
    // dt/dr^2 get. Off by default; set use_implicit_mag=true in the input
    // file to enable.
    bool use_implicit_mag = input.Contains("use_implicit_mag") ? Utils::ConvertStrToBool(input.Get("use_implicit_mag")) : false;
    SetImplicitMagneticDiffusion(use_implicit_mag);

    // Face-flux scheme (face_scheme.cpp): two-point face fluxes, implicit
    // viscous/pressure diffusion, MUSCL (van Leer, ln Sigma) + SSP-RK2 advection,
    // no smoothing. false (default) keeps the original scheme.
    use_face_scheme_ = input.Contains("face_flux_scheme") ? Utils::ConvertStrToBool(input.Get("face_flux_scheme")) : false;
    SetMusclAdvectionB(use_face_scheme_);
    // Time weight of the implicit magnetic diffusion solve (1.0 = backward Euler,
    // 0.5 = Crank-Nicolson). Default 0.55 = previous hard-coded value.
    SetMagTheta(input.Contains("mag_theta") ? std::stod(input.Get("mag_theta")) : 0.55);

    // Semi-implicit coagulation/fragmentation (dust_coagulation_implicit_note.md).
    // Off by default; set use_implicit_coag=true in the input file to enable.
    use_implicit_coag_ = input.Contains("use_implicit_coag") ? Utils::ConvertStrToBool(input.Get("use_implicit_coag")) : false;

    // Fully implicit coagulation/fragmentation (Jacobian + dense linear solve
    // per cell, coagulation_speedup_literature.md Sec.3 / Brauer et al. 2008
    // Appendix C). Off by default; set use_full_implicit_coag=true to enable.
    // Takes priority over use_implicit_coag_ when both are set.
    use_full_implicit_coag_ = input.Contains("use_full_implicit_coag") ? Utils::ConvertStrToBool(input.Get("use_full_implicit_coag")) : false;

    // Elsasser-number(1e-5)/eta_cap(1e20) floor on the magnetic diffusivity
    // used by SetBoundaryConditionForMagneticfield (STS_magnetic_diffusion_note.md
    // Sec.12). On by default (matches historical behavior); a physics choice,
    // independent of use_implicit_mag -- set use_eta_floor=false to use the
    // raw (unfloored) resistivity instead.
    use_eta_floor_ = input.Contains("use_eta_floor") ? Utils::ConvertStrToBool(input.Get("use_eta_floor")) : true;

    // Radial grid spacing, selectable: "log" (default, unchanged behavior),
    // "linear" (uniform dr), "root" (uniform in sqrt(r), i.e. r^0.5).
    radial_grid_type_ = input.Contains("radial_grid_type") ? input.Get("radial_grid_type") : "log";
    // "log_floor" grid: minimum cell width near rin_, given in au in the
    // input file. Deliberately independent of NG (the ghost-cell count, a
    // numerical/stencil detail unrelated to how much resolution should be
    // spent on the inner disk) -- unlike the original coupled rin_/(NG+1)
    // formula, this can be set directly to whatever inner resolution is
    // physically wanted. <=0 (or absent) falls back to the old rin_/(NG+1)
    // behavior for backward compatibility.
    radial_grid_floor_dr_ = input.Contains("radial_grid_floor_dr")
        ? std::stod(input.Get("radial_grid_floor_dr")) * ASTRONOMICAL_UNIT : -1.0;
    // turbulence parameter
    alpha_ = std::stod(input.Get("alpha_t"));

    // central star mass
    star_mass_ = std::stod(input.Get("star_mass")) * SOLAR_MASS;

    // cfl number
    cfl_ = std::stod(input.Get("cfl"));

    // calculation time of simulation
    tend_ = std::stod(input.Get("tend")) * SOLAR_YEAR;
    on_thermal_ = Utils::ConvertStrToBool(input.Get("thermal"));
    rad_diff_fac_ = input.Contains("rad_diff_fac") ? std::stod(input.Get("rad_diff_fac")) : 1.0;
    zeta_CR0_ = input.Contains("zeta_CR0") ? std::stod(input.Get("zeta_CR0")) : 1.0e-17;

    // Initial dust size distribution (Initialization() only).
    dust_init_mode_ = input.Contains("dust_init_mode") ? input.Get("dust_init_mode") : "monomer";
    dust_init_powerlaw_q_ = input.Contains("dust_init_powerlaw_q") ? std::stod(input.Get("dust_init_powerlaw_q")) : 3.5;
    dust_init_amax_mode_ = input.Contains("dust_init_amax_mode") ? input.Get("dust_init_amax_mode") : "fixed";
    dust_init_amax_ = input.Contains("dust_init_amax") ? std::stod(input.Get("dust_init_amax")) : amin_;
    dust_init_St_target_ = input.Contains("dust_init_St_target") ? std::stod(input.Get("dust_init_St_target")) : 0.1;
    dust_init_f_frag_ = input.Contains("dust_init_f_frag") ? std::stod(input.Get("dust_init_f_frag")) : 0.37;
    dust_init_f_drift_ = input.Contains("dust_init_f_drift") ? std::stod(input.Get("dust_init_f_drift")) : 0.55;

    // fragmentation parameter
    on_fragmentation_ = Utils::ConvertStrToBool(input.Get("fragmentation"));

    on_growth_ = Utils::ConvertStrToBool(input.Get("growth"));
    frag_slope_ = std::stod(input.Get("frag_slope"));
    frag_frac_ = std::stod(input.Get("frag_frac"));
    v_col_crit_ = std::stod(input.Get("v_col_crit"));
    v_col_sili_ = std::stod(input.Get("v_col_sili"));
    T_snow_ = input.Contains("T_snow") ? std::stod(input.Get("T_snow")) : 160.0;

    v_col_crit_ *= 1.0e2;
    v_col_sili_ *= 1.0e2;

    v_coag_ = 0.2 * v_col_crit_;
    v_frag_ = 0.2 * v_col_crit_;

    // make directory and copy input file
    String filename;
    mkdir(outdir_.c_str(), 0777);
    mkdir(outdir_mag.c_str(),0777);

    filename = input.GetInputFileName();

    std::string inputf = outdir_ + "/input";
    std::string inputf_mag = outdir_mag + "/input";

    Utils::FileCopy(filename, inputf);
    Utils::FileCopy(filename, inputf_mag);
    
    return;
}

void DustDisk::Initialization()
{
    /*
        Initialization of DustDisk
    */
    std::cout << "DustDisk::Initialization" << std::endl;

    // set grid
    SetRadiusGrid();
    SetVerticalGrid();
    SetDustGrid();

    // set magnetic field
    B_.resize(nr_tot,0.0);
    dB_dr_.resize(nr_tot,0.0);
    beta_disk_.resize(nr_tot, DoubleVector(nz_tot,0.0));
    beta_disk_ini.resize(nr_tot, DoubleVector(nz_tot,0.0));
    // set gas paramter
    sigmag_.resize(nr_tot, 0.0);
    sigmag_bod_.resize(nr_tot+1,0.0);
    dsigmagdt_w.resize(nr_tot,0.0);
    dsigmaddt_w.resize(nr_tot,DoubleVector(nbin_,0.0));
    v_dsigmaddt_w.resize(nr_tot,0.0);

    M_in.resize(nr_tot,0.0);
    dMdt_.resize(nr_tot,0.0);
    T_.resize(nr_tot, 0.0);
    T_irr.resize(nr_tot,0.0);
    T_irr4.resize(nr_tot,0.0);
    cs_.resize(nr_tot, 0.0);
    Omega_.resize(nr_tot, 0.0);
    Hg_.resize(nr_tot, 0.0);
    rhog_mid_.resize(nr_tot, 0.0);
    rhog.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    ng_mid_.resize(nr_tot, 0.0);
    ng_.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    vg_r_.resize(nr_tot, 0.0);
    vgr_bod_.resize(nr_tot+1,0.0);
    vg_phi_.resize(nr_tot,0.0);
    dvg_phi.resize(nr_tot,0.0);
    j_.resize(nr_tot,0.0);
    vg_wind_.resize(nr_tot,0.0);
    vg_rtest_.resize(nr_tot,0.0);
    N_torque.resize(nr_tot,0.0);
    dNdr_torque.resize(nr_tot,0.0);
    dOmegadr_.resize(nr_tot,0.0);

    P_integral.resize(nr_tot, 0.0);
    Ptot_integral.resize(nr_tot,0.0);
    P_mag.resize(nr_tot,0.0);
    P_.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    Ptot_.resize(nr_tot,0.0);

    dP_dr.resize(nr_tot, 0.0);
    dPtot_dr.resize(nr_tot,0.0);
    dPinteg_dr.resize(nr_tot,0.0);
    dPintegtot_dr.resize(nr_tot,0.0);

    eta_P_.resize(nr_tot, 0.0);
    Q_toomre.resize(nr_tot,0.0);
    vK_.resize(nr_tot, 0.0);

    alpha_MRI.resize(nr_tot, 0.0);
    alpha_phiz.resize(nr_tot,0.0);
    alpha_GI.resize(nr_tot,0.0);
    alpha_tot.resize(nr_tot,0.0);

    windDW_zdead_.resize(nr_tot,0.0);
    windDW_zsurf_.resize(nr_tot,0.0);
    windDW_Delta_.resize(nr_tot,0.0);
    windDW_Ameff_.resize(nr_tot,0.0);
    alpha_DW_.resize(nr_tot,0.0);

    vg_T2_.resize(nr_tot, 0.0);
    nu_.resize(nr_tot, 0.0);
    Re_.resize(nr_tot, 0.0);
    inv_sqrt_Re_.resize(nr_tot, 0.0);

    Roe_mean_.resize(nr_tot,DoubleVector(2,0.0));
    Roe_mean_dust.resize(nr_tot,DoubleVector(2,0.0));

    zeta_.resize(nr_tot, DoubleVector(nz_tot + 1, 0.0)); // set ionization rate

    etaO.resize(nr_tot,DoubleVector(nz_tot,0.0));
    etaA.resize(nr_tot,DoubleVector(nz_tot,0.0));
    etaO_R.resize(nr_tot,0.0);
    etaA_R.resize(nr_tot,0.0);
    eta_total.resize(nr_tot,0.0);
    // Warm-start cache for CalcAlphaParameter::GetChargestate()'s phi Newton
    // iteration, per (i,k) grid cell; 1.0e-2 matches the solver's original
    // fixed starting guess, so the first call at each cell is unaffected.
    phi_prev_.resize(nr_tot,DoubleVector(nz_tot,1.0e-2));

    nu_mag.resize(nr_tot,0.0);
    numag_eq.resize(nr_tot,0.0);

    A_.resize(nr_tot,0.0);
    lambda1_.resize(nr_tot,0.0);
    lambda2_.resize(nr_tot,0.0);
    fluxg_.resize(nr_tot,0.0);

    SetDiskGasModel();

    // set dust
    sigmad_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_active.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_vapor.resize(nr_tot,0.0);

    sigmad_floor_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_tot_.resize(nr_tot, 0.0);
    Nd_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Nd_RZ.resize(nz_tot, DoubleVector(nbin_, 0.0));
    St_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    St_mean.resize(nr_tot,DoubleVector(nbin_, 0.0));
    vd_r_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    vd_r_tot.resize(nr_tot,DoubleVector(nbin_,0.0));
    vd_r_diff.resize(nr_tot,DoubleVector(nbin_,0.0));
    vd_phi_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Hd_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Dd_.resize(nr_tot,DoubleVector(nbin_,0.0));
    fdg_ratio.resize(nr_tot,DoubleVector(nbin_,0.0));
    dfdg_dr.resize(nr_tot,DoubleVector(nbin_,0.0));
    fluxd_.resize(nr_tot+1, DoubleVector(nbin_, 0.0));
    fluxd_vapor.resize(nr_tot+1,0.0);
    fluxd_diff_.resize(nr_tot+1, DoubleVector(nbin_, 0.0));
    fluxd_vapor_diff_.resize(nr_tot+1,0.0);

    vd_r_v.resize(nr_tot,0.0);
    vd_r_diff_v.resize(nr_tot,0.0);
    vd_r_tot_v.resize(nr_tot,0.0);
    fdg_ratio_vapor.resize(nr_tot,0.0);
    dfdg_dr_v.resize(nr_tot,0.0);

    // initial dust density
    SetInitialDustDistribution();

    // calculate sigmad_floor ( = one dust particle / disk annual area)
    // double area;
    // for (int i = nr_s; i < nr_e; ++i)
    // {
    //     area = M_PI * (SQR(r_bod_[i + 1]) - SQR(r_bod_[i]));
    //     for (int j = 0; j < nbin_; ++j)
    //     {
    //         sigmad_floor_[i][j] = m_[j] / area;
    //     }
    // }

    // calculate sigmad_floor ( = 1.0e-50)
    double area;
    for (int i = nr_s; i < nr_e; ++i)
    {
        area = M_PI * (SQR(r_bod_[i + 1]) - SQR(r_bod_[i]));
        for (int j = 0; j < nbin_; ++j)
        {
            sigmad_floor_[i][j] = std :: max(1.0e-50, m_[j] / area);
        }
    }

    CalcCoagulationAndFragmentationArrays();
    CalcStokesNumberAndDustScaleHeight();
    CalcDustVelocity();

    return;
}

void DustDisk::Initialization_Restart(std :: string dustfile,std :: string gasfile,std :: string Bfile,double t_restart)
{
    /*
        Initialization of DustDisk
    */
    std::cout << "DustDisk::Initialization_Restart" << std::endl;
   
    // set grid
    SetRadiusGrid();
    SetVerticalGrid();
    SetDustGrid();
    restart_ = true;

    // set magnetic field
    B_.resize(nr_tot,0.0);
    dB_dr_.resize(nr_tot,0.0);
    beta_disk_.resize(nr_tot, DoubleVector(nz_tot,0.0));
    beta_disk_ini.resize(nr_tot, DoubleVector(nz_tot,0.0));
    // set gas paramter
    sigmag_.resize(nr_tot, 0.0);
    sigmag_bod_.resize(nr_tot+1,0.0);
    dsigmagdt_w.resize(nr_tot,0.0);
    dsigmaddt_w.resize(nr_tot,DoubleVector(nbin_,0.0));
    v_dsigmaddt_w.resize(nr_tot,0.0);

    M_in.resize(nr_tot,0.0);
    dMdt_.resize(nr_tot,0.0);
    T_.resize(nr_tot, 0.0);
    T_irr.resize(nr_tot,0.0);
    T_irr4.resize(nr_tot,0.0);
    cs_.resize(nr_tot, 0.0);
    Omega_.resize(nr_tot, 0.0);
    Hg_.resize(nr_tot, 0.0);
    rhog_mid_.resize(nr_tot, 0.0);
    rhog.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    ng_mid_.resize(nr_tot, 0.0);
    ng_.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    vg_r_.resize(nr_tot, 0.0);
    vgr_bod_.resize(nr_tot+1,0.0);
    vg_phi_.resize(nr_tot,0.0);
    dvg_phi.resize(nr_tot,0.0);
    j_.resize(nr_tot,0.0);
    vg_wind_.resize(nr_tot,0.0);
    vg_rtest_.resize(nr_tot,0.0);
    N_torque.resize(nr_tot,0.0);
    dNdr_torque.resize(nr_tot,0.0);
    dOmegadr_.resize(nr_tot,0.0);

    P_integral.resize(nr_tot, 0.0);
    Ptot_integral.resize(nr_tot,0.0);
    P_mag.resize(nr_tot,0.0);
    P_.resize(nr_tot, DoubleVector(nz_tot, 0.0));
    Ptot_.resize(nr_tot,0.0);

    dP_dr.resize(nr_tot, 0.0);
    dPtot_dr.resize(nr_tot,0.0);
    dPinteg_dr.resize(nr_tot,0.0);
    dPintegtot_dr.resize(nr_tot,0.0);

    eta_P_.resize(nr_tot, 0.0);
    Q_toomre.resize(nr_tot,0.0);
    vK_.resize(nr_tot, 0.0);

    alpha_MRI.resize(nr_tot, 0.0);
    alpha_phiz.resize(nr_tot,0.0);
    alpha_GI.resize(nr_tot,0.0);
    alpha_tot.resize(nr_tot,0.0);

    windDW_zdead_.resize(nr_tot,0.0);
    windDW_zsurf_.resize(nr_tot,0.0);
    windDW_Delta_.resize(nr_tot,0.0);
    windDW_Ameff_.resize(nr_tot,0.0);
    alpha_DW_.resize(nr_tot,0.0);

    vg_T2_.resize(nr_tot, 0.0);
    nu_.resize(nr_tot, 0.0);
    Re_.resize(nr_tot, 0.0);
    inv_sqrt_Re_.resize(nr_tot, 0.0);

    Roe_mean_.resize(nr_tot,DoubleVector(2,0.0));
    Roe_mean_dust.resize(nr_tot,DoubleVector(2,0.0));

    zeta_.resize(nr_tot, DoubleVector(nz_tot + 1, 0.0)); // set ionization rate

    etaO.resize(nr_tot,DoubleVector(nz_tot,0.0));
    etaA.resize(nr_tot,DoubleVector(nz_tot,0.0));
    etaO_R.resize(nr_tot,0.0);
    etaA_R.resize(nr_tot,0.0);
    eta_total.resize(nr_tot,0.0);
    // Warm-start cache for CalcAlphaParameter::GetChargestate()'s phi Newton
    // iteration, per (i,k) grid cell; 1.0e-2 matches the solver's original
    // fixed starting guess, so the first call at each cell is unaffected.
    phi_prev_.resize(nr_tot,DoubleVector(nz_tot,1.0e-2));

    nu_mag.resize(nr_tot,0.0);
    numag_eq.resize(nr_tot,0.0);

    A_.resize(nr_tot,0.0);
    lambda1_.resize(nr_tot,0.0);
    lambda2_.resize(nr_tot,0.0);
    fluxg_.resize(nr_tot,0.0);

    // set dust
    sigmad_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_active.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_vapor.resize(nr_tot,0.0);

    sigmad_floor_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    sigmad_tot_.resize(nr_tot, 0.0);
    Nd_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Nd_RZ.resize(nz_tot, DoubleVector(nbin_, 0.0));
    St_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    St_mean.resize(nr_tot,DoubleVector(nbin_, 0.0));
    vd_r_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    vd_r_tot.resize(nr_tot,DoubleVector(nbin_,0.0));
    vd_r_diff.resize(nr_tot,DoubleVector(nbin_,0.0));
    vd_r_v.resize(nr_tot,0.0);
    vd_r_diff_v.resize(nr_tot,0.0);
    vd_r_tot_v.resize(nr_tot,0.0);
    vd_phi_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Hd_.resize(nr_tot, DoubleVector(nbin_, 0.0));
    Dd_.resize(nr_tot,DoubleVector(nbin_,0.0));
    fdg_ratio.resize(nr_tot,DoubleVector(nbin_,0.0));
    fdg_ratio_vapor.resize(nr_tot,0.0);
    dfdg_dr.resize(nr_tot,DoubleVector(nbin_,0.0));
    dfdg_dr_v.resize(nr_tot,0.0);
    fluxd_.resize(nr_tot+1, DoubleVector(nbin_, 0.0));
    fluxd_vapor.resize(nr_tot+1,0.0);
    fluxd_diff_.resize(nr_tot+1, DoubleVector(nbin_, 0.0));
    fluxd_vapor_diff_.resize(nr_tot+1,0.0);

    SetDiskGasModel();

    // // initial deust density
    // for (int i = nr_s; i < nr_e; ++i)
    // {
    //     for (int j = 0; j < nbin_; ++j)
    //     {
    //         if (j == 0)
    //         {
    //             sigmad_[i][j] = fdg_ * sigmag_[i];
    //             Nd_[i][j] = sigmad_[i][j] / m_[j];
    //         }
    //         else
    //         {
    //             sigmad_[i][j] = 0.0;
    //             Nd_[i][j] = 0.0;
    //         }
    //         fdg_ratio [i][j] = sigmad_[i][j]/sigmag_[i];
    //         sigmad_tot_[i] += sigmad_[i][j];
    //     }
    // }

    // calculate sigmad_floor ( = 1.0e-50)
    double area;
    for (int i = nr_s; i < nr_e; ++i)
    {
        area = M_PI * (SQR(r_bod_[i + 1]) - SQR(r_bod_[i]));
        for (int j = 0; j < nbin_; ++j)
        {
            sigmad_floor_[i][j] = std :: max(1.0e-50, m_[j] / area);
        }
    }

    CalcCoagulationAndFragmentationArrays();
    // CalcStokesNumberAndDustScaleHeight();
    // CalcDustVelocity();
    int row_count = 0;
    int col_count = 0;
    int row;
    int col;
    
    t_re = t_restart* SOLAR_YEAR;
    std :: ifstream file_dust(dustfile,std :: ios :: in);
    std :: ifstream file_gas(gasfile,std :: ios :: in);
    std :: ifstream file_B(Bfile,std :: ios :: in);
    std :: string line;
    std :: string item;

    if(!file_dust){
        FAILED_TO_OPEN(file_dust,dustfile);
    };
    if(!file_gas){
        FAILED_TO_OPEN(file_gas,gasfile);
    };
    if(!file_B){
        FAILED_TO_OPEN(file_B,Bfile);
    };
    
    // Dust
    // 5行目までは読み取りをスキップする
    while(std :: getline(file_dust,line)){
        if(row_count >= 5){
            row = row_count-5+nr_s;
            if(row >= nr_e){
                break;
            }
            std :: stringstream ss(line);
            double value;
            ss >> value; // r_cen
            for(int j = 0; j < nbin_ && ss >> value; ++j){
                sigmad_active[row][j] = value;
            }
        }  
        row_count += 1;
    }
    
    // Gas
    row_count = 0;
    col_count = 0;
    // r, sigmag, sigmad_tot, vapor, vg_r, vg_phi, vg_wind,
    // vd_r_tot, vd_r, vd_r_diff, fdg, dsigmagdt_w, T, P, beta,
    // dNdr_torque, alpha_GI, alpha_phiz, alpha_tot, eta_P,
    // [Enomoto2026 diag, optional] windDW_zdead, windDW_zsurf, windDW_Delta, windDW_Ameff, alpha_DW
    while(std :: getline(file_gas,line)){
        row = row_count+nr_s;
        if(row >= nr_e){
            break;
        }

        DoubleVector values;
        double value;
        std :: stringstream ss(line);
        while(ss >> value){
            values.push_back(value);
        }

        if(values.size() >= 19){
            sigmag_[row] = values[1];
            sigmad_tot_[row] = values[2];
            sigmad_vapor[row] = values[3];
            vg_r_[row] = values[4];
            vg_phi_[row] = values[5];
            vg_wind_[row] = values[6];

            const int offset = (values.size() >= 20) ? 1 : 0;
            dsigmagdt_w[row] = values[10 + offset];
            T_[row] = values[11 + offset];
            P_[row][nz_s] = values[12 + offset];
            dNdr_torque[row] = values[14 + offset];
            alpha_GI[row] = values[15 + offset];
            alpha_phiz[row] = values[16 + offset];
            alpha_tot[row] = values[17 + offset];
            eta_P_[row] = values[18 + offset];

            // Enomoto (2026) wind-driven-alpha diagnostics : present only in the extended
            // gas output (5 trailing columns). Older restart files simply skip this block;
            // the values are recomputed by CalcWindAlpha() on the first step after restart.
            if(static_cast<int>(values.size()) >= 24 + offset){
                windDW_zdead_[row] = values[19 + offset];
                windDW_zsurf_[row] = values[20 + offset];
                windDW_Delta_[row] = values[21 + offset];
                windDW_Ameff_[row] = values[22 + offset];
                alpha_DW_[row]   = values[23 + offset];
            }
        }
        row_count += 1;
    }
    
    // Magnetic field
    row_count = 0;
    col_count = 0;
    while(std :: getline(file_B,line)){
        row = row_count+nr_s;
        if(row >= nr_e){
            break;
        }

        DoubleVector values;
        double value;
        std :: stringstream ss(line);
        while(ss >> value){
            values.push_back(value);
        }

        if(values.size() >= 5){
            alpha_MRI[row] = values[1];
            B_[row] = values[2];
            etaO_R[row] = values[3];
            etaA_R[row] = values[4];
        }
        row_count += 1;
    }

    for (int i = nr_s; i < nr_e; ++i)
    {
        sigmad_tot_[i] = sigmad_vapor[i];
        sigmad_[i][0] = sigmad_vapor[i] + sigmad_active[i][0];
        for (int j = 0; j < nbin_; ++j)
        {
            if(j > 0){
                sigmad_[i][j] = sigmad_active[i][j];
            }
            Nd_[i][j] = sigmad_active[i][j] / m_[j];
            fdg_ratio [i][j] = sigmad_active[i][j]/sigmag_[i];
            sigmad_tot_[i] += sigmad_active[i][j];
            // std :: cout << sigmad_[i][j] <<" ";
        }
    }

    UpdatePhysicalParameter();
    return;
}


void DustDisk::SetRadiusGrid()
{
    /*
    set disk radius grid
    disk radius grid grid are scaled in logscale
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2

    std::cout << __func__ << std::endl;
    // std ::cout << "Set Radial Grid" << std ::endl;
    double dlnr;
    double dr;
    NG = 2;
    nr_tot = nr_disk_ + 2.0*NG; 
    nr_s = NG;
    nr_e = nr_disk_ + NG;

    r_cen_.resize(nr_tot, 0.0);
    r_bod_.resize(nr_tot + 1, 0.0);

    dr_bod_.resize(nr_tot , 0.0);
    dr_cen_.resize(nr_tot - 1, 0.0);

    r_vol_.resize(nr_tot, 0.0);
    inv_vol_.resize(nr_tot, 0.0);

    dr = (rout_ - rin_)/double(nr_disk_);
    dlnr = std::log(rout_/rin_) / double(nr_disk_);

    // Radial grid spacing, selected by "radial_grid_type" in the input file
    // ("log" default, "linear", "root" = uniform in sqrt(r), or "log_floor").
    if (radial_grid_type_ == "linear")
    {
        for (int i = 0; i < nr_tot+1; ++i)
        {
            r_bod_[i] = rin_ + dr * double(i - nr_s);
        }
    }
    else if (radial_grid_type_ == "root")
    {
        double p = 0.5;
        double pinv = 1.0 / p;
        double r0p = std::pow(rin_/ASTRONOMICAL_UNIT, p);
        double r1p = std::pow(rout_/ASTRONOMICAL_UNIT, p);
        double drp = (r1p - r0p) / double(nr_disk_);

        for (int i = 0; i < nr_tot+1; ++i)
        {
            r_bod_[i] = r0p + drp * double(i - nr_s);
            r_bod_[i] = std::pow(r_bod_[i], pinv) * ASTRONOMICAL_UNIT;
        }
    }
    else if (radial_grid_type_ == "log_floor")
    {
        // Plain log spacing, except cell widths are floored at floor_dr
        // instead of letting log shrink them arbitrarily small near rin_ --
        // this avoids "wasting" a large fraction of the total cell budget on
        // the innermost decade(s), which is what plain "log" does for a
        // large dynamic range. Continuous by construction: the switch back
        // to plain log happens exactly where log's own candidate step first
        // reaches floor_dr, so there is no discontinuity in dr at the
        // transition (unlike a naive two-segment linear+log hybrid).
        //
        // floor_dr is set from "radial_grid_floor_dr" (an explicit input,
        // independent of NG -- how much resolution to spend on the inner
        // disk is a physical choice, not a numerical-stencil detail) when
        // given; otherwise it falls back to rin_/(NG+1), matching the
        // original ghost-cell spacing below (the previously-coupled
        // behavior, kept for backward compatibility when the new input is
        // absent). Ghost cells (i<=nr_s) themselves always use rin_/(NG+1)
        // regardless -- that is a separate, purely numerical concern (stable,
        // monotonic extrapolation to r=0) unrelated to the inner-disk
        // resolution choice.
        double floor_dr = (radial_grid_floor_dr_ > 0.0) ? radial_grid_floor_dr_ : rin_/double(NG+1.0);
        int g_count = 0;
        bool found_first = false;
        double x, r_first = 0.0;
        for (int i = 0; i < nr_tot+1; ++i)
        {
            if (i <= nr_s)
            {
                r_bod_[i] = rin_ * (1.0 + double(i-nr_s)/double(NG+1.0));
            }
            else if (!found_first)
            {
                r_first = r_bod_[i-1];
                dlnr = (std::log(rout_) - std::log(r_first)) / double(nr_disk_ - g_count);
                x = std::log(r_bod_[i-1]) + double(i-nr_s-g_count)*dlnr;
                r_bod_[i] = std::exp(x);

                if (r_bod_[i] - r_bod_[i-1] < floor_dr)
                {
                    r_bod_[i] = r_bod_[i-1] + floor_dr;
                    g_count += 1;
                }
                else
                {
                    found_first = true;
                }
            }
            else
            {
                x = std::log(r_first) + double(i-nr_s-g_count)*dlnr;
                r_bod_[i] = std::exp(x);
            }
        }
    }
    else // "log" (default)
    {
        for (int i = 0; i < nr_tot+1; ++i)
        {
            r_bod_[i] = rin_ * std::exp(dlnr * (i - nr_s));
        }
    }

    for (int i = 0; i < nr_tot; ++i)
    {
        r_cen_[i] = std :: sqrt(r_bod_[i] * r_bod_[i + 1]);
        dr_bod_[i] = r_bod_[i + 1] - r_bod_[i];
        r_vol_[i] = SQR(r_bod_[i + 1]) - SQR(r_bod_[i]);
        inv_vol_[i] = 1.0 / (SQR(r_bod_[i + 1]) - SQR(r_bod_[i]));
    }
    for (int i = 0;i < nr_tot-1; i++)
    {
        dr_cen_[i] = r_cen_[i+1] - r_cen_[i];
    }

    // Geometry-only weights for the vg_r_ de-checkerboarding round trip
    // (centre -> face -> centre linear interpolation). Grid is fixed, so build
    // these once here rather than per step.
    vg_w_cf_.assign(nr_tot + 1, 0.5);
    vg_w_fc_.assign(nr_tot, 0.5);
    for (int i = 1; i < nr_tot; ++i)
        vg_w_cf_[i] = (r_bod_[i] - r_cen_[i-1]) / (r_cen_[i] - r_cen_[i-1]);
    for (int i = 0; i < nr_tot; ++i)
        vg_w_fc_[i] = (r_cen_[i] - r_bod_[i]) / (r_bod_[i+1] - r_bod_[i]);

    SetDiffrentialCoefficient_R();
    SetDiffrentialCoefficient_R_Half();
    SetDiffrentialCoefficientFor2order_R();
    SetZeroDifferentialCoeffient_R();
    return;
}

void DustDisk::SetVerticalGrid()
{
    /*
    set disk vertical grid
    disk vertical grid are scaled in linear
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2

    double dz;
    std::cout << __func__ << std::endl;

    NG = 1;
    nz_s = NG;
    nz_e = NG + nz_disk_;
    nz_tot = nz_disk_ + 2.0*NG; 

    z_cen_.resize(nz_tot, 0.0);
    z_bod_.resize(nz_tot + 1, 0.0);

    dz_bod_.resize(nz_tot , 0.0);
    dz_cen_.resize(nz_tot - 1, 0.0);
    
    // linear scale
    dz = z_max/(nz_disk_-1);
    for (int i = 0; i < nz_tot; ++i)
    {
        z_cen_[i] = (i-nz_s)*dz;
    }

    for(int i=0;i<nz_tot+1;i++)
    {
        // z_cen_ has nz_tot elements (indices 0..nz_tot-1); reading z_cen_[i]
        // first and overwriting only when i==nz_tot (below) reads one past the
        // end of z_cen_ every time this loop reaches i==nz_tot (heap-buffer-
        // overflow, caught by AddressSanitizer) even though that read's result
        // is discarded. Check the boundary case first instead.
        if(i == nz_tot)
        {
            z_bod_[i] = z_cen_[i-1] + 0.5*dz;
        }
        else
        {
            z_bod_[i] = z_cen_[i] - 0.5*dz;
        }
    }

    for(int i=0;i<nz_tot;i++)
    {
        dz_bod_[i] = z_bod_[i+1] - z_bod_[i];
    }
    
    for (int i = 0;i < nz_tot-1; i++)
    {
        dz_cen_[i] = z_cen_[i+1] - z_cen_[i];
    }

    //SetDiffrentialCoefficient_Z();
    return;
}

void DustDisk :: SetDiffrentialCoefficient_R()
{
    double s,t;
    cff_R.resize(nr_tot,DoubleVector(3,0.0));
    for (int i=1 ; i<nr_tot-1;i++)
    {
        s = r_cen_[i] - r_cen_[i-1];
        t = r_cen_[i+1] - r_cen_[i];
        cff_R[i][0] = s*s / (s*t*(s+t));
        cff_R[i][1] = (t*t - s*s) / (s*t*(s+t));
        cff_R[i][2] = -t*t / (s*t*(s+t)); 
    }
}

void DustDisk :: SetDiffrentialCoefficient_R_Half()
{
    double s,t;
    cff_R_Half.resize(nr_tot,DoubleVector(3,0.0));
    for(int i=1;i < nr_tot-1;i++)
    {
        s = r_bod_[i] - r_cen_[i-1];
        t = r_cen_[i] - r_bod_[i];
        cff_R_Half[i][0] = s*s / (s*t*(s+t));
        cff_R_Half[i][1] = (t*t - s*s) / (s*t*(s+t));
        cff_R_Half[i][2] = -t*t / (s*t*(s+t));
    }
}

void DustDisk :: SetDiffrentialCoefficientFor2order_R()
{
    double s,t;
    cff2_R.resize(nr_tot,DoubleVector(3,0.0));
    for (int i=1 ; i<nr_tot-1;i++)
    {
        s = r_cen_[i] - r_cen_[i-1];
        t = r_cen_[i+1] - r_cen_[i];
        cff2_R[i][0] = 2.0 / (t*(s+t));
        cff2_R[i][1] = -2.0 / (s*t);
        cff2_R[i][2] = 2.0 / (s*(s+t)); 
    }
}

void DustDisk :: SetDiffrentialCoefficient_Z()
{
    double s,t;
    cff_Z.resize(nz_tot,DoubleVector(3,0.0));
    for (int i=1 ; i<nz_tot-1;i++)
    {
        s = z_cen_[i] - z_cen_[i-1];
        t = z_cen_[i+1] - z_cen_[i];
        cff_Z[i][0] = s*s / (s*t*(s+t));
        cff_Z[i][1] = (t*t - s*s) / (s*t*(s+t));
        cff_Z[i][2] = -t*t / (s*t*(s+t)); 
    }
}

void DustDisk :: SetZeroDifferentialCoeffient_R()
{
    zero_cff2_R.resize(nr_tot,DoubleVector(3,0.0));

    for(int i=nr_e;i<nr_tot;i++)
    {
        double C = (dr_cen_[i-1]+dr_cen_[i-2])/(dr_cen_[i-2]+dr_cen_[i-3]);
        double C0 = C*( SQR(dr_cen_[i-2])+dr_cen_[i-2]*dr_cen_[i-3]-dr_cen_[i-1]*dr_cen_[i-2] );
        double C1 = (dr_cen_[i-1]/dr_cen_[i-3])*( SQR(dr_cen_[i-2])-dr_cen_[i-2]*dr_cen_[i-3]+dr_cen_[i-1]*dr_cen_[i-2] );
        double C2 = C*(dr_cen_[i-1]/dr_cen_[i-3])*SQR(dr_cen_[i-2]);
        zero_cff2_R[i][0] = C0/SQR(dr_cen_[i-2]);
        zero_cff2_R[i][1] = C1/SQR(dr_cen_[i-2]);
        zero_cff2_R[i][2] = C2/SQR(dr_cen_[i-2]);
    }
}

void DustDisk ::SetBoundaryConditionForDisk()
{
    for(int i=nr_s-1;i>=0;i--)
    {
        // Zero gradient boundary
        sigmag_[i] = sigmag_[nr_s];
        T_[i] = T_[nr_s];
        // // Zero torque boundary
        // sigmag_[i] = sigmag_[nr_s] * pow(r_cen_[i]/r_cen_[nr_s],-1.5);
        // T_[i] = T_[nr_s] * pow(r_cen_[i]/r_cen_[nr_s],-0.5);

        if (sigmag_[i] < 1.0e-4)
        {
            sigmag_[i] = 1.0e-4;
        }

        M_in[i] = M_in[i+1] - M_PI * sigmag_[i] * r_vol_[i];

        if(T_[i] < 10.0)
        {
            T_[i] = 10.0;
        }

        alpha_MRI[i] = alpha_MRI[nr_s];
        alpha_tot[i] = alpha_tot[nr_s];
        Omega_[i] = std::sqrt(GRAVITATIONAL_CONSTANT * M_in[i] / CUB(r_cen_[i]));
        cs_[i] = std::sqrt(BOLTZMANN_CONSTANT * T_[i] / GAS_MOLECULAR_MASS);
        P_integral[i] = sigmag_[i] * SQR(cs_[i]);
        Hg_[i] = cs_[i]/(Omega_[i]);
        rhog_mid_[i] = sigmag_[i]/(std :: sqrt(2*M_PI)*Hg_[i]);
        P_[i][nz_s] = rhog_mid_[i] * SQR(cs_[i]);
    }

    for(int i=nr_e;i<nr_tot;i++)
    {
        // Zero gradient boundary
        sigmag_[i] = sigmag_[nr_e-1];
        T_[i] = T_[nr_e-1];
        // Zero torque boundary
        // sigmag_[i] = sigmag_[nr_e-1] * pow(r_cen_[i]/r_cen_[nr_e-1],-1.5);
        // T_[i] = T_[nr_e-1] * pow(r_cen_[i]/r_cen_[nr_e-1],-0.5);
        //Zero  boundary
        // sigmag_[i] = zero_cff2_R[i][0]*sigmag_[i-1] + zero_cff2_R[i][1]*sigmag_[i-2] - zero_cff2_R[i][2]*sigmag_[i-3];
        // T_[i] = zero_cff2_R[i][0]*T_[i-1] + zero_cff2_R[i][1]*T_[i-2] - zero_cff2_R[i][2]*T_[i-3];

        if (sigmag_[i] < 1.0e-4)
        {
            sigmag_[i] = 1.0e-4;
        }
        M_in[i] = M_in[i-1] + M_PI * sigmag_[i] * r_vol_[i];

        if(T_[i] < 10.0)
        {
            T_[i] = 10.0;
        }

        alpha_MRI[i] = alpha_MRI[nr_e-1];
        alpha_tot[i] = alpha_tot[nr_e-1];
        Omega_[i] = std::sqrt(GRAVITATIONAL_CONSTANT * M_in[i] / CUB(r_cen_[i]));
        cs_[i] = std::sqrt(BOLTZMANN_CONSTANT * T_[i] / GAS_MOLECULAR_MASS);
        P_integral[i] = sigmag_[i] * SQR(cs_[i]);
        Hg_[i] = cs_[i]/(Omega_[i]);
        rhog_mid_[i] = sigmag_[i]/(std :: sqrt(2*M_PI)*Hg_[i]);
        P_[i][nz_s] = rhog_mid_[i] * SQR(cs_[i]);
    }

    for(int i=nr_s-1;i<nr_e+1;i++)
    {
        dP_dr[i] = cff_R[i][0]*P_[i+1][nz_s] + cff_R[i][1]*P_[i][nz_s] + cff_R[i][2]*P_[i-1][nz_s];
        // dPinteg_dr[i] = cff_R[i][0]*P_integral[i+1] + cff_R[i][1]*P_integral[i] + cff_R[i][2]*P_integral[i-1];
        dOmegadr_[i] = -1.5*Omega_[i]/r_cen_[i];
    }
}

void DustDisk :: SetBoundaryConditionForMagneticfield()
{
    double P_tot,lambda,alf_wv2,e;
    double etaO_P;
    double etaA_P;
    for(int i=nr_s;i<nr_e;i++)
    {
        P_tot = 0.0;
        etaO_P = 0.0;
        etaA_P = 0.0;
        for(int j=nz_s;j<nz_e;j++)
        {
            P_tot += P_[i][j]*dz_bod_[j];
            etaO_P += etaO[i][j]*P_[i][j]*dz_bod_[j];
            etaA_P += etaA[i][j]*P_[i][j]*dz_bod_[j];
        }
        etaO_R[i] = etaO_P/P_tot;
        etaA_R[i] = etaA_P/P_tot;
        eta_total[i] = etaA_R[i]+etaO_R[i];

        // // エルサッサー数Lambdaが10^{-5}以下のときはタイムステップをLambdaが10^{-5}として取るようにする //
        // Historically used to keep the old *explicit* scheme's
        // dt_mag ~ dr^2/eta from collapsing when eta is enormous (dead zone,
        // eta up to ~1e20-1e24); caps eta DOWN, which makes the field diffuse
        // *slower* than the raw physical value would (larger eta = faster
        // diffusion). That original numerical motivation no longer applies
        // once use_implicit_mag=true (dt_mag isn't part of the dt formula,
        // and Calc_Mag_Implicit's accuracy-based sub-cycling with convergence
        // retry -- STS_magnetic_diffusion_note.md Sec.9.1/11 -- resolves
        // however large eta actually is), but removing the floor is a
        // physics change (not just a numerical cleanup), independent of
        // use_implicit_mag, and existing results/calibration were produced
        // with it on -- so it is a separate, explicitly opt-in input flag
        // ("use_eta_floor", default true = unchanged historical behavior),
        // not tied to use_implicit_mag.
        //
        // Split from the hard eta_cap below (Sec.12 postmortem, 2026-09-01):
        // with use_eta_floor=false, eta_total also feeds Ohmic/Joule heating
        // in energy.cpp (F_joule ~ eta_tot*(dB/dr)^2), and unbounded eta there
        // made the Initialization() temperature<->resistivity<->alpha_MRI
        // fixed-point iteration fail to converge. eta_cap is therefore always
        // enforced (unconditionally, independent of use_eta_floor) as a hard
        // numerical safety bound on eta_total wherever it's used downstream;
        // only the (softer, more debatable) Elsasser-number floor is gated by
        // use_eta_floor.
        alf_wv2 = 0.25 * SQR(B_[i]) /(M_PI * rhog_mid_[i]);
        lambda = alf_wv2 /(eta_total[i] * Omega_[i]);
        if (use_eta_floor_ && lambda < 1.0e-5)
        {
            eta_total[i] = alf_wv2 * 1.0e5/Omega_[i];
        }
        if (eta_total[i] > eta_cap)
        {
            eta_total[i] = eta_cap;
        }
    }

    for(int i=0;i<nr_s;i++)
    {
        eta_total[i] = eta_total[nr_s];
        B_[i]        = B_[nr_s]; 
    }

    for(int i=nr_e;i<nr_tot;i++)
    {
        // Outer ghost: zero-gradient copy (B=B, eta=eta) of the last active
        // cell. The 3-point extrapolation used before could drift eta above
        // eta_cap or negative, and drift B. With the diffusive flux forced to
        // zero at both domain faces (Calc_Mag / Calc_Mag_Implicit) the ghost
        // eta is unused anyway, and the ghost B only feeds the outflow-clamped
        // advective reconstruction / magnetic pressure.
        eta_total[i] = eta_total[nr_e-1];
        B_[i]        = B_[nr_e-1];
    }

    CalcMagneticFieldPressure();
}

void DustDisk::SetDustGrid()
{
    std ::cout << "SetDustGrid" << std ::endl;
    int i;
    UnsignedInt n_dust_log;
    double amin_ghost, mmin_ghost;
    double eta, zeta, mb;

    a_.resize(nbin_);
    m_.resize(nbin_);
    ab_.resize(nbin_ + 1);
    mb_.resize(nbin_ + 1);

    if (nbin_lin_ == 0)
    {

        eta = std::pow(mmax_ / mmin_, 1.0 / double(nbin_));
        zeta = std::pow(eta, 1.0 / 3.0);
        for (i = 0; i < nbin_ + 1; ++i)
        {
            mb_[i] = mmin_ * std::pow(eta, double(i));
            ab_[i] = amin_ * std::pow(zeta, double(i));
        }

        for (i = 0; i < nbin_; ++i)
        {
            m_[i] = std::sqrt(mb_[i] * mb_[i + 1]);
            a_[i] = std::pow(3.0 * m_[i] / (4.0 * rhodi_ * M_PI), 1.0 / 3.0);
        }
    }
    else
    {

        /*
         */

        n_dust_log = nbin_ - nbin_lin_;
        mmin_ghost = mmin_ / 2.0;

        for (i = 0; i <= nbin_lin_; ++i)
        {
            mb_[i] = mmin_ghost + mmin_ * double(i);
            ab_[i] = std::pow(3.0 * mb_[i] / (4.0 * M_PI * rhodi_), 1.0 / 3.0);
        }

        mb = mb_[nbin_lin_];
        eta = std::pow(mmax_ / mb, 1.0 / double(n_dust_log));
        for (i = 0; i <= n_dust_log; ++i)
        {
            mb_[nbin_lin_ + i] = mb * std::pow(eta, double(i));
            ab_[nbin_lin_ + i] = std::pow(3.0 * mb_[nbin_lin_ + i] / (4.0 * M_PI * rhodi_), 1.0 / 3.0);
        }

        for (i = 0; i < nbin_; ++i)
        {
            m_[i] = std::sqrt(mb_[i] * mb_[i + 1]);
            a_[i] = std::pow(3.0 * m_[i] / (4.0 * M_PI * rhodi_), 1.0 / 3.0);
        }
    }

    return;
}

double DustDisk::CalcSizeForTargetStokes(int i, double St_target)
{
    /*
        Inverts the Epstein/Stokes Stokes-number formula used by
        CalcStokesNumberAndDustScaleHeight() to find the grain size a such
        that St(a, i) = St_target at cell i. Both regimes are monotonically
        increasing in a, so this is a direct closed-form inversion, not a
        root find. Requires ng_mid_[i] and sigmag_[i] to already be set
        (SetDiskGasModel() runs before this is used in Initialization()).
    */
    const double sigma_h2 = 2.0e-15;
    double mfp = 1.0 / (ng_mid_[i] * sigma_h2);
    double a_transition = 2.25 * mfp;
    double St_at_transition = 0.5 * M_PI * rhodi_ * a_transition / sigmag_[i];

    if (St_target <= St_at_transition)
    {
        // Epstein regime: St = 0.5*pi*rhodi_*a/sigmag_
        return 2.0 * St_target * sigmag_[i] / (M_PI * rhodi_);
    }
    else
    {
        // Stokes regime: St = 2*pi*rhodi_*a^2/(9*mfp*sigmag_)
        return std::sqrt(9.0 * St_target * mfp * sigmag_[i] / (2.0 * M_PI * rhodi_));
    }
}

double DustDisk::CalcInitialStokesLimit(int i)
{
    /*
        Stokes number of the largest grain allowed by fragmentation and radial
        drift (Birnstiel et al. 2012), St_lim = min(St_frag, St_drift), used for
        dust_init_amax_mode_ == "frag_drift".

        St_frag  = f_f * v_frag^2 / (3 alpha cs^2)
        St_drift = f_d * eps * (vK/cs)^2 * |dlnP/dlnr|^-1

        v_frag is the same T-blended collision velocity as in dust_growth.cpp,
        alpha is alpha_MRI[i] (the one entering vg_T2_), eps = fdg_.
        dPtot_dr is not set yet at Initialization(), so dlnP/dlnr is taken from
        the midplane gas pressure rhog_mid_*cs^2 by finite differences (central,
        one-sided at the edges); magnetic pressure is ignored.
        Requires SetDiskGasModel() to have run.
    */
    double w = 0.5 * (1.0 - std::tanh((T_[i] - T_snow_) / (0.1 * T_snow_)));
    double v_frag = w * v_col_crit_ + (1.0 - w) * v_col_sili_;
    double St_frag = dust_init_f_frag_ * SQR(v_frag) / (3.0 * alpha_MRI[i] * SQR(cs_[i]));

    int il = (i > nr_s) ? i - 1 : i;
    int ir = (i < nr_e - 1) ? i + 1 : i;
    double P_l = rhog_mid_[il] * SQR(cs_[il]);
    double P_r = rhog_mid_[ir] * SQR(cs_[ir]);
    double dlnP_dlnr = std::log(P_r / P_l) / std::log(r_cen_[ir] / r_cen_[il]);

    double St_drift = std::numeric_limits<double>::max();
    if (std::fabs(dlnP_dlnr) > 0.0)
        St_drift = dust_init_f_drift_ * fdg_ * SQR(vK_[i] / cs_[i]) / std::fabs(dlnP_dlnr);

    return std::min(St_frag, St_drift);
}

void DustDisk::SetInitialDustDistribution()
{
    /*
        Initial dust size distribution (Initialization() only; restarts read
        sigmad_ from file and none of this runs).

        dust_init_mode_ == "monomer" (default): all dust mass in bin 0 (grain
        size amin_) -- the historical behaviour, unchanged.

        dust_init_mode_ == "powerlaw": n(a) da ~ a^-dust_init_powerlaw_q_ from
        the grid minimum size up to a_max(r), same total dust mass
        fdg_*sigmag_[i] as before, just spread across bins instead of packed
        into bin 0. a_max(r) is either a single fixed size
        (dust_init_amax_mode_ == "fixed", dust_init_amax_) or set per radius
        so that St(a_max(r), r) = dust_init_St_target_
        (dust_init_amax_mode_ == "stokes", via CalcSizeForTargetStokes()).
        or so that St(a_max(r), r) = min(St_frag, St_drift)
        (dust_init_amax_mode_ == "frag_drift", via CalcInitialStokesLimit()).
        Falls back to the monomer case wherever a_max(r) does not clear the
        smallest bin edge.
    */
    bool powerlaw = (dust_init_mode_ == "powerlaw");
    bool amax_stokes = (dust_init_amax_mode_ == "stokes");
    bool amax_frag_drift = (dust_init_amax_mode_ == "frag_drift");
    a_max_init_.assign(sigmad_tot_.size(), 0.0);
    double q = dust_init_powerlaw_q_;
    double exp4mq = 4.0 - q; // dM/da ~ a^3 * a^-q = a^(3-q) -> integral exponent 4-q

    for (int i = nr_s; i < nr_e; ++i)
    {
        // reset so this can be called repeatedly (Iteration())
        sigmad_tot_[i] = 0.0;

        double a_max_i = 0.0;
        if (powerlaw)
        {
            if (amax_frag_drift)
                a_max_i = CalcSizeForTargetStokes(i, CalcInitialStokesLimit(i));
            else if (amax_stokes)
                a_max_i = CalcSizeForTargetStokes(i, dust_init_St_target_);
            else
                a_max_i = dust_init_amax_;
        }
        a_max_init_[i] = a_max_i;

        if (!powerlaw || a_max_i <= ab_[0])
        {
            // monomer: all mass in bin 0
            for (int j = 0; j < nbin_; ++j)
            {
                if (j == 0)
                {
                    sigmad_[i][j] = fdg_ * sigmag_[i];
                    Nd_[i][j] = sigmad_[i][j] / m_[j];
                }
                else
                {
                    sigmad_[i][j] = 0.0;
                    Nd_[i][j] = 0.0;
                }
                fdg_ratio[i][j] = sigmad_[i][j] / sigmag_[i];
                sigmad_tot_[i] += sigmad_[i][j];
            }
            continue;
        }

        // Power-law mass weight per bin: integral of dM/da ~ a^(3-q) over
        // [ab_[j], ab_[j+1]], clipped to a_max_i for the bin it falls in;
        // bins entirely above a_max_i get weight 0.
        DoubleVector w(nbin_, 0.0);
        double w_sum = 0.0;
        for (int j = 0; j < nbin_; ++j)
        {
            double lo = ab_[j];
            double hi = std::min(ab_[j + 1], a_max_i);
            if (hi <= lo)
                continue;
            double wj = (std::fabs(exp4mq) < 1.0e-8)
                ? std::log(hi / lo)
                : (std::pow(hi, exp4mq) - std::pow(lo, exp4mq)) / exp4mq;
            w[j] = wj;
            w_sum += wj;
        }

        for (int j = 0; j < nbin_; ++j)
        {
            sigmad_[i][j] = (w_sum > 0.0) ? fdg_ * sigmag_[i] * (w[j] / w_sum) : 0.0;
            Nd_[i][j] = (sigmad_[i][j] > 0.0) ? sigmad_[i][j] / m_[j] : 0.0;
            fdg_ratio[i][j] = sigmad_[i][j] / sigmag_[i];
            sigmad_tot_[i] += sigmad_[i][j];
        }
    }

    return;
}

void DustDisk::SetDiskGasModel()
{
    std ::cout << "SetDiskGasModel" << std ::endl;
    /*
    set disk vertical grid
    disk vertical grid are scaled in linear
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2

    double cs2, B2;
    double zeta_AL = 7.0e-19;
    double zeta_xr;
    double zeta_cr;

    M_total = star_mass_;

    // MMSN model
    for (int i = nr_s; i < nr_e; ++i)
    {
        // sigmag_[i] = sigmag0_ * std::pow(r_cen_[i] / ASTRONOMICAL_UNIT, -1.5)
        //sigmag_[i] = sigmag0_ * std::pow(r_cen_[i] / ASTRONOMICAL_UNIT, -1.5) * std :: exp(-r_cen_[i]/r_cut);
        sigmag_[i] = sigmag0_ * std::pow(r_cen_[i] / r_cut, -1.0) * std :: exp(-r_cen_[i]/r_cut);

        if(sigmag_[i] < 1.0e-04)
        {
            sigmag_[i] = 1.0e-04;
        }

        if (i==nr_s)
        {
            // // zero gradient
            M_in[i] = star_mass_ + M_PI *sigmag_[i] * SQR(r_bod_[i+1]);
            // zero torque
            //M_in[i] = star_mass_ + 4.0*M_PI *sigmag_[i] * (std :: sqrt(r_bod_[nr_s]) - std :: sqrt(r_bod_[0]) )  * std :: pow(ASTRONOMICAL_UNIT,1.5);
        }
        else
        {
            M_in[i] =  M_in[i-1] + M_PI*sigmag_[i]*r_vol_[i];
        }


        if (i==nr_e-1)
        {
            M_total = M_in[i];
        }
        
        //T_[i] = 2.8e2 * std::pow(r_cen_[i]/ASTRONOMICAL_UNIT, -0.5);
        T_[i] = 2.0e1 * std::pow(r_cen_[i]/r_cut, -0.45);

        if(T_[i] < 10.0)
        {
            T_[i] = 10.0;
        }
        // T_irr (-> T_irr4 -> Q_irr heating in energy.cpp) must respect the same
        // 10 K floor, so copy it AFTER the clamp, not before.
        T_irr[i] = T_[i];
        T_irr4[i] = pow(T_irr[i],4.0);

        cs_[i] = std::sqrt(BOLTZMANN_CONSTANT * T_[i] / GAS_MOLECULAR_MASS);
        cs2 = cs_[i] * cs_[i];
        P_integral[i] = sigmag_[i] * cs2;

        Omega_[i] = std::sqrt(GRAVITATIONAL_CONSTANT * M_in[i] / CUB(r_cen_[i]) );
        Hg_[i] = cs_[i] / Omega_[i];
        rhog_mid_[i] = sigmag_[i] / (std::sqrt(2.0 * M_PI) * Hg_[i]);
        ng_mid_[i] = rhog_mid_[i] / GAS_MOLECULAR_MASS;
        vg_phi_[i] = r_cen_[i]*Omega_[i];
        vK_[i] = vg_phi_[i];
        j_[i] = SQR(r_cen_[i]) * Omega_[i];
        Q_toomre[i] = cs_[i]*Omega_[i]/(M_PI*GRAVITATIONAL_CONSTANT*sigmag_[i]);

        B_[i] = sqrt(8.0 * M_PI * cs2 * rhog_mid_[i]/beta_ini); 

        if(B_[i] < B_out)
        {
            B_[i] = B_out;
        }

        B2 = SQR(B_[i]);
        numag_eq[i] = std :: sqrt(0.5*M_PI) * B2/(sigmag_[i]*cs_[i]*Omega_[i]); //Riols et.al. (2019) 3.5
        nu_mag[i] = std :: sqrt(0.5*M_PI) * B2/(sigmag_[i]*cs_[i]*Omega_[i]);  //Riols et.al. (2019) 2.2
        // gas turbulance
        alpha_MRI[i] = alpha_;
        alpha_phiz[i] = alphaDW_ini;
        if(Q_toomre[i] > 3.0)
        {
            alpha_GI[i] = 0.0;
        }
        else
        {
            alpha_GI[i] = std :: exp(-(SQR(Q_toomre[i])*SQR(Q_toomre[i])));
        }
        alpha_tot[i] = alpha_MRI[i] + alpha_GI[i];
        // calc vertical structure
        for (int j = nz_s; j < nz_e; j++)
        {
            rhog[i][j] = rhog_mid_[i] * exp(-0.5 * SQR(z_cen_[j]) );
            ng_[i][j] = rhog[i][j] / GAS_MOLECULAR_MASS;
            P_[i][j] = cs2 * rhog[i][j];
            beta_disk_[i][j] = 8.0 * M_PI * P_[i][j] / B2;
            beta_disk_ini[i][j] =  beta_disk_[i][j];
            zeta_xr = CalcXrayIonizationrate(sigmag_[i], rhog_mid_[i], Hg_[i], r_cen_[i], z_cen_[j]);
            zeta_cr = CalcCosmicRayIonizationrate(sigmag_[i], rhog_mid_[i], Hg_[i], r_cen_[i], z_cen_[j]);
            zeta_[i][j] = zeta_AL + zeta_cr + zeta_xr;
        }
    }
    
    SetBoundaryConditionForDisk();
    return;
}

double DustDisk ::CalcXrayIonizationrate(double sigma, double rho_mid,double H_g,double R, double Z)
{
    double zeta_XR0 = 2.6e-15; 
    // double luminocity_XR=2.0*10e30;
    double chi_XR = 8.0;
    double t = Z / (sqrt(2.0));
    double chi_D = sqrt(0.5 * M_PI) *H_g * rho_mid * erfc(t);

    // Okuzumi 2009 Sano & Turner 2008
    double Zeta_XR = zeta_XR0 * SQR(ASTRONOMICAL_UNIT/R) * (exp(-chi_D / chi_XR) + exp(-(sigma - chi_D) / chi_XR));

    return Zeta_XR;
}

double DustDisk ::CalcCosmicRayIonizationrate(double sigma, double rho_mid,double H_g,double R, double Z)
{
    double chi_CR = 96;
    double t = Z /(sqrt(2));
    double chi_D = sqrt(0.5 * M_PI) * rho_mid * H_g * erfc(t);
    double a1 = chi_D / chi_CR;
    // T. N. Delage et al. 2022 Umebayashi & Nakano 2009
    double Zeta_CR = 0.5 * zeta_CR0_ * (exp(-chi_D / chi_CR) * pow(1.0 + pow(chi_D / chi_CR, 0.75), -4.0/3.0) + exp(-(sigma - chi_D) / chi_CR) * pow(1.0 + pow((sigma - chi_D) / chi_CR, 0.75), -4.0/3.0));

    return Zeta_CR;
}

void DustDisk :: DetermineVaporDust()
{
    if (no_dust_) return; // fdg_<=0: sigmad_/sigmad_active/sigmad_vapor stay 0, nothing to sublimate/condense
    double a1,a2;
    for(int i=0;i<nr_tot;i++)
    {
        // a1 = 0.5*(1.0 - std :: tanh ((T_[i] - 1.35e3)/50.0) );
        // if(a1 > 9.9999e-1){
        //     a1 = 1.0;
        // }
        // else if(a1 < 1.0e-3){
        //     a1 = 0.0;
        // }
        // a2 = 1.0-a1;

        // sigmad_vapor[i] = sigmad_tot_[i] * a2;
        // for (int j=0;j<nbin_;j++)
        // {
        //     sigmad_active[i][j] = sigmad_[i][j] * a1;
        // }

        if(T_[i] > 1.5e3){
            sigmad_vapor[i] = sigmad_tot_[i];
            for (int j=0;j<nbin_;j++)
            {
                sigmad_active[i][j] = 0.0;
            }
        }
        else{
            sigmad_vapor[i] = 0.0;
            for (int j=0;j<nbin_;j++)
            {
                sigmad_active[i][j] = sigmad_[i][j];
            }
        }
    }
}

// CalcCoagulationAndFragmentationArrays() moved to dust_growth.cpp

void DustDisk::CalcStokesNumberAndDustScaleHeight()
{
    /*
        calculate Stokes number and dust scale height
    */
    if (no_dust_) return; // fdg_<=0: St_/Hd_/Nd_ stay 0 (correct -- no dust to have a Stokes number)

    int i, j;
    double mfp; // mean free path of gas molecules
    double sigma_h2 = 2.0e-15;
    double min_st;

    for (i = nr_s; i < nr_e; ++i)
    { // loop for disk radius
        mfp = 1.0 / (ng_mid_[i] * sigma_h2);
        for (j = 0; j < nbin_; ++j)
        { // loop for dust radius
            // calculate stokes number
            if (a_[j] < 2.25 * mfp)
            {
                // Epstein's law
                St_[i][j] = 0.5 * M_PI * rhodi_ * a_[j] / sigmag_[i];
            }
            else
            {
                // Stokes's law
                St_[i][j] = 2.0 * M_PI * rhodi_ * SQR(a_[j]) / (9.0 * mfp * sigmag_[i]);
            }
            // calculate dust scale height (Brauer 2008, Birnsteil 2010)
            min_st = std::min(St_[i][j], 0.5);
            Hd_[i][j] = Hg_[i] * std::min(1.0, std::sqrt(alpha_MRI[i] / (min_st * (1.0 + SQR(St_[i][j])))));
            Nd_[i][j] = sigmad_active[i][j] / (sqrt(2.0*M_PI) * Hd_[i][j] * m_[j]);
        }
    }

    return;
}

void DustDisk :: CalculateViscousTorque()
{
     /*
    set disk vertical grid
    disk vertical grid are scaled in linear
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2
    for(int i=nr_s-1;i<=nr_e;i++)
    {
        N_torque[i] = CUB(r_cen_[i]) * nu_[i] * sigmag_[i] * dOmegadr_[i];
    }
    N_torque[nr_s-1] = N_torque[nr_s];

    for(int i=nr_s;i<nr_e;i++)
    {
        if( std :: fabs(N_torque[i+1]-N_torque[i])/std :: fabs(N_torque[i]) < 1.0e-15 && std :: fabs(N_torque[i]-N_torque[i-1])/std :: fabs(N_torque[i-1]) < 1.0e-15)
        {
            dNdr_torque[i] = 0.0;
        }
        else
        {
            dNdr_torque[i] = cff_R[i][0]*N_torque[i+1] +  cff_R[i][1]*N_torque[i] +  cff_R[i][2]*N_torque[i-1];
        }
        dNdr_torque[i] =  dNdr_torque[i]/(r_cen_[i]*sigmag_[i]) - 0.75*alpha_DW_[i]*SQR(cs_[i]); // 2026 09/28　修正
    } 
}

void DustDisk :: CalculateDust_to_Gas_massratio()
{
    if (no_dust_) return; // fdg_<=0: fdg_ratio/fdg_ratio_vapor stay 0 (correct)

    for (int j = 0; j < nbin_; ++j)
    {
        sigmad_active[nr_s-1][j] =  sigmad_active[nr_s][j];
        sigmad_active[nr_e][j] =  sigmad_active[nr_e-1][j];
    }

    for(int i=nr_s;i<nr_e;i++)
    {
        fdg_ratio_vapor[i] = sigmad_vapor[i]/sigmag_[i];
        for(int k=0;k<nbin_;k++)
        {
            fdg_ratio[i][k] = sigmad_active[i][k]/sigmag_[i];
            if(fdg_ratio[i][k] < 1.0e-20)
            {
                fdg_ratio[i][k] = 0.0;
            }
        }
    }
    
    for (int j = 0; j < nbin_; j++)
    { // loop for dust radius
        fdg_ratio[nr_s-1][j] = fdg_ratio[nr_s][j];
        fdg_ratio[nr_e+1][j] = fdg_ratio[nr_e][j];
    } 
    fdg_ratio_vapor[nr_s-1] = sigmad_vapor[nr_s-1]/sigmag_[nr_s-1];
    fdg_ratio_vapor[nr_e] = sigmad_vapor[nr_e]/sigmag_[nr_e];

    for(int i=nr_s;i < nr_e;i++)
    {
        for(int j=0;j<nbin_;j++)
        {
            dfdg_dr[i][j] = cff_R[i][0]*fdg_ratio[i+1][j] + cff_R[i][1]*fdg_ratio[i][j] + cff_R[i][2]*fdg_ratio[i-1][j];
        }
        dfdg_dr_v[i] = cff_R[i][0]*fdg_ratio_vapor[i+1] + cff_R[i][1]*fdg_ratio_vapor[i] + cff_R[i][2]*fdg_ratio_vapor[i-1];
    }

    for (int j = 0; j < nbin_; j++)
    { // loop for dust radius
        dfdg_dr[nr_s-1][j] = 0.0;
        dfdg_dr[nr_e][j] = 0.0;
    } 
    dfdg_dr_v[nr_s-1] = 0.0;
    dfdg_dr_v[nr_e] = 0.0;
}

void DustDisk::CalculateDiskCoeffcient()
{
    /*
    set disk vertical grid
    disk vertical grid are scaled in linear
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2
    double C1;
    for(int i=nr_s;i<nr_e;i++)
    {
        A_[i] = 2.0*M_PI*SQR(r_cen_[i])*sigmag_[i]/M_in[i];
        C1 = 1.0 + A_[i];
        if(sigmag_[i] > 1.0e-4)
        {
            vg_wind_[i] = r_cen_[i]*dMdt_[i]/M_in[i];
        }
        else
        {
            vg_wind_[i] = 0.0;
        }
        eta_P_[i] = -0.5 * SQR(cs_[i]/vK_[i]) * r_cen_[i] * (dPtot_dr[i]/Ptot_[i]); 

        lambda1_[i] = 0.0;
        lambda2_[i] = 0.0;
        for(int j=0;j<nbin_;j++)
        {
            lambda1_[i] += St_[i][j]*fdg_ratio[i][j]/(1.0+C1*SQR(St_[i][j]));
            lambda2_[i] += fdg_ratio[i][j]/(1.0+C1*SQR(St_[i][j]));
        }
    }
}

void DustDisk::UpdatePhysicalParameter()
{
    /*
    set disk vertical grid
    disk vertical grid are scaled in linear
    */
    //     Ghost cells  　    |      Real cells    |          Ghost cells
    // -----------------------|--------------------|----------------------------------
    //    　0  　1  ...  NG-1  |   NG   ... nr+NG-1 |   nr+NG      ...    ntot-1
    //  -1/2 1/2 ... NG-2/3 | NG-1/2   .....   nr+NG-1/2  | nr+NG+1/2   ...   ntot-1/2
    //    0   1      NG-1       NG      .....   nr+NG+1   | nr+NG+3/2   ...   ntot+1/2
    double s,t;
    double B2,cs2;
    double zeta_AL = 7.0e-19;
    double zeta_xr;
    double zeta_cr;

    M_total = star_mass_;

    // MMSN model
    // Mdisk_calc_pre = 0.0;
    for (int i = nr_s; i < nr_e; ++i)
    {
        if (i==nr_s)
        {
            //zero gradient
            M_in[i] = star_mass_ +  M_PI *sigmag_[i] * SQR(r_bod_[i+1]);
            // zero torque
            //M_in[i] = star_mass_ + 4.0*M_PI *sigmag_[i] * ( 1.0/std :: sqrt(r_bod_[0]) -  1.0/std :: sqrt(r_bod_[nr_s]))  * std :: pow(ASTRONOMICAL_UNIT,1.5);
        }
        else
        {
            M_in[i] =  M_in[i-1] + M_PI *sigmag_[i] * r_vol_[i];
        }

        if (i == nr_e-1)
        {
            M_total = M_in[i];
        }

        if(T_[i] < 10.0)
        {
            T_[i] = 10.0;
        }
        // Was only recomputed inside the T_[i]<10.0 floor branch above, so
        // for every real cell (T_ essentially never drops to the 10K floor)
        // cs_[i] stayed frozen at whatever Initialization()/restart set it
        // to, never tracking T_'s evolution -- unlike the ghost cells,
        // which SetBoundaryConditionForDisk() below does recompute
        // unconditionally every call. Recompute unconditionally here too.
        cs_[i] = std::sqrt(BOLTZMANN_CONSTANT * T_[i] / GAS_MOLECULAR_MASS);
        if(B_[i] < B_out){
            B_[i] = B_out;
        }
        B2 = SQR(B_[i]);
        cs2 = SQR(cs_[i]);
        P_integral[i] = sigmag_[i] * cs2;
        Omega_[i] = std::sqrt(GRAVITATIONAL_CONSTANT * M_in[i] / CUB(r_cen_[i]));
        Hg_[i] = cs_[i] / Omega_[i];
        rhog_mid_[i] = sigmag_[i] / (std::sqrt(2.0 * M_PI) * Hg_[i]);
        // Moved below the rhog_mid_ update above (was reading last step's
        // rhog_mid_ here before it got refreshed a few lines down).
        Ptot_integral[i] = P_integral[i] + frac1_8pi*SQR(B_[i])*sigmag_[i]/rhog_mid_[i]; // P+B^2/8pi * sigmag/rhog_mid
        ng_mid_[i] = rhog_mid_[i] / GAS_MOLECULAR_MASS;
        vK_[i] = r_cen_[i] * Omega_[i];
        j_[i] = r_cen_[i] * vK_[i];
        // alpha parameter
        Q_toomre[i] = cs_[i]*Omega_[i]/(M_PI*GRAVITATIONAL_CONSTANT*sigmag_[i]);
        if(Q_toomre[i] > 3.0)
        {
            alpha_GI[i] = 0.0;
        }
        else
        {
            alpha_GI[i] = std :: exp(-(SQR(Q_toomre[i])*SQR(Q_toomre[i])));
        }
        // calc vertical structure
        for (int j = nz_s; j < nz_e; j++)
        {
            rhog[i][j] = rhog_mid_[i] * exp(-0.5 * SQR(z_cen_[j]) );
            ng_[i][j] = rhog[i][j] / GAS_MOLECULAR_MASS;
            P_[i][j] = cs2 * rhog[i][j];
            beta_disk_[i][j] = 8.0 * M_PI * P_[i][j] / B2;
            zeta_xr = CalcXrayIonizationrate(sigmag_[i], rhog_mid_[i], Hg_[i], r_cen_[i], z_cen_[j]);
            zeta_cr = CalcCosmicRayIonizationrate(sigmag_[i], rhog_mid_[i], Hg_[i], r_cen_[i], z_cen_[j]);
            zeta_[i][j] = zeta_AL + zeta_cr + zeta_xr;
        }
        Ptot_[i] = P_[i][nz_s] + frac1_8pi * SQR(B_[i]);
    }

    SetBoundaryConditionForDisk();
}

// Gaussian smoothing over cell-centred indices [lo, hi) of f, in place.
// "Nearest" edge handling (matches scipy.ndimage.gaussian_filter1d's default
// mode='nearest'), truncated at radius = ceil(4*sigma) (scipy's default
// truncate=4.0). Values outside [lo, hi) (ghost cells) are read via edge
// clamping but never written.
static void GaussianSmooth1D(DoubleVector &f, int lo, int hi, double sigma)
{
    int n = hi - lo;
    if (n <= 1 || sigma <= 0.0) return;
    int radius = (int)std::ceil(4.0 * sigma);
    if (radius < 1) return;

    std::vector<double> kernel(2 * radius + 1);
    double norm = 0.0;
    for (int k = -radius; k <= radius; ++k)
    {
        double w = std::exp(-0.5 * SQR(double(k) / sigma));
        kernel[k + radius] = w;
        norm += w;
    }
    for (double &w : kernel) w /= norm;

    DoubleVector out(n);
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (int k = -radius; k <= radius; ++k)
        {
            int j = i + k;
            if (j < 0) j = 0;
            if (j >= n) j = n - 1;
            sum += kernel[k + radius] * f[lo + j];
        }
        out[i] = sum;
    }
    for (int i = 0; i < n; ++i) f[lo + i] = out[i];
}

void DustDisk::CalcGasVelocity_TM2018()
{
    if (use_face_scheme_)
    {
        if ((int)lnr_cen_.size() != nr_tot) SetupFaceScheme();
        CalcGasFaceVelocity();
        return;
    }
    // #pragma omp parallel for schedule(static,4)
    for(int i=nr_s;i<nr_e;i++)
    {
        if(sigmag_[i] <= 1.0e-4)
        {
            vg_r_[i] = 0.0;
            dvg_phi[i] = 0.0;
            vg_phi_[i] = vK_[i];
        }
        else
        {
            vg_r_[i] = 2.0*dNdr_torque[i]*(lambda2_[i] + 1.0)/vK_[i] + 2.0*lambda1_[i]*(1.0+A_[i])*vK_[i]*eta_P_[i]; 
            vg_r_[i] = vg_r_[i]/( SQR(1.0+lambda2_[i]) + SQR(lambda1_[i]) * (1.0+A_[i]) );
            vg_r_[i] = vg_r_[i] - vg_wind_[i];

            dvg_phi[i] =  dNdr_torque[i]*lambda1_[i]/vK_[i] - (1.0+lambda2_[i])*vK_[i]*eta_P_[i];
            dvg_phi[i] = dvg_phi[i]/( SQR(1.0+lambda2_[i]) + SQR(lambda1_[i]) * (1.0+A_[i]) );
            vg_phi_[i] = vK_[i] + dvg_phi[i];
        }
        //std :: cout << std :: scientific << i <<" "<<vg_r_[i]<<" "<< 2.0*dNdr_torque[i]*(lambda2_[i] + 1.0)/vK_[i]<<" "<<2.0*lambda1_[i]*(1.0+A_[i])*vK_[i]*eta_P_[i]<<std :: endl;
    }
    //exit(1);

    // De-checkerboard vg_r_. It is a pointwise combination of centred radial
    // derivatives (dNdr_torque, dP_dr via eta_P_) whose centre-stencil weight is
    // ~0 on the near-uniform log grid -> weak odd-even coupling, and each
    // differentiation amplifies the grid-scale mode. Project the cell values to
    // faces (a 2-point interpolation has zero response at the Nyquist mode) and
    // interpolate back to centres; the round trip is a positive-centre 3-point
    // low-pass. Both weights are pure geometry (vg_w_cf_/vg_w_fc_ from
    // SetRadiusGrid). Only the cell-centred vg_r_ used downstream (dust drift,
    // dt_adv_gas, energy) is affected -- CalcGasFlux rebuilds its own face
    // velocity, so Sigma_g transport is unchanged.
    {
        DoubleVector vgf(nr_e + 1, 0.0);
        for (int i = nr_s + 1; i < nr_e; ++i)
            vgf[i] = (1.0 - vg_w_cf_[i]) * vg_r_[i-1] + vg_w_cf_[i] * vg_r_[i];
        vgf[nr_s] = vg_r_[nr_s];
        vgf[nr_e] = vg_r_[nr_e-1];
        for (int i = nr_s; i < nr_e; ++i)
            vg_r_[i] = (1.0 - vg_w_fc_[i]) * vgf[i] + vg_w_fc_[i] * vgf[i+1];
    }

    // Same de-checkerboard for dvg_phi. It is built from the same centred radial
    // derivatives (dNdr_torque, dP_dr via eta_P_) as vg_r_ and carries the same
    // grid-scale mode; it feeds the dust radial/azimuthal drift
    // (2*St*dvg_phi in CalcDustVelocity_TM2018) and vg_phi_. Identical
    // centre->face->centre low-pass with the pure-geometry weights, then refresh
    // vg_phi_ = vK_ + dvg_phi.
    {
        DoubleVector vpf(nr_e + 1, 0.0);
        for (int i = nr_s + 1; i < nr_e; ++i)
            vpf[i] = (1.0 - vg_w_cf_[i]) * dvg_phi[i-1] + vg_w_cf_[i] * dvg_phi[i];
        vpf[nr_s] = dvg_phi[nr_s];
        vpf[nr_e] = dvg_phi[nr_e-1];
        for (int i = nr_s; i < nr_e; ++i)
            dvg_phi[i] = (1.0 - vg_w_fc_[i]) * vpf[i] + vg_w_fc_[i] * vpf[i+1];
        for (int i = nr_s; i < nr_e; ++i)
            vg_phi_[i] = vK_[i] + dvg_phi[i];
    }

    // Additional Gaussian smoothing (sigma=2 cells, matches the convention
    // used for velocity plots throughout thesis.ipynb/compare_beta_ini.ipynb)
    // on top of the round-trip de-checkerboard above: that filter only kills
    // the grid-scale (Nyquist) mode, which is not enough once a few cells'
    // dust grows to extreme Stokes number and its drag back-reaction
    // (lambda1_, lambda2_, A_) makes vg_r_/dvg_phi swing over a few-cell
    // width, not just alternate sign cell-to-cell. This wider Gaussian pass
    // damps that too, before it can feed into CalcDustRelativeVelocity and
    // amplify the coagulation kernel.
    GaussianSmooth1D(vg_r_, nr_s, nr_e, 2.0);
    GaussianSmooth1D(dvg_phi, nr_s, nr_e, 2.0);
    for (int i = nr_s; i < nr_e; ++i)
        vg_phi_[i] = vK_[i] + dvg_phi[i];

    vg_r_[nr_s-1] = vg_r_[nr_s];
    vg_r_[nr_e] = vg_r_[nr_e-1];
}

void DustDisk::CalcDustVelocity_TM2018()
{
    /*
        calculate dust velocity
    */
    // Turbulent dust diffusion is no longer carried as a drift velocity here
    // (dividing the Fick velocity -Dd*d ln(Sigma_d/Sigma_g)/dr by a near-zero
    // fdg_ratio blew up at sharp Sigma_d gradients). It is now a separate
    // conservative face flux built in CalcDustFlux() (fluxd_diff_), and
    // vd_r_diff there is only a back-computed diagnostic (F_diff/Sigma_d).
    // vd_r_tot below is therefore the pure drift velocity.
    if (no_dust_) return; // fdg_<=0: vd_r_/vd_phi_/vd_r_tot stay 0 (correct -- no dust to move)

    // sigmad_tot_ ghost cells are not updated elsewhere; keep them synced for
    // downstream consumers (CalcDustFlux diffusion gate, coagulation).
    sigmad_tot_[nr_s-1] = sigmad_tot_[nr_s];
    sigmad_tot_[nr_e]   = sigmad_tot_[nr_e-1];

    for (int i = nr_s; i < nr_e; ++i)
    { // loop for disk radius
        for (int j = 0; j < nbin_; ++j)
        { // loop for dust radius
            if(sigmad_active[i][j] > sigmad_floor_[i][j])
            {
                vd_r_[i][j] = (1.0+A_[i]*SQR(St_[i][j]))*vg_r_[i] + 2.0*St_[i][j]*dvg_phi[i] - SQR(St_[i][j])*vg_wind_[i];
                vd_r_[i][j] = vd_r_[i][j]/(1.0 + (1.0+A_[i])*SQR(St_[i][j]));

                vd_phi_[i][j] = -St_[i][j]*vg_r_[i] + 2.0*dvg_phi[i] - St_[i][j]*vg_wind_[i];
                vd_phi_[i][j] = 0.5*vd_phi_[i][j]/( 1.0+(1.0+A_[i])*SQR(St_[i][j]) ) + vK_[i];

                vd_r_tot[i][j] = vd_r_[i][j];
            }
            else
            {
                vd_r_[i][j] = 0.0;
                vd_phi_[i][j] = vK_[i];
                vd_r_tot[i][j] = 0.0;
            }
        }

        if(sigmad_vapor[i] > sigmad_floor_[i][0])
        {
            vd_r_v[i] = vg_r_[i];
            vd_r_tot_v[i] = vd_r_v[i];
        }
        else
        {
            vd_r_v[i] = 0.0;
            vd_r_tot_v[i] = 0.0;
        }
    }

    for (int j = 0; j < nbin_; ++j)
    { // loop for dust radius
        vd_r_[nr_s-1][j] = vd_r_[nr_s][j];
        vd_r_[nr_e][j] = vd_r_[nr_e-1][j];
        vd_r_tot[nr_s-1][j] = vd_r_tot[nr_s][j];
        vd_r_tot[nr_e][j] = vd_r_tot[nr_e-1][j];
    } 

    vd_r_v[nr_s-1] = vd_r_v[nr_s];
    vd_r_v[nr_e] = vd_r_v[nr_e-1];
    vd_r_tot_v[nr_s-1] = vd_r_tot_v[nr_s];
    vd_r_tot_v[nr_e] = vd_r_tot_v[nr_e-1];
    if (use_face_scheme_) CalcDustFaceVelocity();
    return;
}

void DustDisk::CalcDustVelocity()
{
    /*
        calculate dust velocity
    */

    for (int i = nr_s; i < nr_e; ++i)
    { // loop for disk radius
        for (int j = 0; j < nbin_; ++j)
        { // loop for dust radius
            // radial velocity
            vd_r_[i][j] = -2.0 * vK_[i] * eta_P_[i] * St_[i][j] / (1.0 + SQR(St_[i][j]));
            // azimuthal velocity
            vd_phi_[i][j] = -2.0 * vK_[i] * eta_P_[i] / (1.0 + SQR(St_[i][j]));
        }
    }

    return;
}

void DustDisk::CalcMagneticFieldPressure()
{
    for(int i=nr_s;i<nr_e;i++)
    {
        // Ptot_[i] = frac1_8pi  * SQR(B_[i]) + P_[i][nz_s];
        dB_dr_[i] = B_[i+1]*cff_R[i][0] + B_[i]*cff_R[i][1] + B_[i-1]*cff_R[i][2];
        dPtot_dr[i] = dP_dr[i] + frac1_4pi  * (B_[i]*dB_dr_[i]);
    }
}

// CalcDustRelativeVelocity() moved to dust_growth.cpp

void DustDisk::CalcGasFlux()
{
    /*
        calculate gas flux at cell boundary
        finite volume method
        dSigma/dt + (1/r)d(rF)/dr = 0, F = Sigma * u
    */

    double u_f;
    double s,t;
    double sigma,vr_;
    // 1. Calculate Roe mean
    for(int i = nr_s; i <= nr_e; i++) {
        s = std::sqrt(sigmag_[i]);
        t = std::sqrt(sigmag_[i-1]);
        Roe_mean_[i][0] = s / (s + t);
        Roe_mean_[i][1] = 1.0 - Roe_mean_[i][0];
        sigmag_bod_[i] = std::sqrt(sigmag_[i] * sigmag_[i-1]);
    }
    // 2. Calculate Roe mean surface density and gas velocity
    for (int i = nr_s+1; i < nr_e; ++i) 
    {
        vgr_bod_[i] = Roe_mean_[i][0]*vg_r_[i] + Roe_mean_[i][1]*vg_r_[i-1];
    }
    vgr_bod_[nr_s] = vg_r_[nr_s];
    vgr_bod_[nr_e] = vg_r_[nr_e-1];        

    // 3. Calculate gas flux (avoid conditional branches)
    for (int i = nr_s+1; i < nr_e; ++i) {
        fluxg_[i] = (vgr_bod_[i] >= 0.0) ? vgr_bod_[i] *sigmag_[i-1] : vgr_bod_[i] *sigmag_[i];
    }
    // 4. Process boundary conditions (sequentially)
    fluxg_[nr_s] = (vgr_bod_[nr_s] < 0.0) ? vgr_bod_[nr_s] * sigmag_[nr_s] : 0.0;
    fluxg_[nr_e] = (vgr_bod_[nr_e] > 0.0) ? vgr_bod_[nr_e] * sigmag_[nr_e - 1] : 0.0;
    // Face-flux scheme: Roe_mean_ above is still used by the magnetic solver,
    // but the face gas velocity is the two-point face value (fluxg_ unused).
    if (use_face_scheme_)
        for (int i = nr_s; i <= nr_e; ++i) vgr_bod_[i] = vgf_[i];
    return;
}

void DustDisk::CalcDustFlux()
{
    /*
        calculate dust flux at cell boundary
        finite volume method
        dSigma/dt + (1/r)d(rF)/dr = 0, F = Sigma * u
    */
    if (no_dust_) return; // fdg_<=0: fluxd_/fluxd_vapor/fluxd_diff_/fluxd_vapor_diff_ stay 0 (correct)
    if (use_face_scheme_) { CalcDustFluxFace(); return; }

    double s,t;
    double roe_dust1,roe_dust2;
    double u_f;

    for (int i = nr_s; i <= nr_e; ++i)
    { // loop for disk radius
        // std::max(...,0.0) guard: fdg_ratio_vapor/sigmad_active (below) can
        // carry roundoff-level negative excursions from the advection update
        // (e.g. a cell already at the vaporization/floor boundary picking up
        // a tiny net outflow). sqrt() of a negative silently returns NaN
        // (no exception), which would otherwise propagate through u_f into
        // fluxd_/fluxd_vapor and from there back into sigmad_active.
        s = std :: sqrt(std :: max(fdg_ratio_vapor[i], 0.0));
        t = std :: sqrt(std :: max(fdg_ratio_vapor[i-1], 0.0));
        if(s+t <= 0.0)
        {
            fluxd_vapor[i] = 0.0;
        }
        else{
            roe_dust1 = s/(s+t);
            roe_dust2 = 1.0 -roe_dust1;
            u_f = roe_dust1*vd_r_tot_v[i] + roe_dust2*vd_r_tot_v[i-1];

            if(i == nr_s)
            {
                fluxd_vapor[nr_s] = (u_f < 0.0) ? u_f * sigmad_vapor[nr_s] : 0.0;
            }
            else if(i==nr_e)
            {
                fluxd_vapor[nr_e] = (u_f > 0.0) ? u_f* sigmad_vapor[nr_e - 1] : 0.0;
            }
            else
            {
                if( (sigmad_vapor[i-1] > 0.0 && sigmad_vapor[i] <= 0.0) || (sigmad_vapor[i-1] <= 0.0 && sigmad_vapor[i] > 0.0) ){
                    u_f = roe_dust1*vd_r_v[i] + roe_dust2*vd_r_v[i-1];
                }
                fluxd_vapor[i] = (u_f >= 0.0) ? u_f *sigmad_vapor[i-1] : u_f*sigmad_vapor[i];
            }
        }

        for (int j = 0; j < nbin_; ++j)
        { // loop for dust bin
            if(sigmad_active[i][j] < sigmad_floor_[i][j] && sigmad_active[i-1][j] < sigmad_floor_[i-1][j])
            {
                fluxd_[i][j] = 0.0;
                continue;
            }

            // 2流体のRoe平均
            // std::max(...,0.0) guard: the floor check above only skips this
            // boundary when BOTH sides are below floor, so a cell that is
            // below floor but has picked up a roundoff-level negative value
            // (see ApplyCoagulationAndCorrectMass_FullImplicit) while its
            // neighbor is populated falls through to sqrt() of a negative,
            // which is NaN (silently, no exception) -- guard it here too.
            s = std :: sqrt(std :: max(sigmad_active[i][j], 0.0));
            t = std :: sqrt(std :: max(sigmad_active[i-1][j], 0.0));

            if(s+t <= 0.0)
            {
                fluxd_[i][j] = 0.0;
                continue;
            }

            roe_dust1 = s/(s+t);
            roe_dust2 = 1.0 -roe_dust1;

            u_f = roe_dust1*vd_r_tot[i][j] + roe_dust2*vd_r_tot[i-1][j];

            // Advect the solid dust bins as sigmad_active only. sigmad_ folds
            // vapor into bin 0 (sigmad_[i][0] = sigmad_vapor[i] + sigmad_active
            // [i][0]) for the phase-transition hand-off / output; using it here
            // would transport the vapor mass a second time (it is already
            // advected via fluxd_vapor at the gas velocity), double-counting it
            // in sigmad_tot_. For j>=1, sigmad_[i][j] == sigmad_active[i][j],
            // so this only changes bin 0.
            if(i == nr_s)
            {
                fluxd_[nr_s][j] = (u_f < 0.0) ? u_f * sigmad_active[nr_s][j] : 0.0;
            }
            else if(i==nr_e)
            {
                fluxd_[nr_e][j] = (u_f > 0.0) ? u_f* sigmad_active[nr_e - 1][j] : 0.0;
            }
            else
            {
                if( (sigmad_vapor[i-1] > 0.0 && sigmad_vapor[i] <= 0.0) || (sigmad_vapor[i-1] <= 0.0 && sigmad_vapor[i] > 0.0) ){
                    u_f = roe_dust1*vd_r_[i][j] + roe_dust2*vd_r_[i-1][j];
                }
                fluxd_[i][j] = (u_f >= 0.0) ? u_f *sigmad_active[i-1][j] : u_f*sigmad_active[i][j];
            }
        }
    }

    return;
}

void DustDisk::ApplyDustDiffusionImplicit(double dt)
{
    /*
        Turbulent radial dust diffusion, operator-split, backward-Euler implicit.

        The concentration f = Sigma_d/Sigma_g diffuses with coefficient Dd_ (dust
        bins) or nu_ (vapor, Schmidt number 1). Conservative finite-volume flux
          F_{i+1/2} = -Dd_face * Sigma_g_face * (f_i - f_{i-1}) / (r_cen_i - r_cen_{i-1})
        with F = 0 at both radial domain boundaries (insulating). Solved once per
        macrostep at the full dt (unconditionally stable -> no CFL / sub-cycling).
        Sigma_g, Dd_, nu_ are frozen over the step.

        Raise ABS_THRESH_VDIFF to switch diffusion off in bins whose surface
        density is below a level considered physically meaningful.
    */
    if (no_dust_) return; // fdg_<=0: sigmad_active/sigmad_vapor stay 0, nothing to diffuse

    const double REL_THRESH_VDIFF = 1.0e-20;
    const double ABS_THRESH_VDIFF = 1.0e-20;

    const int N = nr_e - nr_s;
    if (N <= 1) return;

    DoubleVector sub(N), dia(N), sup(N), rhs(N), cp(N), dp(N), sol(N);
    // Face coefficient A_i = r_bod_[i] * Dd_face_i * Sigma_g_face_i / dr_face_i,
    // built for interior faces i = nr_s+1 .. nr_e-1; A at the two domain-boundary
    // faces (nr_s, nr_e) stays 0 -> zero diffusive flux there.
    DoubleVector Aface(N + 1);

    // ---- dust bins ----
    for (int j = 0; j < nbin_; ++j)
    {
        for (int f = 0; f <= N; ++f) Aface[f] = 0.0;
        for (int i = nr_s + 1; i < nr_e; ++i)
        {
            double thr_i   = std::max(REL_THRESH_VDIFF * sigmad_tot_[i],   ABS_THRESH_VDIFF);
            double thr_im1 = std::max(REL_THRESH_VDIFF * sigmad_tot_[i-1], ABS_THRESH_VDIFF);
            if (sigmad_active[i][j] <= thr_i || sigmad_active[i-1][j] <= thr_im1) continue;
            double Dd_face     = 0.5 * (Dd_[i][j] + Dd_[i-1][j]);
            double sigmag_face = 0.5 * (sigmag_[i] + sigmag_[i-1]);
            double dr_face     = r_cen_[i] - r_cen_[i-1];
            Aface[i - nr_s]    = r_bod_[i] * Dd_face * sigmag_face / dr_face;
        }

        for (int i = nr_s; i < nr_e; ++i)
        {
            int k = i - nr_s;
            double g  = 2.0 * dt * inv_vol_[i];   // FV factor 2: inv_vol_ = 1/(r_{i+1}^2 - r_i^2)
            double AL = Aface[k];       // face i   (left)
            double AR = Aface[k + 1];   // face i+1 (right)
            sub[k] = (k == 0)     ? 0.0 : -g * AL / sigmag_[i-1];
            sup[k] = (k == N - 1) ? 0.0 : -g * AR / sigmag_[i+1];
            dia[k] = 1.0 + g * (AL + AR) / sigmag_[i];
            rhs[k] = sigmad_active[i][j];
        }

        // Thomas algorithm (M-matrix -> stable, positivity-preserving).
        cp[0] = sup[0] / dia[0];
        dp[0] = rhs[0] / dia[0];
        for (int k = 1; k < N; ++k)
        {
            double m = dia[k] - sub[k] * cp[k-1];
            cp[k] = sup[k] / m;
            dp[k] = (rhs[k] - sub[k] * dp[k-1]) / m;
        }
        sol[N-1] = dp[N-1];
        for (int k = N - 2; k >= 0; --k) sol[k] = dp[k] - cp[k] * sol[k+1];

        for (int i = nr_s; i < nr_e; ++i)
        {
            int k = i - nr_s;
            double v = sol[k];
            if (v < 0.0) v = 0.0;               // roundoff guard
            sigmad_active[i][j] = v;
        }

        // Diagnostic: bare face diffusion flux F_i (matches fluxd_ convention,
        // where the r_bod factor is applied separately) from the solved
        // concentration, and the cell-centred effective velocity F_diff/Sigma_d.
        for (int i = nr_s; i < nr_e; ++i)
        {
            int k = i - nr_s;
            fluxd_diff_[i][j] = (Aface[k] > 0.0)
                ? -(Aface[k] / r_bod_[i]) * (sigmad_active[i][j] / sigmag_[i] - sigmad_active[i-1][j] / sigmag_[i-1])
                : 0.0;
        }
        fluxd_diff_[nr_e][j] = 0.0;
        for (int i = nr_s; i < nr_e; ++i)
        {
            vd_r_diff[i][j] = (sigmad_active[i][j] > sigmad_floor_[i][j])
                ? 0.5 * (fluxd_diff_[i][j] + fluxd_diff_[i+1][j]) / sigmad_active[i][j]
                : 0.0;
        }
    }

    // ---- vapor (diffusivity nu_, Schmidt number 1) ----
    for (int f = 0; f <= N; ++f) Aface[f] = 0.0;
    for (int i = nr_s + 1; i < nr_e; ++i)
    {
        if (sigmad_vapor[i] <= ABS_THRESH_VDIFF || sigmad_vapor[i-1] <= ABS_THRESH_VDIFF) continue;
        double nu_face     = 0.5 * (nu_[i] + nu_[i-1]);
        double sigmag_face = 0.5 * (sigmag_[i] + sigmag_[i-1]);
        double dr_face     = r_cen_[i] - r_cen_[i-1];
        Aface[i - nr_s]    = r_bod_[i] * nu_face * sigmag_face / dr_face;
    }
    bool vapor_active = false;
    for (int i = nr_s; i < nr_e; ++i)
    {
        int k = i - nr_s;
        double g  = 2.0 * dt * inv_vol_[i];
        double AL = Aface[k];
        double AR = Aface[k + 1];
        if (AL > 0.0 || AR > 0.0) vapor_active = true;
        sub[k] = (k == 0)     ? 0.0 : -g * AL / sigmag_[i-1];
        sup[k] = (k == N - 1) ? 0.0 : -g * AR / sigmag_[i+1];
        dia[k] = 1.0 + g * (AL + AR) / sigmag_[i];
        rhs[k] = sigmad_vapor[i];
    }
    if (vapor_active)
    {
        cp[0] = sup[0] / dia[0];
        dp[0] = rhs[0] / dia[0];
        for (int k = 1; k < N; ++k)
        {
            double m = dia[k] - sub[k] * cp[k-1];
            cp[k] = sup[k] / m;
            dp[k] = (rhs[k] - sub[k] * dp[k-1]) / m;
        }
        sol[N-1] = dp[N-1];
        for (int k = N - 2; k >= 0; --k) sol[k] = dp[k] - cp[k] * sol[k+1];
        for (int i = nr_s; i < nr_e; ++i)
        {
            double v = sol[i - nr_s];
            sigmad_vapor[i] = (v < 0.0) ? 0.0 : v;
        }
    }
    for (int i = nr_s + 1; i < nr_e; ++i)
    {
        int k = i - nr_s;
        fluxd_vapor_diff_[i] = (Aface[k] > 0.0)
            ? -(Aface[k] / r_bod_[i]) * (sigmad_vapor[i] / sigmag_[i] - sigmad_vapor[i-1] / sigmag_[i-1])
            : 0.0;
    }
    fluxd_vapor_diff_[nr_s] = 0.0;
    fluxd_vapor_diff_[nr_e] = 0.0;
    for (int i = nr_s; i < nr_e; ++i)
    {
        vd_r_diff_v[i] = (sigmad_vapor[i] > sigmad_floor_[i][0])
            ? 0.5 * (fluxd_vapor_diff_[i] + fluxd_vapor_diff_[i+1]) / sigmad_vapor[i]
            : 0.0;
    }

    return;
}

void DustDisk::CalcWindMassLoss_Tabone2022()
{
    double A;
    double Vg_tot_,V_A;
    double V_max,V_min;
    // DoubleVector r_A; // Alfven radius
    // DoubleVector lever_;

    for(int i=nr_s;i<nr_e;i++)
    {
        if(sigmag_[i] > 1.0e-4)
        {
            dsigmagdt_w[i] = alpha_phiz[i]*rhog_mid_[i]*SQR(cs_[i])/(vK_[i])*inv_lambda_w_m1_; // Tabone et al .(2022) , 1/(lambda-1)
            A = dsigmagdt_w[i]/sigmag_[i];
            v_dsigmaddt_w[i] = (sigmad_vapor[i] > 1.0e-50) ? A*sigmad_vapor[i] : 0.0; // Tabone et al .(2022)
            dMdt_[i]       =  dMdt_[i-1] - M_PI*r_vol_[i]*dsigmagdt_w[i];
            for(int j=0;j<nbin_;j++)
            {
                if( sigmad_active[i][j] < sigmad_floor_[i][j] || St_[i][j] > 1.0e-4)
                {
                    dsigmaddt_w[i][j] = 0.0;
                }
                else
                {
                    dsigmaddt_w[i][j] = sigmad_active[i][j] * A;
                }
            }
        }
        else
        {
            dsigmagdt_w[i] = 0.0;
            dMdt_[i]       =  dMdt_[i-1];
            for(int j=0;j<nbin_;j++)
            {
                dsigmaddt_w[i][j] = 0.0;
            }
        }
    }
}

void DustDisk::CalcWindMassLoss_Riols2019()
{
    double A,B;
    for(int i=nr_s;i<nr_e;i++)
    {
        if(sigmag_[i] > 1.0e-4)
        {
            nu_mag[i]      = std :: sqrt(0.5*M_PI) * SQR(B_[i])/(sigmag_[i]*cs_[i]*Omega_[i]);
            B = nu_mag[i]/(numag_eq[i]);
            if (B > inv_Cw_)
            {
                B = inv_Cw_;
            }
            A = Omega_[i]*Cw_*B; // Riols et al. (2019) 3.4
            dsigmagdt_w[i] = sigmag_[i] * A; // Riols et al. (2019) 3.4
            dMdt_[i]       =  dMdt_[i-1] - M_PI*r_vol_[i]*dsigmagdt_w[i];
            for(int j=0;j<nbin_;j++)
            {
                if( sigmad_active[i][j] < sigmad_floor_[i][j] || St_[i][j] > 1.0e-4)
                {
                    dsigmaddt_w[i][j] = 0.0;
                }
                else
                {
                    dsigmaddt_w[i][j] = sigmad_active[i][j] * A;
                }
            }
        }
        else
        {
            dsigmagdt_w[i] = 0.0;
            dMdt_[i]       =  dMdt_[i-1];
            for(int j=0;j<nbin_;j++)
            {
                dsigmaddt_w[i][j] = 0.0;
            }
        }
    }
}

void DustDisk::CalcWindTorque_Suzuki2016()
{
    double C_DW = 0.531923041; // C_DW=4/3 * 1/sqrt(2pi)
    for(int i=nr_s;i<nr_e;i++)
    {
        alpha_phiz[i] = std :: min( 1.0 , alphaDW_ini * std:: pow(beta_disk_[i][nz_s]/beta_disk_ini[i][nz_s] ,-0.66) );
        alpha_DW_[i] = C_DW*(r_cen_[i]/Hg_[i])*alpha_phiz[i];
    }
}

void DustDisk::CalcWindAlpha()
{
    // dispatch : 0 = Suzuki2016 (default) , 1 = Enomoto et al. (2026) self-consistent alpha_DW
    if (wind_alpha_model_ == 1)
    {
        CalcWindDrivenAlpha_Enomoto();
    }
    else
    {
        CalcWindTorque_Suzuki2016();
    }
}

void DustDisk::CalcWindDrivenAlpha_Enomoto()
{
    // Enomoto, Mori & Okuzumi (2026, PASJ) scaling law for wind-driven accretion.
    //   |Bphi/Bz|_surf ~ k beta0^0.27 Ameff^0.48 Delta^1.2            (eq.22, k=0.84)
    //   alpha_DW       ~ 2.1e-2 (beta0/1e4)^-0.73 (eps/0.1)^-1
    //                       * Ameff^0.48 Delta^1.2   , eps = Hg/r     (eq.27)
    // code convention : alpha_phiz = 0.75*sqrt(2pi)*(Hg/r) * alpha_DW
    const double four_pi = 4.0 * M_PI;
    const double pref    = 0.75 * std::sqrt(2.0 * M_PI);   // geometric conversion factor

    for (int i = nr_s; i < nr_e; ++i)
    {
        if (sigmag_[i] <= 1.0e-4)
        {
            alpha_phiz[i]    = 0.0;
            windDW_zsurf_[i] = 0.0;
            windDW_zdead_[i] = 0.0;
            windDW_Delta_[i] = 0.0;
            windDW_Ameff_[i] = 0.0;
            alpha_DW_[i]   = 0.0;
            continue;
        }

        const double beta0 = beta_disk_[i][nz_s];   // midplane plasma beta
        const double B2    = SQR(B_[i]);

        // --- z_surf : beta(z) crosses 1, scanned from the midplane upward (linear interp) ---
        double z_surf = z_max;                       // fallback : no crossing within the grid
        if (beta_disk_[i][nz_s] < 1.0)
        {
            z_surf = 0.0;                            // magnetic pressure dominates at midplane
        }
        else
        {
            for (int k = nz_s; k < nz_e - 1; ++k)
            {
                double b0 = beta_disk_[i][k]     - 1.0;
                double b1 = beta_disk_[i][k + 1] - 1.0;
                if (b0 == 0.0) { z_surf = z_cen_[k]; break; }
                if (b0 * b1 < 0.0)
                {
                    z_surf = z_cen_[k] + (z_cen_[k + 1] - z_cen_[k]) * (-b0) / (b1 - b0);
                    break;
                }
            }
        }

        // --- z_dead : Ohmic Elsasser number Lambda_O(z) crosses 1 ---
        //     Lambda_O = vA^2 / (etaO * Omega) , vA^2 = B^2 / (4 pi rho)
        double z_dead = z_max;                       // fallback : coupled all the way to z_max
        double LO_mid = B2 / (four_pi * rhog[i][nz_s] * etaO[i][nz_s] * Omega_[i]);
        if (LO_mid > 1.0)
        {
            z_dead = 0.0;                            // no Ohmic dead zone at the midplane
        }
        else
        {
            for (int k = nz_s; k < nz_e - 1; ++k)
            {
                double l0 = B2 / (four_pi * rhog[i][k]     * etaO[i][k]     * Omega_[i]) - 1.0;
                double l1 = B2 / (four_pi * rhog[i][k + 1] * etaO[i][k + 1] * Omega_[i]) - 1.0;
                if (l0 * l1 < 0.0)
                {
                    z_dead = z_cen_[k] + (z_cen_[k + 1] - z_cen_[k]) * (-l0) / (l1 - l0);
                    break;
                }
            }
        }

        // --- Delta : normalized active-layer thickness (eq.21) ---
        double Delta = (z_surf > 0.0 && z_dead < z_surf) ? (z_surf - z_dead) / z_surf : 0.0;

        // --- Ameff : <log10 Am> over z_dead < |z| < z_surf, only when Delta > 0 ---
        double Am_eff     = 0.0;
        double alpha_DW = 0.0;
        if (Delta > 0.0)
        {
            double sum = 0.0;
            int    cnt = 0;
            for (int k = nz_s; k < nz_e; ++k)
            {
                if (z_cen_[k] > z_dead && z_cen_[k] < z_surf)
                {
                    double Am = B2 / (four_pi * rhog[i][k] * etaA[i][k] * Omega_[i]);
                    if (Am > 0.0) { sum += std::log10(Am); ++cnt; }
                }
            }
            if (cnt > 0)
            {
                Am_eff = std::pow(10.0, sum / double(cnt));
                alpha_DW = 2.1e-2
                           * std::pow(beta0 / 1.0e4, -0.73)
                           * (0.1 * r_cen_[i] / Hg_[i])       // (eps/0.1)^-1 , eps = Hg/r
                           * std::pow(Am_eff, 0.48)
                           * std::pow(Delta, 1.2);
            }
        }

        double aw = pref * (Hg_[i] / r_cen_[i]) * alpha_DW; // 0.75*sqrt(2pi)*(Hg/r)*alpha_DW
        alpha_phiz[i] = std::min(1.0, std::max(0.0, aw));

        windDW_zsurf_[i] = z_surf;
        windDW_zdead_[i] = z_dead;
        windDW_Delta_[i] = Delta;
        windDW_Ameff_[i] = Am_eff;
        alpha_DW_[i]   = alpha_DW;
    }

    // ghost cells (see CalculationMRIturbulence convention for alpha_)
    alpha_phiz[nr_s - 1] = alpha_phiz[nr_s];
    alpha_phiz[nr_e]     = alpha_phiz[nr_e - 1];
}


void DustDisk::CalculationRun()
{
    if (on_growth_)
    {
        if (on_fragmentation_)
        {
            UpdateDustDiskFragOn();
        }
        else
        {
            UpdateDustDisk();
        }
    }
    else
    {
        UpdateNoGrowthDustDisk();
    }
    return;
}

void DustDisk ::UpdateResistivity()
{
    // Calculate resistivity (all temperatures) : alpha.cpp (CalcAlphaParameter)
    //   ionization equilibrium with grain charging (GetChargestate)
    //   + thermal ionization of potassium by the Saha equation for T > 1000 K
    //     (CalcThermalIonization), eta_O / eta_A from CalcResistivity_ForThermal

    int i, j;

    double alpha_P, P_tot, Hg2, z2;
    double gosa_p, Nd_tot;
    double min_st, gosa;
    int count;

    # pragma omp parallel
    {
        CalcAlphaParameter alpha_solver;

        DoubleMatrix Nd_RZ_local;
        Nd_RZ_local.resize(nr_tot,DoubleVector(nbin_,0.0));
        #pragma omp for schedule(dynamic,8)
        for(int i=nr_s;i<nr_e;i++)
        {
            for (int k = nz_s; k < nz_e; k++)
            {
                // Update Nd_
                double Nd_tot = 0.0;
                for (int j = 0; j < nbin_; j++)
                {
                    Nd_RZ_local[k][j] = Nd_[i][j] * exp(-0.5*SQR(z_cen_[k]*Hg_[i] / Hd_[i][j]));
                    Nd_tot += Nd_RZ_local[k][j];
                }

                if(Nd_tot < 1.0e-50 && sigmad_tot_[i] < sigmag_[i])
                {
                    // No dust here: pass true zero density through. alpha_solver
                    // (CalcAlphaParameter::GetChargestate()) now handles n_d < 1e-50
                    // directly (gas-phase-only ionization/recombination equilibrium),
                    // so faking a trace bin-0 population to dodge its old 0/0 divide
                    // is no longer needed.
                    for (int j = 0; j < nbin_; j++)
                    {
                        Nd_RZ_local[k][j] = 0.0;
                    }
                }
                else if (sigmad_tot_[i] > sigmag_[i])
                {
                    for (int j = 0; j < nbin_; j++)
                    {
                        Nd_RZ_local[k][j] = (j==0) ? rhog[i][k]/m_[0] : 0.0;
                    }
                }

                Resistivity res;
                alpha_solver.SetParameter(zeta_[i][k], T_[i], ng_[i][k], Nd_RZ_local[k], m_, a_);
                alpha_solver.SetParameterForAlpha(beta_disk_[i][k], B_[i], Omega_[i]);
                alpha_solver.GetChargestate(phi_prev_[i][k]);
                phi_prev_[i][k] = alpha_solver.GetPhi();
                alpha_solver.CalcResistivity_ForThermal(&res);
                etaO[i][k] = res._eta_O_;
                etaA[i][k] = res._eta_A_;
            }
        }
    }
    SetBoundaryConditionForMagneticfield();
}

void DustDisk ::UpdateTurbulent()
{
    double nu_turb;
    for (int i = nr_s-1; i <= nr_e; ++i)
    {
        // gas turbulance
        alpha_tot[i] = alpha_MRI[i] + alpha_GI[i];
        nu_[i] = alpha_tot[i] * cs_[i] * Hg_[i];
        vg_T2_[i] = alpha_MRI[i] * SQR(cs_[i]);
        Re_[i] = 0.5 * alpha_MRI[i] * sigmag_[i] * H2_CROSS_SECTION / GAS_MOLECULAR_MASS;
        inv_sqrt_Re_[i] = 1.0 / std::sqrt(Re_[i]);

        nu_turb = alpha_MRI[i] * cs_[i] * Hg_[i];
        for(int j=0;j<nbin_;j++)
        {
            Dd_[i][j] = nu_turb/(1.0 + SQR(St_[i][j]) );
        }
    }
    for(int j=0;j<nbin_;j++)
    {
        Dd_[nr_s-1][j] = Dd_[nr_s][j];
        Dd_[nr_e][j] = Dd_[nr_e-1][j];
    }
}

void DustDisk::UpdateDustDisk()
{
    /*
        main calculate function in the case of only coagulation
        (on_growth = true, on_fragmentation = false).

        *** 工事中 (under construction) ***
        The old explicit-coagulation loop that used to live here was outdated
        (no restart support, T_ never updated after Initialization(), no
        DetermineVaporDust() call, CalcWindMassLoss_Riols2019() instead of
        _Tabone2022(), no implicit magnetic/coagulation solves) and did not
        compute properly, so it has been removed.

        Plan: share the time loop of UpdateDustDiskFragOn(). The collision
        kernels (CoagulateChangeRateFragmentation*() in dust_growth.cpp)
        already treat !on_fragmentation_ as perfect sticking (kij += fij:
        high-velocity collisions also coagulate, the fragmentation channel is
        skipped), but this path has not been validated yet, so for now the
        run stops here with an error.
        Use on_fragmentation = true (UpdateDustDiskFragOn) or
        on_growth = false (UpdateNoGrowthDustDisk) instead.
    */
    std :: cerr << "Error: UpdateDustDisk (coagulation only, fragmentation = false) is under construction (工事中)." << std :: endl;
    std :: cerr << "       Use fragmentation = true or growth = false." << std :: endl;
    std :: exit(1);
}

void DustDisk::UpdateDustDiskFragOn()
{
    /*
        main calculate function in the case of coagulation and fragmentation
    */

    std::cout << "UpdateDustDiskFargOn" << std::endl;

    int i, j, n, tcount, icount, ndisk, nbin , time_count , alpha_count;
    double gosa = 1.0e6;
    double t, time_end,DeltaT,dt_tot ,dt, dt_coag,dt_Esrc, dt_wind ,dt_adv, dt_adv_gas, dt_gain_global;
    double vol_inv,  local_sigma_tot, temp_sigma_coag, prev_sigma_tot, err, delta_sigma, cfl;
    double etaO_P = 0.0;
    double etaA_P = 0.0;
    double Ptotal_ = 0.0;

    DoubleMatrix dNdt_rd(nr_tot, DoubleVector(nbin_, 0.0));
    DoubleMatrix dNdt_coag(nr_tot, DoubleVector(nbin_, 0.0));
    // Only allocated/used when use_implicit_coag_ is true.
    DoubleMatrix Gain_coag(nr_tot, DoubleVector(nbin_, 0.0));
    DoubleMatrix LossRate_coag(nr_tot, DoubleVector(nbin_, 0.0));
    // (use_full_implicit_coag_ scratch is now allocated per-thread inside the
    // parallel coagulation-apply loop, not here.)
    DoubleVector alpha_pre(nr_tot,0.0);
    DoubleVector alpha_invdalpha(nr_tot,0.0);
    DoubleVector dt_alpha_v(nr_tot,0.0);
    DoubleVector gosa_R(nr_tot,0.0);
    DoubleVector outtime(1, 0.0);
    std::string filename;
    std::string filename_mag;

    Mag_field mag;
    Resistivity res;
    Energy_Calc eng;

    // Grid geometry (r_cen_, r_bod_, dr_cen_, dr_bod_, inv_vol_, z_cen_, z_bod_,
    // dz_cen_, dz_bod_, cff_R, cff_R_Half) is fixed for the whole run once
    // Initialization() has set it up, so these only need to be pushed into `mag`
    // once here instead of every substep inside the while(dt_tot < DeltaT) loop below.
    SetRadiusGridParameter(nr_s,nr_e,r_cen_,r_bod_,dr_cen_,dr_bod_,inv_vol_);
    SetVerticalGridParameter(nz_s,nz_e,z_cen_,z_bod_,dz_cen_,dz_bod_);
    Set_Differential_Coefficient(cff_R,cff_R_Half);
    eng.SetGridParameter(nr_s,nr_e,nr_tot);
    eng.SetRadDiffFactor(rad_diff_fac_);

    ndisk = nr_disk_;
    nbin = nbin_;

    outtime[0] = tend_;
    time_end = outtime.back();

    n = 0;
    time_count = 1;
    DeltaT = SOLAR_YEAR;
    tcount = 0;
    icount = 0;
    t = 0.0;
    cfl = cfl_;
    
    if(restart_ == false){
        Iteration();
        filename = outdir_ + "/out0";
        filename_mag = outdir_mag + "/out0";
        DetermineVaporDust();
        OutPut(filename,filename_mag);
        std::cout << icount << ", time (yr) = " << (t / SOLAR_YEAR)  <<" ,M_disk = "<<(M_total - star_mass_)/star_mass_<<" M."<< std::endl;
    }
    else{
        t = t_re;
        time_count = int(t/SOLAR_YEAR);
        DetermineVaporDust();
        UpdatePhysicalParameter();
        CalcStokesNumberAndDustScaleHeight();
        UpdateResistivity();
        CalcWindAlpha();
        UpdateTurbulent();
        CalculateViscousTorque();
        CalculateDust_to_Gas_massratio();
        CalcWindMassLoss_Tabone2022();

        CalculateDiskCoeffcient();
        CalcGasVelocity_TM2018();
        CalcDustVelocity_TM2018();

        CalcGasFlux();
        CalcDustFlux();
        filename = outdir_ + "/restart" + std::to_string(time_count);
        filename_mag = outdir_mag + "/restart" + std::to_string(time_count);
        OutPut(filename,filename_mag);
        std::cout << time_count << ", time (yr) = " << (t / SOLAR_YEAR)  << ": Restart !! "<< std::endl;
    }
    
    // calculation
    while(t < time_end)
    {
        /* code */
        dt_tot = 0.0;
        if(time_count < 1.0e4){DeltaT = 5.0e2*SOLAR_YEAR;}
        else if(time_count < 1.0e5){DeltaT = 5.0e3*SOLAR_YEAR;}
        else{DeltaT = 5.0e4*SOLAR_YEAR;}

        if(t+DeltaT > time_end){DeltaT = time_end-t;}

        while (dt_tot < DeltaT)
        {
            // calculate dSigma/dt
            UpdatePhysicalParameter();
            DetermineVaporDust();
            CalcStokesNumberAndDustScaleHeight();
            UpdateResistivity();
            CalcWindAlpha();
            UpdateTurbulent();
            CalculateViscousTorque();
            CalculateDust_to_Gas_massratio();
            CalcWindMassLoss_Tabone2022();

            CalculateDiskCoeffcient();
            CalcGasVelocity_TM2018();
            CalcDustVelocity_TM2018();

            CalcDustFlux();
            CalcGasFlux();

            // 成長・破砕ソース項の計算
            // use_full_implicit_coag_: nothing to precompute here -- Gain/Loss/
            // Jacobian are built and solved together per cell in the apply step
            // below (ApplyCoagulationAndCorrectMass_FullImplicit), since that
            // path needs no dt-side Gain estimate (dt_coag/dt_gain_global are
            // both unconstrained for it).
            if (no_dust_)
            {
                // fdg_<=0: nothing to coagulate/fragment.
            }
            else if (use_full_implicit_coag_)
            {
            }
            else if (use_implicit_coag_)
            {
                CoagulateChangeRateFragmentation_Implicit(t, sigmad_active, Gain_coag, LossRate_coag);
            }
            else
            {
                CoagulateChangeRateFragmentation(t, sigmad_active, dNdt_coag);
            }

            // determine dt
            t_mag = 1.0e100;
            dt_Esrc = 1.0e100;
            dt_adv_gas = 1.0e100;
            dt_wind = 1.0e100;
            dt_adv = 1.0e100;
            dt_coag = 1.0e100;
            // determine dt for magnetic field
            // (grid geometry already pushed into `mag` once, before the time loop)
            SetDiskParameter(ng_,vgr_bod_,Roe_mean_,P_,Omega_,Hg_); // vgrtestはガスの速度が0として計算
            SetMagneticFieldParameter(etaO,etaA,eta_total,B_,B_out);
            CalculateDeltaTime(&mag);
            // With the implicit magnetic diffusion solve enabled, Calc_Mag
            // advances diffusion stably for any dt (Thomas-algorithm solve),
            // so mag.Delta_time (the single-explicit-stage limit) must not
            // cap the macrostep here.
            t_mag = IsImplicitMagneticDiffusionEnabled() ? 1.0e100 : mag.Delta_time;

            // determine Temepature
            if(on_thermal_ == true){
                eng.SetDiskParameter(r_cen_,sigmag_,T_,T_irr,T_irr4,Omega_,dsigmagdt_w,alpha_MRI,alpha_phiz);
                eng.SetMagneticField(B_,dB_dr_,eta_total);
                dt_Esrc = eng.Calculate_Dt(vg_r_,r_bod_,inv_vol_);
                eng.GetTemperature(T_);
                // eng.GetTemperature_Initialization(T_);
            }

            if (!no_dust_)
            {
                for (int i = nr_s; i < nr_e; ++i)
                { // loop for disk radius
                    vol_inv = inv_vol_[i];
                    for (int j = 0; j < nbin; ++j)
                    { // loop for dust bin
                        // advection term
                        dNdt_rd[i][j] = -2.0 * (r_bod_[i + 1]*fluxd_[i + 1][j] - r_bod_[i] * fluxd_[i][j])*vol_inv - dsigmaddt_w[i][j];
                        if (sigmad_active[i][j] > sigmad_floor_[i][j] && std::abs(vd_r_tot[i][j]) != 0.0)
                        {
                            dt_adv = std::min(dt_adv, std::abs(dr_bod_[i] / vd_r_tot[i][j]));
                        }
                    }
                }
            }

            // for(int i=nr_s;i<nr_e;i++){
            //     std :: cout <<std :: scientific<< sigmad_tot_[i]<<" "<<sigmad_[i][0]<<" "<<sigmad_vapor[i]<<std :: endl;
            // }
            // exit(1);

            for (int i = nr_s; i < nr_e; ++i)
            { // loop for disk radius
                dt_adv_gas = std::min( dt_adv_gas, std::abs(dr_bod_[i]/vg_r_[i]) );
                dt_wind = std::min(dt_wind, std::abs( sigmag_[i]/dsigmagdt_w[i] ));
            }
            if (use_face_scheme_)
            {
                // explicit limits from the face velocities (diffusion is implicit)
                dt_adv_gas = GasAdvectionDtFace();
                if (!no_dust_) dt_adv = DustAdvectionDtFace();
            }

            // The semi-implicit update (y_new = (y+dt*Gain)/(1+dt*LossRate)) is
            // unconditionally non-negative, so the depletion-timescale limit
            // EstimateDtCoagulation exists to enforce is not needed here --
            // same reasoning as t_mag under STS above. The fully implicit
            // solve (use_full_implicit_coag_) needs neither this nor
            // dt_gain_global below: Gain and Loss are linearized and solved
            // together, self-consistently, in one step -- see
            // ApplyCoagulationAndCorrectMass_FullImplicit /
            // coagulation_speedup_literature.md Sec.3.
            dt_coag = (use_implicit_coag_ || use_full_implicit_coag_) ? 1.0e100
                : EstimateDtCoagulation(sigmad_active, dNdt_coag, sigmad_floor_, nr_s, nr_e, nbin);

            // Growth-direction global bound (dust_coagulation_implicit_note.md
            // Sec.12): internal-only n_sub sub-cycling in
            // ApplyCoagulationAndCorrectMass_Implicit proved insufficient for
            // the most fragmentation-active cells, so Gain (unlike Loss) does
            // get a global dt constraint here, restoring what dt_coag used to
            // provide for growth-dominant bins "by accident" (via some other
            // cell's depletion keeping dt small globally) before use_implicit_coag
            // removed dt_coag entirely. Not needed for use_full_implicit_coag_
            // (see comment above).
            dt_gain_global = (use_implicit_coag_ && !use_full_implicit_coag_)
                ? EstimateDtGain(sigmad_tot_, Gain_coag, nr_s, nr_e, nbin)
                : 1.0e100;

            // Operator split: dt_adv (dust radial advection) is NOT allowed to
            // constrain the outer step -- it is handled by sub-cycling the dust
            // advection update below. All the other timescales (gas advection,
            // wind, coagulation, energy source, magnetic) still bound dt here.
            dt = std::min(dt_wind, dt_coag);
            dt = std :: min(dt_adv_gas , dt);
            dt = std :: min(dt_Esrc,dt);
            dt = std :: min(dt_gain_global, dt);
            dt = cfl_ * std :: min(dt,t_mag);

            if (dt_tot + dt > DeltaT)
            {
                dt = DeltaT - dt_tot;
            }

            // Sub-cycle count for the dust radial-advection update: advance it
            // in n_adv_sub explicit steps of dt_adv_sub (<= cfl*dt_adv) while
            // the outer step dt is used for B, gas, thermal, coagulation,
            // resistivity, velocities, and the (gas-timescale) vapor advection.
            double dt_adv_cfl = cfl_ * dt_adv;
            int    n_adv_sub  = 1;
            if (dt_adv_cfl > 0.0 && std::isfinite(dt_adv_cfl) && dt > dt_adv_cfl)
                n_adv_sub = (int)std::ceil(dt / dt_adv_cfl);
            double dt_adv_sub = dt / (double)n_adv_sub;

            // Update Magnetic field
            Calc_Mag(&mag,dt);
            CalculationMRIturbulence(&mag);
            //CalculationMRIturbulence2(&mag);

            B_ = mag.B_;
            beta_disk_ = mag.beta_RZ;
            alpha_MRI = mag.Alpha_;

            // Calcualte advection term
            // Gas surface density
            if (use_face_scheme_)
            {
                AdvanceGasFaceScheme(dt);
            }
            else
            {
            for (int i = nr_s; i < nr_e; ++i)
            { // loop for disk radius
                vol_inv = inv_vol_[i];
                sigmag_[i] = sigmag_[i] - 2.0 *dt*(r_bod_[i + 1]*fluxg_[i + 1] - r_bod_[i]*fluxg_[i]) * vol_inv - dsigmagdt_w[i]*dt;
                if(sigmag_[i] < 1.0e-4)
                {
                    sigmag_[i] = 1.0e-4;
                }
            }
            }

            //determine Temepature
            // if(on_thermal_ == true){
            //     eng.GetTemperature_TimeEvolution(T_,dt);
            // }
            
            // Dust surface density.
            // Radial advection of the size bins is operator-split and
            // sub-cycled: n_adv_sub explicit steps of dt_adv_sub, each
            // recomputing the upwind flux from the (frozen) dust velocities and
            // the current sigma. The wind-removal source (-dsigmaddt_w) stays
            // inside the sub-cycle (frozen rate -> the sub-divided total is
            // identical, and small increments stay non-negative). Vapor rides
            // the gas timescale and is advanced once at the outer dt below.
            if (!no_dust_)
            {
            for (int s = 0; s < n_adv_sub; ++s)
            {
                if (use_face_scheme_)
                {
                    AdvanceDustAdvectionFace(dt_adv_sub);   // MUSCL + SSP-RK2
                    continue;
                }
                // CalcDustFlux advects sigmad_active directly (see the bin-0
                // note there), so nothing needs to be synced into sigmad_ here;
                // sigmad_ is rebuilt once at the end of the outer step.
                CalcDustFlux();
                for (int i = nr_s; i < nr_e; ++i)
                {
                    vol_inv = inv_vol_[i];
                    for (int j = 0; j < nbin; ++j)
                    {
                        sigmad_active[i][j] += ( -2.0*(r_bod_[i + 1]*fluxd_[i + 1][j] - r_bod_[i]*fluxd_[i][j])*vol_inv
                                                 - dsigmaddt_w[i][j] ) * dt_adv_sub;
                    }
                }
            }
            // Vapor advection (gas timescale -> outer dt, once) + dust totals.
            // fluxd_vapor from the last CalcDustFlux above is the start-of-step
            // vapor flux: it depends only on sigmad_vapor / fdg_ratio_vapor,
            // both frozen through the sub-cycle.
            if (use_face_scheme_) AdvanceVaporFace(dt);   // MUSCL + SSP-RK2 at the face gas velocity
            for (int i = nr_s; i < nr_e; ++i)
            {
                vol_inv = inv_vol_[i];
                if (!use_face_scheme_)
                sigmad_vapor[i] = sigmad_vapor[i] - 2.0 *dt*(r_bod_[i + 1]*fluxd_vapor[i + 1] - r_bod_[i]*fluxd_vapor[i])*vol_inv - v_dsigmaddt_w[i]*dt;
                sigmad_tot_[i] = 0.0;
                for (int j = 0; j < nbin; ++j)
                {
                    sigmad_tot_[i] += sigmad_active[i][j];
                }
            }

            // Turbulent radial dust diffusion: operator-split, backward-Euler
            // implicit (unconditionally stable -> full outer dt, no sub-cycle).
            ApplyDustDiffusionImplicit(dt);
            for (int i = nr_s; i < nr_e; ++i)
            {
                sigmad_tot_[i] = 0.0;
                for (int j = 0; j < nbin; ++j)
                {
                    sigmad_tot_[i] += sigmad_active[i][j];
                }
            }

            // Coagulation/fragmentation solve -- per radial cell, no radial
            // coupling in this step, so parallelize over i. The fully-implicit
            // path needs per-thread scratch (Jacobian + linear-solve buffers);
            // the semi-implicit / explicit paths read per-grid Gain/Loss/dNdt
            // rows already indexed by i, so they need none.
            #pragma omp parallel
            {
                DoubleVector Gain_t, LossRate_t, rhs_t, sigma_o_t;
                DoubleMatrix J_t, A_t;
                if (use_full_implicit_coag_)
                {
                    Gain_t.assign(nbin, 0.0);
                    LossRate_t.assign(nbin, 0.0);
                    rhs_t.assign(nbin, 0.0);
                    sigma_o_t.assign(nbin, 0.0);
                    J_t.assign(nbin, DoubleVector(nbin, 0.0));
                    A_t.assign(nbin, DoubleVector(nbin, 0.0));
                }

                #pragma omp for schedule(dynamic, 4)
                for (int i = nr_s; i < nr_e; ++i)
                {
                    if (use_full_implicit_coag_)
                    {
                        ApplyCoagulationAndCorrectMass_FullImplicit(sigmad_active, sigmad_tot_,
                            Gain_t, LossRate_t, J_t, A_t, rhs_t, sigma_o_t, dt, i, nbin);
                    }
                    else if (use_implicit_coag_)
                    {
                        ApplyCoagulationAndCorrectMass_Implicit(sigmad_active, sigmad_tot_, Gain_coag, LossRate_coag, dt, i, nbin);
                    }
                    else
                    {
                        ApplyCoagulationAndCorrectMass(sigmad_active, sigmad_tot_, dNdt_coag, dt, i, nbin);
                    }
                }
            }

            double sig_tot , a1;
            for (int i = nr_s; i < nr_e; ++i)
            {
                sigmad_tot_[i] = sigmad_tot_[i]+sigmad_vapor[i];
                sigmad_[i][0] = sigmad_vapor[i]+ sigmad_active[i][0]; //凝結して０.１ミクロンのダストに
                for(int j=1;j<nbin;j++)
                {
                    sigmad_[i][j] = sigmad_active[i][j];
                }
            }
            } // !no_dust_

            // if(icount < 10 && !restart_){
            //     // time_count = int(t/SOLAR_YEAR);
            //     filename = outdir_ + "/try" + std::to_string(icount);
            //     filename_mag = outdir_mag + "/try" + std::to_string(icount);
            //     OutPut(filename,filename_mag);
            //     std::cout <<"try "<< icount << ", time (yr) = " << (dt_tot / SOLAR_YEAR) << ", dt = " << (dt / SOLAR_YEAR) <<" ,M_disk = "<<(M_total - star_mass_)/star_mass_<<" M."<< std::endl;
            // }
            
            dt_tot += dt;
            icount++;
        }

        // bool was_early = (t < 1.0*SOLAR_YEAR);
        t += DeltaT;
        // if(was_early) {
        //     time_count = (int)round(t * 10.0 / SOLAR_YEAR); // 0.1yr→1, 0.2yr→2, ..., 1.0yr→10
        // } else {
        //     time_count = int(t/SOLAR_YEAR); // 11yr→11, 20yr→20, ...
        // }
        time_count = int(t/SOLAR_YEAR);
        filename = outdir_ + "/out" + std::to_string(time_count);
        filename_mag = outdir_mag + "/out" + std::to_string(time_count);
        OutPut(filename,filename_mag);
        std::cout << time_count << ", time (yr) = " << (t / SOLAR_YEAR) << ", dt = " << (dt / SOLAR_YEAR) <<" ,M_disk = "<<(M_total - star_mass_)/star_mass_<<" M."<< std::endl;
    }

    return;
}

void DustDisk::OutPut(std::string filename,std::string filename_mag)
{
    std :: string filename_gas;
    std :: string filename_magZ;
    filename_magZ = filename_mag + "_B";
    filename_gas = filename + "_gas";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    std::ofstream file_gas(filename_gas, std::ios::trunc | std::ios::out);
    std :: ofstream file_Bz(filename_magZ, std::ios::trunc | std::ios::out);
  
    FAILED_TO_OPEN(file, filename);
    FAILED_TO_OPEN(file_gas, filename);
    FAILED_TO_OPEN(file_Bz,filename_mag);

    file << std::scientific;
    for (int i = 0; i < nbin_; ++i)
    {
        file << m_[i] << " ";
    }
    file << std::endl;

    for (int i = 0; i < nbin_ + 1; ++i)
    {
        file << mb_[i] << " ";
    }
    file << std::endl;

    for (int i = 0; i < nbin_; ++i)
    {
        file << a_[i] << " ";
    }
    file << std::endl;

    for (int i = 0; i < nbin_ + 1; ++i)
    {
        file << ab_[i] << " ";
    }
    file << std::endl;

    file << std::endl;
    for (int i = nr_s; i < nr_e; ++i)
    {
        file << r_cen_[i] << " ";
        for (int j = 0; j < nbin_; ++j)
        {
            file << sigmad_active[i][j] << " ";
        }
        file << std::endl;
    }

    for (int i = nr_s; i < nr_e; ++i)
    {
        file_gas << std :: scientific << r_cen_[i] <<" " <<sigmag_[i]<<" "<<sigmad_tot_[i]<<" "<<sigmad_vapor[i]<<" "<<vg_r_[i]<<" "<<vg_phi_[i]<<" "<<vg_wind_[i]<<" "<<vd_r_tot[i][0]<<" "<<vd_r_[i][0]<<" "<<vd_r_diff[i][0]<<" "<<fdg_ratio[i][0]<<" "<<dsigmagdt_w[i]<<" "<<T_[i]<<" "<<P_[i][nz_s]<<" "<<P_[i][nz_s]/beta_disk_[i][nz_s]<< " "<<dNdr_torque[i]<<" "<<alpha_GI[i]<<" "<<alpha_phiz[i]<<" "<<alpha_tot[i]<<" "<<eta_P_[i]<<" "<<windDW_zdead_[i]<<" "<<windDW_zsurf_[i]<<" "<<windDW_Delta_[i]<<" "<<windDW_Ameff_[i]<<" "<<alpha_DW_[i]<<" ";
        // [face_flux_scheme] 25:r_face(inner face of cell i) 26:vg_r_face 27:dvg_phi_face 28:v_Bdiff_face
        if (use_face_scheme_ && (int)vgf_.size() > i)
            file_gas << r_bod_[i] << " " << vgf_[i] << " " << dvphif_[i] << " " << vBdifff_[i] << " ";
        file_gas << std :: endl; 
        file_Bz << std :: scientific << r_cen_[i] << " "<< alpha_MRI[i] <<" "<<B_[i]<<" "<<etaO_R[i]<<" "<<etaA_R[i]<<" "<<beta_disk_[i][nz_s]<<std::endl;
    }

    file.close();
    file_gas.close();
    file_Bz.close();
    return;
}

void DustDisk::UpdateNoGrowthDustDisk()
{
    /*
        main calculate function in the case of no dust growth (on_growth = false):
        gas / magnetic field evolve, dust is only advected, removed by the wind
        and diffused radially (no coagulation / fragmentation).
        Time stepping, output and restart follow UpdateDustDiskFragOn(): the
        outer loop advances one output interval DeltaT (5e2 / 5e3 / 5e4 yr for
        time_count < 1e4 / < 1e5 / otherwise), the inner loop sub-steps it with
        the CFL-limited dt, and out<int(t/yr)> is written after each DeltaT.
    */

    std::cout << "NoGrowth-UpdateDisk" << std::endl;

    int icount, time_count, nbin;
    double t, time_end, DeltaT, dt_tot, dt, dt_adv, dt_adv_gas, dt_wind, vol_inv;
    DoubleMatrix dNdt_rd(nr_tot, DoubleVector(nbin_, 0.0));
    std::string filename;
    std::string filename_mag;
    Mag_field mag;

    // Grid geometry is fixed for the whole run -> push it into `mag` once here
    // (same as UpdateDustDiskFragOn).
    SetRadiusGridParameter(nr_s,nr_e,r_cen_,r_bod_,dr_cen_,dr_bod_,inv_vol_);
    SetVerticalGridParameter(nz_s,nz_e,z_cen_,z_bod_,dz_cen_,dz_bod_);
    Set_Differential_Coefficient(cff_R,cff_R_Half);

    nbin = nbin_;
    time_end = tend_;
    time_count = 1;
    DeltaT = SOLAR_YEAR;
    icount = 0;
    t = 0.0;

    if(restart_ == false){
        Iteration();
        filename = outdir_ + "/out0";
        filename_mag = outdir_mag + "/out0";
        OutPut(filename,filename_mag);
        std::cout << icount << ", time (yr) = " << (t / SOLAR_YEAR)  <<" ,M_disk = "<<(M_total - star_mass_)/star_mass_<<" M."<< std::endl;
    }
    else{
        t = t_re;
        time_count = int(t/SOLAR_YEAR);
        UpdatePhysicalParameter();
        CalcStokesNumberAndDustScaleHeight();
        UpdateResistivity();
        CalcWindAlpha();
        UpdateTurbulent();
        CalculateViscousTorque();
        CalculateDust_to_Gas_massratio();
        CalcWindMassLoss_Tabone2022();

        CalculateDiskCoeffcient();
        CalcGasVelocity_TM2018();
        CalcDustVelocity_TM2018();

        CalcGasFlux();
        CalcDustFlux();
        filename = outdir_ + "/restart" + std::to_string(time_count);
        filename_mag = outdir_mag + "/restart" + std::to_string(time_count);
        OutPut(filename,filename_mag);
        std::cout << time_count << ", time (yr) = " << (t / SOLAR_YEAR)  << ": Restart !! "<< std::endl;
    }

    // calculation
    while (t < time_end)
    {
        dt_tot = 0.0;
        if(time_count < 1.0e4){DeltaT = 5.0e2*SOLAR_YEAR;}
        else if(time_count < 1.0e5){DeltaT = 5.0e3*SOLAR_YEAR;}
        else{DeltaT = 5.0e4*SOLAR_YEAR;}

        if(t+DeltaT > time_end){DeltaT = time_end-t;}

        while (dt_tot < DeltaT)
        {
            // calculate dSigma/dt
            UpdatePhysicalParameter();
            CalcStokesNumberAndDustScaleHeight();
            UpdateResistivity();
            CalcWindAlpha();
            UpdateTurbulent();
            CalculateViscousTorque();
            CalculateDust_to_Gas_massratio();
            CalcWindMassLoss_Tabone2022();

            CalculateDiskCoeffcient();
            CalcGasVelocity_TM2018();
            CalcDustVelocity_TM2018();

            CalcDustFlux();
            CalcGasFlux();

            // determine dt
            t_mag = 1.0e100;
            dt_adv_gas = 1.0e100;
            dt_wind = 1.0e100;
            dt_adv = 1.0e100;
            // determine dt for magnetic field
            // (grid geometry already pushed into `mag` once, before the time loop)
            SetDiskParameter(ng_,vgr_bod_,Roe_mean_,P_,Omega_,Hg_); // vgrtestはガスの速度が0として計算
            SetMagneticFieldParameter(etaO,etaA,eta_total,B_,B_out);
            CalculateDeltaTime(&mag);
            // implicit magnetic diffusion is stable for any dt -> no t_mag cap
            t_mag = IsImplicitMagneticDiffusionEnabled() ? 1.0e100 : mag.Delta_time;

            if (!no_dust_)
            {
                for (int i = nr_s; i < nr_e; ++i)
                { // loop for disk radius
                    vol_inv = inv_vol_[i];
                    for (int j = 0; j < nbin; ++j)
                    { // loop for dust bin
                        // advection term
                        dNdt_rd[i][j] = -2.0 * (r_bod_[i + 1]*fluxd_[i + 1][j] - r_bod_[i] * fluxd_[i][j])*vol_inv - dsigmaddt_w[i][j];
                        if (sigmad_active[i][j] > sigmad_floor_[i][j] && std::abs(vd_r_tot[i][j]) != 0.0)
                        {
                            dt_adv = std::min(dt_adv, std::abs(dr_bod_[i] / vd_r_tot[i][j]));
                        }
                    }
                }
            }

            for (int i = nr_s; i < nr_e; ++i)
            { // loop for disk radius
                dt_adv_gas = std::min( dt_adv_gas, std::abs(dr_bod_[i]/vg_r_[i]) );
                dt_wind = std::min(dt_wind, std::abs( sigmag_[i]/dsigmagdt_w[i] ));
            }

            // dust advection is explicit (no sub-cycling here) -> dt_adv bounds dt
            dt = std::min(dt_adv, dt_wind);
            dt = std :: min(dt_adv_gas , dt);
            dt = cfl_ * std :: min(dt,t_mag);

            if (dt_tot + dt > DeltaT)
            {
                dt = DeltaT - dt_tot;
            }

            // Update Magnetic field
            Calc_Mag(&mag,dt);
            CalculationMRIturbulence(&mag);

            B_ = mag.B_;
            beta_disk_ = mag.beta_RZ;
            alpha_MRI = mag.Alpha_;

            // Calcualte advection term
            // Gas surface density
            for (int i = nr_s; i < nr_e; ++i)
            { // loop for disk radius
                vol_inv = inv_vol_[i];
                sigmag_[i] = sigmag_[i] - 2.0 *dt*(r_bod_[i + 1]*fluxg_[i + 1] - r_bod_[i]*fluxg_[i]) * vol_inv - dsigmagdt_w[i]*dt;
                if(sigmag_[i] < 1.0e-10)
                {
                    sigmag_[i] = 1.0e-10;
                }
            }

            // Dust surface density
            if (!no_dust_)
            {
                for (int i = nr_s; i < nr_e; ++i)
                {
                    sigmad_tot_[i] = 0.0;
                    for (int j = 0; j < nbin; ++j)
                    {
                        sigmad_active[i][j] += dNdt_rd[i][j] * dt;
                        sigmad_tot_[i] += sigmad_active[i][j];
                    }
                }

                // Turbulent radial dust diffusion (operator-split, implicit).
                ApplyDustDiffusionImplicit(dt);
                for (int i = nr_s; i < nr_e; ++i)
                {
                    sigmad_tot_[i] = 0.0;
                    for (int j = 0; j < nbin; ++j) sigmad_tot_[i] += sigmad_active[i][j];
                }
            }

            dt_tot += dt;
            icount++;
        }

        t += DeltaT;
        time_count = int(t/SOLAR_YEAR);
        filename = outdir_ + "/out" + std::to_string(time_count);
        filename_mag = outdir_mag + "/out" + std::to_string(time_count);
        OutPut(filename,filename_mag);
        std::cout << time_count << ", time (yr) = " << (t / SOLAR_YEAR) << ", dt = " << (dt / SOLAR_YEAR) <<" ,M_disk = "<<(M_total - star_mass_)/star_mass_<<" M."<< std::endl;
    }

    return;
}

void DustDisk :: Iteration()
{
    std :: cout << "Iteration" << std :: endl;
    Energy_Calc eng;
    Mag_field mag;
    Resistivity res;
    DoubleVector gosa_R;
    int n = 0;
    double gosa = 1.0e6;

    sigmad_active = sigmad_;

    DetermineVaporDust();
    UpdatePhysicalParameter();
    CalcStokesNumberAndDustScaleHeight();
    UpdateResistivity();
    CalcWindAlpha();
    UpdateTurbulent();
    CalculateViscousTorque();
    CalculateDust_to_Gas_massratio();
    CalcWindMassLoss_Tabone2022();

    CalculateDiskCoeffcient();
    CalcGasVelocity_TM2018();
    CalcDustVelocity_TM2018();
    
    mag.B_ = B_;
    gosa_R.resize(nr_tot,0.0);
    // Rebuild the initial dust distribution every pass when a_max(r) depends
    // on T / cs / alpha_MRI (amax_mode = stokes or frag_drift), so dust,
    // resistivity, alpha_MRI and T converge together.
    bool redo_dust_init = (dust_init_mode_ == "powerlaw") && (dust_init_amax_mode_ != "fixed");
    DoubleVector a_max_old;
    // determine Temepature
    
    while (gosa > 5.0e-2)
    {
        // // determine Temepature
        if(on_thermal_ == true){
            DoubleVector T_old = T_;
            eng.SetGridParameter(nr_s,nr_e,nr_tot);
            eng.SetRadDiffFactor(rad_diff_fac_);
            eng.SetWindLeverArm(lambda_w_);
            eng.SetDiskParameter(r_cen_,sigmag_,T_,T_irr,T_irr4,Omega_,dsigmagdt_w,alpha_MRI,alpha_phiz);
            eng.SetMagneticField(B_,dB_dr_,eta_total);
            eng.GetTemperature(T_);
            //eng.GetTemperature_Initialization(T_);
            // Relax the update: T_ used to be overwritten outright every
            // pass. The K/Na charge-exchange fix (2026-09-03) made
            // eta_total's T-dependence steeper at the thermal-ionization
            // front, turning this previously-convergent, undamped update
            // into a T<->eta_total<->alpha_MRI oscillation there. A 50/50
            // geometric mean of T_old and T_new (matching alpha_MRI's own
            // damping below) settles into a stable period-4 limit cycle
            // instead of decaying -- average in 1/T space instead (harmonic
            // mean), which weights toward the smaller of the two.
            for (int i = nr_s; i < nr_e; ++i)
            {
                T_[i] = 2.0*T_[i]*T_old[i]/(T_[i]+T_old[i]);
            }
        }
        
        UpdatePhysicalParameter();

        if (redo_dust_init)
        {
            a_max_old = a_max_init_;
            SetInitialDustDistribution();
        }

        DetermineVaporDust();
        
        CalcStokesNumberAndDustScaleHeight();
        UpdateResistivity();
        DoubleVector alpha_phiz_old = alpha_phiz;
        CalcWindAlpha();
        double gosa_wind = 0.0;
        for (int i = nr_s; i < nr_e; ++i)
        {
            if (alpha_phiz_old[i] > 0.0 && i < nr_e - 1)
                gosa_wind = std :: max(gosa_wind, fabs(alpha_phiz[i] - alpha_phiz_old[i])/alpha_phiz_old[i]);
            if (alpha_phiz_old[i] > 0.0 && alpha_phiz[i] > 0.0)
                alpha_phiz[i] = std :: sqrt(alpha_phiz[i]*alpha_phiz_old[i]);
        }
        alpha_phiz[nr_s - 1] = alpha_phiz[nr_s];
        alpha_phiz[nr_e]     = alpha_phiz[nr_e - 1];
        UpdateTurbulent();
        CalculateViscousTorque();

        // determine dt for magnetic field
        SetRadiusGridParameter(nr_s,nr_e,r_cen_,r_bod_,dr_cen_,dr_bod_,inv_vol_);
        SetVerticalGridParameter(nz_s,nz_e,z_cen_,z_bod_,dz_cen_,dz_bod_);
        Set_Differential_Coefficient(cff_R,cff_R_Half);
        SetDiskParameter(ng_,vgr_bod_,Roe_mean_,P_,Omega_,Hg_); // vgrtestはガスの速度が0として計算
        SetMagneticFieldParameter(etaO,etaA,eta_total,B_,B_out); 
        CalculateDeltaTime(&mag);
        mag.beta_RZ = beta_disk_ini;
        
        CalculationMRIturbulence(&mag);
        //CalculationMRIturbulence2(&mag);
        for (int i=nr_s;i<nr_e;i++)
        {
            gosa_R[i] = fabs(mag.Alpha_[i] - alpha_MRI[i])/alpha_MRI[i];
        }
        gosa = *std :: max_element(gosa_R.begin()+nr_s,gosa_R.begin()+nr_e-1);
        //std :: cout << "Iteration " << n << " : gosa(alpha) = " << gosa << " , gosa(alpha_phiz) = " << gosa_wind << std :: endl;
        gosa = std :: max(gosa, gosa_wind);
        if (redo_dust_init)
        {
            double gosa_amax = 0.0;
            for (int i=nr_s;i<nr_e;i++)
            {
                if (a_max_old[i] > 0.0)
                    gosa_amax = std :: max(gosa_amax, fabs(a_max_init_[i] - a_max_old[i])/a_max_old[i]);
            }
            std :: cout << "Iteration " << n << " : gosa(alpha) = " << gosa << " , gosa(a_max) = " << gosa_amax << std :: endl;
            gosa = std :: max(gosa, gosa_amax);
        }
        for (int i=nr_s;i<nr_e;i++)
        {
            alpha_MRI[i] = std :: sqrt(mag.Alpha_[i]*alpha_MRI[i]);
        }

        n += 1;
        if (n>30)
        {
            std :: cerr << "Iteration error !! : Initialization : "<<gosa<< std :: endl;
            exit(1);
        }
    }
}