#ifndef DUST_DISK_HPP_
#define DUST_DISK_HPP_

/*
DustDisk

Created by kawasaki 

version 1.0
*/


// include header file
#include "physical_constants.hpp"
#include "alpha.hpp"
#include "./magnetic_field1D.hpp"
#include "typedef.hpp"
#include "utils.hpp"

// Shared constants, used by both dust_disk.cpp and dust_growth.cpp.
// inline constexpr (C++17) gives these external linkage with no ODR risk
// across translation units, unlike a plain constexpr at namespace scope.
inline constexpr double sqrt2 = M_SQRT2;                // sqrt(2)
inline constexpr double sqrt1_2 = M_SQRT1_2;            // 1 / sqrt(2)
inline constexpr double sqrt1_2pi = 0.3989422804014327; // 1.0 / sqrt(2*pi)
inline constexpr double sqrtpi_8 = 0.6266570686577501;  // sqrt(pi/8)
inline constexpr double sqrt8_pi = 1.5957691216057308;  // sqrt(8/pi)
inline constexpr double pi_8 = M_PI / 8.0;              // pi / 8
inline constexpr double frac1_pi = 1.0/M_PI;            // 1/pi
inline constexpr double frac1_4pi = 0.25/M_PI;          // 1/4pi
inline constexpr double frac1_8pi = 0.125/M_PI;         // 1/8pi
inline constexpr double num_limit = __DBL_EPSILON__;

inline constexpr double coe_vbr = 8.0 * BOLTZMANN_CONSTANT / M_PI;

class DustDisk : public CalcAlphaParameter,Magnetic_Field2D
{
public:

    DustDisk();
    ~DustDisk();

    void InputParameter(Utils::InputConfigure& input);

    void Initialization();
    void Initialization_Restart(std :: string dustfile,std :: string gasfile,std :: string Bfile,double t_restart);

    void CalculationRun();

private:

    /* menber functions (private) */
    void SetRadiusGrid();
    void SetVerticalGrid();
    void SetDiffrentialCoefficient_R();
    void SetDiffrentialCoefficient_R_Half();
    void SetDiffrentialCoefficientFor2order_R();
    void SetDiffrentialCoefficient_Z();
    void SetZeroDifferentialCoeffient_R();
    void SetBoundaryConditionForDisk();
    void SetBoundaryConditionForMagneticfield();
    void SetDustGrid();
    void SetDiskGasModel();
    void SetInitialDustDistribution();
    double CalcSizeForTargetStokes(int i, double St_target);
    double CalcInitialStokesLimit(int i);

    double CalcXrayIonizationrate(double sigma,double rho_mid,double H_g,double R,double Z);
    double CalcCosmicRayIonizationrate(double sigma,double rho_mid,double H_g,double R,double Z);

    void CalcCoagulationAndFragmentationArrays();
    void CalcStokesNumberAndDustScaleHeight();
    void CalculateViscousTorque();
    void CalculateDust_to_Gas_massratio();
    void CalculateDiskCoeffcient();
    void CalcGasVelocity_TM2018();
    void CalcDustVelocity_TM2018();
    void CalcDustVelocity();
    void CalcDustRelativeVelocity(int ir, int id, int jd, double &vd_ij, double &sigma_ij, double &vrel);
    void CalcGasFlux();
    void CalcDustFlux();
    void ApplyDustDiffusionImplicit(double dt);

    // Face-flux scheme (face_scheme.cpp, input flag face_flux_scheme)
    void SetupFaceScheme();
    void CalcGasFaceVelocity();
    void AdvanceGasFaceScheme(double dt);
    void CalcDustFaceVelocity();
    void CalcDustFluxFace();
    void CalcVaporFluxFace(const DoubleVector &Sv);
    double DustAdvectionDtFace() const;
    double GasAdvectionDtFace() const;
    void AdvanceDustAdvectionFace(double h);
    void AdvanceVaporFace(double dt);
    void CalcMagneticFieldPressure();
    void CalcWindMassLoss_Riols2019(); 
    void CalcWindMassLoss_Tabone2022(); 
    void CalcWindTorque_Suzuki2016();// Wind torqueあり
    void CalcWindAlpha();            // wind_alpha_model_ で Suzuki2016 / Enomoto2026 を分岐
    void CalcWindDrivenAlpha_Enomoto(); // Enomoto et al. (2026) 自己無撞着 alpha_DW
    void UpdatePhysicalParameter();

    void DetermineVaporDust();

    void CoagulateChangeRate(double t, DoubleMatrix &y, DoubleMatrix &dydt);
    void CoagulateChangeRateFragmentation(double t, DoubleMatrix &y, DoubleMatrix &dydt);

    double EstimateDtCoagulation(const std::vector<std::vector<double>>& sigmad,
        const std::vector<std::vector<double>>& dNdt_coag,const std::vector<std::vector<double>>& sigmad_floor,
        int nr_s, int nr_e, int nbin);

    // Global (not per-cell-internal) Gain-timescale bound for use_implicit_coag,
    // added 2026-08-30 (dust_coagulation_implicit_note.md Sec.12): internal
    // sub-cycling alone proved insufficient/fragile for the most fragmentation-
    // active cells, so this restores a *growth*-direction analog of
    // EstimateDtCoagulation's depletion-direction dt_coag -- Loss stays
    // unbounded (protected by the implicit division), only Gain feeds this.
    // Uses the cell-total sigmad_tot[i] as the reference (Sec.6's original
    // criterion), NOT the per-bin sigmad_floor_ reference used by the
    // *internal, N_SUB_MAX-capped* n_sub criterion in Sec.9 -- that floor is
    // fine when a hard cap bounds the fallout, but with no such cap here a
    // freshly-empty bin anywhere in the whole grid collapsed dt to ~1e-45yr
    // (Sec.12 postmortem).
    double EstimateDtGain(const std::vector<double>& sigmad_tot,
        const std::vector<std::vector<double>>& Gain,
        int nr_s, int nr_e, int nbin);

    void ApplyCoagulationAndCorrectMass(std::vector<std::vector<double>>& sigmad,
        std::vector<double>& sigmad_tot,const std::vector<std::vector<double>>& dNdt_coag,
        double dt, int i, int nbin);

    // Semi-implicit coagulation/fragmentation (see dust_coagulation_implicit_note.md):
    // Loss is linearized and solved implicitly per bin (no matrix), Gain stays
    // explicit. Selected via the "use_implicit_coag" input flag.
    void CoagulateChangeRateFragmentation_Implicit(double t, DoubleMatrix &y, DoubleMatrix &Gain, DoubleMatrix &LossRate);
    // Per-cell body of the above, re-invoked mid-sub-cycle by
    // ApplyCoagulationAndCorrectMass_Implicit() to refresh Gain[ir]/LossRate[ir]
    // from the current (not the stale sigmad^n) distribution.
    void CoagulateChangeRateFragmentation_Implicit_Cell(int ir, const DoubleMatrix &y, DoubleMatrix &Gain, DoubleMatrix &LossRate);
    void ApplyCoagulationAndCorrectMass_Implicit(std::vector<std::vector<double>>& sigmad,
        std::vector<double>& sigmad_tot,
        std::vector<std::vector<double>>& Gain,
        std::vector<std::vector<double>>& LossRate,
        double dt, int i, int nbin);

    // Fully implicit coagulation/fragmentation (Jacobian + dense linear solve
    // per cell; see coagulation_speedup_literature.md Sec.3 / Brauer et al.
    // 2008 Appendix C). Selected via the "use_full_implicit_coag" input flag
    // (takes priority over "use_implicit_coag" if both are set).
    void CoagulateChangeRateFragmentation_FullImplicit_Cell(int ir, const DoubleMatrix &y,
        DoubleVector &Gain, DoubleVector &LossRate, DoubleMatrix &J);
    void ApplyCoagulationAndCorrectMass_FullImplicit(std::vector<std::vector<double>>& sigmad,
        std::vector<double>& sigmad_tot,
        DoubleVector &Gain_scratch, DoubleVector &LossRate_scratch, DoubleMatrix &J_scratch,
        DoubleMatrix &A_scratch, DoubleVector &rhs_scratch, DoubleVector &sigma_o_scratch,
        double dt, int ir, int nbin);

    void UpdateResistivity();
    void UpdateTurbulent();
    void UpdateDustDisk();

    void Iteration();
    void UpdateNoGrowthDustDisk();
    void UpdateDustDiskFragOn();

    void OutPut(std::string filename,std::string filename_mag);

    /* menber variables */
    
    // disk radius
    int NG;        // number of ghost cell grid
    int nr_tot;
    int nr_s;
    int nr_e;
    int nr_disk_;  // number of disk radius grid

    double rin_;           // inner disk radius [au]
    double rout_;          // outer disk radius [au]

    DoubleVector r_cen_;   // cell center [au]
    DoubleVector r_bod_;   // cell boudary [au]
    DoubleVector r_vol_;   // cell volume [au]
    DoubleVector dr_bod_;  // specing between cell boundaries [au]
    DoubleVector dr_cen_;  // specing between cell centers [au]
    DoubleVector inv_vol_; //
    // Geometry-only interpolation weights (precomputed once in SetRadiusGrid),
    // used to de-checkerboard vg_r_ via a centre->face->centre round trip.
    DoubleVector vg_w_cf_; // centres (i-1,i) -> face i   : (r_bod[i]-r_cen[i-1])/(r_cen[i]-r_cen[i-1])
    DoubleVector vg_w_fc_; // faces (i,i+1)   -> centre i : (r_cen[i]-r_bod[i])/(r_bod[i+1]-r_bod[i])
    DoubleMatrix cff_R;
    DoubleMatrix cff_R_Half;
    DoubleMatrix cff2_R;
    DoubleMatrix zero_cff2_R;

    // disk height
    double z_max;
    int nz_disk_; // disk gas height grid
    int nz_tot; 
    int nz_s;
    int nz_e;

    DoubleVector z_cen_;   // cell center [au]
    DoubleVector z_bod_;   // cell boudary [au]
    DoubleVector dz_bod_;  // 
    DoubleVector dz_cen_;  // specing between cell centers [au]
    DoubleMatrix cff_Z;

    double sigmag0_;

    // disk magnetic field
    double B_out;           // magnetic field in outer boundary
    double beta_ini;        // initial plasma beta
    DoubleVector B_ ;       // disk magnetic field at radius
    DoubleVector dB_dr_;
    DoubleMatrix beta_disk_;// disk plasma beta
    DoubleMatrix beta_disk_ini;// initial disk plasma beta
    // ionization 
    DoubleMatrix zeta_; // ionization in disk[s^-1]

    // dust size
    double amin_;           // minimum dust size [cm]            
    double amax_;           // maximum dust size [cm]
    double mmax_;           // mimimum dust mass [g]
    double mmin_;           // maximum dust mass [g]
    double rhodi_;          // dust internal(material) density [g/cm^3]
    
    UnsignedInt nbin_;      // number of dust bin
    UnsignedInt nbin_lin_;  // number of liner bin
    DoubleVector a_;        // dust size bin [cm]
    DoubleVector ab_;       // dust size bin at boundary [cm]
    DoubleVector m_;        // dust mass bin [g]
    DoubleVector mb_;       // dust mass bin at boundary [g]

    // dust disk
    DoubleMatrix sigmad_;        // dust surface density [g/cm^2]
    DoubleMatrix sigmad_floor_;  // floor value of dust surface density [g/cm^2]
    DoubleVector sigmad_tot_;    // total dust surface density at each radius

    DoubleMatrix sigmad_active;  // dust surface density [g/cm^2]
    DoubleVector sigmad_vapor;   // dust surface density [g/cm^2]

    DoubleMatrix Nd_;            // dust surface number density [/cm^2]
    DoubleMatrix Nd_RZ;          // dust number density at each area
    DoubleMatrix St_;            // Stokes number 
    DoubleMatrix St_mean;        // Mean Stokes number 
    DoubleMatrix vd_r_;          // dust radial velocity [cm/s]
    DoubleMatrix vd_r_tot;       // dust radial total velocity [cm/s]
    DoubleMatrix vd_r_diff;      // diagnostic only: effective turbulent-diffusion radial velocity, back-computed as F_diff/Sigma_d [cm/s] (NOT used for transport; see fluxd_diff_)
    DoubleMatrix vd_phi_;        // dust azimuthal velocity [cm/s]
    DoubleMatrix Hd_;            // dust scale height [cm]
    DoubleMatrix Dd_;            // turbulent diffusion coefficient for dust
    DoubleMatrix fdg_ratio;      // dust to gas mass ratio
    DoubleMatrix dfdg_dr;        // difd
    double fdg_;                 // dust to gas mass ratio
    DoubleMatrix fluxd_;
    DoubleVector fluxd_vapor;

    // gas
    double r_cut;
    double M_total;            // total Mass
    double alphaDW_ini;        // initial disk wind parameter
    int    wind_alpha_model_;  // 0 = Suzuki2016 (default) , 1 = Enomoto2026 self-consistent alpha_DW
    double lambda_w_;          // magnetic lever arm of the disk wind (Tabone et al. 2022), must be > 1
    double inv_lambda_w_m1_;   // 1/(lambda_w_ - 1)
    double Cw_;
    double inv_Cw_;
    DoubleVector gas_dustrate; // gas dust massrate
    DoubleVector sigmag_;       // gas surface density [g/cm^2]
    DoubleVector sigmag_bod_;      // 
    DoubleVector dsigmagdt_w;   // wind mass loss rate for gas
    DoubleMatrix dsigmaddt_w;   // wind mass loss rate for dust
    DoubleVector v_dsigmaddt_w; // wind mass loss rate for vapor dust
    DoubleVector M_in;          // inner mass
    
    // Vapor
    
    DoubleVector vd_r_v;
    DoubleVector vd_r_diff_v;
    DoubleVector vd_r_tot_v;
    DoubleVector fdg_ratio_vapor;
    DoubleVector dfdg_dr_v;

    DoubleVector dMdt_;         // total mass loss rate
    DoubleVector T_;            // gas temperature [K]
    DoubleVector T_irr;         // irradiation gas temperature [K]
    DoubleVector T_irr4;         // irradiation gas temperature [K]
    DoubleVector cs_;           // sound speed [cm/s]
    DoubleVector Omega_;        // angular speed [/s]
    DoubleVector Hg_;           // gas scale height [cm]
    DoubleVector rhog_mid_;     // gas mass density at disk midplane [g/cm^3]
    DoubleMatrix rhog;          // gas mass density at disk[g/cm^3]
    DoubleVector ng_mid_;       // gas number density at disk midplane [/cm^3]
    DoubleMatrix ng_;           // gas number density in disk[/cm^3]
    DoubleVector vg_r_;         // gas radial velocity [cm/s]
    DoubleVector vgr_bod_;         // gas radial velocity at boundary
    DoubleVector vg_wind_;         // gas radial velocity [cm/s]
    DoubleVector vg_phi_;       // gas radial velocity [cm/s]
    DoubleVector j_;
    DoubleVector dvg_phi;
    DoubleVector vg_rtest_;
    DoubleVector N_torque;      //
    DoubleVector dNdr_torque;      //
    DoubleVector dOmegadr_;     
    DoubleMatrix P_;            // gas pressure
    DoubleVector Ptot_;        // total gas pressure
    DoubleVector dP_dr;            // gas pressure gradient
    DoubleVector dPinteg_dr;            // gas pressure graditent
    DoubleVector dPtot_dr;            // gas pressure
    DoubleVector dPintegtot_dr;      // total gas pressure gradient
    DoubleVector P_mag;         // Magnetic field pressure
    DoubleVector P_integral;    // gas pressure
    DoubleVector Ptot_integral;    // gas pressure
    DoubleVector eta_P_;        // power index of Pressure gradient
    DoubleVector vK_;           // keplerian velocity [cm/s]
    DoubleVector Q_toomre;      // Toomre Q
    double       alpha_;        // turbulent alpha parameter (value)
    DoubleVector alpha_MRI;     // turbulent alpha parameter for MRI (array)
    DoubleVector alpha_GI;      //  turbulent alpha parameter for GI (array)
    DoubleVector alpha_phiz;    // turbulent alpha parameter for disk wind (array)
    DoubleVector alpha_tot;      // turbulent alpha parameter  (array)
    DoubleVector alpha_MRIite;    // turbulent alpha parameter for iteration
    // Enomoto et al. (2026) wind-driven alpha diagnostics (wind_alpha_model_ == 1)
    DoubleVector windDW_zdead_;  // Ohmic dead-zone boundary height [Hg]
    DoubleVector windDW_zsurf_;  // disk surface height beta=1 [Hg]
    DoubleVector windDW_Delta_;  // normalized active-layer thickness Delta (eq.21)
    DoubleVector windDW_Ameff_;  // effective ambipolar Elsasser number (eq. sec4.2)
    DoubleVector alpha_DW_;    // scaled alpha_DW from eq.(27) (before geometric conversion)
    DoubleVector vg_T2_;        // square of gas turbulence velocity [(cm/s)^2]
    DoubleVector nu_;           // turbulent viscosity [cm^2/s]
    DoubleVector Re_;           // Reynolds number
    DoubleVector inv_sqrt_Re_;  // = 1 / sqrt(Re_)

    DoubleMatrix Roe_mean_;     // Roe mean
    DoubleMatrix Roe_mean_dust; // Roe mean for dust
    
    double       eta_cap;       // cap resistivity
    DoubleMatrix etaO;          // Ohmic resistivity   
    DoubleMatrix etaA;          // Ambipolar resistivity 
    DoubleVector etaO_R;        // radial distribution of the Ohmic diffusivity
    DoubleVector etaA_R;        // radial distribution of the Ambipolar diffusivity
    DoubleVector eta_total;     //  radial distribution of the total diffusivity
    DoubleMatrix phi_prev_;     // per-(i,k) warm-start cache for CalcAlphaParameter's phi Newton solve

    DoubleVector nu_mag;        // magnetisation
    DoubleVector numag_eq;      // 

    DoubleVector A_;
    DoubleVector lambda1_;
    DoubleVector lambda2_;
    DoubleVector fluxg_;
    // arrays for coagulation 
    DoubleMatrix mass_transfer_fraction_;  // array used for coagulation
    IntMatrix    mass_transfer_index_c1_;  // array used for coagulation 
    IntMatrix    mass_transfer_index_c2_;  // array used for coagulation

    double zeta_CR0_;               // cosmic-ray ionization rate at r=1au, z=0 [s^-1] (CalcCosmicRayIonizationrate; input "zeta_CR0", default 1.0e-17)

    // fragmentation parameter
    bool on_fragmentation_;        // true : fragmentation on, false : fragmentation off
    bool on_growth_;              // true : dust growth on  , false : dust growth off
    double v_col_crit_;            // critical collision velocity [cm/s]
    double v_col_sili_;            // silicate collision velocity [cm/s]
    double T_snow_;                // snowline temperature [K] (icy/silicate v_frag switch center; width = 0.1*T_snow_). input "T_snow", default 160
    double v_coag_;                // upper velocity for coagulation [cm/s]
    double v_frag_;                // lower velocity for fragmentation [cm/s]
    double frag_slope_;            // power index of mass distribution of fragments
    double frag_frac_;             // frag_frac * mj is the upper mass of mass distribution of fragments
    DoubleVector frag_mass_dist_;  // array used for fragmentation mass distribution
    double m0_gamma_p2_;
    DoubleVector mf_gamma_p2_;
    double _p_beta_; // plasma beta in midplane

    // other variables
    double t_re;
    bool restart_;
    bool on_thermal_;            // true : thermal calculation on, false : thermal calculation off
    bool use_implicit_coag_;     // true : semi-implicit coagulation/fragmentation (dust_coagulation_implicit_note.md), false : explicit (default)
    bool use_full_implicit_coag_; // true : fully implicit (Jacobian + linear solve, coagulation_speedup_literature.md Sec.3 / Brauer et al. 2008), false : explicit (default). Takes priority over use_implicit_coag_ if both are true.
    // Face-flux scheme (face_scheme.cpp). false (default) : original cell-centred
    // velocities + smoothing + Roe-averaged fluxes.
    bool use_face_scheme_ = false;
    DoubleVector lnr_cen_, lnr_bod_;   // ln r at cell centres / faces
    DoubleVector fs_Gv_, fs_Gp_, fs_u_; // face coefficients: viscous, pressure, advective velocity
    DoubleVector vgf_;                 // face gas radial velocity [cm/s]
    DoubleVector dvphif_;              // face gas azimuthal deviation from Kepler [cm/s]
    DoubleVector vwindf_, Af_;         // face v_wind and A
    DoubleVector vBdifff_;             // face magnetic diffusion velocity -eta dB/dr / B [cm/s] (diagnostic)
    DoubleMatrix vdf_;                 // face dust radial drift velocity per bin [cm/s]
    bool use_eta_floor_; // true (default) : historical Elsasser-number(1e-5)/eta_cap(1e20) floor on eta_total in SetBoundaryConditionForMagneticfield, false : use the raw (unfloored) resistivity. Independent of use_implicit_mag (STS_magnetic_diffusion_note.md Sec.12) -- a physics choice, not a performance one.
    std::string radial_grid_type_; // "log" (default), "linear", "root", or "log_floor" -- radial grid spacing used by SetRadiusGrid()
    double radial_grid_floor_dr_; // "log_floor" only: minimum cell width [cm] near rin_, in cm. <=0 means "unset" -> falls back to rin_/(NG+1) (matches ghost-cell spacing, the original coupled behavior).
    double tend_;               // calculation time of simulation [solar year]
    double star_mass_;          // central star mass [solar mass]
    double cfl_;                // cfl (Courant-Friedrichs-Lewy) number
    double rad_diff_fac_;       // radial radiative diffusion strength in the energy eq (input "rad_diff_fac"; 0 = off, 1 = physical)
    std::string outdir_;        // output dust directory
    std::string outdir_mag;        // output magnetic field directory
    double t_mag; // magnetic calculation time

    // turbulent-diffusion face fluxes (declared last to keep the layout of all
    // pre-existing members unchanged): -Dd*Sigma_g*d(Sigma_d/Sigma_g)/dr
    DoubleMatrix fluxd_diff_;
    DoubleVector fluxd_vapor_diff_;

    bool no_dust_; // true : fdg_ <= 0 -- skip all dust motion/growth (see InputParameter())

    // Initial dust size distribution (Initialization() only; restarts read
    // sigmad_ from file and ignore all of this).
    std::string dust_init_mode_;       // "monomer" (default, all mass in bin 0) or "powerlaw"
    double dust_init_powerlaw_q_;      // powerlaw only: n(a) da ~ a^-q  (MRN = 3.5)
    std::string dust_init_amax_mode_;  // powerlaw only: "fixed", "stokes" or "frag_drift"
    double dust_init_amax_;            // powerlaw + amax_mode=fixed: upper cutoff size [cm], same at all r
    double dust_init_St_target_;       // powerlaw + amax_mode=stokes: a_max(r) set so St(a_max(r), r) = this value
    double dust_init_f_frag_;          // powerlaw + amax_mode=frag_drift: f_f in St_frag = f_f v_frag^2/(3 alpha cs^2) (Birnstiel+2012: 0.37)
    double dust_init_f_drift_;         // powerlaw + amax_mode=frag_drift: f_d in St_drift = f_d eps (vK/cs)^2 |dlnP/dlnr|^-1 (Birnstiel+2012: 0.55)
    DoubleVector a_max_init_;          // powerlaw only: a_max(r) [cm] used by the last SetInitialDustDistribution() call (Iteration() convergence check)
};


#endif /* DUST_DISK_HPP_  */