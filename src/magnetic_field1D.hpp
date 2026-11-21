#ifndef MAGNETIC_FIELD1D_HPP
#define MAGNETIC_FIELD1D_HPP

#include "physical_constants.hpp"
#include "alpha.hpp"
#include "typedef.hpp"

class Magnetic_Field2D
{
public:
    Magnetic_Field2D();
    ~Magnetic_Field2D();
    void SetRadiusGridParameter(int ns_r,int ne_r,const DoubleVector &radi,const DoubleVector &radii_b,const DoubleVector &dr_cen,const DoubleVector &dr_bnd,const DoubleVector &inv_R_volume);
    void SetVerticalGridParameter(int ns_z,int ne_z,const DoubleVector &Z_,const DoubleVector &Z_b,const DoubleVector &dz_cen,const DoubleVector &dz_bnd);
    void Set_Differential_Coefficient(const DoubleMatrix &coeff_R,const DoubleMatrix &coeff2_R);
    void SetDiskParameter(const DoubleMatrix &n_g,const DoubleVector &vgr,const DoubleMatrix &Roe_,const DoubleMatrix &P,const DoubleVector &Omega,const DoubleVector &Hg);
    void SetMagneticFieldParameter(const DoubleMatrix &eta_Ohm,const DoubleMatrix &eta_Am,const DoubleVector &eta_,const DoubleVector &B_,double B_out);
    void CalculateDeltaTime(Mag_field *mag);
    void Calc_Mag(Mag_field *mag,double dt);
    void CalculationMRIturbulence(Mag_field *mag);
    void CalculationMRIturbulence2(Mag_field *mag);

    // --- Super-time-stepping (SUPERSEDED -- see note/STS_magnetic_diffusion_note.md) ---
    // Tried two STS variants (Alexiades-Amiez-Gremaud, then Runge-Kutta-
    // Gegenbauer/arXiv:2606.13398) for the resistive-diffusion part of the
    // induction equation. Both are only stable up to some finite stage count,
    // and eta varies by ~18 orders of magnitude (combined with the log grid)
    // across this disk, which pushes the required stage count past even
    // RKG's much larger stability limit -- both diverged in practice. Left
    // here (commented out) for reference; replaced by an implicit
    // (unconditionally stable) tridiagonal solve below.
    // void SetSuperTimeStepping(bool use_sts, double alpha = 10.0);
    // bool IsSuperTimeSteppingEnabled() const { return use_sts_; }

    // Implicit (backward Euler, Thomas-algorithm tridiagonal solve) treatment
    // of Ohmic/ambipolar diffusion in the induction equation. Unconditionally
    // stable regardless of how large eta or dt/dr^2 get, unlike any STS
    // variant -- the right tool given eta's ~18-orders-of-magnitude spatial
    // range here. Off by default so Calc_Mag's behavior is unchanged unless
    // explicitly enabled.
    void SetImplicitMagneticDiffusion(bool use_implicit) { use_implicit_mag_ = use_implicit; }
    // When implicit diffusion is enabled, t_mag (=dt_mag_, the single-
    // explicit-stage diffusion limit) must NOT be folded into the overall
    // macrostep dt selection -- doing so would cap dt right back down to the
    // value the implicit treatment exists to avoid.
    bool IsImplicitMagneticDiffusionEnabled() const { return use_implicit_mag_; }
    // Face-flux scheme: advective flux -v_f B_f with an upwind MUSCL B_f (van Leer,
    // ln B in ln r) and SSP-RK2 in the implicit path. Off by default.
    void SetMusclAdvectionB(bool on) { use_muscl_B_ = on; }
    // Time weight of the implicit diffusion solve (1 = backward Euler, 0.5 = CN).
    void SetMagTheta(double th) { mag_theta_ = th; }
private:
    void Calculate_Differential(DoubleVector &dfdx , DoubleVector f_ ,DoubleMatrix coeff_);
    void Calculate_Advection_B();
    void Calculate_Diffuse_B();

    // // Sets B_field = Y, refreshes B_bnd/B_diff from it, and returns the
    // // (dt-independent) diffusion flux-divergence D(Y) in D_out.
    // void EvaluateDiffusionOperator(const DoubleVector &Y, DoubleVector &D_out);
    // void Calc_Mag_STS(Mag_field *mag, double dt);

    void Calc_Mag_Implicit(Mag_field *mag, double dt);

    double Func_Am(double beta,double Am);
    double PrimeFunc_Am(double beta,double Am);
    double Newton_Am(double beta,double Am_init);

    void Apply_Gaussianfilter_MRI(int radii,double sigma_Gaus);

    int nr_start;
    int nr_end;
    int nz_start;
    int nz_end;

    int nradius;
    int nzvertical;
    double nu_time;
    double dt_;

    DoubleVector dt_mag_;

    DoubleVector r_;
    DoubleVector r_bnd_;
    DoubleVector dr_cen_;
    DoubleVector dr_bnd_;
    DoubleVector inv_vol_;
    DoubleMatrix coeff_R_;
    DoubleMatrix coeff2_R_;

    DoubleMatrix Roe_mean;

    DoubleVector z_;
    DoubleVector z_bnd_;
    DoubleVector dz_cen_;
    DoubleVector dz_bnd_;

    DoubleVector omega_K;
    DoubleVector hg;
    DoubleVector vgr_;
    DoubleVector sigma_gas;
    DoubleMatrix n_g_;
    DoubleMatrix P_;
    DoubleVector P_tot;

    DoubleVector dfdr;
    DoubleVector B_field;
    DoubleVector B_bnd;
    DoubleVector dBdr_cen;
    DoubleVector dBdr_;
    DoubleVector B_adv;
    DoubleVector B_diff;
    DoubleVector B_flux;

    DoubleMatrix beta_;
    DoubleMatrix eta_O;
    DoubleMatrix eta_A;
    DoubleVector eta_tot_;

    DoubleVector alpha_;
    DoubleMatrix alpha_RZ_;

    // Warm-start cache for Newton_Am: last converged Am per (i,j) cell, so the
    // Newton iteration in CalculationMRIturbulence starts near the root instead
    // of always from the same fixed initial guess (Am_prev_[i][j] <= 0 means "no
    // previous value yet", handled by falling back to the default initial guess).
    DoubleMatrix Am_prev_;

    // bool use_sts_ = false;
    // double sts_alpha_ = 10.0; // Gegenbauer parameter (0.5 = classical RKL)
    bool use_implicit_mag_ = false;
    bool use_muscl_B_ = false;
    double mag_theta_ = 0.55;
    double UpwindB(const DoubleVector &B, int f, double u) const;

    double Bini_;
    double alpha_hydro = 1.0e-4;
    double B_boundary = 1.0e-4;
};

#endif