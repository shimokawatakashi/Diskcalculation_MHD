#include "magnetic_field1D.hpp"
#include "alpha.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <stdio.h>
#include <limits>

Magnetic_Field2D :: Magnetic_Field2D()
{

}

Magnetic_Field2D :: ~Magnetic_Field2D()
{

}

void Magnetic_Field2D :: SetRadiusGridParameter(int ns_r,int ne_r,const DoubleVector &radi,const DoubleVector &radii_b,const DoubleVector &dr_cen,const DoubleVector &dr_bnd,const DoubleVector &Inv_R_volume)
{
    nradius = radi.size();
    nr_start = ns_r;
    nr_end = ne_r;

    r_ = radi;
    r_bnd_ = radii_b;
    inv_vol_ = Inv_R_volume;
    dr_cen_ = dr_cen;
    dr_bnd_ = dr_bnd;
}

void Magnetic_Field2D :: SetVerticalGridParameter(int ns_z,int ne_z,const DoubleVector &Z_,const DoubleVector &Z_b,const DoubleVector &dz_cen,const DoubleVector &dz_bnd)
{
    nzvertical = Z_.size();
    nz_start = ns_z;
    nz_end = ne_z;

    z_ = Z_;
    z_bnd_ = Z_b;
    dz_cen_ = dz_cen;
    dz_bnd_ = dz_bnd;
}

void Magnetic_Field2D :: Set_Differential_Coefficient(const DoubleMatrix &coeff_R,const DoubleMatrix &coeff2_R)
{
    coeff_R_ = coeff_R;
    coeff2_R_ = coeff2_R;
}

void Magnetic_Field2D :: SetDiskParameter(const DoubleMatrix &n_g,const DoubleVector &vgr,const DoubleMatrix &Roe_,const DoubleMatrix &P,const DoubleVector &Omega,const DoubleVector &Hg)
{
    P_tot.resize(nradius,0.0);

    n_g_ = n_g;
    vgr_ = vgr;
    Roe_mean = Roe_;
    P_ = P;
    omega_K = Omega;
    hg = Hg;
}

void Magnetic_Field2D :: SetMagneticFieldParameter(const DoubleMatrix &eta_Ohm,const DoubleMatrix &eta_Am,const DoubleVector &eta_,const DoubleVector &B_,double B_out)
{
    dfdr.resize(nradius,0.0);

    dBdr_cen.resize(nradius,0.0);
    dBdr_.resize(nradius+1,0.0);
    B_field.resize(nradius,0.0);
    B_bnd.resize(nradius+1,0.0);
    beta_.resize(nradius,DoubleVector(nzvertical,0.0));

    B_adv.resize(nradius,0.0);
    B_diff.resize(nradius,0.0);
    dt_mag_.resize(nradius,0.0);

    eta_tot_.resize(nradius,0.0);

    B_boundary = B_out;
    B_field = B_;
    eta_O = eta_Ohm;
    eta_A = eta_Am;
    eta_tot_ = eta_;
}

void Magnetic_Field2D :: CalculateDeltaTime(Mag_field *mag)
{
    for(int i=nr_start;i<nr_end;i++)
    {
        dt_mag_[i] =  SQR(dr_bnd_[i])/eta_tot_[i];
    }

    dt_ = *std :: min_element(dt_mag_.begin()+nr_start,dt_mag_.end()-nr_start);
    mag ->Delta_time = dt_;
}

void Magnetic_Field2D :: Calculate_Differential(DoubleVector &dfdx , DoubleVector f_ , DoubleMatrix coeff_)
{
    int N = dfdx.size();
    for(int i=1;i<nradius-1;i++)
    {
        dfdx[i] = coeff_[i][0] * f_[i+1];
        dfdx[i] += coeff_[i][1] * f_[i];
        dfdx[i] += coeff_[i][2] * f_[i-1];
    }
}

// Upwind MUSCL value of B at face f (between cells f-1 and f) for face velocity u.
// van Leer limiter on ln B vs ln r; first order next to the domain edges or if
// B <= 0 in any of the three cells.
double Magnetic_Field2D :: UpwindB(const DoubleVector &B, int f, double u) const
{
    const int c = (u >= 0.0) ? f - 1 : f;
    const int cc = std::min(std::max(c, nr_start), nr_end - 1);
    if (cc - 1 < nr_start || cc + 1 >= nr_end) return B[cc];
    const double Bm = B[cc-1], Bc = B[cc], Bp = B[cc+1];
    if (Bm <= 0.0 || Bc <= 0.0 || Bp <= 0.0) return Bc;
    const double xc = std::log(r_[cc]);
    const double dL = (std::log(Bc) - std::log(Bm)) / (xc - std::log(r_[cc-1]));
    const double dR = (std::log(Bp) - std::log(Bc)) / (std::log(r_[cc+1]) - xc);
    const double sl = (dL*dR <= 0.0) ? 0.0 : 2.0*dL*dR/(dL+dR);
    return std::exp(std::log(Bc) + sl * (std::log(r_bnd_[f]) - xc));
}

void Magnetic_Field2D :: Calculate_Advection_B()
{
    for(int i=nr_start;i <= nr_end;i++)
    {
        // B_bnd (Roe average) is also used by the diffusion flux -- keep it.
        B_bnd[i] =  Roe_mean[i][0]*B_field[i] + Roe_mean[i][1]*B_field[i-1];
        B_adv[i] =  use_muscl_B_ ? -vgr_[i]*UpwindB(B_field, i, vgr_[i]) : -vgr_[i]*B_bnd[i];
    }
}
void Magnetic_Field2D :: Calculate_Diffuse_B()
{
    double eta_mean;

    for(int i=nr_start;i <= nr_end;i++)
    {
        eta_mean = Roe_mean[i][0]*eta_tot_[i] + Roe_mean[i][1]*eta_tot_[i-1];
        if(B_field[i] == B_field[i-1])
        {
            dBdr_[i] = 0.0;
        }
        else
        {
            dBdr_[i] = coeff2_R_[i][0]*B_field[i] + coeff2_R_[i][1]*B_bnd[i] + coeff2_R_[i][2]*B_field[i-1];
        }
        B_diff[i] = eta_mean*dBdr_[i];
    }
}

void Magnetic_Field2D :: Calc_Mag(Mag_field *mag,double dt)
{
    if (use_implicit_mag_)
    {
        Calc_Mag_Implicit(mag, dt);
        return;
    }

    int j;
    double vB_;
    double B2;

    B_flux.resize(nradius+1,0.0);

    Calculate_Differential(dBdr_cen,B_field,coeff_R_);
    Calculate_Advection_B();
    Calculate_Diffuse_B();

    // Insulating (zero diffusive flux) condition at both domain faces: the
    // field can only leave by advection (outflow-clamped below), and diffusion
    // conserves the total flux exactly.
    B_diff[nr_start] = 0.0;
    B_diff[nr_end]   = 0.0;

    // Calculate magnetic field flux
    for(int i= nr_start;i<=nr_end;i++)
    {
        B_flux[i] = B_adv[i] + B_diff[i];
    }
    // 境界での磁場のfluxの計算
    // inner boundary
    B_flux[nr_start] = ( B_flux[nr_start] > 0.0) ? B_flux[nr_start] : 0.0;
    // outer boundary
    B_flux[nr_end] = ( B_flux[nr_end] < 0.0) ? B_flux[nr_end] : 0.0;

    // 磁場の計算
    for(int i=nr_start;i<nr_end;i++)
    {
        B_field[i] = B_field[i] + 2.0*dt*inv_vol_[i] * ( r_bnd_[i+1]*B_flux[i+1] -  r_bnd_[i]*B_flux[i] );
        B2 = SQR(B_field[i]);
        for(int j=nz_start;j<nz_end;j++)
        {
            beta_[i][j] = 8.0*M_PI*P_[i][j]/B2;
        }
    }

    mag -> B_ = B_field;
    mag -> beta_RZ = beta_;
}

// --- Super-time-stepping (SUPERSEDED -- see note/STS_magnetic_diffusion_note.md) ---
// void Magnetic_Field2D :: SetSuperTimeStepping(bool use_sts, double alpha)
// {
//     use_sts_ = use_sts;
//     sts_alpha_ = alpha;
// }

// --- Superseded Alexiades-Amiez-Gremaud (AAG) STS implementation ---
// AAG-STS is only numerically stable for a *small* stage count -- it breaks
// down once N reaches "a few tens" -- so growing N to cover a large dt in
// one cycle (an earlier version of this function) diverges. The fix tried
// next -- a fixed, safe N repeated over many cycles until dt is covered --
// was stable but could need a huge number of cycles when eta spans many
// orders of magnitude, defeating the point of STS. Replaced below by the
// Runge-Kutta-Gegenbauer (RKG) scheme, which stays stable at far larger
// stage counts and reaches dt in a single pass.
//
// int Magnetic_Field2D :: Choose_STS_Nstages(double nu)
// {
//     int N = int(std::lround(0.5/std::sqrt(nu)));
//     return std::max(2, std::min(N, 20));
// }
//
// double Magnetic_Field2D :: STS_SuperstepTime(int N, double dt_exp, double nu)
// {
//     double sqrt_nu = std::sqrt(nu);
//     double growth = (1.0+sqrt_nu)/(1.0-sqrt_nu);
//     double ratio = std::pow(growth, 2.0*N);
//     double frac = std::isfinite(ratio) ? (ratio-1.0)/(ratio+1.0) : 1.0;
//     return dt_exp * (double(N)/(2.0*sqrt_nu)) * frac;
// }
//
// (Calc_Mag_STS's cycling AAG-STS body omitted here; see version control
// history / STS_magnetic_diffusion_note.md for the full previous attempts.)

// void Magnetic_Field2D :: EvaluateDiffusionOperator(const DoubleVector &Y, DoubleVector &D_out)
// {
//     B_field = Y;
//     Calculate_Advection_B();
//     Calculate_Diffuse_B();
//     DoubleVector diff_flux = B_diff;
//     diff_flux[nr_start] = (diff_flux[nr_start] > 0.0) ? diff_flux[nr_start] : 0.0;
//     diff_flux[nr_end]   = (diff_flux[nr_end]   < 0.0) ? diff_flux[nr_end]   : 0.0;
//     D_out.assign(nradius, 0.0);
//     for (int i = nr_start; i < nr_end; i++)
//     {
//         D_out[i] = 2.0*inv_vol_[i] * ( r_bnd_[i+1]*diff_flux[i+1] - r_bnd_[i]*diff_flux[i] );
//     }
// }

// void Magnetic_Field2D :: Calc_Mag_STS(Mag_field *mag, double dt)
// {
//     ... RKG2 stage-recursion body (arXiv:2606.13398 Appendix B) omitted here;
//     see STS_magnetic_diffusion_note.md and version control history. Diverged
//     in practice (grid-scale checkerboard blow-up at the innermost cells)
//     because eta's ~18-orders-of-magnitude spatial range demands a stage
//     count beyond what even RKG's enlarged stability region covers.
// }

// Implicit (theta-method: Crank-Nicolson by default, THETA below) treatment of
// Ohmic/ambipolar diffusion, solved via the Thomas algorithm (tridiagonal,
// O(N)). Unconditionally stable for any dt/eta/dr combination, which is what
// this disk's ~18-orders-of-magnitude spatial range in eta needs -- no STS
// variant has a stage count large enough to remain stable here (see
// STS_magnetic_diffusion_note.md). Crank-Nicolson is used instead of backward
// Euler because BE is numerically dissipative and over-smooths the radial B
// profile at the large dt this solver is run with.
//
// Derivation: Calculate_Diffuse_B's dBdr_[i] = coeff2_R_[i][0]*B[i] +
// coeff2_R_[i][1]*B_bnd[i] + coeff2_R_[i][2]*B[i-1], and B_bnd[i] is itself
// linear in B[i], B[i-1] (Calculate_Advection_B), so dBdr_[i] reduces to
// P[i]*B[i] + Q[i]*B[i-1] with
//   P[i] = coeff2_R_[i][0] + coeff2_R_[i][1]*Roe_mean[i][0]
//   Q[i] = coeff2_R_[i][1]*Roe_mean[i][1] + coeff2_R_[i][2]
// The diffusion flux at face i is then eta_face[i]*(P[i]*B[i]+Q[i]*B[i-1]),
// and the same flux-divergence D(B)_i = 2*inv_vol[i]*(r_bnd[i+1]*flux[i+1] -
// r_bnd[i]*flux[i]) used by the explicit scheme becomes a plain linear
// combination of B[i-1], B[i], B[i+1] -- i.e. exactly tridiagonal. Backward
// Euler (B^{n+1}_i - dt*D(B^{n+1})_i = B^n_i) is then one Thomas-algorithm
// solve. The domain-edge flux clamp (same as the explicit scheme's) is
// evaluated once from the pre-diffusion field to keep the system linear.
void Magnetic_Field2D :: Calc_Mag_Implicit(Mag_field *mag, double dt)
{
    B_flux.resize(nradius+1,0.0);

    // --- Advection: one explicit step over the full dt (its CFL ~dr is far
    // less restrictive than diffusion's ~dr^2, so it doesn't need an
    // implicit treatment) ---
    Calculate_Differential(dBdr_cen,B_field,coeff_R_);
    Calculate_Advection_B();

    DoubleVector adv_flux = B_adv;
    adv_flux[nr_start] = (adv_flux[nr_start] > 0.0) ? adv_flux[nr_start] : 0.0;
    adv_flux[nr_end]   = (adv_flux[nr_end]   < 0.0) ? adv_flux[nr_end]   : 0.0;

    DoubleVector B_after_adv = B_field; // preserves ghost values outside [nr_start,nr_end)
    for (int i = nr_start; i < nr_end; i++)
    {
        B_after_adv[i] = B_field[i] + 2.0*dt*inv_vol_[i] * ( r_bnd_[i+1]*adv_flux[i+1] - r_bnd_[i]*adv_flux[i] );
    }

    if (use_muscl_B_)
    {
        // SSP-RK2 second stage: B = 0.5 B^n + 0.5 (B1 + dt L(B1))
        DoubleVector B0 = B_field;
        B_field = B_after_adv;
        Calculate_Advection_B();
        DoubleVector adv2 = B_adv;
        adv2[nr_start] = (adv2[nr_start] > 0.0) ? adv2[nr_start] : 0.0;
        adv2[nr_end]   = (adv2[nr_end]   < 0.0) ? adv2[nr_end]   : 0.0;
        for (int i = nr_start; i < nr_end; i++)
        {
            B_after_adv[i] = 0.5*B0[i] + 0.5*( B_field[i] + 2.0*dt*inv_vol_[i] * ( r_bnd_[i+1]*adv2[i+1] - r_bnd_[i]*adv2[i] ) );
        }
    }

    B_field = B_after_adv;

    // --- Diffusion: implicit tridiagonal solve, accuracy-limited sub-cycling ---
    // Refresh B_bnd/B_diff for the post-advection field, only to read off
    // the boundary-flux direction (same clamp rule as the explicit scheme).
    // Fixed for the whole macrostep (all sub-cycles below), same as before.
    Calculate_Advection_B();
    Calculate_Diffuse_B();
    // Insulating (zero diffusive flux) condition at both domain faces -- both
    // boundary faces are dropped from the diffusion operator (same effect as
    // B_diff[nr_start]=B_diff[nr_end]=0 in the explicit scheme). Advection
    // keeps its outflow clamp above.
    bool inner_active = false;
    bool outer_active = false;

    DoubleVector P(nradius+1,0.0), Q(nradius+1,0.0), eta_face(nradius+1,0.0);
    for (int i = nr_start; i <= nr_end; i++)
    {
        eta_face[i] = Roe_mean[i][0]*eta_tot_[i] + Roe_mean[i][1]*eta_tot_[i-1];
        P[i] = coeff2_R_[i][0] + coeff2_R_[i][1]*Roe_mean[i][0];
        Q[i] = coeff2_R_[i][1]*Roe_mean[i][1] + coeff2_R_[i][2];
    }

    int n = nr_end - nr_start;

    // Per-cell (dt-independent) operator coefficients: D(B)_i = coefL[k]*B[i-1]
    // + coefD[k]*B[i] + coefU[k]*B[i+1] -- same linear combination as before,
    // just built once here so it can be reused both to estimate the
    // accuracy-limited sub-step count below and to assemble the dt_sub-scaled
    // tridiagonal matrix (eta is fixed for the whole macrostep, so these
    // coefficients don't change between sub-cycles).
    DoubleVector coefL(n,0.0), coefD(n,0.0), coefU(n,0.0);
    for (int i = nr_start; i < nr_end; i++)
    {
        int k = i - nr_start;
        bool face_i_active   = (i   == nr_start) ? inner_active : true;
        bool face_ip1_active = (i+1 == nr_end)   ? outer_active : true;

        double cL       = face_i_active   ? -2.0*inv_vol_[i]*r_bnd_[i]  *eta_face[i]  *Q[i]   : 0.0;
        double cD_left  = face_i_active   ?  2.0*inv_vol_[i]*r_bnd_[i]  *eta_face[i]  *P[i]   : 0.0;
        double cD_right = face_ip1_active ?  2.0*inv_vol_[i]*r_bnd_[i+1]*eta_face[i+1]*Q[i+1] : 0.0;
        double cU       = face_ip1_active ?  2.0*inv_vol_[i]*r_bnd_[i+1]*eta_face[i+1]*P[i+1] : 0.0;

        coefL[k] = cL;
        coefD[k] = cD_right - cD_left;
        coefU[k] = cU;
    }

    // Explicit stability step of this operator: for the grid-scale (checkerboard)
    // mode the eigenvalue of L is ~2*coefD (most negative), so forward Euler is
    // stable only for dt <= 1/|coefD|. Take the tightest over all cells.
    double dt_expl_cfl = std::numeric_limits<double>::max();
    for (int k = 0; k < n; k++)
    {
        if (coefD[k] != 0.0)
            dt_expl_cfl = std::min(dt_expl_cfl, 1.0/std::fabs(coefD[k]));
    }

    // Accuracy-based sub-step count (same idea as the dust-coagulation Gain-
    // based sub-cycling -- see dust_coagulation_implicit_note.md Sec.6 and
    // STS_magnetic_diffusion_note.md Sec.9): how fast would B change under
    // one *explicit* application of this operator, evaluated at the current
    // (pre-diffusion) field? The implicit solve is unconditionally stable
    // regardless of dt -- this bound is purely about not under-resolving a
    // fast-changing cell's transient, not about preventing blow-up. Used only
    // as the *initial guess* for n_sub below -- coagulation_speedup_literature.md
    // Sec.7 found that a pre-estimated n_sub can still be wrong (the estimate
    // is evaluated once, at the start, and can go stale), so this alone is
    // not trusted; see the convergence check in the loop below.
    const double ACC_TOL = 0.05;    // max fractional change in B allowed per sub-step

    // Sub-step cap = dt_macro / dt_explicit_CFL: the point past which sub-cycling
    // this implicit solve is no finer than an explicit scheme would be, so there
    // is nothing more to gain by going further. Computed from the operator so it
    // tracks eta -- at start-up (field far from equilibrium, huge eta) it is
    // large enough not to bail prematurely; in the quasi-steady phase it settles
    // to ~1e3. The 8e6 ceiling is only a guard against a degenerate/NaN operator.
    int N_SUB_MAX = 10000;
    if (std::isfinite(dt_expl_cfl) && dt_expl_cfl > 0.0)
    {
        long est = (long)std::ceil(dt / dt_expl_cfl);
        N_SUB_MAX = (int)std::max(1L, std::min(est, 8000000L));
    }

    // Time discretisation weight for the diffusion solve:
    //   (I - THETA*dt*L) B^{n+1} = (I + (1-THETA)*dt*L) B^n
    // THETA = 1.0 -> backward Euler (L-stable, 1st order, numerically
    // dissipative -> over-smooths the profile at large dt).
    // THETA = 0.5 -> Crank-Nicolson (2nd order, non-dissipative). The stiffest
    // modes get amplification factor -> -1 (they ring rather than decay); the
    // ACC_TOL / CONV_TOL sub-step limiter keeps per-sub-step |lambda*dt_sub|
    // bounded, which keeps those modes away from the ringing limit.
    // Raise slightly (e.g. 0.55) if high-wavenumber ringing appears.
    const double THETA = mag_theta_;   // input mag_theta (default 0.55)

    double dt_acc_min = std::numeric_limits<double>::max();
    for (int i = nr_start; i < nr_end; i++)
    {
        int k = i - nr_start;
        double rate = coefL[k]*B_field[i-1] + coefD[k]*B_field[i] + coefU[k]*B_field[i+1];
        double Bref = std::fabs(B_field[i]);
        if (Bref > 0.0 && rate != 0.0)
        {
            dt_acc_min = std::min(dt_acc_min, Bref/std::fabs(rate));
        }
    }

    int n_sub = 1;
    if (std::isfinite(dt_acc_min) && dt_acc_min > 0.0)
    {
        n_sub = std::max(1, (int)std::ceil(dt / (ACC_TOL*dt_acc_min)));
        n_sub = std::min(n_sub, N_SUB_MAX);
    }

    // Same principle as coagulation_speedup_literature.md Sec.7 (applied here
    // per the user's request to use the same fix for the magnetic solve):
    // don't just trust the pre-estimated n_sub, and don't silently paper over
    // a bad sub-step -- verify, after actually solving each sub-step, that no
    // cell's B changed by more than the intended tolerance, and if it did,
    // redo the *whole* macro dt from the pre-advection-diffusion B with n_sub
    // doubled. Unlike the (nonlinear) coagulation Newton solve, this theta-
    // method diffusion step is linear and (theta >= 0.5) unconditionally stable
    // -- it cannot blow up -- so this check is purely about accuracy
    // (under-resolution, and for CN also keeping the stiffest modes off the
    // g -> -1 ringing limit),
    // not about catching a diverging solve, but the retry mechanism is the
    // same: shrink the step internally to this cell/operator until the
    // observed behavior actually satisfies the tolerance the estimate was
    // supposed to guarantee.
    const double CONV_TOL = 2.0*ACC_TOL; // allow some slack over the estimate's own target before rejecting
    DoubleVector B_field_start = B_field;
    DoubleVector rhs(n,0.0), cp(n,0.0), dp(n,0.0), Bsol(n,0.0);
    DoubleVector lower(n,0.0), diag(n,0.0), upper(n,0.0);
    bool converged = false;
    while (!converged)
    {
        double dt_sub = dt/double(n_sub);
        B_field = B_field_start;

        // Tridiagonal LHS matrix (I - THETA*dt_sub*L), assembled once per n_sub
        // attempt (coefL/coefD/coefU are dt-independent, identical every
        // sub-cycle within this attempt); only the RHS is rebuilt each sub-step.
        for (int k = 0; k < n; k++)
        {
            lower[k] = (k == 0)   ? 0.0 : -THETA*dt_sub*coefL[k];
            diag[k]  = 1.0 - THETA*dt_sub*coefD[k];
            upper[k] = (k == n-1) ? 0.0 : -THETA*dt_sub*coefU[k];
        }

        converged = true;
        for (int s = 0; s < n_sub && converged; s++)
        {
            bool rhs_ok = true;
            for (int i = nr_start; i < nr_end; i++)
            {
                int k = i - nr_start;
                // RHS = (I + (1-THETA)*dt_sub*L) B^n
                double LBn = coefL[k]*B_field[i-1] + coefD[k]*B_field[i] + coefU[k]*B_field[i+1];
                rhs[k] = B_field[i] + (1.0-THETA)*dt_sub*LBn;
                // Implicit part of a known (Dirichlet) ghost -> moved to RHS.
                if (i == nr_start && inner_active) rhs[k] += THETA*dt_sub*coefL[k]*B_field[i-1];
                if (i+1 == nr_end && outer_active) rhs[k] += THETA*dt_sub*coefU[k]*B_field[i+1];
                // coefL/coefD/coefU carry eta_face, which with use_eta_floor=false
                // has no ceiling -- LBn (and hence rhs) can overflow to Inf/NaN
                // directly here, before the Thomas algorithm even starts.
                if (!std::isfinite(rhs[k])) rhs_ok = false;
            }

            // Thomas algorithm. No partial pivoting (tridiagonal), so a
            // degenerate/near-zero pivot (diag[0], or an m below, going
            // through zero) is possible when the operator is spatially
            // ill-conditioned -- e.g. eta_tot_ spanning many orders of
            // magnitude between adjacent cells with use_eta_floor=false,
            // which can happen early on before the resistivity profile has
            // settled. Detect it explicitly (same remedy as everywhere else
            // in this session for a degenerate solve: treat as non-converged,
            // shrink dt_sub, retry) rather than silently dividing by ~0 or
            // letting an overflowed Inf/NaN through -- "NaN/Inf > x" is false
            // under IEEE754, so a bare magnitude-only check would silently
            // accept a broken solve as "converged".
            const double PIVOT_FLOOR = 1.0e-300;
            bool pivot_ok = rhs_ok
                && std::isfinite(diag[0]) && std::fabs(diag[0]) > PIVOT_FLOOR;
            if (pivot_ok)
            {
                cp[0] = upper[0]/diag[0];
                dp[0] = rhs[0]/diag[0];
                pivot_ok = std::isfinite(cp[0]) && std::isfinite(dp[0]);
                for (int k = 1; k < n && pivot_ok; k++)
                {
                    double m = diag[k] - lower[k]*cp[k-1];
                    if (!std::isfinite(m) || std::fabs(m) <= PIVOT_FLOOR) { pivot_ok = false; break; }
                    cp[k] = upper[k]/m;
                    dp[k] = (rhs[k] - lower[k]*dp[k-1])/m;
                    if (!std::isfinite(cp[k]) || !std::isfinite(dp[k])) { pivot_ok = false; break; }
                }
            }
            if (!pivot_ok)
            {
                converged = false;
                break;
            }
            Bsol[n-1] = dp[n-1];
            for (int k = n-2; k >= 0; k--)
            {
                Bsol[k] = dp[k] - cp[k]*Bsol[k+1];
            }

            for (int i = nr_start; i < nr_end; i++)
            {
                int k = i - nr_start;
                // isfinite check first: SolveLinearSystemGaussianElim-style
                // pivot guard above catches an exact-zero divide, but a
                // near-zero (not caught by PIVOT_FLOOR) pivot can still blow
                // Bsol up to Inf/NaN via catastrophic cancellation. Without
                // this, "NaN > x" is false under IEEE754, so the next check
                // would silently accept a NaN Bsol as "converged".
                if (!std::isfinite(Bsol[k]))
                {
                    converged = false;
                    break;
                }
                double Bref = std::max(std::fabs(B_field[i]), B_boundary);
                if (std::fabs(Bsol[k] - B_field[i]) > CONV_TOL*Bref)
                {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;

            for (int i = nr_start; i < nr_end; i++)
            {
                B_field[i] = Bsol[i-nr_start];
            }
        }

        if (!converged)
        {
            if (n_sub >= N_SUB_MAX)
            {
                // Exhausted the sub-cycling budget without meeting the
                // tolerance -- restore the pre-diffusion field rather than
                // accept a known-under-resolved result.
                B_field = B_field_start;
                break;
            }
            n_sub *= 2;
        }
    }

    for (int i = nr_start; i < nr_end; i++)
    {
        double B2 = SQR(B_field[i]);
        for(int j=nz_start;j<nz_end;j++)
        {
            beta_[i][j] = 8.0*M_PI*P_[i][j]/B2;
        }
    }

    mag -> B_ = B_field;
    mag -> beta_RZ = beta_;
}

void Magnetic_Field2D :: CalculationMRIturbulence(Mag_field *mag)
{
    B_field = mag->B_;
    beta_ = mag->beta_RZ;
    
    alpha_.resize(nradius,0.0);
    if ((int)Am_prev_.size() != nradius || (nradius > 0 && (int)Am_prev_[0].size() != nzvertical))
    {
        Am_prev_.assign(nradius, DoubleVector(nzvertical, 0.0));
    }
    # pragma omp parallel
    {
        #pragma omp for
        for(int i=nr_start;i<nr_end;i++)
        {
            P_tot[i] = 0.0;
            double alpha_P_ = 0.0;

            // 論文 Sec.2.7 :
            // (1) 高さごとの判定
            //   alpha(z) = alpha_loc(z) (式52) : lambda_MRI(z) < H_g (式62)
            //            = alpha_hydro          : それ以外
            //   Lambda_O > 1 と beta > beta_min(Am) の効果は式(52)の tanh (err2, err1) に連続的に入っているので
            //   別の閾値は課さない.
            // (2) MRI の層 L_MRI = {z | alpha(z) > alpha_hydro} のうち最も厚い連続領域の厚さ h_MRI と,
            //   その領域の中だけの lambda_MRI の最大値 lambda_max を比べ, lambda_max < h_MRI でなければ
            //   その半径は全高度で alpha_hydro (式(64). lambda < H_g は (1) で課している)
            DoubleVector alpha_z(nzvertical, alpha_hydro);
            DoubleVector lambda_z(nzvertical, 0.0);
            for(int j=nz_start;j<nz_end;j++)
            {
                double alpha = 0.0;
                P_tot[i] += P_[i][j]*dz_bnd_[j];
                double V_A = B_field[i]/sqrt(4.0*M_PI*n_g_[i][j]*GAS_MOLECULAR_MASS); // Alfven wave
                double A   = V_A*V_A/(eta_O[i][j]*omega_K[i]);          //Elsasser数(オーム散逸)
                double A_M = V_A*V_A/(eta_A[i][j]*omega_K[i]);         //Elsasser数(両極性拡散)
                double lambda_MRI = 2.0*M_PI*V_A/(omega_K[i]);        // 式(62) Sano & Miyama 1999
                lambda_z[j] = lambda_MRI;

                // beta <= 1 : 磁気圧優勢で beta_min(Am) > 1 を満たす Am が存在しない
                // (Newton_Am は必ず発散する) ため、計算せず alpha_hydro とする
                // lambda_MRI >= H_g : MRI の波長が円盤に収まらないので alpha_hydro
                // 1/(3 beta) < alpha_hydro : 式(52) は Am, Lambda_O によらず alpha <= alpha_hydro となるので
                //                            計算せず alpha_hydro (Newton_Am も呼ばない)
                if (beta_[i][j] <= 1.0 || !(lambda_MRI < hg[i]) || 1.0/(3.0*beta_[i][j]) < alpha_hydro)
                {
                    alpha_z[j] = alpha_hydro;
                    continue;
                }

                // T.N.Delage et.al 2022(39)を連続にした
                alpha = 1.0/(3.0*beta_[i][j]);
                double alpha_max = alpha;
                // Warm start from the previous substep's converged Am at this
                // cell so Newton_Am typically needs only 1-2 iterations instead
                // of always restarting from a fixed guess.
                double Am_beta = Newton_Am(beta_[i][j], Am_prev_[i][j]);
                Am_prev_[i][j] = Am_beta;
                
                double Am_log = std :: log10(A_M);
                double Am_logmax = std :: log10(Am_beta);

                double alpha_logmax = std :: log10(alpha_max);
                double alpha_logmin = std :: log10(alpha_hydro);

                double err1 = 0.4 * std :: min ((Am_log-Am_logmax) , 0.0);
                double err2 = 0.4 * std :: min ( log10(A), 0.0);

                double x = alpha_logmin + (alpha_logmax - alpha_logmin) * (1.0 - std :: tanh( -err1 )) * (1.0 - std :: tanh(  -err2 ) ); 
                alpha = std :: pow(10.0,x);                       // 式(52)
                alpha_z[j] = std :: max(alpha, alpha_hydro);      // 床値 alpha_hydro (式51)
            }

            // (2) 最も厚い連続領域の厚さ h_MRI と, その領域内の lambda_max
            //   z は H_g 単位, j=nz_start が赤道面 (z=0). 赤道面を含む領域は上下対称なので -z..+z の厚さ (= 2 * 上端)
            double h_MRI = 0.0;
            double lambda_max = 0.0;
            int j = nz_start;
            while (j < nz_end)
            {
                if (alpha_z[j] > alpha_hydro)
                {
                    int j0 = j;
                    double lam = 0.0;
                    while (j < nz_end && alpha_z[j] > alpha_hydro)
                    {
                        lam = std :: max(lam, lambda_z[j]);
                        j++;
                    }
                    double thick = (j0 == nz_start) ? 2.0*z_bnd_[j]
                                                    : (z_bnd_[j] - z_bnd_[j0]);
                    if (thick > h_MRI)
                    {
                        h_MRI = thick;
                        lambda_max = lam;
                    }
                }
                else
                {
                    j++;
                }
            }
            h_MRI *= hg[i]; // [cm]

            bool mri_active = (h_MRI > 0.0) && (lambda_max < h_MRI);

            for (int jj = nz_start; jj < nz_end; jj++)
            {
                double alpha = mri_active ? alpha_z[jj] : alpha_hydro;
                alpha_P_ += alpha * P_[i][jj] * dz_bnd_[jj];      // 式(49)
            }

            alpha_[i] = alpha_P_ / P_tot[i];
        }
    }

    // inner boundary
    for(int i=0;i<nr_start;i++)
    {
        alpha_[i] = alpha_[nr_start];
    }
    // outer boundary
    for(int i=nr_end;i<nradius;i++)
    {
        alpha_[i] = alpha_[nr_end-1];
    }

    // Gaussian filterを設定して計算する
    // Apply_Gaussianfilter_MRI(6,2.0); // 前後3つの数字を用いてfittingする(通常は標準偏差の3倍) / 標準偏差は1.0

    mag -> Alpha_ = alpha_;
    mag -> Alpha_RZ_ = alpha_RZ_; 
    mag -> beta_RZ = beta_;
    return;
}

void Magnetic_Field2D :: CalculationMRIturbulence2(Mag_field *mag)
{
    B_field = mag->B_;
    beta_ = mag->beta_RZ;
    
    alpha_.resize(nradius,0.0);
    
    for(int i=nr_start;i<nr_end;i++)
    {
        P_tot[i] = 0.0;
        double alpha_P_ = 0.0;
        double h_MRI = 0.0;
        double lambda_max = 0.0;
        for(int j=nz_start;j<nz_end;j++)
        {
            // T.N.Delage et.al 2022
            double alpha = 0.0;
            P_tot[i] += P_[i][j]*dz_bnd_[j];
            double V_A = B_field[i]/sqrt(4.0*M_PI*n_g_[i][j]*GAS_MOLECULAR_MASS); // Alfven wave
            double A   = V_A*V_A/(eta_O[i][j]*omega_K[i]);          //Elsasser数(オーム散逸)
            double A_M = V_A*V_A/(eta_A[i][j]*omega_K[i]);         //Elsasser数(両極性拡散)
            
            alpha = 1.0/(3.0*beta_[i][j]);
            double alpha_1 = 50.0/pow(A_M,1.2); 
            double alpha_2 = 8.0/pow(A_M,0.3)+1;
            double beta_min=sqrt(alpha_1*alpha_1 + alpha_2*alpha_2);//Delage et al.(34)
            double lambda = 2.0*M_PI*V_A/(omega_K[i]); // Sano & Miyama 1999;

            // MRI-active requires BOTH the Ohmic (A>1) and ambipolar
            // (beta>beta_min) Elsasser conditions (Delage et al. 2022 A&A
            // 658, A97, eq. 30-31) -- the A>1 check was missing here (it is
            // present in the equivalent per-cell criterion in alpha.cpp),
            // letting tenuous, negligible-mass upper-atmosphere cells with
            // A<<1 register as "active" and dominate lambda_max below.
            if(A > 1.0 && beta_[i][j] > beta_min)
            {
                alpha = (lambda < hg[i] && alpha>alpha_hydro) ?  alpha : alpha_hydro;
                lambda_max = std :: max(lambda,lambda_max);
                h_MRI += dz_bnd_[j];
            }
            else{
                alpha = alpha_hydro;
            }
            alpha_P_ += alpha*P_[i][j]*dz_bnd_[j];
        }
        alpha_[i] = alpha_P_/P_tot[i];

        if(alpha_[i] > alpha_hydro)
        {
            h_MRI = h_MRI*hg[i];
            if(h_MRI < lambda_max){alpha_[i] = alpha_hydro;};
        }
    }
    

    // inner boundary
    for(int i=0;i<nr_start;i++)
    {
        alpha_[i] = alpha_[nr_start];
    }
    // outer boundary
    for(int i=nr_end;i<nradius;i++)
    {
        alpha_[i] = alpha_[nr_end-1];
    }

    // Gaussian filterを設定して計算する
    // Apply_Gaussianfilter_MRI(6,2.0); // 前後3つの数字を用いてfittingする(通常は標準偏差の3倍) / 標準偏差は1.0

    mag -> Alpha_ = alpha_;
    mag -> Alpha_RZ_ = alpha_RZ_; 
    mag -> beta_RZ = beta_;
    return;
}


double Magnetic_Field2D :: Func_Am(double beta,double Am)
{
    double A1,A2,B2;
    A1 = 50.0 * std :: pow(Am,-1.2);
    A2 = 8.0*std :: pow(Am,-0.3)+1.0;
    B2 = std :: sqrt(A1*A1 + A2*A2);

    return B2;
}

double Magnetic_Field2D :: PrimeFunc_Am(double beta,double Am)
{
    double A1,A2,B2;
    A1 = -6.0e3 * std :: pow(Am,-3.4);
    A2 = -4.8*( 8.0*std :: pow(Am,-0.3)+1.0 )*std :: pow(Am,-1.3);
    B2 = A1+A2;

    return B2;
}

double Magnetic_Field2D :: Newton_Am(double beta,double Am_init)
{
    double Am,Am_next;
    double f0,f1;
    double gosa = 1.0e6;
    int n =0;

    // Am_init is normally the previous substep's converged value at this cell
    // (warm start); fall back to the old fixed guess if none is available yet.
    const double AM_COLD = 1.0e-4;
    bool warm = (Am_init > 0.0 && std::isfinite(Am_init));
    Am = warm ? Am_init : AM_COLD;
    // NaN-safe loop condition: a NaN gosa must NOT count as converged
    // (with "gosa > tol" a NaN made the loop exit and return NaN as Am).
    while (!(gosa <= 1.0e-5))
    {
        f0 = Func_Am(beta,Am);
        f1 = PrimeFunc_Am(beta,Am);

        // Am_next = Am - f0/f1;
        Am_next = Am - 2.0 * (SQR(f0) - beta*f0)/f1;

        // Newton diverged (Am -> inf / 0 / NaN). This happens when the warm
        // start Am_prev_ is far from the new root, e.g. when beta jumps by a
        // large factor in one step. Retry once from the cold guess.
        if (!std::isfinite(Am_next) || Am_next == 0.0)
        {
            if (warm)
            {
                warm = false;
                Am = AM_COLD;
                gosa = 1.0e6;
                n = 0;
                continue;
            }
            std :: cout << "Calculation ERROR :: Newton_Am diverged : beta = " << beta << " , Am_init = " << Am_init << std :: endl;
            exit(1);
        }

        gosa = fabs(Am_next - Am)/ fabs(Am);
        Am = Am_next;
        n += 1;

        if (Am < 0)
        {
            Am = -Am;
        }

        if (n>100)
        {
            std :: cout << "Calculation ERROR :: I can't find Am :" <<gosa << std :: endl;
            exit(1);
        }
    }
    return Am;
}

void Magnetic_Field2D :: Apply_Gaussianfilter_MRI(int radii,double sigma_Gaus)
{
    double sum = 0.0;
    DoubleVector Kernel_Gaus;
    DoubleVector alpha_smooth;
    Kernel_Gaus.resize(2*radii+1,0.0);
    alpha_smooth.resize(nradius,0.0);
    // カーネルの設定
    for(int i = -radii;i <= radii;i++)
    {
        Kernel_Gaus[i+radii] = std :: exp( -0.5 * SQR(i/sigma_Gaus) );
        sum += Kernel_Gaus[i+radii];
    }
    // 規格化
    for(double &v : Kernel_Gaus) v = v/sum;
    // 畳み込み
    for(int i=nr_start;i<nr_end;i++)
    {
        alpha_smooth[i] = 0.0;
        if(i < radii)
        {
            for(int j=-radii;j<=radii;j++)
            {
                if(i+j < 0)
                {
                    alpha_smooth[i] += log10(alpha_[0])*Kernel_Gaus[j+radii];
                }
                else
                {
                    alpha_smooth[i] += log10(alpha_[i+j])*Kernel_Gaus[j+radii];
                }
            }
        }
        else if(i > nr_end-radii)
        {
            for(int j=-radii;j<=radii;j++)
            {
                if(i+j >= nradius)
                {
                    alpha_smooth[i] += log10(alpha_[nradius-1])*Kernel_Gaus[j+radii];
                }
                else
                {
                    alpha_smooth[i] += log10(alpha_[i+j])*Kernel_Gaus[j+radii];
                }
            }
        }
        else
        {
            for(int j=-radii;j<=radii;j++)
            {
                alpha_smooth[i] += log10(alpha_[i+j])*Kernel_Gaus[j+radii];
            }
        }
        alpha_[i] = std :: pow(10.0,alpha_smooth[i]);
    }
}
