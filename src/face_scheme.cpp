// ============================================================================
// Face-flux scheme (input flag: face_flux_scheme = true)
//
// Replaces the cell-centred velocity + smoothing + Roe-average approach with
// fluxes built directly on cell faces from the two adjacent cells only:
//
//   gas  : Sigma_g flux at face f (between cells f-1 and f) is split into
//          (a) viscous torque   F = Gv_f (N_f - N_{f-1})          , N = r^3 nu Sigma dOmega/dr
//          (b) pressure (dust back-reaction)
//                               F = Gp_f (Ptot_f - Ptot_{f-1})
//          (c) advection        F = u_f Sigma_upwind (MUSCL)       , u = wind torque + v_wind
//          (a)+(b) are diffusive -> backward Euler (Thomas algorithm),
//          (c) and the wind mass-loss sink are explicit (SSP-RK2).
//   dust : face velocity v_d,f from face gas quantities, upwind MUSCL, SSP-RK2.
//   B    : advective flux uses an upwind MUSCL B_f with the face gas velocity
//          (see magnetic_field1D.cpp), diffusion unchanged except theta.
//
// MUSCL: piecewise-linear reconstruction of ln(Sigma) in ln(r) with the van
// Leer limiter. Works on any (non-uniform) radial grid; on the log grid the
// ln r spacing is uniform. Falls back to first-order upwind next to the domain
// edges and where any of the three cells is at/below the floor.
//
// No smoothing (no Gaussian filter, no centre->face->centre round trip).
// ============================================================================

#include "dust_disk.hpp"
#include <cmath>
#include <algorithm>

namespace
{
inline double VanLeerSlope(double dL, double dR)
{
    if (dL * dR <= 0.0) return 0.0;
    return 2.0 * dL * dR / (dL + dR);
}

// Upwind MUSCL value at face f for a positive field S (accessor S(i)).
// c = upwind cell (f-1 if u >= 0, else f). Real cells are [lo, hi).
template <class Get>
inline double MusclUpwindValue(Get S, int c, int f, int lo, int hi, double floor_c,
                               const DoubleVector &lnrc, const DoubleVector &lnrb)
{
    const double Sc = S(c);
    if (c - 1 < lo || c + 1 >= hi) return Sc;
    const double Sm = S(c - 1);
    const double Sp = S(c + 1);
    if (Sc <= floor_c || Sm <= floor_c || Sp <= floor_c) return Sc;
    const double qc = std::log(Sc);
    const double dL = (qc - std::log(Sm)) / (lnrc[c] - lnrc[c - 1]);
    const double dR = (std::log(Sp) - qc) / (lnrc[c + 1] - lnrc[c]);
    const double s = VanLeerSlope(dL, dR);
    return std::exp(qc + s * (lnrb[f] - lnrc[c]));
}

// Thomas algorithm for a tridiagonal system (sub[0], sup[n-1] unused).
inline void SolveTridiagonal(const DoubleVector &sub, const DoubleVector &dia, const DoubleVector &sup,
                             const DoubleVector &rhs, DoubleVector &sol)
{
    const int n = (int)dia.size();
    DoubleVector cp(n), dp(n);
    cp[0] = sup[0] / dia[0];
    dp[0] = rhs[0] / dia[0];
    for (int k = 1; k < n; ++k)
    {
        const double m = dia[k] - sub[k] * cp[k - 1];
        cp[k] = sup[k] / m;
        dp[k] = (rhs[k] - sub[k] * dp[k - 1]) / m;
    }
    sol[n - 1] = dp[n - 1];
    for (int k = n - 2; k >= 0; --k) sol[k] = dp[k] - cp[k] * sol[k + 1];
}
} // namespace

void DustDisk::SetupFaceScheme()
{
    lnr_cen_.assign(nr_tot, 0.0);
    lnr_bod_.assign(nr_tot + 1, 0.0);
    for (int i = 0; i < nr_tot; ++i) lnr_cen_[i] = std::log(r_cen_[i]);
    for (int i = 0; i <= nr_tot; ++i) lnr_bod_[i] = std::log(r_bod_[i]);

    fs_Gv_.assign(nr_tot + 1, 0.0);
    fs_Gp_.assign(nr_tot + 1, 0.0);
    fs_u_.assign(nr_tot + 1, 0.0);
    vgf_.assign(nr_tot + 1, 0.0);
    dvphif_.assign(nr_tot + 1, 0.0);
    vwindf_.assign(nr_tot + 1, 0.0);
    Af_.assign(nr_tot + 1, 0.0);
    vBdifff_.assign(nr_tot + 1, 0.0);
    vdf_.assign(nr_tot + 1, DoubleVector(nbin_, 0.0));
}

// ---------------------------------------------------------------------------
// Gas: face coefficients, face velocities, and cell-centred velocities
// (cell value = geometric-weight interpolation of its two faces; no smoothing).
// Requires (same as CalcGasVelocity_TM2018): N_torque, Ptot_, P_, lambda1_,
// lambda2_, A_, vK_, cs_, alpha_DW_, vg_wind_, sigmag_.
// ---------------------------------------------------------------------------
void DustDisk::CalcGasFaceVelocity()
{
    const double SIG_FLOOR = 1.0e-4;

    for (int f = nr_s + 1; f < nr_e; ++f)
    {
        const int i = f, im = f - 1;
        const double lam1 = 0.5 * (lambda1_[i] + lambda1_[im]);
        const double lam2 = 0.5 * (lambda2_[i] + lambda2_[im]);
        const double Af   = 0.5 * (A_[i] + A_[im]);
        const double vKf  = std::sqrt(vK_[i] * vK_[im]);
        const double cs2f = cs_[i] * cs_[im];
        const double Sigf = std::sqrt(sigmag_[i] * sigmag_[im]);
        const double Ptf  = std::sqrt(Ptot_[i] * Ptot_[im]);
        const double aDWf = 0.5 * (alpha_DW_[i] + alpha_DW_[im]);
        const double vwf  = 0.5 * (vg_wind_[i] + vg_wind_[im]);
        const double Df   = SQR(1.0 + lam2) + SQR(lam1) * (1.0 + Af);
        const double dr   = r_cen_[i] - r_cen_[im];
        const double rf   = r_bod_[f];

        Af_[f] = Af;
        vwindf_[f] = vwf;

        if (sigmag_[i] <= SIG_FLOOR && sigmag_[im] <= SIG_FLOOR)
        {
            fs_Gv_[f] = fs_Gp_[f] = fs_u_[f] = 0.0;
            vgf_[f] = 0.0;
            dvphif_[f] = 0.0;
            continue;
        }

        const double cv = 2.0 * (1.0 + lam2) / (vKf * Df);                 // multiplies (1/(r Sigma)) dN/dr
        fs_Gv_[f] = cv / (rf * dr);                                          // F_visc = Gv (N_i - N_im)
        const double etaP_per_dP = -0.5 * cs2f / SQR(vKf) * rf / Ptf;        // eta_P = etaP_per_dP * dPtot/dr
        fs_Gp_[f] = 2.0 * lam1 * (1.0 + Af) * vKf / Df * etaP_per_dP * Sigf / dr;   // F_P = Gp (Pt_i - Pt_im)
        fs_u_[f] = -cv * 0.75 * aDWf * cs2f - vwf;                           // advective velocity (wind)

        const double dNdr = (N_torque[i] - N_torque[im]) / dr;
        const double dPdr = (Ptot_[i] - Ptot_[im]) / dr;
        const double Fdiff = fs_Gv_[f] * (N_torque[i] - N_torque[im]) + fs_Gp_[f] * (Ptot_[i] - Ptot_[im]);
        vgf_[f] = Fdiff / Sigf + fs_u_[f];

        const double dNdr_p = dNdr / (rf * Sigf) - 0.75 * aDWf * cs2f;
        const double etaPf = etaP_per_dP * dPdr;
        dvphif_[f] = (dNdr_p * lam1 / vKf - (1.0 + lam2) * vKf * etaPf) / Df;
    }
    // Domain-edge faces: zero-gradient copy of the adjacent interior face.
    // (Their flux is outflow-only and explicit; see AdvanceGasFaceScheme.)
    vgf_[nr_s] = vgf_[nr_s + 1];
    vgf_[nr_e] = vgf_[nr_e - 1];
    dvphif_[nr_s] = dvphif_[nr_s + 1];
    dvphif_[nr_e] = dvphif_[nr_e - 1];
    vwindf_[nr_s] = vg_wind_[nr_s];
    vwindf_[nr_e] = vg_wind_[nr_e - 1];
    Af_[nr_s] = A_[nr_s];
    Af_[nr_e] = A_[nr_e - 1];
    fs_Gv_[nr_s] = fs_Gp_[nr_s] = fs_u_[nr_s] = 0.0;
    fs_Gv_[nr_e] = fs_Gp_[nr_e] = fs_u_[nr_e] = 0.0;

    // Cell-centred values (used by dust drift/coagulation, energy, output).
    for (int i = nr_s; i < nr_e; ++i)
    {
        if (sigmag_[i] <= SIG_FLOOR)
        {
            vg_r_[i] = 0.0;
            dvg_phi[i] = 0.0;
        }
        else
        {
            const double w = vg_w_fc_[i];
            vg_r_[i]   = (1.0 - w) * vgf_[i] + w * vgf_[i + 1];
            dvg_phi[i] = (1.0 - w) * dvphif_[i] + w * dvphif_[i + 1];
        }
        vg_phi_[i] = vK_[i] + dvg_phi[i];
    }
    vg_r_[nr_s - 1] = vg_r_[nr_s];
    vg_r_[nr_e] = vg_r_[nr_e - 1];

    // Face gas velocity handed to the magnetic-field advection.
    for (int f = nr_s; f <= nr_e; ++f) vgr_bod_[f] = vgf_[f];

    // Diagnostic: magnetic diffusion velocity at faces, v_B,diff = -eta_f (dB/dr)_f / B_f
    // (two-point difference, Roe-weighted eta_f and B_f as in the code's diffusion flux).
    for (int f = nr_s + 1; f < nr_e; ++f)
    {
        const double s = std::sqrt(sigmag_[f]), t = std::sqrt(sigmag_[f - 1]);
        const double w0 = s / (s + t), w1 = 1.0 - w0;
        const double etaf = w0 * eta_total[f] + w1 * eta_total[f - 1];
        const double Bf = w0 * B_[f] + w1 * B_[f - 1];
        const double dBdr = (B_[f] - B_[f - 1]) / (r_cen_[f] - r_cen_[f - 1]);
        vBdifff_[f] = (Bf != 0.0) ? -etaf * dBdr / Bf : 0.0;
    }
    vBdifff_[nr_s] = vBdifff_[nr_s + 1];
    vBdifff_[nr_e] = vBdifff_[nr_e - 1];
}

// ---------------------------------------------------------------------------
// Gas update over dt:
//   1) explicit SSP-RK2 : advection (u_f, MUSCL upwind), domain-edge outflow
//      (total face velocity, first-order upwind), wind mass-loss sink
//   2) implicit backward Euler : viscous + pressure diffusion (tridiagonal)
// Coefficients (Gv, Gp, u, k = N/Sigma, q = P_mid/Sigma, P_B) frozen at the
// start of the step.
// ---------------------------------------------------------------------------
void DustDisk::AdvanceGasFaceScheme(double dt)
{
    const double SIG_FLOOR = 1.0e-4;
    const int N = nr_e - nr_s;

    auto explicit_rate = [&](const DoubleVector &S, DoubleVector &L)
    {
        DoubleVector F(nr_tot + 1, 0.0);
        auto get = [&](int i) { return S[i]; };
        for (int f = nr_s + 1; f < nr_e; ++f)
        {
            const double u = fs_u_[f];
            const int c = (u >= 0.0) ? f - 1 : f;
            F[f] = u * MusclUpwindValue(get, c, f, nr_s, nr_e, 0.0, lnr_cen_, lnr_bod_);
        }
        F[nr_s] = (vgf_[nr_s] < 0.0) ? vgf_[nr_s] * S[nr_s] : 0.0;
        F[nr_e] = (vgf_[nr_e] > 0.0) ? vgf_[nr_e] * S[nr_e - 1] : 0.0;
        for (int i = nr_s; i < nr_e; ++i)
            L[i] = -2.0 * (r_bod_[i + 1] * F[i + 1] - r_bod_[i] * F[i]) * inv_vol_[i] - dsigmagdt_w[i];
    };

    // ---- 1) explicit SSP-RK2 ----
    DoubleVector S0 = sigmag_, S1 = sigmag_, L(nr_tot, 0.0);
    explicit_rate(S0, L);
    for (int i = nr_s; i < nr_e; ++i) S1[i] = std::max(S0[i] + dt * L[i], SIG_FLOOR);
    explicit_rate(S1, L);
    for (int i = nr_s; i < nr_e; ++i)
        sigmag_[i] = std::max(0.5 * S0[i] + 0.5 * (S1[i] + dt * L[i]), SIG_FLOOR);

    // ---- 2) implicit diffusion: S + g [r_{i+1} F_{i+1} - r_i F_i] = S* ----
    // F_f = alpha_f S_f - beta_f S_{f-1} + gamma_f (interior faces only)
    DoubleVector k(nr_tot, 0.0), q(nr_tot, 0.0), pB(nr_tot, 0.0);
    for (int i = nr_s; i < nr_e; ++i)
    {
        k[i]  = (S0[i] > 0.0) ? N_torque[i] / S0[i] : 0.0;
        q[i]  = (S0[i] > 0.0) ? P_[i][nz_s] / S0[i] : 0.0;
        pB[i] = Ptot_[i] - P_[i][nz_s];
    }
    DoubleVector al(nr_tot + 1, 0.0), be(nr_tot + 1, 0.0), ga(nr_tot + 1, 0.0);
    for (int f = nr_s + 1; f < nr_e; ++f)
    {
        al[f] = fs_Gv_[f] * k[f]     + fs_Gp_[f] * q[f];
        be[f] = fs_Gv_[f] * k[f - 1] + fs_Gp_[f] * q[f - 1];
        ga[f] = fs_Gp_[f] * (pB[f] - pB[f - 1]);
    }
    DoubleVector sub(N), dia(N), sup(N), rhs(N), sol(N);
    for (int i = nr_s; i < nr_e; ++i)
    {
        const int kk = i - nr_s;
        const double g = 2.0 * dt * inv_vol_[i];
        sub[kk] = (kk == 0)     ? 0.0 : g * r_bod_[i] * be[i];
        sup[kk] = (kk == N - 1) ? 0.0 : g * r_bod_[i + 1] * al[i + 1];
        dia[kk] = 1.0 + g * (-r_bod_[i + 1] * be[i + 1] - r_bod_[i] * al[i]);
        rhs[kk] = sigmag_[i] - g * (r_bod_[i + 1] * ga[i + 1] - r_bod_[i] * ga[i]);
    }
    SolveTridiagonal(sub, dia, sup, rhs, sol);
    for (int i = nr_s; i < nr_e; ++i) sigmag_[i] = std::max(sol[i - nr_s], SIG_FLOOR);
}

// ---------------------------------------------------------------------------
// Dust: face drift velocity per bin from face gas quantities (St at the face =
// geometric mean). Called at the end of CalcDustVelocity_TM2018.
// ---------------------------------------------------------------------------
void DustDisk::CalcDustFaceVelocity()
{
    for (int f = nr_s; f <= nr_e; ++f)
    {
        const int i  = std::min(f, nr_e - 1);
        const int im = std::max(f - 1, nr_s);
        const double A  = Af_[f];
        const double vg = vgf_[f];
        const double dvp = dvphif_[f];
        const double vw = vwindf_[f];
        for (int j = 0; j < nbin_; ++j)
        {
            const double St = std::sqrt(St_[i][j] * St_[im][j]);
            const double St2 = SQR(St);
            vdf_[f][j] = ((1.0 + A * St2) * vg + 2.0 * St * dvp - St2 * vw) / (1.0 + (1.0 + A) * St2);
        }
    }
}

// ---------------------------------------------------------------------------
// Dust (and vapor) advective face fluxes: upwind MUSCL with face velocities.
// Same conventions as CalcDustFlux (fluxd_[f][j], fluxd_vapor[f]; outflow-only
// at the domain edges; zero between two below-floor cells).
// ---------------------------------------------------------------------------
void DustDisk::CalcDustFluxFace()
{
    if (no_dust_) return;

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < nbin_; ++j)
    {
        auto get = [&](int i) { return sigmad_active[i][j]; };
        for (int f = nr_s; f <= nr_e; ++f)
        {
            const int i = std::min(f, nr_e - 1), im = std::max(f - 1, nr_s);
            if (sigmad_active[i][j] < sigmad_floor_[i][j] && sigmad_active[im][j] < sigmad_floor_[im][j])
            {
                fluxd_[f][j] = 0.0;
                continue;
            }
            const double u = vdf_[f][j];
            if (f == nr_s)
                fluxd_[f][j] = (u < 0.0) ? u * sigmad_active[nr_s][j] : 0.0;
            else if (f == nr_e)
                fluxd_[f][j] = (u > 0.0) ? u * sigmad_active[nr_e - 1][j] : 0.0;
            else
            {
                const int c = (u >= 0.0) ? f - 1 : f;
                fluxd_[f][j] = u * MusclUpwindValue(get, c, f, nr_s, nr_e, sigmad_floor_[c][j], lnr_cen_, lnr_bod_);
            }
        }
    }
    CalcVaporFluxFace(sigmad_vapor);
}

void DustDisk::CalcVaporFluxFace(const DoubleVector &Sv)
{
    auto getv = [&](int i) { return Sv[i]; };
    for (int f = nr_s; f <= nr_e; ++f)
    {
        const double u = vgf_[f];
        if (f == nr_s)
            fluxd_vapor[f] = (u < 0.0) ? u * Sv[nr_s] : 0.0;
        else if (f == nr_e)
            fluxd_vapor[f] = (u > 0.0) ? u * Sv[nr_e - 1] : 0.0;
        else
        {
            const int c = (u >= 0.0) ? f - 1 : f;
            fluxd_vapor[f] = u * MusclUpwindValue(getv, c, f, nr_s, nr_e, 0.0, lnr_cen_, lnr_bod_);
        }
    }
}

// Explicit dust-advection time-step limit from the face velocities.
double DustDisk::DustAdvectionDtFace() const
{
    double dt = 1.0e100;
    for (int f = nr_s + 1; f < nr_e; ++f)
    {
        const double dr = std::min(dr_bod_[f - 1], dr_bod_[f]);
        for (int j = 0; j < nbin_; ++j)
        {
            if (sigmad_active[f][j] <= sigmad_floor_[f][j] && sigmad_active[f - 1][j] <= sigmad_floor_[f - 1][j]) continue;
            const double u = std::fabs(vdf_[f][j]);
            if (u > 0.0) dt = std::min(dt, dr / u);
        }
    }
    return dt;
}

// Explicit gas time-step limit (advection of Sigma_g and of B by the face velocity).
double DustDisk::GasAdvectionDtFace() const
{
    double dt = 1.0e100;
    for (int f = nr_s; f <= nr_e; ++f)
    {
        const int i = std::min(f, nr_e - 1), im = std::max(f - 1, nr_s);
        const double dr = std::min(dr_bod_[i], dr_bod_[im]);
        const double u = std::max(std::fabs(vgf_[f]), std::fabs(fs_u_[f]));
        if (u > 0.0) dt = std::min(dt, dr / u);
    }
    return dt;
}

// One dust-advection sub-step of length h with SSP-RK2 (fluxes recomputed at
// each stage; face velocities and the wind sink frozen).
void DustDisk::AdvanceDustAdvectionFace(double h)
{
    if (no_dust_) return;
    DoubleMatrix S0 = sigmad_active;
    for (int stage = 0; stage < 2; ++stage)
    {
        CalcDustFluxFace();
        #pragma omp parallel for schedule(static)
        for (int i = nr_s; i < nr_e; ++i)
        {
            const double vol_inv = inv_vol_[i];
            for (int j = 0; j < nbin_; ++j)
            {
                const double Lij = -2.0 * (r_bod_[i + 1] * fluxd_[i + 1][j] - r_bod_[i] * fluxd_[i][j]) * vol_inv
                                   - dsigmaddt_w[i][j];
                if (stage == 0)
                    sigmad_active[i][j] = S0[i][j] + h * Lij;
                else
                    sigmad_active[i][j] = 0.5 * S0[i][j] + 0.5 * (sigmad_active[i][j] + h * Lij);
            }
        }
    }
}

// Vapor advection over dt with SSP-RK2 at the face gas velocity.
void DustDisk::AdvanceVaporFace(double dt)
{
    if (no_dust_) return;
    DoubleVector S0 = sigmad_vapor, S1 = sigmad_vapor;
    CalcVaporFluxFace(S0);
    for (int i = nr_s; i < nr_e; ++i)
        S1[i] = S0[i] - 2.0 * dt * (r_bod_[i + 1] * fluxd_vapor[i + 1] - r_bod_[i] * fluxd_vapor[i]) * inv_vol_[i]
                - v_dsigmaddt_w[i] * dt;
    CalcVaporFluxFace(S1);
    for (int i = nr_s; i < nr_e; ++i)
        sigmad_vapor[i] = 0.5 * S0[i] + 0.5 * (S1[i] - 2.0 * dt * (r_bod_[i + 1] * fluxd_vapor[i + 1] - r_bod_[i] * fluxd_vapor[i]) * inv_vol_[i]
                                               - v_dsigmaddt_w[i] * dt);
}
