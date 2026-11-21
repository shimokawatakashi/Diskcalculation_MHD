/*
    dust_growth.cpp

    Dust coagulation and fragmentation (DustDisk's "dust growth" methods),
    split out of dust_disk.cpp for easier review. These are still DustDisk
    member functions -- only their definitions live here; declarations stay
    in dust_disk.hpp.

    Moved from dust_disk.cpp (no behavior change):
      - CalcCoagulationAndFragmentationArrays()
      - CalcDustRelativeVelocity()
      - CoagulateChangeRate()
      - CoagulateChangeRateFragmentation()
      - EstimateDtCoagulation()
      - ApplyCoagulationAndCorrectMass()
*/

#include "dust_disk.hpp"
#include "alpha.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

void DustDisk::CalcCoagulationAndFragmentationArrays()
{
    /*
        set arrays for coagulation
        see Kawasaki et al. 2022 appendix A
    */

    int nbin = nbin_;
    int nbinm = nbin - 1;
    int c1, c2;
    double mi, mj, mij, mij_m, f;

    mass_transfer_fraction_.resize(nbin_, DoubleVector(nbin_));
    mass_transfer_index_c1_.resize(nbin_, IntVector(nbin_));
    mass_transfer_index_c2_.resize(nbin_, IntVector(nbin_));

    for (int i = 0; i < nbin; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {

            mi = m_[i];
            mj = m_[j];
            mij = mi + mj;
            mij_m = mb_[i] + mb_[j];

            for (int k = i; k < nbin; ++k)
            {
                if (mb_[k] <= mij_m && mij_m < mb_[k + 1])
                {
                    c1 = k;
                    c2 = c1 + 1;
                    break;
                }
            }
            c1 = std::min(c1, nbinm);
            c2 = std::min(c2, nbinm);
            mass_transfer_index_c1_[i][j] = c1;
            mass_transfer_index_c2_[i][j] = c2;

            if (c1 == c2)
            {
                mass_transfer_fraction_[i][j] = 1.0;
            }
            else
            {
                if (mj + mb_[i] > mb_[c1 + 1])
                {
                    mass_transfer_fraction_[i][j] = 0.0;
                }
                else if (mj + mb_[i + 1] < mb_[c1 + 1])
                {
                    mass_transfer_fraction_[i][j] = 1.0;
                }
                else
                {
                    mass_transfer_fraction_[i][j] = (m_[c2] - mij) / (m_[c2] - m_[c1]);
                }
            }

        } // end for j
    }     // end for i

    // fragmentation
    double gamma_p2 = frag_slope_ + 2.0;

    frag_mass_dist_.resize(nbin, 0.0);
    mf_gamma_p2_.resize(nbin, 0.0);

    m0_gamma_p2_ = std::pow(mb_[0], gamma_p2);
    for (int i = 0; i < nbin; ++i)
    {
        frag_mass_dist_[i] = std::pow(mb_[i + 1], gamma_p2) - std::pow(mb_[i], gamma_p2);
        mf_gamma_p2_[i] = 1.0 / (std::pow(mb_[i + 1], gamma_p2) - m0_gamma_p2_);
    }

    return;
}

void DustDisk::CalcDustRelativeVelocity(int ir, int id, int jd, double &vd_ij, double &sigma_ij, double &vrel)
{
    /*
    calculate dust relative velocity

    ir : index of disk radius
    id : index of i-th dust bin
    jd : index of j-th dust bin

    vd_ij : relative velocity betawenn dust,i and dust,j due to drift (mean) motion.
    sigma_ij : relative velocity betawenn dust,i and dust,j due to random (Brownian and turbulence) motion.
    vrel : total relative velocity between dust,i and dust,j
    */

    double vbr, delta_vr_drift, vn, delta_vazim_drift, delta_vz_drift,vam_drift;
    double v_turb2, v_turb_I2, v_turb_II2, turb_coeff, StT, StP, Stmin, StTPast;

    // Brownian motion
    vbr = coe_vbr * T_[ir] * (m_[id] + m_[jd]) / (m_[id] * m_[jd]);

    // turbulnce (Ormel & Cuzzi 2007)
    StT = St_[ir][id];
    StP = St_[ir][jd];
    if (StT < StP)
    {
        StT = St_[ir][jd];
        StP = St_[ir][id];
    }
    Stmin = inv_sqrt_Re_[ir];
    StTPast = std::max(Stmin, std::min(1.6 * StT, 1.0));

    v_turb_I2 = (StT - StP)/(StT + StP) * ( SQR(StT)/(StTPast + StT) - SQR(StT)/(1 + StT) - SQR(StP)/(StTPast + StP) + SQR(StP)/(1 + StP) );
    v_turb_II2 = (StTPast - Stmin) * (((StTPast + Stmin) * StT + StTPast * Stmin) / ((StT + StTPast) * (StT + Stmin)) + ( (StTPast + Stmin)*StP + StTPast*Stmin ) /( (StP+StTPast) * (StP+Stmin) ));
    turb_coeff = std :: abs(v_turb_I2)+std :: abs(v_turb_II2);
    v_turb2 = vg_T2_[ir]*turb_coeff; // v_turb2 = v_rurb**2
    // radial drift
    delta_vr_drift = std::abs(vd_r_[ir][id] - vd_r_[ir][jd]);

    // azimuthal drift Okuzumi 2012
    delta_vazim_drift = std::abs(vd_phi_[ir][id] - vd_phi_[ir][jd]);

    // vertical drift Birnstiel 2010
    delta_vz_drift = std::abs(Hd_[ir][id] * std::min(St_[ir][id], 0.5) - Hd_[ir][jd] * std::min(St_[ir][jd], 0.5)) * Omega_[ir];

    // ralative velocity due to drift motion
    vd_ij = std::sqrt(SQR(delta_vr_drift) + SQR(delta_vazim_drift) + SQR(delta_vz_drift));

    // relative velocity due to random motion
    sigma_ij = std::sqrt(pi_8 * (vbr + v_turb2));

    // total relative velocity
    vrel = std::sqrt(vbr + v_turb2 + SQR(vd_ij));

    return;
}

void DustDisk::CoagulateChangeRate(double t, DoubleMatrix &y, DoubleMatrix &dydt)
{
    /*
        t : time
        y : previous dust mass density
        dydt : change rate

        only coagulation
    */

    int ir, i, j, nbin, nbinm, c1, c2;
    double f, kij, dndt, mi, mj, ai, aj, mij, a_sij, y_min;
    double vel, vd_ij, sigma_ij, Av, Bv, Af, Bf, coe1, coe2, C, x, z;

    nbin = nbin_;
    nbinm = nbin - 1;

    // reset dydt
    for (i = nr_s; i < nr_e; ++i)
    {
        std::fill(dydt[i].begin(), dydt[i].end(), 0.0);
    }

    // calculate dydt
    for (ir = nr_s; ir < nr_e; ++ir)
    { // loop for disk radius

        for (i = 0; i < nbin; ++i)
        { // loop for dust bin

            if (y[ir][i] < sigmad_floor_[ir][i])
                continue;
            mi = m_[i];
            ai = a_[i];

            for (j = 0; j <= i; ++j)
            { // loop for dust bin

                if (y[ir][j] < sigmad_floor_[ir][j])
                    continue;

                mj = m_[j];
                aj = a_[j];
                mij = mi + mj;

                // coagulation
                c1 = mass_transfer_index_c1_[i][j];
                c2 = mass_transfer_index_c2_[i][j];
                f = mass_transfer_fraction_[i][j];

                CalcDustRelativeVelocity(ir, i, j, vd_ij, sigma_ij, vel);
                // vd_ij    : relative velocity due to drift motion
                // sigma_ij : relative velocity due to random motion (turbulence and Brownian motion)
                // vel      : total relative velocity
                a_sij = M_PI * SQR(ai + aj); // corss section
                C = 1.0 / std::sqrt(2.0 * M_PI * (SQR(Hd_[ir][i]) + SQR(Hd_[ir][j])));

                if (vd_ij < 5.0e-2 * sigma_ij)
                {
                    // Garaud et al. 2013 eq.(35)
                    x = SQR(v_coag_) / (2.0 * SQR(sigma_ij));
                    kij = C * a_sij * sqrt8_pi * sigma_ij * (1.0 - (1.0 + x) * std::exp(-x));
                }
                else
                {
                    // Garaud et al. 2013 eq.(33), (34)
                    x = 1.0 / (sqrt2 * sigma_ij);
                    z = 1.0 / vd_ij;
                    Av = std::erf((v_coag_ + vd_ij) * x) - std::erf((v_coag_ - vd_ij) * x);
                    Bv = (v_coag_ + vd_ij) * std::exp(-SQR((v_coag_ - vd_ij) * x)) - (v_coag_ - vd_ij) * std::exp(-SQR((v_coag_ + vd_ij) * x));
                    kij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * (2.0 * std::erf(vd_ij * x) - Av) + sqrt1_2pi * sigma_ij * z * (2.0 * vd_ij * std::exp(-SQR(vd_ij * x)) - Bv));
                }

                kij = kij / (mi * mj);
                if (i == j)
                    kij = kij * 0.5;
                dndt = kij * y[ir][i] * y[ir][j];

                dydt[ir][i] -= mi * dndt;
                dydt[ir][j] -= mj * dndt;
                dydt[ir][c1] += f * mij * dndt;
                dydt[ir][c2] += (1.0 - f) * mij * dndt;

            } // end for j
        }     // dne for i
    }         // end for ir

    return;
}

void DustDisk::CoagulateChangeRateFragmentation(double t, DoubleMatrix &y, DoubleMatrix &dydt)
{
    /*
        t : time
        y : previous dust mass density
        dydt : change rate

        This function calculate the change rate of dust mass density due to not only coagulation but also fragmentation.
        Fragmentation model is based on the mothod of Kawasaki et al.2022
    */

    int ir, i, j, k, nbin, nbinm, c1, c2;
    double f, kij, dndt, mi, mj, ai, aj, mij, a_sij, vel_coag;
    double vel, vd_ij, sigma_ij, Av, Bv, Af, Bf, coe1, coe2, C, x, z;
    int kf;
    double mej, mrm, mf, inv_mfrag, inv_vfrag, fij, gamma_p2, m0_gamma_p2;
    // Per-iteration locals (not the v_coag_/v_frag_ members) so each ir is
    // independent under the #pragma omp parallel for below.
    double v_coag, v_frag;

    nbin = nbin_;
    nbinm = nbin - 1;

    gamma_p2 = frag_slope_ + 2.0;
    m0_gamma_p2 = m0_gamma_p2_; // std::pow(mb_[0], gamma_p2);

    // sigmad_floor_ (~1e-50) never actually excludes anything; bins that are only
    // populated at the numerical-noise level of the coagulation front still get
    // collision rates computed, which can crater EstimateDtCoagulation's dt_coag.
    // Skip collisions for bins that are negligible relative to the local total.
    const double REL_THRESH_COAG = 1.0e-20;
    // Absolute floor (g/cm^2): bins below this are skipped regardless of the
    // cell's total, so truly negligible populations (far below any physical
    // significance) don't enter the O(nbin^2) kernel loop / Jacobian.
    const double ABS_THRESH_COAG = 1.0e-20;

    // Each radial cell ir is fully independent (dydt[ir] is only ever written by
    // the iteration that owns it), so this loop parallelizes cleanly over ir --
    // as long as v_coag/v_frag are per-iteration locals rather than the v_coag_/
    // v_frag_ member variables (those would otherwise be a data race between
    // threads processing different ir at once, which is why this was previously
    // left serial).
    #pragma omp parallel for private(ir, i, j, k, c1, c2, f, kij, dndt, mi, mj, ai, aj, mij, a_sij, vel_coag, vel, vd_ij, sigma_ij, Av, Bv, Af, Bf, coe1, coe2, C, x, z, kf, mej, mrm, mf, inv_mfrag, inv_vfrag, fij, v_coag, v_frag) firstprivate(nbin, nbinm, gamma_p2, m0_gamma_p2) schedule(static,4)
    for (ir = nr_s; ir < nr_e; ++ir)
    {
        // loop for disk radius
        std::fill(dydt[ir].begin(), dydt[ir].end(), 0.0);
        if(T_[ir] > 1.5e3)
        {
            continue;
        }

        // reset dydt

        // if (T_[ir] < 1.5e2)
        // {
        //     inv_vfrag = 1.0/v_col_crit_;
        //     v_coag = 0.2 * v_col_crit_;
        //     v_frag = 0.2 * v_col_crit_;
        // }
        // else
        // {
        //     inv_vfrag =  1.0/v_col_sili_;
        //     v_coag = 0.2 * v_col_sili_;
        //     v_frag = 0.2 * v_col_sili_;
        // }
        // w_icy = 1 (T≪Tf, 氷) … 0 (T≫Tf, シリケイト)
        double w   = 0.5*(1.0 - std::tanh((T_[ir] - T_snow_)/(0.1*T_snow_)));
        double vc  = w*v_col_crit_ + (1.0-w)*v_col_sili_;   // 衝突速度を線形ブレンド
        inv_vfrag  = 1.0/vc;
        v_coag     = 0.2*vc;
        v_frag     = 0.2*vc;

        double coag_min_level = std::max(REL_THRESH_COAG * sigmad_tot_[ir], ABS_THRESH_COAG);

        for (i = 0; i < nbin; ++i)
        {
            // loop for dust bin
            if (y[ir][i] < sigmad_floor_[ir][i])
                continue;
            if (y[ir][i] < coag_min_level)
                continue;

            mi = m_[i];
            ai = a_[i];

            for (j = 0; j <= i; ++j)
            {
                // loop for dust bin
                if (y[ir][j] < sigmad_floor_[ir][j])
                    continue;
                if (y[ir][j] < coag_min_level)
                    continue;

                mj = m_[j];
                aj = a_[j];
                mij = mi + mj;

                /* coagulation */
                c1 = mass_transfer_index_c1_[i][j];
                c2 = mass_transfer_index_c2_[i][j];
                f = mass_transfer_fraction_[i][j];

                CalcDustRelativeVelocity(ir, i, j, vd_ij, sigma_ij, vel);
                a_sij = M_PI * SQR(ai + aj);
                C = 1.0 / std::sqrt(2.0 * M_PI * (SQR(Hd_[ir][i]) + SQR(Hd_[ir][j])));

                if (vd_ij < 5.0e-2 * sigma_ij)
                {
                    // Garaud et al. 2013 eq.(35)
                    x = SQR(v_coag) / (2.0 * SQR(sigma_ij));
                    kij = C * a_sij * sqrt8_pi * sigma_ij * (1.0 - (1.0 + x) * std::exp(-x));
                    fij = std::max(C * a_sij * sqrt8_pi * sigma_ij - kij, 0.0);

                    if (fij / kij < 1.0e-3)
                    {
                        kij = C * a_sij * sqrt8_pi * sigma_ij;
                        fij = 0.0;
                    }
                }
                else
                {
                    // Garaud et al. 2013 eq.(33), (34)
                    x = 1.0 / (sqrt2 * sigma_ij);
                    z = 1.0 / vd_ij;
                    Av = std::erf((v_coag + vd_ij) * x) - std::erf((v_coag - vd_ij) * x);
                    Bv = (v_coag + vd_ij) * std::exp(-SQR((v_coag - vd_ij) * x)) - (v_coag - vd_ij) * std::exp(-SQR((v_coag + vd_ij) * x));
                    kij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * (2.0 * std::erf(vd_ij * x) - Av) + sqrt1_2pi * sigma_ij * z * (2.0 * vd_ij * std::exp(-SQR(vd_ij * x)) - Bv));
                    fij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * Av + sqrt1_2pi * sigma_ij * z * Bv);
                }

                // fragmentation off -> perfect sticking: high-velocity collisions (fij) also coagulate
                if (!on_fragmentation_) { kij += fij; fij = 0.0; }
                kij = kij / (mi * mj);
                if (i == j)
                    kij = kij * 0.5;
                dndt = kij * y[ir][i] * y[ir][j];

                dydt[ir][i] -= mi * dndt;
                dydt[ir][j] -= mj * dndt;
                dydt[ir][c1] += f * mij * dndt;
                dydt[ir][c2] += (1.0 - f) * mij * dndt;

                /* fragmentation */
                if (!on_fragmentation_)
                    continue;
                if (i == 0 && j == 0)
                    continue;
                if (fij == 0.0)
                    continue;

                mej = mj * vel * inv_vfrag;
                if (mej < m_[0])
                    continue;
                mrm = mij - mej;
                if (mrm < m_[0])
                {
                    mej = mij;
                    mrm = 0.0;
                }

                mf = std::max(frag_frac_ * mej, mb_[1]);

                // mb_ is sorted ascending, so the bin containing mf can be found
                // with a binary search (O(log nbin)) instead of a linear scan
                // (O(nbin)) -- this loop sits inside the O(nbin^2) (i,j) loop, so
                // the linear version could cost O(nbin^3) per radial cell.
                {
                    auto it = std::upper_bound(mb_.begin(), mb_.begin() + (nbin + 1), mf);
                    int k_found = int(it - mb_.begin()) - 1;
                    if (k_found >= 0 && k_found < nbin)
                    {
                        kf = k_found;
                        mf = mb_[kf + 1];
                    }
                }

                fij = fij / (mi * mj);
                if (i == j)
                    fij *= 0.5;
                dndt = fij * y[ir][i] * y[ir][j];
                dydt[ir][i] -= mi * dndt;
                dydt[ir][j] -= mj * dndt;

                if (mf > mrm)
                {
                    inv_mfrag = mij * dndt * mf_gamma_p2_[kf];
                }
                else
                {
                    // Same binary-search replacement as above, restricted to [kf, nbin)
                    // since mrm can only fall in a bin at or above kf.
                    {
                        auto it2 = std::upper_bound(mb_.begin() + kf, mb_.begin() + (nbin + 1), mrm);
                        int k2_found = int(it2 - mb_.begin()) - 1;
                        if (k2_found >= kf && k2_found < nbin)
                        {
                            dydt[ir][k2_found] += mrm * dndt;
                        }
                    }

                    inv_mfrag = mej * dndt * mf_gamma_p2_[kf];
                }
                if( std :: isnan(vd_r_[ir][j]) == true || std :: isnan(vd_phi_[ir][j]) == true || std :: isnan(vg_r_[ir]) == true || std :: isnan(vg_phi_[ir]) == true || std :: isnan(fdg_ratio[ir][j]) == true || std :: isnan(lambda1_[ir]) == true  || std :: isnan(mf_gamma_p2_[kf]) == true || std :: isnan(vd_r_[ir][i]) == true || std :: isnan(dndt) == true ||  std :: isnan(vd_ij) == true)
                {
                    std :: cerr <<std :: scientific<<ir<<" "<< i <<" "<<j <<" "<<sigmad_[ir][j]<<" "<<a_[i]<<" "<< a_[j] <<" "<<vd_r_[ir][j] <<" "<<vd_phi_[ir][j]<<" "<<sigmag_[ir]<<" "<<vg_r_[ir]<<" "<<vg_phi_[ir]<<" "<<dNdr_torque[ir]<<" "<<lambda1_[ir]<<" "<<lambda2_[ir]<<" "<<vg_wind_[ir]<<" "<<St_[ir][j]<<" "<<fdg_ratio[ir][j]<<" "<<A_[ir]<<" "<<vd_ij<<" "<<dndt<<" "<<sigma_ij
                        <<" Hd_i="<<Hd_[ir][i]<<" Hd_j="<<Hd_[ir][j]<<" C="<<C<<" kij="<<kij<<" fij="<<fij<<" a_sij="<<a_sij<<" St_i="<<St_[ir][i]<< std :: endl;
                    exit(1);
                }

                for (k = 0; k <= kf; ++k)
                {
                    dydt[ir][k] += frag_mass_dist_[k] * inv_mfrag;
                }
            }
        }
    }
    return;
}

// Semi-implicit counterpart of CoagulateChangeRateFragmentation(): same
// pairwise (i,j) loop and physics, but instead of writing one combined rate
// dydt[i] = Gain_i - Loss_i, it separates Gain_i (explicit, evaluated at the
// old y) from the *rate coefficient* LossRate_i such that Loss_i =
// LossRate_i * y_i -- i.e. every "dydt[k] -= m*dndt" below becomes
// "LossRate[k] += (dndt/y[k])" instead. This lets ApplyCoagulationAndCorrectMass_Implicit
// solve y_i implicitly against its own (self-proportional) loss term without
// building any matrix: y_i^{n+1} = (y_i^n + dt*Gain_i) / (1 + dt*LossRate_i),
// which is unconditionally non-negative regardless of dt (see
// dust_coagulation_implicit_note.md). Gain terms are untouched (they aren't
// the source of the depletion instability that dt_coag guards against).
void DustDisk::CoagulateChangeRateFragmentation_Implicit(double t, DoubleMatrix &y, DoubleMatrix &Gain, DoubleMatrix &LossRate)
{
    #pragma omp parallel for schedule(static,4)
    for (int ir = nr_s; ir < nr_e; ++ir)
    {
        CoagulateChangeRateFragmentation_Implicit_Cell(ir, y, Gain, LossRate);
    }
    return;
}

// Per-cell body of CoagulateChangeRateFragmentation_Implicit(), split out so
// that ApplyCoagulationAndCorrectMass_Implicit() can re-invoke it for a
// single cell mid-sub-cycle (refreshing Gain[ir]/LossRate[ir] from the
// current, evolving y instead of the stale y^n snapshot -- see the "residual
// of the semi-implicit method itself" discussion in
// dust_coagulation_implicit_note.md Sec.8) without touching any other cell.
void DustDisk::CoagulateChangeRateFragmentation_Implicit_Cell(int ir, const DoubleMatrix &y, DoubleMatrix &Gain, DoubleMatrix &LossRate)
{
    int i, j, k, nbin, nbinm, c1, c2;
    double f, kij, dndt, mi, mj, ai, aj, mij, a_sij, vel_coag;
    double vel, vd_ij, sigma_ij, Av, Bv, Af, Bf, coe1, coe2, C, x, z;
    int kf;
    double mej, mrm, mf, inv_mfrag, inv_vfrag, fij, gamma_p2, m0_gamma_p2;
    double v_coag, v_frag;

    nbin = nbin_;
    nbinm = nbin - 1;

    gamma_p2 = frag_slope_ + 2.0;
    m0_gamma_p2 = m0_gamma_p2_;

    const double REL_THRESH_COAG = 1.0e-20;
    // Absolute floor (g/cm^2): bins below this are skipped regardless of the
    // cell's total, so truly negligible populations (far below any physical
    // significance) don't enter the O(nbin^2) kernel loop / Jacobian.
    const double ABS_THRESH_COAG = 1.0e-20;

    std::fill(Gain[ir].begin(), Gain[ir].end(), 0.0);
    std::fill(LossRate[ir].begin(), LossRate[ir].end(), 0.0);
    if(T_[ir] > 1.5e3)
    {
        return;
    }

        // if (T_[ir] < 1.5e2)
        // {
        //     inv_vfrag = 1.0/v_col_crit_;
        //     v_coag = 0.2 * v_col_crit_;
        //     v_frag = 0.2 * v_col_crit_;
        // }
        // else
        // {
        //     inv_vfrag =  1.0/v_col_sili_;
        //     v_coag = 0.2 * v_col_sili_;
        //     v_frag = 0.2 * v_col_sili_;
        // }
        double w   = 0.5*(1.0 - std::tanh((T_[ir] - T_snow_)/(0.1*T_snow_)));
        double vc  = w*v_col_crit_ + (1.0-w)*v_col_sili_;
        inv_vfrag  = 1.0/vc;
        v_coag     = 0.2*vc;
        v_frag     = 0.2*vc;

        double coag_min_level = std::max(REL_THRESH_COAG * sigmad_tot_[ir], ABS_THRESH_COAG);

        for (i = 0; i < nbin; ++i)
        {
            if (y[ir][i] < sigmad_floor_[ir][i])
                continue;
            if (y[ir][i] < coag_min_level)
                continue;

            mi = m_[i];
            ai = a_[i];

            for (j = 0; j <= i; ++j)
            {
                if (y[ir][j] < sigmad_floor_[ir][j])
                    continue;
                if (y[ir][j] < coag_min_level)
                    continue;

                mj = m_[j];
                aj = a_[j];
                mij = mi + mj;

                /* coagulation */
                c1 = mass_transfer_index_c1_[i][j];
                c2 = mass_transfer_index_c2_[i][j];
                f = mass_transfer_fraction_[i][j];

                CalcDustRelativeVelocity(ir, i, j, vd_ij, sigma_ij, vel);
                a_sij = M_PI * SQR(ai + aj);
                C = 1.0 / std::sqrt(2.0 * M_PI * (SQR(Hd_[ir][i]) + SQR(Hd_[ir][j])));

                if (vd_ij < 5.0e-2 * sigma_ij)
                {
                    x = SQR(v_coag) / (2.0 * SQR(sigma_ij));
                    kij = C * a_sij * sqrt8_pi * sigma_ij * (1.0 - (1.0 + x) * std::exp(-x));
                    fij = std::max(C * a_sij * sqrt8_pi * sigma_ij - kij, 0.0);

                    if (fij / kij < 1.0e-3)
                    {
                        kij = C * a_sij * sqrt8_pi * sigma_ij;
                        fij = 0.0;
                    }
                }
                else
                {
                    x = 1.0 / (sqrt2 * sigma_ij);
                    z = 1.0 / vd_ij;
                    Av = std::erf((v_coag + vd_ij) * x) - std::erf((v_coag - vd_ij) * x);
                    Bv = (v_coag + vd_ij) * std::exp(-SQR((v_coag - vd_ij) * x)) - (v_coag - vd_ij) * std::exp(-SQR((v_coag + vd_ij) * x));
                    kij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * (2.0 * std::erf(vd_ij * x) - Av) + sqrt1_2pi * sigma_ij * z * (2.0 * vd_ij * std::exp(-SQR(vd_ij * x)) - Bv));
                    fij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * Av + sqrt1_2pi * sigma_ij * z * Bv);
                }

                // fragmentation off -> perfect sticking: high-velocity collisions (fij) also coagulate
                if (!on_fragmentation_) { kij += fij; fij = 0.0; }
                kij = kij / (mi * mj);
                if (i == j)
                    kij = kij * 0.5;
                dndt = kij * y[ir][i] * y[ir][j];

                // Loss (self-proportional): factor out y[i] / y[j] into the rate
                // coefficient instead of writing -mi*dndt / -mj*dndt directly.
                LossRate[ir][i] += mi * kij * y[ir][j];
                LossRate[ir][j] += mj * kij * y[ir][i];
                // Gain stays explicit (evaluated at the old y, same as before).
                Gain[ir][c1] += f * mij * dndt;
                Gain[ir][c2] += (1.0 - f) * mij * dndt;

                /* fragmentation */
                if (!on_fragmentation_)
                    continue;
                if (i == 0 && j == 0)
                    continue;
                if (fij == 0.0)
                    continue;

                mej = mj * vel * inv_vfrag;
                if (mej < m_[0])
                    continue;
                mrm = mij - mej;
                if (mrm < m_[0])
                {
                    mej = mij;
                    mrm = 0.0;
                }

                mf = std::max(frag_frac_ * mej, mb_[1]);

                {
                    auto it = std::upper_bound(mb_.begin(), mb_.begin() + (nbin + 1), mf);
                    int k_found = int(it - mb_.begin()) - 1;
                    if (k_found >= 0 && k_found < nbin)
                    {
                        kf = k_found;
                        mf = mb_[kf + 1];
                    }
                }

                fij = fij / (mi * mj);
                if (i == j)
                    fij *= 0.5;
                dndt = fij * y[ir][i] * y[ir][j];
                LossRate[ir][i] += mi * fij * y[ir][j];
                LossRate[ir][j] += mj * fij * y[ir][i];

                if (mf > mrm)
                {
                    inv_mfrag = mij * dndt * mf_gamma_p2_[kf];
                }
                else
                {
                    {
                        auto it2 = std::upper_bound(mb_.begin() + kf, mb_.begin() + (nbin + 1), mrm);
                        int k2_found = int(it2 - mb_.begin()) - 1;
                        if (k2_found >= kf && k2_found < nbin)
                        {
                            Gain[ir][k2_found] += mrm * dndt;
                        }
                    }

                    inv_mfrag = mej * dndt * mf_gamma_p2_[kf];
                }
                if( std :: isnan(vd_r_[ir][j]) == true || std :: isnan(vd_phi_[ir][j]) == true || std :: isnan(vg_r_[ir]) == true || std :: isnan(vg_phi_[ir]) == true || std :: isnan(fdg_ratio[ir][j]) == true || std :: isnan(lambda1_[ir]) == true  || std :: isnan(mf_gamma_p2_[kf]) == true || std :: isnan(vd_r_[ir][i]) == true || std :: isnan(dndt) == true ||  std :: isnan(vd_ij) == true)
                {
                    std :: cerr <<std :: scientific<<ir<<" "<< i <<" "<<j <<" "<<sigmad_[ir][j]<<" "<<a_[i]<<" "<< a_[j] <<" "<<vd_r_[ir][j] <<" "<<vd_phi_[ir][j]<<" "<<sigmag_[ir]<<" "<<vg_r_[ir]<<" "<<vg_phi_[ir]<<" "<<dNdr_torque[ir]<<" "<<lambda1_[ir]<<" "<<lambda2_[ir]<<" "<<vg_wind_[ir]<<" "<<St_[ir][j]<<" "<<fdg_ratio[ir][j]<<" "<<A_[ir]<<" "<<vd_ij<<" "<<dndt<<" "<<sigma_ij
                        <<" Hd_i="<<Hd_[ir][i]<<" Hd_j="<<Hd_[ir][j]<<" C="<<C<<" kij="<<kij<<" fij="<<fij<<" a_sij="<<a_sij<<" St_i="<<St_[ir][i]<< std :: endl;
                    exit(1);
                }

                for (k = 0; k <= kf; ++k)
                {
                    Gain[ir][k] += frag_mass_dist_[k] * inv_mfrag;
                }
            }
        }
    return;
}

// Fully implicit coagulation/fragmentation (coagulation_speedup_literature.md
// Sec.3, Brauer et al. 2008 Appendix C): same (i,j) pairwise physics as
// CoagulateChangeRateFragmentation_Implicit_Cell() -- same kij/fij, same
// Gain/LossRate accumulation -- but *also* builds the full Jacobian
// J_kl = dF_k/dsigma_l of the net rate F_k(sigma) = Gain_k(sigma) -
// sigma_k*LossRate_k(sigma), so the whole cell can be solved as one linear
// system instead of dividing each bin by its own (self only) loss term.
//
// Differentiation rule used throughout below: every accumulation
// "LossRate[k] += COEF*y[l]" (COEF is y-independent) contributes to
// F_k = ... - y_k*LossRate_k, so its two partial-derivative pieces are added
// alongside it:
//   J[k][k] -= COEF*y[l];   // the "-delta_kl*LossRate_k" piece (always the
//                           // diagonal, regardless of which l the term uses)
//   J[k][l] -= y[k]*COEF;   // the "-y_k * dLossRate_k/dy_l" piece
// and every accumulation "Gain[k] += COEF*y[a]*y[b]" contributes
//   J[k][a] += COEF*y[b];
//   J[k][b] += COEF*y[a];
// (when a==b or k==l these collapse onto the same matrix entry and correctly
// sum to the product-rule result for a squared term, e.g. d(C*y_a^2)/dy_a =
// 2*C*y_a instead of C*y_a).
void DustDisk::CoagulateChangeRateFragmentation_FullImplicit_Cell(int ir, const DoubleMatrix &y,
    DoubleVector &Gain, DoubleVector &LossRate, DoubleMatrix &J)
{
    int i, j, k, nbin, nbinm, c1, c2;
    double f, kij, dndt, mi, mj, ai, aj, mij, a_sij, vel_coag;
    double vel, vd_ij, sigma_ij, Av, Bv, Af, Bf, coe1, coe2, C, x, z;
    int kf;
    double mej, mrm, mf, inv_mfrag, inv_mfrag_coef, inv_vfrag, fij, gamma_p2, m0_gamma_p2;
    double v_coag, v_frag;

    nbin = nbin_;
    nbinm = nbin - 1;

    gamma_p2 = frag_slope_ + 2.0;
    m0_gamma_p2 = m0_gamma_p2_;

    const double REL_THRESH_COAG = 1.0e-20;
    // Absolute floor (g/cm^2): bins below this are skipped regardless of the
    // cell's total, so truly negligible populations (far below any physical
    // significance) don't enter the O(nbin^2) kernel loop / Jacobian.
    const double ABS_THRESH_COAG = 1.0e-20;

    std::fill(Gain.begin(), Gain.end(), 0.0);
    std::fill(LossRate.begin(), LossRate.end(), 0.0);
    for (i = 0; i < nbin; ++i)
    {
        std::fill(J[i].begin(), J[i].end(), 0.0);
    }
    if(T_[ir] > 1.5e3)
    {
        return;
    }

    // if (T_[ir] < 1.6e2)
    // {
    //     inv_vfrag = 1.0/v_col_crit_;
    //     v_coag = 0.2 * v_col_crit_;
    //     v_frag = 0.2 * v_col_crit_;
    // }
    // else
    // {
    //     inv_vfrag =  1.0/v_col_sili_;
    //     v_coag = 0.2 * v_col_sili_;
    //     v_frag = 0.2 * v_col_sili_;
    // }
    double w   = 0.5*(1.0 - std::tanh((T_[ir] - T_snow_)/(0.1*T_snow_)));
    double vc  = w*v_col_crit_ + (1.0-w)*v_col_sili_;
    inv_vfrag  = 1.0/vc;
    v_coag     = 0.2*vc;
    v_frag     = 0.2*vc;

    double coag_min_level = std::max(REL_THRESH_COAG * sigmad_tot_[ir], ABS_THRESH_COAG);

    for (i = 0; i < nbin; ++i)
    {
        if (y[ir][i] < sigmad_floor_[ir][i])
            continue;
        if (y[ir][i] < coag_min_level)
            continue;

        mi = m_[i];
        ai = a_[i];

        for (j = 0; j <= i; ++j)
        {
            if (y[ir][j] < sigmad_floor_[ir][j])
                continue;
            if (y[ir][j] < coag_min_level)
                continue;

            mj = m_[j];
            aj = a_[j];
            mij = mi + mj;

            /* coagulation */
            c1 = mass_transfer_index_c1_[i][j];
            c2 = mass_transfer_index_c2_[i][j];
            f = mass_transfer_fraction_[i][j];

            CalcDustRelativeVelocity(ir, i, j, vd_ij, sigma_ij, vel);
            a_sij = M_PI * SQR(ai + aj);
            C = 1.0 / std::sqrt(2.0 * M_PI * (SQR(Hd_[ir][i]) + SQR(Hd_[ir][j])));

            if (vd_ij < 5.0e-2 * sigma_ij)
            {
                x = SQR(v_coag) / (2.0 * SQR(sigma_ij));
                kij = C * a_sij * sqrt8_pi * sigma_ij * (1.0 - (1.0 + x) * std::exp(-x));
                fij = std::max(C * a_sij * sqrt8_pi * sigma_ij - kij, 0.0);

                if (fij / kij < 1.0e-3)
                {
                    kij = C * a_sij * sqrt8_pi * sigma_ij;
                    fij = 0.0;
                }
            }
            else
            {
                x = 1.0 / (sqrt2 * sigma_ij);
                z = 1.0 / vd_ij;
                Av = std::erf((v_coag + vd_ij) * x) - std::erf((v_coag - vd_ij) * x);
                Bv = (v_coag + vd_ij) * std::exp(-SQR((v_coag - vd_ij) * x)) - (v_coag - vd_ij) * std::exp(-SQR((v_coag + vd_ij) * x));
                kij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * (2.0 * std::erf(vd_ij * x) - Av) + sqrt1_2pi * sigma_ij * z * (2.0 * vd_ij * std::exp(-SQR(vd_ij * x)) - Bv));
                fij = C * a_sij * ((SQR(vd_ij) + SQR(sigma_ij)) * 0.5 * z * Av + sqrt1_2pi * sigma_ij * z * Bv);
            }

            // fragmentation off -> perfect sticking: high-velocity collisions (fij) also coagulate
            if (!on_fragmentation_) { kij += fij; fij = 0.0; }
            kij = kij / (mi * mj);
            if (i == j)
                kij = kij * 0.5;
            dndt = kij * y[ir][i] * y[ir][j];

            // Loss/Gain (see differentiation rule in the function comment above).
            // When the product lands back in bin i (c1 == i), the i-loss and
            // the c1-gain are merged into one net coefficient for dsigma_i/dt
            // = net*y_i*y_j, net = kij*(f*m_j - (1-f)*m_i) (algebraically equal
            // to kij*(f*m_ij - m_i), but avoids cancelling m_i against m_ij).
            // The sign decides where it goes, keeping LossRate/Gain >= 0.
            double coef_j = mj * kij; // LossRate[j] += coef_j * y[i]
            LossRate[j] += coef_j * y[ir][i];
            J[j][j] -= coef_j * y[ir][i];
            J[j][i] -= y[ir][j] * coef_j;

            if (c1 == i)
            {
                double net_i = kij * (f * mj - (1.0 - f) * mi);
                if (net_i < 0.0)
                {
                    double loss_i = -net_i; // LossRate[i] += loss_i * y[j]
                    LossRate[i] += loss_i * y[ir][j];
                    J[i][i] -= loss_i * y[ir][j];
                    J[i][j] -= y[ir][i] * loss_i;
                }
                else
                {
                    Gain[i] += net_i * y[ir][i] * y[ir][j];
                    J[i][i] += net_i * y[ir][j];
                    J[i][j] += net_i * y[ir][i];
                }
            }
            else
            {
                double coef_i = mi * kij; // LossRate[i] += coef_i * y[j]
                LossRate[i] += coef_i * y[ir][j];
                J[i][i] -= coef_i * y[ir][j];
                J[i][j] -= y[ir][i] * coef_i;

                double coef_c1 = f * mij * kij;
                Gain[c1] += coef_c1 * y[ir][i] * y[ir][j];
                J[c1][i] += coef_c1 * y[ir][j];
                J[c1][j] += coef_c1 * y[ir][i];
            }
            {
                double coef_c2 = (1.0 - f) * mij * kij;
                Gain[c2] += coef_c2 * y[ir][i] * y[ir][j];
                J[c2][i] += coef_c2 * y[ir][j];
                J[c2][j] += coef_c2 * y[ir][i];
            }

            /* fragmentation */
            if (!on_fragmentation_)
                continue;
            if (i == 0 && j == 0)
                continue;
            if (fij == 0.0)
                continue;

            mej = mj * vel * inv_vfrag;
            if (mej < m_[0])
                continue;
            mrm = mij - mej;
            if (mrm < m_[0])
            {
                mej = mij;
                mrm = 0.0;
            }

            mf = std::max(frag_frac_ * mej, mb_[1]);

            {
                auto it = std::upper_bound(mb_.begin(), mb_.begin() + (nbin + 1), mf);
                int k_found = int(it - mb_.begin()) - 1;
                if (k_found >= 0 && k_found < nbin)
                {
                    kf = k_found;
                    mf = mb_[kf + 1];
                }
            }

            fij = fij / (mi * mj);
            if (i == j)
                fij *= 0.5;
            dndt = fij * y[ir][i] * y[ir][j];

            int k2_found = -1;
            if (!(mf > mrm))
            {
                auto it2 = std::upper_bound(mb_.begin() + kf, mb_.begin() + (nbin + 1), mrm);
                int k2_cand = int(it2 - mb_.begin()) - 1;
                if (k2_cand >= kf && k2_cand < nbin)
                    k2_found = k2_cand;
            }
            bool remnant_self = (k2_found == i);

            if (remnant_self)
            {
                // Remnant returns to bin i: net dsigma_i/dt = fij*(m_rm - m_i) = fij*(m_j - m_ej)
                // (also holds when m_rm was clipped to 0, since then m_ej = m_ij).
                double net_i = fij * (mj - mej);
                if (net_i < 0.0)
                {
                    double loss_i = -net_i;
                    LossRate[i] += loss_i * y[ir][j];
                    J[i][i] -= loss_i * y[ir][j];
                    J[i][j] -= y[ir][i] * loss_i;
                }
                else
                {
                    Gain[i] += net_i * y[ir][i] * y[ir][j];
                    J[i][i] += net_i * y[ir][j];
                    J[i][j] += net_i * y[ir][i];
                }
            }
            else
            {
                double coef_i = mi * fij;
                LossRate[i] += coef_i * y[ir][j];
                J[i][i] -= coef_i * y[ir][j];
                J[i][j] -= y[ir][i] * coef_i;
            }
            {
                double coef_j = mj * fij;
                LossRate[j] += coef_j * y[ir][i];
                J[j][j] -= coef_j * y[ir][i];
                J[j][i] -= y[ir][j] * coef_j;
            }

            if (mf > mrm)
            {
                inv_mfrag_coef = mij * mf_gamma_p2_[kf] * fij;
                inv_mfrag = inv_mfrag_coef * y[ir][i] * y[ir][j];
            }
            else
            {
                if (k2_found >= 0 && !remnant_self)
                {
                    double coef_k2 = mrm * fij;
                    Gain[k2_found] += coef_k2 * y[ir][i] * y[ir][j];
                    J[k2_found][i] += coef_k2 * y[ir][j];
                    J[k2_found][j] += coef_k2 * y[ir][i];
                }

                inv_mfrag_coef = mej * mf_gamma_p2_[kf] * fij;
                inv_mfrag = inv_mfrag_coef * y[ir][i] * y[ir][j];
            }
            if( std :: isnan(vd_r_[ir][j]) == true || std :: isnan(vd_phi_[ir][j]) == true || std :: isnan(vg_r_[ir]) == true || std :: isnan(vg_phi_[ir]) == true || std :: isnan(fdg_ratio[ir][j]) == true || std :: isnan(lambda1_[ir]) == true  || std :: isnan(mf_gamma_p2_[kf]) == true || std :: isnan(vd_r_[ir][i]) == true || std :: isnan(dndt) == true ||  std :: isnan(vd_ij) == true)
            {
                std :: cerr <<std :: scientific<<ir<<" "<< i <<" "<<j <<" "<<sigmad_[ir][j]<<" "<<a_[i]<<" "<< a_[j] <<" "<<vd_r_[ir][j] <<" "<<vd_phi_[ir][j]<<" "<<sigmag_[ir]<<" "<<vg_r_[ir]<<" "<<vg_phi_[ir]<<" "<<dNdr_torque[ir]<<" "<<lambda1_[ir]<<" "<<lambda2_[ir]<<" "<<vg_wind_[ir]<<" "<<St_[ir][j]<<" "<<fdg_ratio[ir][j]<<" "<<A_[ir]<<" "<<vd_ij<<" "<<dndt<<" "<<sigma_ij
                    <<" Hd_i="<<Hd_[ir][i]<<" Hd_j="<<Hd_[ir][j]<<" C="<<C<<" kij="<<kij<<" fij="<<fij<<" a_sij="<<a_sij<<" St_i="<<St_[ir][i]
                    <<" y_i="<<y[ir][i]<<" y_j="<<y[ir][j]<<" mej="<<mej<<" mrm="<<mrm<<" mf="<<mf<<" kf="<<kf<<" vel="<<vel<<" mij="<<mij<<" mi="<<mi<<" mj="<<mj<< std :: endl;
                exit(1);
            }

            for (k = 0; k <= kf; ++k)
            {
                double coef_k = frag_mass_dist_[k] * inv_mfrag_coef;
                Gain[k] += coef_k * y[ir][i] * y[ir][j];
                J[k][i] += coef_k * y[ir][j];
                J[k][j] += coef_k * y[ir][i];
            }
        }
    }
    return;
}

// Solves A*x = b in place via Gaussian elimination with partial pivoting
// (A and b are overwritten; the solution is left in b). Returns false if a
// pivot is degenerate (singular/near-singular matrix) -- the caller must
// then fall back to leaving sigmad unchanged for that cell rather than
// trusting garbage output. O(n^3), fine for the dense n=nbin_ (~300) system
// built by CoagulateChangeRateFragmentation_FullImplicit_Cell().
static bool SolveLinearSystemGaussianElim(DoubleMatrix &A, DoubleVector &b, int n)
{
    for (int col = 0; col < n; ++col)
    {
        int pivot_row = col;
        double pivot_val = std::fabs(A[col][col]);
        for (int r = col + 1; r < n; ++r)
        {
            if (std::fabs(A[r][col]) > pivot_val)
            {
                pivot_val = std::fabs(A[r][col]);
                pivot_row = r;
            }
        }
        if (pivot_val < 1.0e-300)
        {
            return false;
        }
        if (pivot_row != col)
        {
            std::swap(A[col], A[pivot_row]);
            std::swap(b[col], b[pivot_row]);
        }

        double inv_pivot = 1.0 / A[col][col];
        for (int r = col + 1; r < n; ++r)
        {
            double factor = A[r][col] * inv_pivot;
            if (factor == 0.0) continue;
            for (int c = col; c < n; ++c)
            {
                A[r][c] -= factor * A[col][c];
            }
            b[r] -= factor * b[col];
        }
    }

    for (int row = n - 1; row >= 0; --row)
    {
        double sum = b[row];
        for (int c = row + 1; c < n; ++c)
        {
            sum -= A[row][c] * b[c];
        }
        b[row] = sum / A[row][row];
    }
    return true;
}

// Fully implicit counterpart of ApplyCoagulationAndCorrectMass(): solves the
// nonlinear backward-Euler equation G(sigma) = sigma - sigma_prev - dt_sub*F(sigma) = 0
// (F_k(sigma) = Gain_k(sigma) - sigma_k*LossRate_k(sigma)) via a fixed number
// of Newton iterations per sub-step (Brauer et al. 2008 Appendix C used a
// single such step; coagulation_speedup_literature.md Sec.5 postmortem: for
// the most fragmentation-active cells, a *single* linearization noticeably
// over-populated the smallest bins even though total mass came out right --
// one linear solve isn't enough to resolve the nonlinearity when dt is this
// large, so this re-linearizes and re-solves a few times per sub-step,
// converging toward the true nonlinear backward-Euler solution):
//   G'(sigma) = I - dt_sub*J(sigma)
//   delta = -G'(sigma^(m))^-1 * G(sigma^(m)),   sigma^(m+1) = sigma^(m) + delta
//
// Sec.7 postmortem (frg10_beta3 out5400: nearly all mass at one cell
// collapsed into the single largest bin for one snapshot, then recovered):
// a Newton solve can occasionally fail to converge for a given dt (producing
// wildly negative bins that get clamped, with the subsequent mass-correction
// then silently papering over the resulting wrong *shape* by rescaling the
// total back to sigma_tot_entry). Rescaling a broken result to *look*
// conserved is not an acceptable fix, but neither is silently skipping the
// physics for that cell/step -- the correct fix is the same principle used
// throughout this codebase for other stiff terms: shrink the step and
// sub-cycle *internally* (this cell only) until the nonlinear solve actually
// converges, exactly like Sec.9/10's n_sub did for the semi-implicit scheme,
// but here the trigger is Newton convergence itself rather than a pre-estimated
// timescale. A sub-step is judged converged only if, after NEWTON_ITERS
// iterations, no bin needed a physically-significant negative-value clamp
// (small clamps from floating-point-level roundoff are fine and expected).
// If not converged, the *entire* macro dt is retried from the original
// sigma_o with n_sub doubled (dt_sub halved) -- not just the failed
// sub-step, since later sub-steps built on a wrong intermediate state would
// themselves be suspect.
//
// Gain_scratch/LossRate_scratch/J_scratch/A_scratch/rhs_scratch/sigma_o_scratch
// are caller-owned buffers (sized nbin / nbin*nbin), reused across the
// serial per-cell loop to avoid reallocating a 300x300 matrix every
// cell/iteration.
void DustDisk :: ApplyCoagulationAndCorrectMass_FullImplicit(std::vector<std::vector<double>>& sigmad,
    std::vector<double>& sigmad_tot,
    DoubleVector &Gain_scratch, DoubleVector &LossRate_scratch, DoubleMatrix &J_scratch,
    DoubleMatrix &A_scratch, DoubleVector &rhs_scratch, DoubleVector &sigma_o_scratch,
    double dt, int ir, int nbin)
{
    const int NEWTON_ITERS = 3;
    const int N_SUB_MAX = 1024;
    // A sub-step is rejected as non-converged if any bin would need clamping
    // by more than this fraction of the cell's total mass -- distinguishes a
    // genuine Newton blow-up (Sec.7: one bin absorbing ~100% of the mass)
    // from ordinary floating-point-level negative roundoff.
    const double NEGATIVE_CLAMP_TOL = 1.0e-8;
    // DustPy-style step-acceptance (Stammler & Birnstiel 2022): the checks
    // below only catch a solution going negative/non-finite, which is a
    // necessary but not sufficient condition for the fixed 3-iteration
    // Newton solve having actually converged -- a stiff cell (e.g. bins near
    // the top of the mass grid, where the collision kernel changes sharply)
    // can land on a "shape" that is algebraically valid (non-negative,
    // finite, and conserves the cell total via the rescaling below) but
    // physically wrong, e.g. most of the cell's mass jumping into a single
    // near-empty bin in one sub-step. Reject a sub-step where any single
    // bin's implicit update moves by more than this fraction of the cell's
    // total dust mass -- a real collision cascade redistributes mass among
    // *adjacent* bins each step, so a jump of this size in one bin is a
    // Newton-convergence failure, not fast physics.
    const double REL_CHANGE_TOL = 0.3;
    // Newton residual tolerance: the loop below used to run exactly
    // NEWTON_ITERS passes and trust the result as long as it was
    // non-negative/finite/within REL_CHANGE_TOL, without ever checking
    // whether the linearization had actually converged. A step that still
    // wants to move by more than this fraction of the cell's total mass
    // after NEWTON_ITERS iterations has not converged -- accepting it
    // anyway (as before) meant the "solution" could be an artifact of
    // stopping mid-iteration on a stiff cell, not a genuine fixed point of
    // F(sigma)=0. Checked every iteration so a cell that converges in fewer
    // than NEWTON_ITERS passes can stop early too.
    const double NEWTON_CONV_TOL = 1.0e-3;

    double sigma_tot_entry = sigmad_tot[ir];
    DoubleVector sigma_start(nbin);
    for (int k = 0; k < nbin; ++k)
    {
        sigma_start[k] = sigmad[ir][k];
    }

    // No a priori dt_coag_cell estimate here (tried and reverted -- it
    // slowed down ordinary, non-stiff cells across the board, apparently by
    // being systematically over-conservative). n_sub always starts at 1;
    // the a posteriori REL_CHANGE_TOL/positivity check below is what
    // actually prevents a bad (wrong-shape) solution from being accepted --
    // it re-tries with a smaller dt_sub regardless of where n_sub started,
    // so starting at 1 costs at most a few extra doubling rounds for a
    // genuinely stiff cell, not a correctness regression.
    int n_sub = 1;
    bool converged = false;
    while (!converged)
    {
        double dt_sub = dt / double(n_sub);
        for (int k = 0; k < nbin; ++k)
        {
            sigmad[ir][k] = sigma_start[k];
        }

        converged = true;
        for (int s = 0; s < n_sub && converged; ++s)
        {
            for (int k = 0; k < nbin; ++k)
            {
                sigma_o_scratch[k] = sigmad[ir][k];
            }

            double newton_resid = std::numeric_limits<double>::max();
            for (int newton_iter = 0; newton_iter < NEWTON_ITERS; ++newton_iter)
            {
                // F/J evaluated at the CURRENT iterate sigma^(m) (sigmad[ir],
                // updated in place below each pass), not the stale
                // sigma_o_scratch -- this is what makes it a genuine Newton
                // iteration instead of repeating the same single-linearization
                // solve.
                CoagulateChangeRateFragmentation_FullImplicit_Cell(ir, sigmad, Gain_scratch, LossRate_scratch, J_scratch);

                // A = I - dt_sub*J(sigma^(m)) ; rhs = -G(sigma^(m))
                //   = sigma_o - sigma^(m) + dt_sub*F(sigma^(m))
                for (int k = 0; k < nbin; ++k)
                {
                    for (int l = 0; l < nbin; ++l)
                    {
                        A_scratch[k][l] = -dt_sub * J_scratch[k][l];
                    }
                    A_scratch[k][k] += 1.0;
                    double F_k = Gain_scratch[k] - sigmad[ir][k] * LossRate_scratch[k];
                    rhs_scratch[k] = sigma_o_scratch[k] - sigmad[ir][k] + dt_sub * F_k;
                }

                bool ok = SolveLinearSystemGaussianElim(A_scratch, rhs_scratch, nbin);
                if (!ok)
                {
                    // Degenerate system: treat exactly like a non-converged
                    // sub-step (shrink and retry) rather than accepting or
                    // silently discarding whatever partial result exists.
                    converged = false;
                    break;
                }

                // Newton step size (before damping), relative to the cell's
                // total mass -- the actual measure of whether the
                // linearization has converged, not just whether the result
                // looks non-negative. Small enough -> stop iterating now
                // rather than spending NEWTON_ITERS-1-newton_iter more
                // passes recomputing a step that would barely move sigmad.
                newton_resid = 0.0;
                double resid_denom = std::max(std::abs(sigma_tot_entry), 1.0e-20);
                for (int k = 0; k < nbin; ++k)
                {
                    double r = std::fabs(rhs_scratch[k]) / resid_denom;
                    if (r > newton_resid) newton_resid = r;
                }

                // Positivity-preserving step-length damping: the raw Newton
                // step can wildly overshoot for this stiff, densely-coupled
                // (nbin x nbin) system -- observed driving bins to negative
                // values many orders of magnitude past the cell's own total
                // mass (negative-value failures in the trace). Shrink the step by
                // the largest lambda in (0,1] that keeps every bin within the
                // same tolerance the convergence check below already accepts
                // (-NEGATIVE_CLAMP_TOL*|sigma_tot_entry|), so a step that would
                // have passed anyway is never damped, and one that would have
                // failed is pulled back to just barely satisfy it.
                double lambda = 1.0;
                double neg_tol_abs = NEGATIVE_CLAMP_TOL * std::abs(sigma_tot_entry);
                for (int k = 0; k < nbin; ++k)
                {
                    if (rhs_scratch[k] >= 0.0) continue;
                    double next = sigmad[ir][k] + rhs_scratch[k];
                    if (next < -neg_tol_abs)
                    {
                        double lambda_k = (-neg_tol_abs - sigmad[ir][k]) / rhs_scratch[k];
                        if (lambda_k < 0.0) lambda_k = 0.0;
                        if (lambda_k < lambda) lambda = lambda_k;
                    }
                }
                for (int k = 0; k < nbin; ++k)
                {
                    sigmad[ir][k] += lambda * rhs_scratch[k];
                }

                if (newton_resid <= NEWTON_CONV_TOL) break; // converged in fewer than NEWTON_ITERS passes
            }
            if (!converged) break;

            if (newton_resid > NEWTON_CONV_TOL)
            {
                // Ran out of Newton iterations without the step size
                // actually shrinking below tolerance -- the result is
                // whatever the last linearization happened to land on, not
                // a converged fixed point. Same remedy as the other
                // convergence failures below: shrink dt_sub and retry.
                converged = false;
                break;
            }

            // Convergence check: a bin going meaningfully negative here is
            // the signature of the Newton solve blowing up for this dt_sub
            // (Sec.7), not ordinary roundoff. A NaN/Inf entry is caught
            // separately: SolveLinearSystemGaussianElim can return ok=true
            // (a pivot was found) yet still produce a non-finite value from
            // overflow on a near-singular system, and IEEE754 makes
            // "NaN < x" false, so the negative check above would silently
            // let it through as "converged" and propagate NaN into the next
            // Newton iteration/cell. Treat non-finite exactly like the
            // negative-overshoot case: not converged, shrink dt_sub, retry.
            // sigma_tot_entry (the cell's total, summed upstream from
            // per-bin values that can themselves carry roundoff-level
            // negative excursions from advection) can itself be a small
            // negative number. Without abs(), a negative sigma_tot_entry
            // flips -NEGATIVE_CLAMP_TOL*sigma_tot_entry positive, turning
            // the "tolerate tiny negative" tolerance into "reject almost
            // everything" -- forcing sub-cycling to exhaust N_SUB_MAX and
            // bail out via the restore-original-state path below, which
            // skips the clamp-to-zero code entirely and lets the negative
            // (or, once contaminated, NaN) value survive into the next
            // step. abs() keeps the tolerance a small window around zero
            // regardless of the total's sign.
            double rel_change_denom = std::max(std::abs(sigma_tot_entry), 1.0e-20);
            for (int k = 0; k < nbin; ++k)
            {
                if (!std::isfinite(sigmad[ir][k]))
                {
                    converged = false;
                    break;
                }
                if (sigmad[ir][k] < -NEGATIVE_CLAMP_TOL * std::abs(sigma_tot_entry))
                {
                    converged = false;
                    break;
                }
                // DustPy-style relative-change step-acceptance (see REL_CHANGE_TOL
                // comment above): sigma_o_scratch[k] is this bin's value at the
                // start of the current sub-step.
                if (std::fabs(sigmad[ir][k] - sigma_o_scratch[k]) > REL_CHANGE_TOL * rel_change_denom)
                {
                    converged = false;
                    break;
                }
            }
        }

        if (!converged)
        {
            if (n_sub >= N_SUB_MAX)
            {
                // Exhausted the sub-cycling budget without the full-implicit
                // Newton solve converging (should now be rare with the
                // positivity-preserving step damping above). Rather than
                // freezing coagulation for this cell (the old "restore
                // sigma_start and return" fallback -- silently skips growth/
                // fragmentation here every time this triggers), fall back to
                // the semi-implicit update for just this cell:
                //   sigma_k_new = (sigma_k_old + dt_sub*Gain_k) / (1 + dt_sub*LossRate_k)
                // Gain_k >= 0 and 1+dt_sub*LossRate_k > 0 always, so this is
                // unconditionally non-negative *in exact arithmetic* for any
                // dt_sub -- no Newton, cannot fail to produce a non-negative
                // shape. But it is NOT part of Brauer, Dullemond & Henning
                // (2008) Appendix C (their implicit scheme is a single
                // linearization around f_o with no iteration or sub-cycling
                // at all) -- this fallback, and the sub-cycling below, are
                // this codebase's own addition, so "non-negative" is the only
                // guarantee inherited from the formula's structure, and even
                // that guarantee is only in real arithmetic: applied in one
                // shot with the full (un-sub-cycled) dt, sigma_start/Gain
                // that have already drifted to extreme (but finite)
                // magnitudes over repeated triggerings of this same fallback
                // can make dt*Gain overflow to a literal `inf` in double
                // precision (observed directly: sigma_start~1e+171, Gain=inf).
                // Sub-cycle at dt/N_SUB_MAX
                // (the finest resolution Newton itself was already trying and
                // failing to reach) instead of one shot with the full dt, and
                // re-evaluate Gain/LossRate from the current state each
                // sub-step (not frozen at sigma_start, per the "coefficients
                // frozen across dt" issue in dust_coagulation_implicit_note.md
                // Sec.8) -- this keeps each step's dt_sub*Gain far smaller and
                // makes overflow far less likely. If it still overflows
                // (isfinite check per sub-step), stop and freeze this cell at
                // sigma_start for this macro step, matching the originally
                // documented safe fallback (coagulation_speedup_literature.md
                // Sec.7: "leave sigma at sigma_o") -- never write inf/NaN into
                // sigmad_.
                for (int k = 0; k < nbin; ++k)
                {
                    sigmad[ir][k] = sigma_start[k];
                }
                double dt_fallback_sub = dt / double(N_SUB_MAX);
                bool fallback_ok = true;
                for (int s = 0; s < N_SUB_MAX && fallback_ok; ++s)
                {
                    CoagulateChangeRateFragmentation_FullImplicit_Cell(ir, sigmad, Gain_scratch, LossRate_scratch, J_scratch);
                    for (int k = 0; k < nbin; ++k)
                    {
                        double val = (sigmad[ir][k] + dt_fallback_sub * Gain_scratch[k]) / (1.0 + dt_fallback_sub * LossRate_scratch[k]);
                        if (!std::isfinite(val))
                        {
                            fallback_ok = false;
                            break;
                        }
                        if (val < 0.0) val = 0.0; // defensive only; formula is non-negative by construction
                        sigmad[ir][k] = val;
                    }
                }
                if (!fallback_ok)
                {
                    for (int k = 0; k < nbin; ++k)
                    {
                        sigmad[ir][k] = sigma_start[k];
                    }
                }
                double temp_sigma_fallback = 0.0;
                for (int k = 0; k < nbin; ++k)
                {
                    temp_sigma_fallback += sigmad[ir][k];
                }
                sigmad_tot[ir] = temp_sigma_fallback;

                // Mass-conservation rescale (same as the converged path's
                // Sec.11 correction, dust_coagulation_implicit_note.md): the
                // semi-implicit formula above (Gain explicit-ish, Loss alone
                // implicit) is only approximately mass-conserving per
                // sub-step -- fine as a one-off, but this fallback cell is,
                // by construction, one where Newton could not converge, so
                // it tends to hit this same fallback again on later macro
                // steps too. Left uncorrected, the per-substep mass error
                // compounds every time this triggers (observed directly:
                // beta4_new/beta3_new dust-to-gas ratio in the inner disk
                // grew from the fdg=0.01 initial value to 0.1-0.3, spiking
                // over 1 at some cells, by t=160-200yr). Guard the same way
                // as Sec.11: skip if the fallback's own total is zero/tiny
                // (nothing to rescale) or non-finite (never divide toward
                // the kind of overflow Tier 1 just eliminated).
                if (sigma_tot_entry > 0.0 && std::isfinite(sigmad_tot[ir]) && sigmad_tot[ir] > 1.0e-250)
                {
                    double err = (sigmad_tot[ir] - sigma_tot_entry) / sigma_tot_entry;
                    if (std::abs(err) > 1.0e-6)
                    {
                        double delta_sigma = sigma_tot_entry / sigmad_tot[ir];
                        double corrected = 0.0;
                        for (int k = 0; k < nbin; ++k)
                        {
                            sigmad[ir][k] *= delta_sigma;
                            corrected += sigmad[ir][k];
                        }
                        sigmad_tot[ir] = corrected;
                    }
                }
                return;
            }
            n_sub *= 2;
        }
    }

    // Converged: clamp only the floating-point-level negative residue that
    // NEGATIVE_CLAMP_TOL deliberately tolerates (this is not "papering over"
    // a wrong shape -- a real shape error would have been caught and
    // retried with a smaller dt_sub above).
    double temp_sigma_coag = 0.0;
    for (int k = 0; k < nbin; ++k)
    {
        if (sigmad[ir][k] < 0.0)
        {
            sigmad[ir][k] = 0.0;
        }
        temp_sigma_coag += sigmad[ir][k];
    }
    sigmad_tot[ir] = temp_sigma_coag;

    // Residual (roundoff-level, not shape-fixing) mass correction over the
    // whole step -- see dust_coagulation_implicit_note.md Sec.11.
    if (sigma_tot_entry > 0.0)
    {
        double err = (sigmad_tot[ir] - sigma_tot_entry) / sigma_tot_entry;
        if (std::abs(err) > 1.0e-6 && sigmad_tot[ir] > 0.0)
        {
            double delta_sigma = sigma_tot_entry / sigmad_tot[ir];
            double corrected = 0.0;
            for (int k = 0; k < nbin; ++k)
            {
                sigmad[ir][k] *= delta_sigma;
                corrected += sigmad[ir][k];
            }
            sigmad_tot[ir] = corrected;
        }
    }
}

// Semi-implicit counterpart of ApplyCoagulationAndCorrectMass(): solves
// y_i^{n+1} = (y_i^n + dt*Gain_i) / (1 + dt*LossRate_i) instead of the
// explicit y_i^{n+1} = y_i^n + dNdt*dt. The Loss side is unconditionally
// non-negative regardless of dt (dividing by 1+dt*LossRate can never blow
// up) -- that is the whole point of skipping/relaxing EstimateDtCoagulation's
// dt_coag, and must NOT be undone by adding a Loss-based dt limit back in.
//
// Gain, however, is fully explicit (evaluated at the old sigmad, never
// divided by anything), so it has no such protection. A bin that starts
// near-empty and is suddenly fed by fragmentation debris can jump by an
// arbitrarily large amount if dt is large -- exactly the "dt_coag disabled"
// case this note warns about (dust_coagulation_implicit_note.md Sec.6). So
// only Gain drives sub-cycling here: dt_gain[j] = sigmad_tot[i]/Gain[i][j]
// (using the *cell's total* dust mass, not the bin's own sigmad, so a
// near-zero starting bin doesn't spuriously force an ever-smaller dt_gain),
// and the macro dt is internally split into enough equal sub-steps to keep
// each sub-step's fractional gain within GAIN_TOL. The macro dt used by the
// rest of the simulation (advection, magnetic diffusion, ...) is untouched.
//
// A second, separate residual (dust_coagulation_implicit_note.md Sec.8):
// Gain[i]/LossRate[i] are themselves nonlinear functions of the *whole*
// current size distribution (sums of k_ij*sigma_i*sigma_j over all bin
// pairs), computed once from sigmad^n by the caller. Reusing that one
// snapshot for every sub-step below only refines the integration of
// dsigma/dt with Gain/LossRate held fixed -- it does nothing about those
// coefficients going stale as sigmad evolves within the (now much larger,
// since dt_coag no longer bounds the global dt) macro dt. So each sub-step
// after the first re-derives Gain[i]/LossRate[i] from the current sigmad via
// CoagulateChangeRateFragmentation_Implicit_Cell(), re-linearizing around the
// latest state instead of the stale sigmad^n.
//
// A third residual, and empirically the dominant one (dust_coagulation_implicit_note.md
// Sec.9): the sigmad_tot[i]-based criterion below only bounds each bin's gain
// to GAIN_TOL of the *cell's total* mass -- it says nothing about a bin's
// gain relative to *its own* mass. A bin that starts near-empty (e.g. freshly
// exposed by fragmentation) can therefore be force-filled, in a single
// sub-step, by an amount that is tiny next to sigmad_tot[i] but enormous
// (many orders of magnitude) next to its own prior value, with n_sub never
// reacting to that at all. So a second, per-bin criterion is added: no bin
// may gain more than GAIN_TOL of max(its own current value, the mass of one
// representative grain in that bin, sigmad_floor_[i][j]) per sub-step. Using
// sigmad_floor_ (~m_[j]/area, already a real physical scale, not an ad-hoc
// epsilon) as the floor -- rather than the bin's raw value, which can be
// exactly 0 -- avoids n_sub blowing up while still forcing fine resolution
// for bins that are genuinely being built up from ~nothing.
void DustDisk :: ApplyCoagulationAndCorrectMass_Implicit(std::vector<std::vector<double>>& sigmad,
    std::vector<double>& sigmad_tot,
    std::vector<std::vector<double>>& Gain,
    std::vector<std::vector<double>>& LossRate,
    double dt, int i, int nbin)
{
    const double GAIN_TOL = 0.2;      // max fractional gain (of sigmad_tot, or of a bin's own mass) allowed per sub-step
    const int N_SUB_MAX = 1000;     // safety net, not expected to bind in practice

    double max_gain = 0.0;
    double min_dt_gain_own = std::numeric_limits<double>::max();
    double min_dt_loss = std::numeric_limits<double>::max();
    for (int j = 0; j < nbin; ++j)
    {
        if (Gain[i][j] > max_gain) max_gain = Gain[i][j];
        if (Gain[i][j] > 0.0)
        {
            double ref = std::max(sigmad[i][j], sigmad_floor_[i][j]);
            double dt_gain_j = ref / Gain[i][j];
            if (dt_gain_j < min_dt_gain_own) min_dt_gain_own = dt_gain_j;
        }
        // A donor bin with a large LossRate[i][j] can empty a large fraction
        // of itself within a single sub-step (the implicit division makes
        // that unconditionally stable, so it's deliberately left unbounded
        // above). But Gain fed to *other* bins (e.g. fragmentation debris
        // into small bins) is computed from this donor's y at the *start* of
        // the sub-step and applied at that fixed rate for the whole dt_sub --
        // if the donor itself collapses a lot within that same dt_sub, that
        // Gain is stale/too-large for the back half of the sub-step. Bound
        // dt_sub against each bin's own Loss relaxation time (purely local to
        // this cell's n_sub, not the shared macro dt -- see
        // dust_coagulation_implicit_note.md Sec.10) so the donor's decay is
        // itself resolved finely enough that the Gain it exports stays
        // accurate.
        if (LossRate[i][j] > 0.0)
        {
            double dt_loss_j = 1.0 / LossRate[i][j];
            if (dt_loss_j < min_dt_loss) min_dt_loss = dt_loss_j;
        }
    }

    int n_sub = 1;
    if (max_gain > 0.0 && sigmad_tot[i] > 0.0)
    {
        double dt_gain = sigmad_tot[i] / max_gain;
        n_sub = std::max(n_sub, (int)std::ceil(dt / (GAIN_TOL * dt_gain)));
    }
    if (std::isfinite(min_dt_gain_own) && min_dt_gain_own > 0.0)
    {
        n_sub = std::max(n_sub, (int)std::ceil(dt / (GAIN_TOL * min_dt_gain_own)));
    }
    if (std::isfinite(min_dt_loss) && min_dt_loss > 0.0)
    {
        n_sub = std::max(n_sub, (int)std::ceil(dt / (GAIN_TOL * min_dt_loss)));
    }
    n_sub = std::min(n_sub, N_SUB_MAX);
    double dt_sub = dt / double(n_sub);

    // Mass-conservation correction: applied ONCE, over the whole macro dt,
    // exactly like the explicit ApplyCoagulationAndCorrectMass() does --
    // *not* every sub-step. Gain (explicit, old y) and Loss (implicit, new
    // y) don't exactly cancel at the discrete level even for a single pair,
    // which is what this correction exists to paper over (see the comment
    // above CoagulateChangeRateFragmentation_Implicit). If that per-step
    // mismatch doesn't shrink proportionally with dt_sub (e.g. a structural
    // issue rather than a vanishing truncation error), applying the same
    // uniform proportional rescale n_sub times compounds it multiplicatively
    // -- and a near-empty bin freshly fed by fragmentation is exactly the
    // kind of value a repeated multiplicative rescale can crush toward zero.
    // Matching the explicit scheme's once-per-macrostep cadence removes that
    // compounding regardless of which of the two explanations is right.
    double sigma_tot_entry = sigmad_tot[i];

    for (int s = 0; s < n_sub; ++s)
    {
    // Gain[i]/LossRate[i] passed in already reflect sigmad at the top of
    // this sub-cycle (computed by the caller for s==0; refreshed at the end
    // of the previous iteration below for s>0), so no recompute needed here.
    double temp_sigma_coag = 0.0;

    for (int j = 0; j < nbin; ++j)
    {
        sigmad[i][j] = (sigmad[i][j] + dt_sub*Gain[i][j]) / (1.0 + dt_sub*LossRate[i][j]);
        if(sigmad[i][j] < 0.0)
        {
            sigmad[i][j] = 0.0;
        }
        temp_sigma_coag += sigmad[i][j];
    }

    sigmad_tot[i] = temp_sigma_coag;

    // Re-linearize for the next sub-step: rebuild Gain[i]/LossRate[i] from
    // the sigmad we just produced, instead of reusing the sigmad^n snapshot
    // for the rest of the macro dt (see comment above the function).
    if (s + 1 < n_sub)
    {
        CoagulateChangeRateFragmentation_Implicit_Cell(i, sigmad, Gain, LossRate);
    }
    } // end for (s) sub-cycle

    // Guard against sigmad_tot[i] having decayed to exactly 0 (possible now
    // that nothing rescales it back mid-cycle, unlike the old per-substep
    // correction): sigma_tot_entry/0 = +inf, and inf times an already-zero
    // bin below is 0*inf = NaN, which would silently poison this whole cell
    // (and, via sigmad being shared grid-wide, everything downstream that
    // reads it). If the cell is genuinely empty there is nothing to correct.
    double err = (sigma_tot_entry > 0.0) ? (sigmad_tot[i] - sigma_tot_entry) / sigma_tot_entry : 0.0;
    if (std::abs(err) > 1.0e-6 && sigmad_tot[i] > 0.0)
    {
        double delta_sigma = sigma_tot_entry / sigmad_tot[i];
        double temp_sigma_coag = 0.0;
        for (int j = 0; j < nbin; ++j)
        {
            sigmad[i][j] *= delta_sigma;
            temp_sigma_coag += sigmad[i][j];
        }
        sigmad_tot[i] = temp_sigma_coag;
    }
}

double :: DustDisk :: EstimateDtCoagulation(const std::vector<std::vector<double>>& sigmad,
    const std::vector<std::vector<double>>& dNdt_coag,const std::vector<std::vector<double>>& sigmad_floor,
    int nr_s, int nr_e, int nbin)
{
    double dt_coag = std::numeric_limits<double>::max();

    for (int i = nr_s; i < nr_e; ++i)
    {
        if(T_[i] > 1.5e3)
        {
            continue;
        }

        for (int j = 0; j < nbin; ++j)
        {
            if (sigmad[i][j] > sigmad_floor[i][j] && dNdt_coag[i][j] < 0.0)
            {
                dt_coag = std::min(dt_coag, std::abs(sigmad[i][j] / dNdt_coag[i][j]));
            }
        }
    }
    return dt_coag;
}

// Growth-direction analog of EstimateDtCoagulation, for use_implicit_coag
// (dust_coagulation_implicit_note.md Sec.12): three different internal-only
// n_sub criteria (cell-total-relative Gain, per-bin-floor-relative Gain,
// donor Loss timescale) all failed to resolve the most fragmentation-active
// (innermost) cells' small-dust accuracy, and removing per-substep mass
// renormalization even exposed a division-by-zero/NaN path once the internal
// sub-cycling was left unchecked for long enough. So this restores a GLOBAL
// bound, like EstimateDtCoagulation's dt_coag but for Gain instead of
// depletion -- Loss/depletion still isn't bounded here (the implicit
// division already makes it unconditionally stable, so re-adding that would
// defeat the point of going implicit; only Gain, which stays fully explicit,
// needs this).
//
// Uses sigmad_tot[i] (the cell's *total* dust mass, Sec.6's original
// criterion) as the reference, NOT the per-bin sigmad_floor_ reference used
// by the internal, N_SUB_MAX-capped n_sub criterion in Sec.9. That per-bin
// floor is fine internally because N_SUB_MAX bounds the fallout; used here,
// uncapped, a single freshly-empty bin anywhere in the whole grid collapsed
// the global dt to ~1e-45 yr (Sec.12 postmortem) -- sigmad_tot[i] is a much
// less easily-vanishing scale.
double DustDisk :: EstimateDtGain(const std::vector<double>& sigmad_tot,
    const std::vector<std::vector<double>>& Gain,
    int nr_s, int nr_e, int nbin)
{
    double dt_gain = std::numeric_limits<double>::max();

    for (int i = nr_s; i < nr_e; ++i)
    {
        if(T_[i] > 1.5e3)
        {
            continue;
        }
        if (sigmad_tot[i] <= 0.0)
        {
            continue;
        }

        double max_gain = 0.0;
        for (int j = 0; j < nbin; ++j)
        {
            if (Gain[i][j] > max_gain) max_gain = Gain[i][j];
        }
        if (max_gain > 0.0)
        {
            dt_gain = std::min(dt_gain, sigmad_tot[i] / max_gain);
        }
    }
    return dt_gain;
}

void DustDisk :: ApplyCoagulationAndCorrectMass(std::vector<std::vector<double>>& sigmad,
    std::vector<double>& sigmad_tot,
    const std::vector<std::vector<double>>& dNdt_coag,
    double dt, int i, int nbin)
{
    double prev_sigma_tot = sigmad_tot[i];
    double temp_sigma_coag = 0.0;

    for (int j = 0; j < nbin; ++j)
    {
        sigmad[i][j] += dNdt_coag[i][j] * dt;
        if(sigmad[i][j] < 0.0)
        {
            sigmad[i][j] = 0.0;
        }
        temp_sigma_coag += sigmad[i][j];
    }

    double err = (temp_sigma_coag - prev_sigma_tot) / prev_sigma_tot;

    if (std::abs(err) > 1.0e-6)
    {
        double delta_sigma = prev_sigma_tot / temp_sigma_coag;
        temp_sigma_coag = 0.0;  // 修正：補正後に再計算
        for (int j = 0; j < nbin; ++j)
        {
            sigmad[i][j] *= delta_sigma;
            temp_sigma_coag += sigmad[i][j];
        }
    }

    sigmad_tot[i] = temp_sigma_coag;
}
