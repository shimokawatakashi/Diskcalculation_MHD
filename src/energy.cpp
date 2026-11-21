#include "./energy.hpp" 
#include <iostream>
#include <algorithm>
#include <cmath>

Energy_Calc :: Energy_Calc()
{

}
Energy_Calc :: ~Energy_Calc()
{

}

void Energy_Calc :: SetGridParameter(int ns,int ne,int n_tot)
{
    nstart = ns;
    nend = ne;
    nr = n_tot;

}

void Energy_Calc :: SetDiskParameter(DoubleVector R_,DoubleVector Sigma_,DoubleVector T_,DoubleVector T_IRR,DoubleVector T_IRR4,DoubleVector Omega_,DoubleVector dSwdt_,DoubleVector Alpha_,DoubleVector Alpha_PZ)
{
    sigmag.resize(nr,0.0);
    omega.resize(nr,0.0);
    alpha.resize(nr,0.0);

    cs_in.resize(nr,0.0);
    cs.resize(nr,0.0);
    cs_out.resize(nr,0.0);

    dSigmaw_dt.resize(nr,0.0);

    rhomid.resize(nr,0.0);
    rhomid_in.resize(nr,0.0);
    rhomid_out.resize(nr,0.0);

    vz_.resize(nr,0.0);
    vz_in.resize(nr,0.0);
    vz_out.resize(nr,0.0);

    T_mid.resize(nr,0.0);
    T_vis.resize(nr,0.0);
    T_irr.resize(nr,0.0);

    T_inner.resize(nr,0.0);
    T_outer.resize(nr,0.0);

    T_irr4.resize(nr,0.0);

    Tin_vis4.resize(nr,0.0);
    T_vis4.resize(nr,0.0);
    Tout_vis4.resize(nr,0.0);

    KR_.resize(nr,0.0);
    KR_out.resize(nr,0.0);
    KR_in.resize(nr,0.0);

    tauR.resize(nr,0.0);
    tauP.resize(nr,0.0);
    tauR_in.resize(nr,0.0);
    tauP_in.resize(nr,0.0);
    tauR_out.resize(nr,0.0);
    tauP_out.resize(nr,0.0);

    Ez_.resize(nr,0.0);
    Ez_in.resize(nr,0.0);
    Ez_out.resize(nr,0.0);

    F_joule.resize(nr,0.0);
    Fin_joule.resize(nr,0.0);
    Fout_joule.resize(nr,0.0);

    F_wind.resize(nr,0.0);
    Fin_wind.resize(nr,0.0);
    Fout_wind.resize(nr,0.0);

    Fin_rad.resize(nr,0.0);
    F_rad.resize(nr,0.0);
    Fout_rad.resize(nr,0.0);
    
    gosa_in.resize(nr,0.0);
    gosa_mid.resize(nr,0.0);
    gosa_out.resize(nr,0.0);

    r_ = R_;
    sigmag = Sigma_;
    dSigmaw_dt = dSwdt_;
    omega = Omega_;
    alpha = Alpha_;
    alpha_pz = Alpha_PZ;

    // Only the raw disk state (T guesses) is set up here. cs_*/rhomid_* are
    // derived quantities that both the equilibrium (bisection) solver and the
    // time-evolution solver recompute for themselves from T_inner/T_mid/T_outer
    // -- via CalcViscousHeating(), which is always called before either solver
    // reads cs_*/rhomid_* -- so they must not be precomputed here too. Doing so
    // used to leave stale values (with a rhomid formula that didn't even match
    // CalcViscousHeating's, missing the sqrt(2*pi) normalization) that were
    // always overwritten before use; keeping the two solvers' setup and derived
    // state cleanly separated avoids that trap.
    for(int i=nstart;i<nend;i++)
    {
        T_outer[i] = 5.0e3;
        T_mid[i]   = T_[i];
        T_irr[i]   = T_IRR[i];
        T_inner[i] = T_irr[i];

        T_irr4[i]  = T_IRR4[i];
    }
}
void Energy_Calc :: SetMagneticField(DoubleVector Bz_,DoubleVector DBDR_,DoubleVector Eta_tot)
{
    B_ = Bz_;
    dBdr_ = DBDR_;
    eta_tot = Eta_tot;
}

void Energy_Calc :: CalculateVerticalVelocity()
{
    for(int i=nstart;i<nend;i++)
    {
        vz_in[i] = dSigmaw_dt[i]/rhomid_in[i];
        vz_[i] = dSigmaw_dt[i]/rhomid[i];
        vz_out[i] =  dSigmaw_dt[i]/rhomid_out[i];
    }
}
void Energy_Calc :: RosselandOpacity()
{
     // According to Kunitomo etal. 2020 (9)
     for (int i=nstart;i<nend;i++)
     {
        KR_[i] = 2.25 * (1.0 - std :: tanh( (T_mid[i]-1.5e3)/1.5e2 ) );
        KR_[i] *= std :: min( 1.0 , SQR(T_mid[i]/1.5e2) );
 
        KR_in[i] = 2.25 * (1.0 - std :: tanh( (T_inner[i]-1.5e3)/1.5e2 ) );
        KR_in[i] *= std :: min( 1.0 , SQR(T_inner[i]/1.5e2) );
 
        KR_out[i] = 2.25 * (1.0 - std :: tanh( (T_outer[i]-1.5e3)/1.5e2 ) );
        KR_out[i] *= std :: min( 1.0 , SQR(T_outer[i]/1.5e2) );
     }
}
void Energy_Calc :: CalcViscousHeating()
{
    for(int i=nstart;i<nend;i++)
    {
        cs_in[i]    = std :: sqrt(BOLTZMANN_CONSTANT*T_inner[i]/GAS_MOLECULAR_MASS);
        cs[i]       = std :: sqrt(BOLTZMANN_CONSTANT*T_mid[i]/GAS_MOLECULAR_MASS);
        cs_out[i]   = std :: sqrt(BOLTZMANN_CONSTANT*T_outer[i]/GAS_MOLECULAR_MASS);

        rhomid_in[i] = sigmag[i]*omega[i]/(std :: sqrt(2.0*M_PI)*cs_in[i]);
        rhomid[i]    = sigmag[i]*omega[i]/(std :: sqrt(2.0*M_PI)*cs[i]);
        rhomid_out[i] = sigmag[i]*omega[i]/(std :: sqrt(2.0*M_PI)*cs_out[i]);

        Fin_rad[i] = 2.25*std :: sqrt(2.0*M_PI)*rhomid_in[i]*CUB(cs_in[i])*alpha[i];
        F_rad[i] = 2.25*std :: sqrt(2.0*M_PI)*rhomid[i]*CUB(cs[i])*alpha[i];
        Fout_rad[i] = 2.25*std :: sqrt(2.0*M_PI)*rhomid_out[i]*CUB(cs_out[i])*alpha[i];
    }
    CalculateVerticalVelocity();
}

void Energy_Calc :: CalcJouleheating()
{
    // Ohmic/ambipolar (Joule) heating of the gas column disabled: the resistive
    // dissipation is concentrated in the low-eta surface active layer, spatially
    // separated from the optically-thick dead midplane, so depositing the full
    // column-integrated eta_tot*(dB/dr)^2 flux into T_mid (with the buried
    // 3*tauR/8 weighting, and a midplane-dominated P-weighted eta_tot ~ eta_cap)
    // over-heats the midplane. F_joule / Fin_joule / Fout_joule stay 0.
    // double a = 0.25/M_PI;
    // for(int i=nstart;i<nend;i++)
    // {
    //     F_joule[i] = a*eta_tot[i] *  SQR(dBdr_[i]) * sigmag[i]/rhomid[i];
    //     Fin_joule[i] =  a*eta_tot[i] *  SQR(dBdr_[i]) * sigmag[i]/rhomid_in[i];
    //     Fout_joule[i] =  a*eta_tot[i] * SQR(dBdr_[i]) * sigmag[i]/rhomid_out[i];
    // }
}

void Energy_Calc :: CalcWindHeating()
{
    for(int i=nstart;i<nend;i++)
    {
        F_wind[i] = wind_heat_fac_*r_[i]*omega[i]*alpha_pz[i]*rhomid[i]*SQR(cs[i]);
        Fin_wind[i] = wind_heat_fac_*r_[i]*omega[i]*alpha_pz[i]*rhomid_in[i]*SQR(cs_in[i]);
        Fout_wind[i] = wind_heat_fac_*r_[i]*omega[i]*alpha_pz[i]*rhomid_out[i]*SQR(cs_out[i]);
    }
}

void Energy_Calc :: CalculateVerticalEnergyFlux()
{
    for(int i=nstart;i<nend;i++)
    {
        Ez_[i] = rhomid[i]*vz_[i]*(0.5*SQR(r_[i]*omega[i]));
        Ez_in[i] = rhomid_in[i]*vz_in[i]*(0.5*SQR(r_[i]*omega[i]) );
        Ez_out[i] = rhomid_out[i]*vz_out[i]*(0.5*SQR(r_[i]*omega[i]) );
    }
}

void Energy_Calc :: CalcTemperature()
{
    // According to Kunitomo etal. 2020 (4)~(8)
    for(int i=nstart;i<nend;i++)
    {
        tauR[i] = KR_[i]*sigmag[i];
        tauP[i] = std :: max(2.4*tauR[i],0.5);

        tauR_in[i] = KR_in[i]*sigmag[i];
        tauP_in[i] = std :: max(2.4*tauR_in[i],0.5);

        tauR_out[i] = KR_out[i]*sigmag[i];
        tauP_out[i] = std :: max(2.4*tauR_out[i],0.5);

        Tin_vis4[i] = 0.5/STEFAN_BOLTZMANN_CONSTANT * (3.0*tauR_in[i]/8.0 + 0.5/tauP_in[i]) * (Fin_rad[i]/*+Fin_joule[i]*/+Fin_wind[i]);
        T_vis4[i] = 0.5/STEFAN_BOLTZMANN_CONSTANT * (3.0*tauR[i]/8.0 + 0.5/tauP[i]) * (F_rad[i]/*+F_joule[i]*/ +F_wind[i]);
        Tout_vis4[i] = 0.5/STEFAN_BOLTZMANN_CONSTANT * (3.0*tauR_out[i]/8.0 + 0.5/tauP_out[i]) * (Fout_rad[i]/*+Fout_joule[i]*/+Fout_wind[i]);

        gosa_in[i] = 1.0 - (T_irr4[i] + Tin_vis4[i])/ std :: pow(T_inner[i],4.0);
        gosa_mid[i] = 1.0 - (T_irr4[i] + T_vis4[i])/ std :: pow(T_mid[i],4.0) ;
        gosa_out[i] = 1.0 - (T_irr4[i] + Tout_vis4[i])/std :: pow(T_outer[i],4.0);
    }
}

void Energy_Calc :: GetTemperature_Initialization(DoubleVector &T_)
{
    int N = 0;
    double gosa = 1.0e6;

    while (gosa > 1.0e-6)
    {
        RosselandOpacity();
        CalcViscousHeating();
        CalcTemperature(); 

        double gosa_local = 0.0;

        #pragma omp parallel for reduction(max: gosa_local)
        for (int i = nstart; i < nend; i++) {
            double corr = fabs(gosa_mid[i]);
            gosa_local = std::max(gosa_local, corr);

            // Converged cells keep their T_mid -- see the note in GetTemperature().
            if (corr <= 1.0e-6) {
                continue;
            }

            if (gosa_in[i] <= 0 && gosa_out[i] >= 0) {
                if (gosa_mid[i] < 0) {
                    T_inner[i] = T_mid[i];
                } else {
                    T_outer[i] = T_mid[i];
                }
                T_mid[i] = 0.5 * (T_inner[i] + T_outer[i]);
            } else {
                std::cerr << "#ERROR Thermal_Structure :: GetTemperature : I cannot solve bisection method" << std::endl;
                std::cout << std::scientific << N << " " << i << " " << gosa_in[i] << " " << gosa_mid[i] << " " << gosa_out[i] << std::endl;
                exit(1);
            }
        }

        gosa = gosa_local;
       
        N += 1;
        if (N>100)
        {
            std :: cerr << "#ERROR Thermal_Structure :: GetTemperature :" <<gosa << std :: endl;
            exit(1);
            break;
        }
    }
    T_ = T_mid;
}

void Energy_Calc :: GetTemperature(DoubleVector &T_)
{
    int N = 0;
    double gosa = 1.0e6;
    DoubleVector rel_corr;

    rel_corr.resize(nr,0.0);

    while (gosa > 1.0e-6)
    {
        RosselandOpacity();
        CalcViscousHeating();
        // CalcJouleheating();   // Joule heating of the midplane disabled (see CalcJouleheating)
        CalcWindHeating();
        CalcTemperature();
        for (int i=nstart;i<nend;i++)
        {
            rel_corr[i] = fabs(gosa_mid[i]);

            // Check convergence BEFORE moving T_mid. Without this, a cell whose
            // incoming T_mid is already the equilibrium (gosa_mid ~ 0 on the
            // first pass, before the [T_irr, 5000] bracket is ever narrowed)
            // still gets one bracket update -> T_mid jumps to 0.5*(T_irr+5000)
            // ~ 2500 K and the loop then exits on the pre-update residual,
            // returning that spurious midpoint (see DetermineVaporDust
            // sublimation cascade). Converged cells keep their T_mid.
            if (rel_corr[i] <= 1.0e-6)
            {
                continue;
            }

            if(gosa_in[i] <= 0 && gosa_out[i] >= 0)
            {
                if(gosa_mid[i]<0)
                {
                    T_inner[i] = T_mid[i];
                }
                else
                {
                    T_outer[i] = T_mid[i];
                }
                T_mid[i] = 0.5*(T_inner[i]+T_outer[i]);
            }
            else
            {
                std :: cerr << "#ERROR Thermal_Structure :: GetTemperature : I cannnot solve bisection method" << std :: endl;
                std :: cout <<std :: scientific <<N<<" "<<i << " "<<gosa_in[i] <<" "<<gosa_mid[i]<<" "<<gosa_out[i]<< std :: endl;
                exit(1);
            }
        }
        gosa = *std :: max_element(rel_corr.begin(),rel_corr.end());
       
        N += 1;
        if (N>100)
        {
            std :: cerr << "#ERROR Thermal_Structure :: GetTemperature :" <<gosa << std :: endl;
            exit(1);
            break;
        }
    }
   
    for(int i=nstart;i<nend;i++)
    {
        if(T_mid[i] < 10.0)
        {
            T_mid[i] = 10.0;
        }
        //std :: cout << std :: scientific << i <<" "<<T_[i]<<" "<<F_rad[i]<<" "<<F_joule[i]<< std :: endl;
    }
    T_ = T_mid;
    // exit(1);
}

double Energy_Calc :: Calculate_Dt(DoubleVector &vgr_,DoubleVector &r_bod_,DoubleVector &vol_inv)
{
    double heat;
    double dt_src,C;
    double T4;
    DoubleVector dt_Esrc;

    vg_ = vgr_;
    rbod_ = r_bod_;
    inv_vol_ = vol_inv;

    Q_irr.resize(nr,0.0);
    dEdt_src.resize(nr,0.0);
    E_tot.resize(nr,0.0);
    E_mag.resize(nr,0.0);
    Lambda_rad.resize(nr,0.0);
    dt_Esrc.resize(nr,0.0);
    P_int.resize(nr,0.0);
    Etot_flux.resize(nr+1,0.0);
    v_bod_.resize(nr+1,0.0);

    //　運動エネルギーと重力エネルギーについては無視する
    RosselandOpacity();
    CalcViscousHeating();
    // CalcJouleheating();   // Joule heating of the midplane disabled (see CalcJouleheating)
    CalcWindHeating();

    for(int i=nstart;i<nend;i++){
        T4 = SQR( SQR(T_mid[i]) );
        tauR[i] = KR_[i]*sigmag[i];
        C = 8.0*tauR[i]*STEFAN_BOLTZMANN_CONSTANT/(3.0*SQR(tauR[i])+1.0);

        P_int[i] = sigmag[i]*SQR(cs[i]);
        E_tot[i] = P_int[i]*frac_1gamma_ad;

        Q_irr[i] = C*T_irr4[i];
        Lambda_rad[i] = C*T4;

        heat = F_rad[i]/*+F_joule[i]*/+Q_irr[i]+F_wind[i];

        dEdt_src[i] =  F_rad[i]/*+F_joule[i]*/+Q_irr[i]+F_wind[i]-Lambda_rad[i];
        // Lambda_rad (T^4 radiative cooling) is handled implicitly below in
        // GetTemperature_TimeEvolution, so it is unconditionally stable and must
        // not drive the CFL-like dt estimate here. Near radiative equilibrium
        // heat and cool are individually large but nearly cancel (dEdt_src~0);
        // using max(heat,cool) as before kept dt pinned to that large-but-
        // cancelling scale forever. heat alone (the genuinely explicit terms)
        // is a much less stiff function of T, so it is the right dt bound now.
        dt_Esrc[i] = std :: fabs(E_tot[i]/ heat );
    };

    E_tot[nstart-1] = E_tot[nstart];
    E_tot[nend] = E_tot[nend-1];
    dEdt_src[nstart-1] = dEdt_src[nstart];
    dEdt_src[nend] = dEdt_src[nend-1];
    P_int[nstart-1] = P_int[nstart];
    P_int[nend] = P_int[nend-1]; 

    double s,t;
    double vbod,Ebod;

    // E_tot is internal energy only (kinetic/potential energy are ignored, see
    // above), so the flux it is advected with is just v*E -- not v*(E+P). The
    // v*(E+P) "enthalpy-flux" form is only correct for a genuine total-energy
    // equation (or if a compensating v.grad(P) source is added back); using it
    // here without that correction silently introduces a spurious v.grad(P)
    // term into the internal-energy budget. Compression/expansion (PdV) work is
    // instead added explicitly below as a -P*div(v) source term in dEdt_src,
    // which is the standard (non-conservative) form of the internal-energy
    // equation: dE/dt + div(E*v) + P*div(v) = Q.
    for(int i=nstart;i<=nend;i++){
        s = std :: sqrt(sigmag[i])/(std :: sqrt(sigmag[i])+std :: sqrt(sigmag[i-1]));
        t = 1.0 - s;

        vbod = s*vg_[i]+t*vg_[i-1];
        Ebod = s*E_tot[i]+t*E_tot[i-1];
        v_bod_[i] = vbod;
        Etot_flux[i] = vbod*Ebod;
    }

    // Compressive heating/cooling: -P*div(v), differenced the same way as the
    // advective fluxes above (2*(r_bod[i+1]*v[i+1]-r_bod[i]*v[i])*inv_vol[i]).
    for(int i=nstart;i<nend;i++){
        double div_v = 2.0*(rbod_[i+1]*v_bod_[i+1]-rbod_[i]*v_bod_[i])*inv_vol_[i];
        dEdt_src[i] -= P_int[i]*div_v;
    }

    // --- Radial radiative diffusion : frozen face conductivity for the
    // implicit solve in GetTemperature_TimeEvolution ---
    // The original model had no radial energy transport except advection, so a
    // grid-scale T kink (e.g. at the v_frag T=160K fragmentation switch, via the
    // dust traffic-jam -> drag back-reaction -> gas bump chain) had nothing to
    // damp it. Column-integrated Rosseland radial conductivity
    //   chi_col = (16/3) sigma_SB T^3 H_g / (kappa_R rho_mid)   [erg/(s K)]
    // giving a radial radiative flux  F = chi_col dT/dr  [erg/(s cm)], frozen
    // here at the current T (kappa/rho too) so the backward-Euler solve stays
    // linear, and flux-limited to the free-streaming value c a T^4 H so it
    // degrades where the disk is optically thin. Whenever it is on it is solved
    // implicitly (tridiagonal), never explicit -- so no dt constraint. Strength
    // is the input parameter "rad_diff_fac" (default 1; 0 = off).
    chi_face_.assign(nr+1, 0.0);
    if (rad_diff_fac_ > 0.0)
    {
        for (int i = nstart+1; i < nend; ++i)        // interior faces; edges stay 0 (insulating)
        {
            double Hl   = cs[i-1]/std::max(omega[i-1],1e-30);
            double Hr   = cs[i]  /std::max(omega[i],  1e-30);
            double chil = (16.0/3.0)*STEFAN_BOLTZMANN_CONSTANT*CUB(T_mid[i-1])*Hl
                          / (std::max(KR_[i-1],1e-30)*std::max(rhomid[i-1],1e-30));
            double chir = (16.0/3.0)*STEFAN_BOLTZMANN_CONSTANT*CUB(T_mid[i])*Hr
                          / (std::max(KR_[i],1e-30)*std::max(rhomid[i],1e-30));
            double chif = rad_diff_fac_ * 2.0*chil*chir/(chil+chir+1e-300);  // harmonic (series) mean
            double drc  = r_[i]-r_[i-1];
            double Fun  = chif*(T_mid[i]-T_mid[i-1])/drc;                    // unlimited [erg/(s cm)]

            double Tf   = 0.5*(T_mid[i-1]+T_mid[i]);
            double Hf   = 0.5*(Hl+Hr);
            double Ffs  = LIGHT_SPEED*RADIATION_CONSTANT*SQR(SQR(Tf))*Hf;    // free-streaming cap
            double lam  = 1.0/(1.0 + std::fabs(Fun)/std::max(Ffs,1e-300));
            chi_face_[i] = lam*chif;
        }
    }

    dt_src = *std :: min_element(dt_Esrc.begin()+nstart,dt_Esrc.begin()+nend);
    return dt_src;
}

void Energy_Calc :: GetTemperature_TimeEvolution(DoubleVector &T_,double dt)
{
    // T = C_gas * E / Sigma  (ideal-gas internal energy <-> temperature).
    // Operator split, 1st order:
    //   (1) explicit  : advection + all non-radiative sources        -> E_star
    //   (2) implicit  : radial radiative diffusion (tridiagonal)     -> E1
    //   (3) implicit  : per-cell Lambda_rad ~ T^4 radiative cooling  -> E_new
    // Steps 2 and 3 are both unconditionally stable, so dt is unconstrained by
    // either the stiff cooling or the diffusion.
    double C_gas = GAS_MOLECULAR_MASS*gamma_ad_1/(BOLTZMANN_CONSTANT);
    int n = nend - nstart;

    // --- (1) explicit: advection + non-radiative sources ---
    DoubleVector E_start(nr,0.0), E_star(nr,0.0);
    for(int i=nstart;i<nend;i++){
        E_start[i] = E_tot[i];   // start-of-step E (linearization ref for step 3)
        E_star[i]  = E_tot[i]
                    - 2.0*dt*(rbod_[i+1]*Etot_flux[i+1]-rbod_[i]*Etot_flux[i])*inv_vol_[i]
                    + (dEdt_src[i] + Lambda_rad[i])*dt;   // dEdt_src has -Lambda_rad folded in; add it back
    }

    // --- (2) implicit radial radiative diffusion: (I - dt*L_diff) E1 = E_star ---
    // L_diff acts on T = C_gas*E/Sigma with the frozen face conductivity
    // chi_face_ from Calculate_Dt. Insulating domain edges (chi_face_=0 there),
    // so the tridiagonal system is naturally closed. M-matrix -> Thomas stable.
    DoubleVector E1 = E_star;
    if (rad_diff_fac_ > 0.0 && (int)chi_face_.size() >= nend+1 && n > 1)
    {
        DoubleVector lo(n,0.0), di(n,1.0), up(n,0.0), rhs(n,0.0), cp(n,0.0), dp(n,0.0);
        for(int i=nstart;i<nend;i++){
            int k = i-nstart;
            double g   = dt*2.0*inv_vol_[i]*C_gas;
            double aL  = rbod_[i]  *chi_face_[i]  /std::max(r_[i]-r_[i-1],1e-30);   // face i
            double aR  = rbod_[i+1]*chi_face_[i+1]/std::max(r_[i+1]-r_[i],1e-30);   // face i+1
            lo[k]  = -g*aL/std::max(sigmag[i-1],1e-30);
            up[k]  = -g*aR/std::max(sigmag[i+1],1e-30);
            di[k]  =  1.0 + g*(aL+aR)/std::max(sigmag[i],1e-30);
            rhs[k] =  E_star[i];
        }
        cp[0]=up[0]/di[0]; dp[0]=rhs[0]/di[0];
        for(int k=1;k<n;k++){
            double m = di[k]-lo[k]*cp[k-1];
            cp[k]=up[k]/m;
            dp[k]=(rhs[k]-lo[k]*dp[k-1])/m;
        }
        E1[nend-1]=dp[n-1];
        for(int k=n-2;k>=0;k--) E1[nstart+k]=dp[k]-cp[k]*E1[nstart+k+1];
    }

    // --- (3) implicit per-cell Lambda_rad ~ T^4 cooling on E1 ---
    // Linearize Lambda_rad(T_new) about T_old = T_mid[i] (opacity/tauR frozen),
    // dLambda/dT = 4*Lambda_rad/T ; backward-Euler in E.
    for(int i=nstart;i<nend;i++){
        double T_old = T_mid[i];
        double dLambda_dT = 4.0*Lambda_rad[i]/std::max(T_old,1.0);
        double alpha = C_gas/sigmag[i];           // dT/dE
        double beta  = dt*dLambda_dT*alpha;       // implicit stiffness factor
        E_tot[i] = (E1[i] - dt*Lambda_rad[i] + beta*E_start[i]) / (1.0 + beta);

        rhomid[i] = sigmag[i]*omega[i]/(std :: sqrt(2.0*M_PI)*cs[i]);
        T_mid[i]  = C_gas*E_tot[i]/sigmag[i];
        if(T_mid[i] < 10.0){
            T_mid[i] = 10.0;
        }
    }
    T_ = T_mid;
}