#include "./calc_resistivity.hpp"

Calc_Resistivity :: Calc_Resistivity()
{

}

Calc_Resistivity :: ~Calc_Resistivity()
{

}

void Calc_Resistivity :: SetActiveSpecies(bool doH, bool doNa, bool doK)
{
    // Marchand 2022 treats a single thermally ionized species s; independent solves per species
    // cannot simply be summed (n_i, n_e and grain charges would be counted several times)
    if (int(doH) + int(doNa) + int(doK) != 1) {
        std::cerr << "SetActiveSpecies: exactly one species (H, Na or K) must be active" << std::endl;
        exit(1);
    }
    active_H  = doH;
    active_Na = doNa;
    active_K  = doK;
    // Reset warm-start when species selection changes
    n_i_warmstart  = 0.0;
    n_i_H_ws_frac  = 0.0;
    n_i_Na_ws_frac = 0.0;
    n_i_K_ws_frac  = 0.0;
    n_H_ws_frac  = 0.0;
    n_Na_ws_frac = 0.0;
    n_K_ws_frac  = 0.0;
    phi0 = 1.0e-2;
    epsilon_ws = 0.9999;
    phi_ws = 0.0;
    has_phi_ws = false;
}

void Calc_Resistivity :: SetGasParameter_For_Resistivity(double ionizationrate_,double temp,double ng_)
{
    zeta = ionizationrate_;
    T = temp;
    n_gas = ng_;
    n_g = NH_PER_NG * ng_;
}

void Calc_Resistivity :: SetDustParameter_For_Resistivity(DoubleVector nd_,DoubleVector m_,DoubleVector a_)
{
    n_d = 0;
    a   = 0;
    s_d = 0;
    m_d = 0;

    nd_bin = nd_.size();

    // dustのパラメーターの設定
    // dust information

    d_n.resize(nd_bin,0.0);
    d_size.resize(nd_bin,0.0);
    d_s.resize(nd_bin,0.0);
    d_mass.resize(nd_bin,0.0);

    Z_d.resize(nd_bin,0.0);
    tau_k.resize(nd_bin,0.0);
    alphak_.resize(nd_bin,0.0);

    if (nd_.size() != a_.size() || nd_.size() != m_.size()) {
        std::cerr << "Size mismatch: nd=" << nd_.size() << ", a=" << a_.size() << ", m=" << m_.size() << std::endl;
        exit(1);
    }

    for(int i=0;i<nd_bin;i++)
    {
        d_n[i] = nd_[i];
        d_size[i] = a_[i];
        d_s[i] = M_PI*a_[i]*a_[i];
        d_mass[i] = m_[i];
    }

    // Calulate mean value for dust
    n_d = 0;
    a = 0;
    s_d = 0;
    m_d = 0;

    for (int i=0 ; i<nd_bin ; i++)
    {
        n_d += nd_[i];
        a += a_[i] * nd_[i];
        s_d += M_PI*a_[i]*a_[i] * nd_[i];
        m_d += m_[i]*nd_[i];
    }

    a = double (a)/n_d;
    total_cross_section = s_d;
    s_d = double(s_d)/n_d;
    m_d = double(m_d)/n_d;

    // 
    n_Na = Na_H_rate * n_g;
    n_K = K_H_rate * n_g;

    // Thermal velocity
    u_i =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_i));
    u_e =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_e));
    u_H =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_g));
    u_Na =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_Na));
    u_K =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_K));

    // Tau
    tau = double(a*BOLTZMANN_CONSTANT*T)/(CHARGE_UNIT*CHARGE_UNIT); // (2)
    alpha_k = std :: sqrt(8.0/(M_PI*tau)); // Marchand (2021) (39)-(40)
    tau_tot = 0.0;

    for(int i=0;i<nd_bin;i++)
    {
        tau_k[i] = (a_[i]*BOLTZMANN_CONSTANT*T)/(CHARGE_UNIT*CHARGE_UNIT); // (2)
        alphak_[i] = std :: sqrt(8.0/(M_PI*tau_k[i])); // Marchand (2021) (39)-(40)
        tau_tot += nd_[i]*tau_k[i]; // (A.5)
    }

    mui_ = 25.0; // = m_i/PROTON_MASS
    Theta = s_e * std :: sqrt(mui_*HYDROGEN_MASS/ELECTRON_MASS);

    // < Sigma v >  (sigma_v_ie: Marchand 2021 Sect. 3)
    sigma_v_ie = 2.0e-7 * std :: pow(T/3.0e2, -0.5);
    sigma_v_se_H = 3.5e-12 *std :: pow(T/3.0e2,-0.75);
    sigma_v_se_Na =  2.78e-12 * std :: pow(T/3.0e2,-0.68);
    sigma_v_se_K = 2.78e-12 *std :: pow(T/3.0e2,-0.68);
}

void Calc_Resistivity :: SetMagneticField_For_Resistivity(double magneticfield)
{
    mag = magneticfield;
}

void Calc_Resistivity :: SetDebugContext(int i, int k, double r_au, double z, double T, double ng)
{
    prev_dbg_i_ = dbg_i_;
    prev_dbg_k_ = dbg_k_;
    dbg_i_  = i;
    dbg_k_  = k;
    dbg_r_au_ = r_au;
    dbg_z_  = z;
    dbg_T_  = T;
    dbg_ng_ = ng;
}

void Calc_Resistivity :: Calc_ChargeState()
{
    // H,Na,Kに対してそれぞれが独立して電離平衡を計算している
    // そのため、電離平衡が成り立っていても総粒子数が異なる可能性が高い
    // ただ、NaやKは存在量が非常に少ないため総数にはあまり影響を与えない

    double ns_tot;
    double log_ni,log_ns;
    double ni_T,ne_T,Z_T;
    double gosa = 1.0e6;
    DoubleVector Zs_;
    DoubleVector gosa_;

    ns_tot = 0.0;
    ni_T = 0.0;
    ne_T = 0.0;
    Z_T = 0.0;

    std::fill(Z_d.begin(), Z_d.end(), 0.0);
    Zs_.resize(nd_bin,0.0);

    Calculate_Initial_Phi0();

    for(int i=0;i<3;i++)
    {
        gosa = 1.0e6;
        dndt_thermal = 0.0;
        if(i==0)
        {
            if (!active_H) { n_H = 0.0; continue; }
            // Hydrogen
            mus_= 1.0;
            n_s0 = n_g;
            dndt_thermal = beta_H*n_g*std :: sqrt(T)*std :: exp(-T0_H/T); //(23)
            k_i_s = k_i_H;
            k_s_i = k_H_i;
            u_s = u_H;
            sigma_v_se = sigma_v_se_H;
        }
        else if (i==1)
        {
            if (!active_Na) { n_Na = 0.0; continue; }
            // Sodium
            mus_ = 22.99;
            n_s0 = n_Na;
            dndt_thermal = beta_Na*n_g*std :: sqrt(T)*std :: exp(-T0_Na/T); //(23)
            k_i_s = k_i_Na;
            k_s_i = k_Na_i;
            u_s = u_Na;
            sigma_v_se = sigma_v_se_Na;
        }
        else
        {
            if (!active_K) { n_K = 0.0; continue; }
            //Potasium
            mus_ = 39.09; // Table 1
            n_s0 = n_K;
            dndt_thermal = beta_K*n_g*std :: sqrt(T)*std :: exp(-T0_K/T); //(23)
            k_i_s = k_i_K;
            k_s_i = k_K_i;
            u_s = u_K;
            sigma_v_se = sigma_v_se_K;
        }

        if (phi0 < -30.0 || phi0 > 5.0) {
            phi0 = 1.0e-2;
        }
        // first point: Eqs. 37-40, later points: previous solution (Sect. 3.3)
        phi = has_phi_ws ? phi_ws : phi0;
        epsilon = (epsilon_ws > 0.0 && epsilon_ws <= 1.0) ? epsilon_ws : 0.9999;
        epsilon0 = 1.0/Theta;
        delta_epsilon = epsilon - epsilon0; //(33)

        double ni_ws_frac_i = (i==0) ? n_i_H_ws_frac : (i==1) ? n_i_Na_ws_frac : n_i_K_ws_frac;
        n_i = (ni_ws_frac_i > 0.0) ? ni_ws_frac_i*n_g : 1.0e-7*n_g;

        double ns_ws_frac_i = (i==0) ? n_H_ws_frac : (i==1) ? n_Na_ws_frac : n_K_ws_frac;
        n_s = (ns_ws_frac_i > 0.0) ? ns_ws_frac_i*n_s0 : 1.0e-4*n_i;
        n_s = std::min(n_s, n_s0);

        n_e = epsilon*(n_i + qis*n_s);

        dbg_species_ = i;
        Calculate_Newton(gosa,Zs_);

        for(int j=0;j<nd_bin;j++)
        {
            Z_d[j] += Zs_[j];
        }

        ne_T += n_e;
        ni_T += n_i;
        Z_T += Z_k_tot;

        // warm start for the next call
        if      (i==0) { n_H   = n_s; n_i_H_ws_frac  = n_i/n_g; n_H_ws_frac  = (n_s0>0.0) ? n_s/n_s0 : 0.0; }
        else if (i==1) { n_Na  = n_s; n_i_Na_ws_frac = n_i/n_g; n_Na_ws_frac = (n_s0>0.0) ? n_s/n_s0 : 0.0; }
        else           { n_K   = n_s; n_i_K_ws_frac  = n_i/n_g; n_K_ws_frac  = (n_s0>0.0) ? n_s/n_s0 : 0.0; }
    }

    ns_tot = n_H+n_Na+n_K;
    n_i = ni_T;
    n_e = ne_T;
    Z_k_tot = Z_T;

    double C = n_i - n_e + Z_k_tot;
    //std :: cout <<std :: scientific<<n_i/n_g<<" "<<n_e/n_g<<" "<<n_H/n_g<<" "<<n_Na/n_g<<" "<<n_K/n_g<<" "<<Z_d[nd_bin-1]<<" "<<C/n_g<< std :: endl;
    // // exit(1);

    // warm start for the next call
    phi_ws        = phi;
    has_phi_ws    = true;
    epsilon_ws    = epsilon;
    n_i_warmstart = n_i;    // n_i here is ni_T (total HCO+)
}

void Calc_Resistivity :: Calc_MagneticResistivity()
{
    // Collision rates (Marchand et al. 2016 eq. 8; Pinto & Galli 2008), n(H2) = n_H/2.
    // nu_j = a_He n(H2) <sigma v>_j m_H2/(m_j + m_H2)
    //   e    : 3.16e-11 v^1.3,  ions (i, Na+, K+): 2.40e-9 v^0.6  [v in km/s]
    //   grain: Epstein drag 4/3 * 1.3 * pi a^2 v_th(H2) (Tsukamoto & Okuzumi 2022)
    //   H+   : Pinto & Galli 2008 Table 2, H2 and He (n(He) = 0.085 n_H)
    const double M_4pi = 12.56637;
    const double mH2   = 2.0 * PROTON_MASS;
    const double n_H2  = 0.5 * n_g;
    const double a_He_e = 1.16, a_He_i = 1.14;

    double mHCO = m_i;
    double mNa  = m_Na;
    double mK   = m_K;
    double mHp  = HYDROGEN_MASS;

    auto v_kms = [&](double mu) { return std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mu)) * 1.0e-5; };
    auto ion_rate = [&](double m) { return 2.40e-9 * std::pow(v_kms(m*mH2/(m + mH2)), 0.6); };
    auto nu_of = [&](double aHe, double rate, double m) { return aHe * n_H2 * rate * mH2 / (m + mH2); };

    double nu_e  = nu_of(a_He_e, 3.16e-11 * std::pow(v_kms(ELECTRON_MASS), 1.3), ELECTRON_MASS);
    double nu_i  = nu_of(a_He_i, ion_rate(mHCO), mHCO);

    double nu_Na = nu_of(a_He_i, ion_rate(mNa), mNa);
    double nu_K  = nu_of(a_He_i, ion_rate(mK),  mK);

    const double mHe  = 4.0 * HYDROGEN_MASS;
    const double n_He = 0.085 * n_g;
    auto mu_red  = [](double m1, double m2) { return m1*m2/(m1 + m2); };
    auto nu_s = [&](double rate_H2, double rate_He, double m) {
        return n_H2*mH2*rate_H2/(m + mH2) + n_He*mHe*rate_He/(m + mHe); };
    double nu_Hp = nu_s(3.89e-9 * std::pow(v_kms(mu_red(HYDROGEN_MASS, mH2)), -0.02),
                        1.48e-9 * std::pow(v_kms(mu_red(HYDROGEN_MASS, mHe)), -0.02), HYDROGEN_MASS);

    DoubleVector nu_d(nd_bin);
    const double vth_H2 = std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mH2));
    for (int i = 0; i < nd_bin; i++) {
        nu_d[i] = nu_of(1.0, 4.0/3.0*1.3*M_PI*SQR(d_size[i])*vth_H2, d_mass[i]);
    }

    // Use quad precision where available because eta_A cancels strongly in
    // the coupled regime. Fall back to the widest standard floating type;
    // on Apple Silicon, long double has the same precision as double.
#if defined(__SIZEOF_FLOAT128__)
    typedef __float128 quad;
#else
    typedef long double quad;
#endif
    quad beta_e  = -(quad)CHARGE_UNIT * mag / ((quad)ELECTRON_MASS * LIGHT_SPEED) / (quad)nu_e;
    quad beta_i  =  (quad)CHARGE_UNIT * mag / ((quad)mHCO * LIGHT_SPEED) / (quad)nu_i;
    quad beta_Na =  (quad)CHARGE_UNIT * mag / ((quad)mNa  * LIGHT_SPEED) / (quad)nu_Na;
    quad beta_K  =  (quad)CHARGE_UNIT * mag / ((quad)mK   * LIGHT_SPEED) / (quad)nu_K;
    quad beta_Hion = (quad)CHARGE_UNIT * mag / ((quad)mHp * LIGHT_SPEED) / (quad)nu_Hp;

    // Charged grains resolved by charge (Tsukamoto & Okuzumi 2022 Sect. 2.5): per size bin the
    // small-grain populations Z = -1, +1 (fractions f(-1), f(1), Marchand 2022 Eqs. 10-11) plus the
    // large-grain population Z = psi*tau_k, summed. A single mean charge per bin would cancel the
    // Z = -1 and +1 contributions and underestimate the grain conductivities by orders of magnitude.
    std::vector<double> g_n, g_Z, g_m, g_nu;
    for (int i = 0; i < nd_bin; i++) {
        const double x = epsilon*Theta, D = 1.0 + alphak_[i]*x + x*x;
        const double n_pop[3] = { d_n[i]*x*x/D, d_n[i]/D, d_n[i] };   // Z = -1, +1, psi*tau_k
        const double Z_pop[3] = { -1.0, 1.0, phi*tau_k[i] };
        for (int p = 0; p < 3; p++) {
            g_n.push_back(n_pop[p]); g_Z.push_back(Z_pop[p]); g_m.push_back(d_mass[i]); g_nu.push_back(nu_d[i]);
        }
    }
    const int ng_sp = g_n.size();
    std::vector<quad> beta_d(ng_sp);
    for (int i = 0; i < ng_sp; i++) {
        beta_d[i] = (quad)g_Z[i] * CHARGE_UNIT * mag / ((quad)g_m[i] * LIGHT_SPEED) / (quad)g_nu[i];
    }

    // Ohmic conductivity sigma_O
    quad term_e  = -(quad)n_e*CHARGE_UNIT*beta_e;
    quad term_i  =  (quad)n_i*CHARGE_UNIT*beta_i;
    quad term_Na =  (quad)n_Na*CHARGE_UNIT*beta_Na;
    quad term_K  =  (quad)n_K*CHARGE_UNIT*beta_K;
    quad term_H  =  (quad)n_H*CHARGE_UNIT*beta_Hion;
    auto qabs = [](quad v) -> quad { return (v < (quad)0.0) ? -v : v; };
    quad sigma_O = term_e + term_i + term_Na + term_K + term_H;
    quad abs_sum = qabs(term_e) + qabs(term_i) + qabs(term_Na) + qabs(term_K) + qabs(term_H);
    for (int i = 0; i < ng_sp; i++) {
        quad term_d = (quad)g_n[i] * g_Z[i]*CHARGE_UNIT * beta_d[i];
        sigma_O += term_d;
        abs_sum += qabs(term_d);
    }
    // warn on strong cancellation in sigma_O
    if (abs_sum > (quad)0.0 && qabs(sigma_O)/abs_sum < (quad)1.0e-3)
    {
        std::cerr << "sigma_O cancellation: nH=" << n_g << " T=" << T
                  << " i=" << dbg_i_ << " k=" << dbg_k_ << " r=" << dbg_r_au_ << " z=" << dbg_z_
                  << " n_e=" << n_e << " n_i=" << n_i << " n_H=" << n_H
                  << " n_Na=" << n_Na << " n_K=" << n_K
                  << " sigma_O=" << (double)sigma_O << " abs_sum=" << (double)abs_sum
                  << " cancel_ratio=" << (double)(qabs(sigma_O)/abs_sum) << std::endl;
    }
    sigma_O *= (quad)LIGHT_SPEED / (quad)mag;

    // Hall conductivity sigma_H
    quad sigma_H = (quad)n_e*CHARGE_UNIT  * SQR(beta_e)  / ((quad)1.0 + SQR(beta_e))
                 - (quad)n_i*CHARGE_UNIT  * SQR(beta_i)  / ((quad)1.0 + SQR(beta_i))
                 - (quad)n_Na*CHARGE_UNIT * SQR(beta_Na) / ((quad)1.0 + SQR(beta_Na))
                 - (quad)n_K*CHARGE_UNIT  * SQR(beta_K)  / ((quad)1.0 + SQR(beta_K))
                 - (quad)n_H*CHARGE_UNIT  * SQR(beta_Hion) / ((quad)1.0 + SQR(beta_Hion));
    for (int i = 0; i < ng_sp; i++) {
        sigma_H -= (quad)g_n[i] * g_Z[i]*CHARGE_UNIT * SQR(beta_d[i]) / ((quad)1.0 + SQR(beta_d[i]));
    }
    sigma_H *= (quad)LIGHT_SPEED / (quad)mag;

    // Pedersen conductivity sigma_P
    quad sigma_P = -(quad)n_e*CHARGE_UNIT  * beta_e  / ((quad)1.0 + SQR(beta_e))
                 +  (quad)n_i*CHARGE_UNIT  * beta_i  / ((quad)1.0 + SQR(beta_i))
                 +  (quad)n_Na*CHARGE_UNIT * beta_Na / ((quad)1.0 + SQR(beta_Na))
                 +  (quad)n_K*CHARGE_UNIT  * beta_K  / ((quad)1.0 + SQR(beta_K))
                 +  (quad)n_H*CHARGE_UNIT  * beta_Hion / ((quad)1.0 + SQR(beta_Hion));
    for (int i = 0; i < ng_sp; i++) {
        sigma_P += (quad)g_n[i] * g_Z[i]*CHARGE_UNIT * beta_d[i] / ((quad)1.0 + SQR(beta_d[i]));
    }
    sigma_P *= (quad)LIGHT_SPEED / (quad)mag;

    _sigmaO_ = (double)sigma_O;
    _sigmaH_ = (double)sigma_H;
    _sigmaP_ = (double)sigma_P;

    quad c2_4pi = (quad)LIGHT_SPEED*(quad)LIGHT_SPEED / (quad)M_4pi;
    quad etaO_l = c2_4pi / sigma_O;
    quad etaH_l = c2_4pi * sigma_H / (SQR(sigma_H) + SQR(sigma_P));
    // eta_A without the cancellation of sigma_P/(sigma_H^2+sigma_P^2) - 1/sigma_O:
    // with p_j = (c/B) n_j q_j beta_j/(1+beta_j^2) (>= 0 for every species),
    //   sigma_O sigma_P - sigma_P^2 - sigma_H^2 = sum_{j<k} p_j p_k (beta_j - beta_k)^2 >= 0
    std::vector<quad> p_s, b_s;
    auto add_species = [&](quad nq, quad beta) {
        p_s.push_back((quad)LIGHT_SPEED / (quad)mag * nq * beta / ((quad)1.0 + SQR(beta)));
        b_s.push_back(beta);
    };
    add_species(-(quad)n_e*CHARGE_UNIT, beta_e);
    add_species( (quad)n_i*CHARGE_UNIT, beta_i);
    add_species( (quad)n_Na*CHARGE_UNIT, beta_Na);
    add_species( (quad)n_K*CHARGE_UNIT, beta_K);
    add_species( (quad)n_H*CHARGE_UNIT, beta_Hion);
    for (int i = 0; i < ng_sp; i++) add_species((quad)g_n[i] * g_Z[i]*CHARGE_UNIT, beta_d[i]);
    quad pair_sum = 0.0;
    for (size_t j = 0; j < p_s.size(); j++)
        for (size_t k = j + 1; k < p_s.size(); k++)
            pair_sum += p_s[j] * p_s[k] * SQR(b_s[j] - b_s[k]);
    quad etaA_l = c2_4pi * pair_sum / (sigma_O * (SQR(sigma_H) + SQR(sigma_P)));

    _etaO_ = (double)etaO_l;
    _etaH_ = (double)etaH_l;
    _etaA_ = (double)etaA_l;
}

void Calc_Resistivity :: Caluclate_F()
{
//  Linear
    double C = n_i+qis*n_s;
    double C1 = 1.0/((zeta+k_s_i*n_s)*n_g);
    double C2 = 1.0/(zeta + dndt_thermal);
    double C3 = C2 * 1.0/n_s0;

    // neutral reservoir of species i (Eqs. 22/29); +1e-20 n_H keeps F[2] dependent on n_i for H
    double neutral_i_for_creation = (n_g - n_s0 - n_i) + 1.0e-20 * n_g;
    double neutral_i_for_backreaction = std::max(0.0, n_g - n_s0 - n_i);

    F[0] = (1.0-phi)/(epsilon*Theta*expo_phi) - 1.0; //(27)
    F[1] = 1.0 - (epsilon*C - Z_k_tot)/ni_ns_tot; // (28)
    F[2] = neutral_i_for_creation/n_g - (sigma_v_ie*epsilon*n_i*C + n_i*u_i*total_jk + k_i_s*n_i*(n_s0-n_s))*C1; // (29)
    F[3] = (1.0 - n_s/n_s0)*(1.0 + k_i_s*n_i*C2) - ((sigma_v_se*epsilon*n_s*C) + n_s*u_s*total_jk + k_s_i*n_s*neutral_i_for_backreaction)*C3; //(30)
}

void Calc_Resistivity :: J1_Component()
{
    J[0][0] = (phi-2.0)/(epsilon*Theta*expo_phi); //(A.1)
    J[0][1] = -(1.0-phi)/( SQR(epsilon)*Theta*expo_phi ); //(A.2)
    J[0][2] = 0.0; // (A.3)
    J[0][3] = 0.0; // (A.4)
}

void Calc_Resistivity :: J2_Component()
{
    double C = 1.0/ni_ns_tot;
    J[1][0] = C * tau_tot; //(A.5)
    J[1][1] = (-(n_i+qis*n_s) + total_dZdEpsilon) * C; //(A.6)
    J[1][2] = -(epsilon*n_s*(1.0-qis) + Z_k_tot) * SQR(C); // (A.7)
    J[1][3] = -(epsilon*n_i*(qis-1.0) + Z_k_tot) * SQR(C); // (A.8)
}

void Calc_Resistivity :: J3_Component()
{
    double C1 = 1.0/((zeta + k_s_i*n_s) * n_g);
    double C2 = 1.0/(SQR(zeta + k_s_i*n_s) * n_g);

    J[2][0] = n_i*u_i*total_cross_section*C1; //(A.9)
    J[2][1] = (-sigma_v_ie*n_i*(n_i+qis*n_s) - n_i*u_i*total_dJdEpsilon)*C1; //(A.10)
    J[2][2] = -1.0/n_g - (sigma_v_ie*epsilon*(2.0*n_i + qis*n_s) + u_i*total_jk + k_i_s*(n_s0-n_s))*C1; // (A.11)
    J[2][3] = -( sigma_v_ie*epsilon*n_i*(qis*zeta - k_s_i*n_i) - k_s_i*n_i*u_i*total_jk - k_i_s*n_i*(zeta + k_s_i*n_s0) )*C2; // (A.12)
}

void Calc_Resistivity :: J4_Component()
{
    double C1 = 1.0/((zeta + dndt_thermal)*n_s0);
    // C2 = 1/(zeta+thermal), used for the (1-ns/ns0)*ki,s term in (A.15)
    double C2 = C1 * n_s0;
    J[3][0] = n_s*u_s*total_cross_section*C1; //(A.13)
    J[3][1] = -(sigma_v_se*n_s*(n_i+qis*n_s) + n_s*u_s*total_dJdEpsilon)*C1; //(A.14)
    // (A.15): first term divides by (zeta+thermal), second term by (zeta+thermal)*ns0
    J[3][2] = (1.0 - n_s/n_s0)*k_i_s*C2 + (-sigma_v_se*epsilon*n_s + k_s_i*n_s)*C1; // (A.15)
    double neutral_i_for_J = std::max(0.0, n_g - n_s0 - n_i);
    J[3][3] = -1.0/n_s0 - (k_i_s*n_i + sigma_v_se*epsilon*(n_i + 2.0*qis*n_s) + u_s*total_jk + k_s_i*neutral_i_for_J)*C1; // (A.16)
}

void Calc_Resistivity :: Calculate_Initial_Phi0()
{
    int count = 0;
    double phi_n;
    double f,df_dphi;
    double gosa;

    gosa = 1.0e6;

    // // Newton Rapson Method
    phi = 0.5;
    while (gosa > 1.0e-6)
    {
        phi_n = phi;
        expo_phi = std :: exp(phi);
        f = Theta*expo_phi + phi -1.0;
        df_dphi = Theta*expo_phi + 1.0;

        phi -= f/df_dphi;
        gosa = std :: abs((phi - phi_n)/phi);
        if(count > 100)
        {
            std :: cerr << "Charge Calculation ERROR !! Could not get initial phi"
                        << " (i=" << dbg_i_ << " k=" << dbg_k_ << " r=" << dbg_r_au_
                        << " z=" << dbg_z_ << " T=" << dbg_T_ << " ng=" << dbg_ng_ << ")" << std :: endl;
            exit(1);
        }
        count += 1;
    }
    phi0 = phi;

    // if( std :: fabs(phi0*tau) > 10.0)
    // {
    //     std::cerr << "Zk became too Large! Zk = " << phi0*tau << std::endl;
    //     exit(1);
    // }
    // std :: cout <<std :: scientific << phi0 <<" "<<phi<<" "<<phi_n<<" "<<Theta*expo_phi + phi-1.0<<" "<<phi0*tau<<" "<<count<< std :: endl;
    // exit(1);

}

void Calc_Resistivity :: Calculate_Newton(double gosa,DoubleVector &Zd_)
{
    int count = 0;
    double  log_ni,log_ns,u_eps;
    DoubleVector gosa_;

    J.resize(4,DoubleVector(4,0.0));
    F.resize(4,0.0);
    X.resize(4,0.0);

    gosa_.resize(4,0.0);
    qis = std :: sqrt(mui_/mus_); // (11)~(12)
    log_ni = std :: log(n_i);
    log_ns = std :: log(n_s);
    // Newton variables X = (phi, u = eps'/eps0, log n_i, log n_s); eps = eps' + eps0 (Sect. 3.1)
    u_eps = delta_epsilon/epsilon0;

    // DIAG: record entry state (warm start) and reset the history buffer
    diag_n_ = 0;
    diag_init_[0] = phi;
    diag_init_[1] = epsilon;
    diag_init_[2] = n_i/n_g;
    diag_init_[3] = (n_s0 > 0.0) ? n_s/n_s0 : 0.0;

    while (gosa > 1.0e-6)
        {
            expo_phi = std :: exp(phi);
            ni_ns_tot = n_i + n_s;
    
            Z_k_tot = 0.0;
            total_jk = 0.0;
            total_dJdEpsilon = 0.0;
            total_dZdEpsilon = 0.0;
            
            for(int j=0;j<nd_bin;j++)
            {
                Z_k = phi*tau_k[j] - (2.0*delta_epsilon*Theta + SQR(delta_epsilon*Theta))/(1.0 + alphak_[j]*epsilon*Theta + SQR(epsilon*Theta) ); // (31)
                j_k = (1.0-phi) + (2.0/tau_k[j] * ( SQR(epsilon*Theta) + epsilon*Theta ))/( SQR(Theta*epsilon) + alphak_[j]*epsilon*Theta + 1.0); // (32)
    
                dZdEpsilon = -(SQR(epsilon)*CUB(Theta)*alphak_[j] + 4.0*epsilon*SQR(Theta) + alphak_[j]*Theta)/SQR(1.0 + alphak_[j]*epsilon*Theta + SQR(epsilon*Theta)); // (A.17)
                dJdEpsilon = 2.0/tau_k[j] * ((alphak_[j]-1.0)*SQR(epsilon)*CUB(Theta) + Theta + 2.0*epsilon*SQR(Theta) )/SQR(1.0 + alphak_[j]*epsilon*Theta + SQR(epsilon*Theta)); // (A.18)
    
                Z_k_tot += d_n[j]*Z_k; // (A.7)
                total_jk += M_PI*SQR(d_size[j])*d_n[j]*j_k; // (29)
    
                total_dZdEpsilon += d_n[j]*dZdEpsilon; // (A.6)
                total_dJdEpsilon += M_PI*SQR(d_size[j])*d_n[j]*dJdEpsilon; // (A.10)

                Zd_[j] = Z_k;
            }

            // (27) - (30)
            // F(phi,delta_epsilon,a,b) = (f1,f2,f3,f4)
            // X = (phi,epsilon,e^a,e^b) a=log(n_i) , b=log(n_s)
            // dy/da = dy/dn_i dn_i/da = n_i dy/dn_i
            // F_{i+1} = F_{i} + df/dx dx 
            // dx = (F_{i+1} -F_{i} )/dfdx
            // dx = -F_{i}/(df/dx) (F_{i+1} = 0 ) → J * dx = -F

            Caluclate_F();


            J1_Component();
            J2_Component();
            J3_Component();
            J4_Component();

            for(int i=0;i<4;i++)
            {
                J[i][1] = epsilon0*J[i][1]; // d/du
                J[i][2] = n_i*J[i][2];
                J[i][3] = n_s*J[i][3];
            }

            // X = (phi,u=eps'/eps0,log n_i,log n_s)
            std::vector<double> F_neg(F.size());
            for (size_t j = 0; j < F.size(); ++j) {
                F_neg[j] = -F[j];
            }

            X = CALC_TOOL :: Calc_LinearEquation_Pivot(F_neg,J);

            // Guard against NaN/Inf from singular Jacobian (can occur for H where F[2] is ill-posed)
            for (size_t kk = 0; kk < X.size(); ++kk)
                if (!std::isfinite(X[kk])) X[kk] = 0.0;

            NewtonDiagRec &rec = diag_buf_[diag_n_ % kDiagN]; // DIAG
            rec.it = count;
            for (int kk = 0; kk < 4; ++kk) { rec.F[kk] = F[kk]; rec.Xraw[kk] = X[kk]; }

            // Clamp Newton steps to prevent ill-conditioned Jacobian from causing huge jumps
            const double max_log_step = 10.0;
            const double max_phi_step = 2.0;
            if (std::fabs(X[0]) > max_phi_step) X[0] = std::copysign(max_phi_step, X[0]);
            if (std::fabs(X[2]) > max_log_step) X[2] = std::copysign(max_log_step, X[2]);
            if (std::fabs(X[3]) > max_log_step) X[3] = std::copysign(max_log_step, X[3]);
            {   // keep eps within a factor e^{+-max_log_step} per step
                double e_now = 1.0 + u_eps;   // eps/eps0
                double du_lo = e_now*(std::exp(-max_log_step) - 1.0);
                double du_hi = e_now*(std::exp( max_log_step) - 1.0);
                if (X[1] < du_lo) X[1] = du_lo;
                if (X[1] > du_hi) X[1] = du_hi;
            }

            // convergence: relative step (psi, u with a 1e-9 floor; log n_i, log n_s)
            double phi_scale  = std::fabs(phi) + 1.0e-9;
            double u_scale    = std::fabs(u_eps) + 1.0e-9;
            gosa = std :: max( std :: fabs(X[0])/phi_scale, std :: fabs(X[1])/u_scale );
            gosa = std :: max( gosa,std :: fabs(X[2]) );
            gosa = std :: max( gosa,std :: fabs(X[3]) );

            { // DIAG
                double g[4] = { std::fabs(X[0])/phi_scale, std::fabs(X[1])/u_scale, std::fabs(X[2]), std::fabs(X[3]) };
                rec.gosa_arg = 0;
                for (int kk = 1; kk < 4; ++kk) if (g[kk] > g[rec.gosa_arg]) rec.gosa_arg = kk;
                for (int kk = 0; kk < 4; ++kk) rec.X[kk] = X[kk];
                rec.gosa = gosa;
                rec.ns_clamp = 0;
            }

            double alpha=1.0;
            bool accepted = false;
            double phi_trial,u_trial,delta_eps_trial,eps_trial,logni_trial,logns_trial,ni_trial;

            for(int i=0;i<20;i++)
            {
                phi_trial = phi+alpha*X[0];
                u_trial = u_eps+alpha*X[1];
                delta_eps_trial = u_trial*epsilon0;
                eps_trial = delta_eps_trial + epsilon0;
                logni_trial = log_ni+alpha*X[2];
                logns_trial = log_ns+alpha*X[3];

                double ns_trial_val = std::exp(logns_trial);
                if( std :: fabs(logns_trial) > 300.0 || std :: fabs(logni_trial) > 300.0
                    || !(eps_trial > 0.0) || eps_trial > 1.0e3
                    || !std::isfinite(logns_trial) || !std::isfinite(logni_trial) || !std::isfinite(u_trial) )
                {
                    alpha *= 0.5;
                    continue;
                }

                accepted = true;
                break;
            }

            // n_s cannot exceed n_s0
            if (accepted) {
                double ns_trial_val = std::exp(logns_trial);
                if (ns_trial_val > n_s0) {
                    logns_trial = std::log(n_s0 * (1.0 - 1.0e-12));
                    rec.ns_clamp = 1; // DIAG
                }
            }
            rec.alpha = alpha; // DIAG

            if (!accepted) {
                std::cerr << "Newton ERROR: alpha backtracking failed after 20 steps (nH="
                          << n_g << " T=" << T << " i=" << dbg_i_ << " k=" << dbg_k_
                          << " r=" << dbg_r_au_ << " z=" << dbg_z_ << ")" << std::endl;
                diag_n_++;
                DumpNewtonDiag("alpha backtracking failed", Zd_);
            }

            phi = phi_trial;
            u_eps = u_trial;
            epsilon = eps_trial;
            delta_epsilon = delta_eps_trial;
            log_ni = logni_trial;
            log_ns = logns_trial;

            n_i = std :: exp(log_ni);
            // Floor n_i to prevent Jacobian column becoming exactly zero (e.g. for H where F[2] forces ni→0)
            if (n_i < 1.0e-40 * n_g) { n_i = 1.0e-40 * n_g; log_ni = std::log(n_i); }
            n_s = std :: exp(log_ns);
            // n_e = epsilon*(ni+qis*ns) from charge neutrality (f2=0 => eps*(ni+qis*ns)=ne)
            n_e = epsilon*(n_i + qis*n_s);

            // DIAG: state after this iteration's update
            rec.phi = phi; rec.eps = epsilon; rec.ni_ng = n_i/n_g;
            rec.ns_ns0 = (n_s0 > 0.0) ? n_s/n_s0 : 0.0;
            diag_n_++;

            // std :: cout << std :: scientific <<" "<<phi<<" "<<epsilon<<" "<<n_e<<" "<<n_i<<" "<<n_s<<" "<<X[0]<<" "<<X[1]<<" "<<X[2]<<" "<<X[3]<<  std :: endl;
            // exit(1);

            if (count > 5000 || epsilon < 0.0 )
            {
                std :: cerr << std :: scientific << count
                            << " Charge Calculation ERROR: did not converge (nH="
                            << n_g << " T=" << T << " i=" << dbg_i_ << " k=" << dbg_k_
                            << " r=" << dbg_r_au_ << " z=" << dbg_z_ << ")" << std :: endl;
                DumpNewtonDiag("did not converge (count > 5000)", Zd_);
            }
            count += 1;
        }
}

// Print the initial state and the last kDiagN Newton iterations, re-solve once
// from the Eqs. 37-40 guess, then exit(1).
void Calc_Resistivity :: DumpNewtonDiag(const char *reason, DoubleVector &Zd_)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(6);
    const char *sp = (dbg_species_==0) ? "H" : (dbg_species_==1) ? "Na" : (dbg_species_==2) ? "K" : "?";
    os << "===== DIAG " << (diag_retry_ ? "[retry from default guess] " : "") << reason
       << " species=" << sp << " i=" << dbg_i_ << " k=" << dbg_k_ << " r=" << dbg_r_au_
       << " z=" << dbg_z_ << " T=" << T << " nH=" << n_g << " zeta=" << zeta
       << " n_s0=" << n_s0 << " Theta=" << Theta << " phi0=" << phi0 << " eps0=" << epsilon0 << "\n";
    os << "  warm start from previous point (i,k)=(" << prev_dbg_i_ << "," << prev_dbg_k_ << ")"
       << "  init: phi=" << diag_init_[0] << " eps=" << diag_init_[1]
       << " ni/nH=" << diag_init_[2] << " ns/ns0=" << diag_init_[3] << "\n";
    os << "  it  F0 F1 F2 F3 | Xraw0..3 | X0..3 | gosa arg alpha nsclamp | phi eps ni/nH ns/ns0\n";
    int n0 = std::max(0, diag_n_ - kDiagN);
    for (int m = n0; m < diag_n_; ++m) {
        const NewtonDiagRec &r = diag_buf_[m % kDiagN];
        os << "  " << r.it;
        for (int kk = 0; kk < 4; ++kk) os << " " << r.F[kk];
        os << " |";
        for (int kk = 0; kk < 4; ++kk) os << " " << r.Xraw[kk];
        os << " |";
        for (int kk = 0; kk < 4; ++kk) os << " " << r.X[kk];
        os << " | " << r.gosa << " " << r.gosa_arg << " " << r.alpha << " " << r.ns_clamp
           << " | " << r.phi << " " << r.eps << " " << r.ni_ng << " " << r.ns_ns0 << "\n";
    }
    std::cerr << os.str() << std::flush;

    if (!diag_retry_) {
        diag_retry_ = true;
        phi = phi0;
        epsilon = 0.9999;
        delta_epsilon = epsilon - epsilon0;
        n_i = 1.0e-7*n_g;
        n_s = std::min(1.0e-4*n_i, n_s0);
        n_e = epsilon*(n_i + qis*n_s);
        Calculate_Newton(1.0e6, Zd_); // exits via DumpNewtonDiag again if it fails
        std::ostringstream os2;
        os2 << std::scientific << std::setprecision(6)
            << "===== DIAG [retry from default guess] CONVERGED: phi=" << phi << " eps=" << epsilon
            << " ni/nH=" << n_i/n_g << " ns/ns0=" << n_s/n_s0 << " Zk_tot=" << Z_k_tot << "\n";
        std::cerr << os2.str() << std::flush;
    }
    exit(1);
}
