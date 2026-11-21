#include "typedef.hpp"
#include "physical_constants.hpp"
#include "alpha.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>

// constexpr double SQRT_1_2PI = 0.398942;
// constexpr double n_H2fra    = 0.840336;
// constexpr double n_Hefra    = 0.159664;
// constexpr double M_4pi       = 12.56637;


CalcAlphaParameter :: CalcAlphaParameter()
{

}

CalcAlphaParameter :: ~CalcAlphaParameter()
{

}

/*
    Parameter Set
*/
void CalcAlphaParameter :: SetParameter(double ionizationrate_,double Temp,double ng_, DoubleVector nd_,DoubleVector m_,DoubleVector a_)
{
    n_g = ng_;
    zeta = ionizationrate_;
    T = Temp;

    n_d = 0;
    a   = 0;
    s_d = 0;
    m_d = 0;

    nd_bin = nd_.size();

    // dust_ = dust;
    // int nb = dust_._number_bin_;

    // dustのパラメーターの設定
    // dust information

    d_n.resize(nd_bin,0.0);
    d_size.resize(nd_bin,0.0);
    d_s.resize(nd_bin,0.0);
    d_mass.resize(nd_bin,0.0);

    for(int i=0;i<nd_bin;i++)
    {
        d_n[i] = nd_[i];
        d_size[i] = a_[i];
        d_s[i] = M_PI*a_[i]*a_[i];
        d_mass[i] = m_[i];
        //std :: cout <<std :: scientific << d_n[i] <<" " <<  d_size[i] << " " <<  d_s[i] <<" " << d_mass[i] <<" " <<n_g<< std :: endl;
    }

    // dust total_numberdensity
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
    s_d = double(s_d)/n_d;
    m_d = double(m_d)/n_d;

    // Temeratureの設定
    // Thermal velocity
    u_i =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_i));
    u_e =  sqrt(8.0*BOLTZMANN_CONSTANT*T/double(M_PI*m_e));

    // Recombination Rate
    beta = 2.4e-7*pow(double(T)/300.0,-0.69);
    // Tau
    tau = double(a*BOLTZMANN_CONSTANT*T)/(CHARGE_UNIT*CHARGE_UNIT);

    //std :: cout << std :: scientific << n_g <<" " << n_d << " " << m_d << " " << s_d << " " << beta << " " <<tau << std :: endl;
    return;
}

void CalcAlphaParameter :: SetParameterForAlpha(double plasmabeta,double magneticfield,double Omega_Kepler)
{
    p_beta = plasmabeta;
    mag = magneticfield;
    omega_k = Omega_Kepler;
}

/*
    Chemical Reaction
*/

double CalcAlphaParameter :: Omega_cross()
{
    double a_1;
    double a_2;

    a_1=j_0;                                  //j_0
    a_2=(1+1.0/tau)*(1+sqrt(2.0/(tau+2.0)));  //j_1

    return  double(a_1)/a_2 ;
}
// Cross Section 

double CalcAlphaParameter :: Ion_cross_section()
{
    double g_1 = exp_phi_/double(1-phi)+1.0;
    double g_2 = func_theta;

    return  j_0*double(g_1)/g_2 + 1.0-phi;
}

double CalcAlphaParameter :: Electron_cross_section()
{
    double g_1=exp_phi_/double(1-phi) * (exp_phi_/double(1-phi) + 1);
    double g_2=func_theta;

    return j_0*double(g_1)/g_2 + exp_phi_;
}

// dustの平均電荷の計算
double CalcAlphaParameter :: total_dustcharge()
{
    double g_1=exp_2phi_/((1-phi)*(1-phi)) - 1.0;
    double g_2=func_theta;

    return omega*g_1/double(g_2)+phi*tau;;
}

// Ion Electoronの数密度の計算

double CalcAlphaParameter :: total_ion_numberdensity()
{
    double F1=sqrt(F)-1.0;
    return C1_i*j_e*F1;
}

double CalcAlphaParameter :: total_electron_numberdensity()
{
    double F1=sqrt(F)-1.0;

    return  C1_e*j_i*F1;
}

// Newton法で使う微分の計算
// CrossSectionの微分
double CalcAlphaParameter :: diff_cross_section_i()
{
    //分子
    double func_g2 = exp_2phi_ /SQR((1.0-phi)) + 1.0;
    double diff_g2 = 2.0 * (exp_2phi_/CUB((1.0-phi)) + exp_2phi_/SQR((1.0-phi)));

    //分母
    double func_g = func_theta;
    double diff_g = diff_theta;

    return j_0*(-diff_g*func_g2 + diff_g2*func_g)/SQR(func_g) - 1.0;
}

double CalcAlphaParameter :: diff_cross_section_e()
{
    //分子
    double func_g3 = exp_2phi_ / SQR(1.0 - phi)+ exp_phi_ / (1.0-phi);
    double diff_g3 = 2.0*exp_2phi_/CUB(1.0-phi) + (2.0*exp_2phi_ + exp_phi_) / SQR((1.0-phi)) + exp_phi_/(1.0-phi);

    //分母
    double func_g = func_theta;
    double diff_g = diff_theta;

    return j_0 * (-diff_g*func_g3 + diff_g3*func_g)/SQR(func_g) + exp_phi_;
}

double CalcAlphaParameter :: diff_function_f()
{
    double F_1 = 0.5*C2/sqrt(F);
    double F_2 = diff_j_e/double(j_e*j_e*j_i) + diff_j_i/double(j_e*j_i*j_i);

    return  -F_1*F_2; 
}

// 数密度の微分

double CalcAlphaParameter :: diff_total_dustcharge()
{
    //分子
    double func_g1 = exp_2phi_/((1.0-phi)*(1.0-phi)) - 1.0;
    double diff_g1 = 2.0*exp_2phi_*pow(1.0-phi,-2.0) + 2.0*exp_2phi_*pow(1.0-phi,-3.0);
    //分母
    double func_g = func_theta;
    double diff_g = diff_theta;
   
    return omega*(-diff_g*func_g1+diff_g1*func_g)/(func_g*func_g)+tau;
}

double CalcAlphaParameter :: diff_total_ion_numberdensity()
{
    double F1=sqrt(F)-1.0;
    
    return C1_i*(diff_j_e*F1 + j_e*diff_f);
}

double CalcAlphaParameter :: diff_total_electron_numberdensity()
{
    double F1=sqrt(F)-1.0;

    return C1_e * (diff_j_i*F1 + j_i*diff_f);
}

// 熱電離の計算
void CalcAlphaParameter :: CalcThermalIonization() // Saha's function
{
    // Saha's function
    double S_k;
    double g_Kion,g_K,gion_g;
    double K_eV;
    g_Kion = 1.0;
    g_K = 2.0;
    gion_g = g_Kion/g_K;
    n_K = n_g*2.2e-10; // カリウムの数密度
    // n_Na = n_g*3.10e-9; // ナトリウムの数密度

    K_eV = 4.34*ELECTRON_VOLT;
    S_k = std :: pow(2.0*gion_g*(2.0*M_PI*m_e*BOLTZMANN_CONSTANT*T/ SQR(PLANCK_CONSTANT) ) ,1.5)*std :: exp(-K_eV/(BOLTZMANN_CONSTANT*T) );

    ne_th = 0.5*( -S_k+std :: sqrt( SQR(S_k)+4.0*S_k*n_K ) );
    n_e = n_e + ne_th;
    nK_th = ne_th;
}

// Newton法の計算
void CalcAlphaParameter :: GetChargestate(double phi_init)
{
    // Without dust the grain-charging Newton system below is singular: it is
    // built around the mean-grain quantities n_d/s_d/tau (n_d = total dust
    // number density), so n_d=0 makes C1_i, C1_e, C2 divide by zero (and
    // a/s_d/m_d in SetParameter() are already 0/0 before that). Charge
    // neutrality without charged grains reduces to zeta*n_g = beta*n_i*n_e
    // with n_i = n_e, i.e. n_i = n_e = sqrt(zeta*n_g/beta) -- pure gas-phase
    // ionization/recombination equilibrium. Skip the Newton solve entirely
    // and use that directly; every dust bin is empty here anyway (n_d is a
    // sum of non-negative per-bin densities), so zeroing the low/high-tau
    // charge distributions below is exact, not an approximation -- and
    // CalcResistivity()'s per-bin dust loop already treats a zero-density
    // bin as a zero contribution.
    if (n_d < 1.0e-50)
    {
        phi = 0.0;
        q_d = 0.0;
        n_i = std::sqrt(zeta * n_g / beta);
        n_e = n_i;
        // Thermal (Saha) ionization of potassium, T > 1000K: an independent
        // equilibrium channel on top of the cosmic-ray-sourced n_i/n_e pair
        // above, same additive convention as CalcThermalIonization()'s own
        // "n_e = n_e + ne_th" (K+/e- treated as a separate species, not
        // folded into n_i, so CalcResistivity_ForThermal() can weight it by
        // its own mobility beta_K). nK_th must be reset to 0 below 1000K --
        // it is a member, so a stale value from a hotter cell/step would
        // otherwise leak into this step's sigma_O/H/P.
        if (T > 1000.0)
        {
            CalcThermalIonization();
        }
        else
        {
            nK_th = 0.0;
        }

        const int N = 101;
        _dnumber_low_.assign(nd_bin, std::vector<double>(3, 0.0));
        _dnumber_high_.assign(nd_bin, std::vector<double>(N - 1, 0.0));
        _dust_charge_low_.assign(nd_bin, std::vector<double>(3, 0.0));
        _dust_charge_high_.assign(nd_bin, std::vector<double>(N - 1, 0.0));
        return;
    }

    // Size-resolved charge neutrality (no mean grain size).
    // phi (normalized grain potential, eps = e^phi/(1-phi)) is common to all
    // bins -- the same assumption the per-bin charge distribution below
    // already makes -- while tau_k, j_0k, omega_k, and the grain cross section
    // s_k are evaluated per bin. Mean-grain quantities are replaced by sums:
    //   A_i = u_i s_i sum_k n_k s_k J_i(tau_k,phi)   (ion loss rate onto grains / n_i)
    //   A_e = u_e s_e sum_k n_k s_k J_e(tau_k,phi)   (electron loss rate / n_e)
    //   Q   = sum_k n_k Z_k(tau_k,phi)               (total grain charge / e)
    // Ionization balance  zeta n_g = beta n_i n_e + A_i n_i = beta n_i n_e + A_e n_e
    // gives n_i = A_e/(2beta)(sqrt(F)-1), n_e = A_i/(2beta)(sqrt(F)-1),
    // F = 1 + 4 beta zeta n_g/(A_i A_e). With one bin this reduces exactly to
    // the former mean-grain formulas (C1_i, C1_e, C2).
    std::vector<double> j0_k(nd_bin), omega_k(nd_bin), tau_k(nd_bin), ns_k(nd_bin);
    for (int k = 0; k < nd_bin; k++)
    {
        tau_k[k]   = d_size[k]*BOLTZMANN_CONSTANT*T / (CHARGE_UNIT*CHARGE_UNIT);
        j0_k[k]    = 1.0 + sqrt(M_PI/(2.0*tau_k[k]));
        omega_k[k] = j0_k[k] / ((1.0+1.0/tau_k[k])*(1.0+sqrt(2.0/(tau_k[k]+2.0))));
        ns_k[k]    = d_n[k]*d_s[k];
    }
    const double coef_i = u_i*s_i;
    const double coef_e = u_e*s_e;
    const double C_F    = 4.0*beta*zeta*n_g;

    int count=0;
    int max_count=1.0e2;
    double gosa=1.0;

    double diff_n_e,diff_n_i,diff_Z;
    double charge_neutrality,diff_charge_neutrality;
    double phi_next;

    phi = phi_init;

    while (gosa >1.0e-5)
    {
        exp_phi_ = exp(phi);
        exp_2phi_ = exp_phi_*exp_phi_;

        const double om1   = 1.0 - phi;
        const double eps   = exp_phi_/om1;                     // e^phi/(1-phi)
        const double eps2  = eps*eps;                          // e^{2phi}/(1-phi)^2
        const double deps  = exp_phi_/om1 + exp_phi_/(om1*om1); // d eps/dphi
        const double deps2 = 2.0*eps*deps;                     // d eps^2/dphi

        // phi-only parts of the cross sections (high-tau terms) and their numerators
        const double gi  = eps + 1.0,        dgi  = deps;          // ion   low-tau numerator
        const double ge  = eps2 + eps,       dge  = deps2 + deps;  // electron
        const double gz  = eps2 - 1.0,       dgz  = deps2;         // charge

        double A_i = 0.0, A_e = 0.0, Q = 0.0;
        double dA_i = 0.0, dA_e = 0.0, dQ = 0.0;

        for (int k = 0; k < nd_bin; k++)
        {
            const double om  = omega_k[k];
            const double th  = om*eps2 + eps + om;    // theta_k
            const double dth = om*deps2 + deps;
            const double ith = 1.0/th;
            const double ith2 = ith*ith;

            const double Ji  = j0_k[k]*gi*ith + om1;
            const double dJi = j0_k[k]*(dgi*th - gi*dth)*ith2 - 1.0;
            const double Je  = j0_k[k]*ge*ith + exp_phi_;
            const double dJe = j0_k[k]*(dge*th - ge*dth)*ith2 + exp_phi_;
            const double Z   = om*gz*ith + phi*tau_k[k];
            const double dZ  = om*(dgz*th - gz*dth)*ith2 + tau_k[k];

            A_i  += ns_k[k]*Ji;   dA_i += ns_k[k]*dJi;
            A_e  += ns_k[k]*Je;   dA_e += ns_k[k]*dJe;
            Q    += d_n[k]*Z;     dQ   += d_n[k]*dZ;
        }
        A_i *= coef_i;  dA_i *= coef_i;
        A_e *= coef_e;  dA_e *= coef_e;

        F = 1.0 + C_F/(A_i*A_e);
        const double sF  = sqrt(F);
        const double dsF = -0.5*(F-1.0)/sF * (dA_i/A_i + dA_e/A_e);   // d sqrt(F)/dphi

        n_i = A_e/(2.0*beta)*(sF-1.0);
        n_e = A_i/(2.0*beta)*(sF-1.0);
        q_d = Q/n_d;   // number-weighted mean grain charge (diagnostic only)

        diff_n_i = (dA_e*(sF-1.0) + A_e*dsF)/(2.0*beta);
        diff_n_e = (dA_i*(sF-1.0) + A_i*dsF)/(2.0*beta);
        diff_Z   = dQ;

        charge_neutrality = n_i - n_e + Q;
        diff_charge_neutrality = diff_n_i - diff_n_e + diff_Z;

        phi_next = phi - double(charge_neutrality)/diff_charge_neutrality;
       
        gosa = fabs((phi_next - phi)/phi_next);

        if (count == max_count)
        {
            std :: cout << std :: scientific << n_g <<" "<<n_e <<" " << n_i <<" " <<n_d<<" " <<gosa <<" "<<F<< " Charge Calculation: Error" << std :: endl;
            exit(1);
        }

        phi=phi_next;
        count=count+1;
    }

    // Thermal (Saha) ionization of potassium, T > 1000K: computed as an
    // independent equilibrium (separate K+/e- channel, own mobility beta_K
    // in CalcResistivity_ForThermal()) on top of the grain-charging
    // Newton solve above, rather than folding it into the coupled
    // ion/electron/grain system -- same additive convention as the n_d<1e-50
    // branch. nK_th is a member, so it must be reset to 0 below 1000K or a
    // stale value from a hotter cell/step would leak into this step.
    if (T > 1000.0)
    {
        CalcThermalIonization();
    }
    else
    {
        nK_th = 0.0;
    }

    //Low tau,high tau の数密度
    int N = 101;

    _dnumber_low_.resize(nd_bin,std:: vector<double>(3));
    _dnumber_high_.resize(nd_bin,std:: vector<double>(N-1));

   //Low tau,high tau の電荷量
    _dust_charge_low_.resize(nd_bin,std:: vector<double>(3)); 
    _dust_charge_high_.resize(nd_bin,std:: vector<double>(N-1)); // N is divided number of Z

    double sum = 0.0;
    double epsilon = exp(phi)/(1.0-phi);
    double tau_i,j_0i,j_1i,omega_i,theta_i;
    double n0;

    double dhcharge;
    double d2Z,Zmax,Zmin;
    double dZ;

    double Z_dhigh[N];
    double dn_dZ[N];

    // exp(-(Z_dhigh[j]-dhcharge)^2/(2*d2Z)) reduces to exp(-x_j^2/2) with
    // x_j = -6+12j/(N-1) *exactly*, independent of d2Z/dhcharge (and hence of
    // dust bin, grid cell, or temperature): Z_dhigh[j] is always dhcharge +/-
    // 6*sqrt(d2Z) sampled at N equally spaced points, so the d2Z cancels in the
    // ratio. This is eq.(23) of Tsukamoto & Okuzumi (2022) sampled at fixed
    // z-scores -- not an approximation -- so the shape only needs computing
    // once instead of every dust bin / grid cell / substep.
    static const std::vector<double> gaussian_kernel = [N]() {
        std::vector<double> k(N);
        for (int j = 0; j < N; j++) {
            double x = -6.0 + 12.0*j/double(N-1);
            k[j] = std::exp(-0.5*x*x);
        }
        return k;
    }();

    for (int i = 0; i < nd_bin ;i++)
    {
        // このbinでのtauの値
        tau_i =d_size[i]*BOLTZMANN_CONSTANT*T / (CHARGE_UNIT*CHARGE_UNIT);
        
        // j_0i = 1+pow(M_PI/double(2.0*tau_i),0.5);
        // j_1i = (1.0+1.0/tau_i)*(1.0+pow(2.0/(tau_i+2.0),0.5));
        j_0i = 1+sqrt(M_PI/double(2.0*tau_i));                      // pow(x,0.5) == sqrt(x)
        j_1i = (1.0+1.0/tau_i)*(1.0+sqrt(2.0/(tau_i+2.0)));         // this loop runs nbin times per grid cell
                                                                     // per substep, so the pow->sqrt swap matters here
        omega_i = j_0i/j_1i;

        theta_i = omega_i/epsilon + 1.0 + omega_i*epsilon;

        n0 = d_n[i]/theta_i;

        // 式(14)
        for(int j=0 ; j<3 ; j++)
        {
            if(j == 0)
            {
                _dnumber_low_[i][j] = omega_i/epsilon * n0;
                _dust_charge_low_[i][j] = -1.0*CHARGE_UNIT;
            }
            else if(j == 1)
            {
                _dnumber_low_[i][j] = n0;
                _dust_charge_low_[i][j] = 0.0;
            }
            else
            {
                _dnumber_low_[i][j] = omega_i*epsilon * n0;
                _dust_charge_low_[i][j] = 1.0*CHARGE_UNIT;
            }
        }

        //high tauの電荷量
        dhcharge = phi*tau_i;
        d2Z = tau_i * (1.0-phi)/(2.0-phi);
        Zmax = dhcharge + 6.0*sqrt(d2Z);
        Zmin = dhcharge - 6.0*sqrt(d2Z);
        dZ = (Zmax-Zmin)/(N-1.0);

        // sqrt(2*pi*d2Z) doesn't depend on j (d2Z is fixed for this dust bin), so
        // it belongs outside the j-loop -- it was being recomputed all N=101
        // times per dust bin per grid cell per substep.
        double norm_high = d_n[i]/sqrt(2.0*M_PI*d2Z);
        for(int j=0 ; j<N ; j++)
        {
            Z_dhigh[j] = Zmin + j*dZ;
            // dn_dZ[j] = d_n[i]/(sqrt(2.0*M_PI*d2Z)) * exp(- (( Z_dhigh[j]-dhcharge)*( Z_dhigh[j]-dhcharge))/(2.0*d2Z));
            dn_dZ[j] = norm_high * gaussian_kernel[j];
        }

        for(int j=0;j<N-1;j++)
        {
            _dnumber_high_[i][j] = 0.5*(dn_dZ[j]+dn_dZ[j+1]) * dZ;
            if (_dnumber_high_[i][j] == 0)
            {
                _dust_charge_high_[i][j] = 0.0;
            }
            else
            {
                _dust_charge_high_[i][j]= 0.5*(dn_dZ[j]*Z_dhigh[j] + dn_dZ[j+1]*Z_dhigh[j+1])*dZ/(_dnumber_high_[i][j]);
            }
        }
    }

}


/*
    Resistivity
*/

void CalcAlphaParameter:: CalcResistivity(Resistivity *res)
{
    // Pinto Galli 2009
    // momentum transfer
    //n_g = n_g*GAS_MOLECULAR_MASS/(1.4*HYDROGEN_MASS);
    double SQRT_1_2PI = 0.398942;
    double n_H2fra    = 0.840336;
    double n_Hefra    = 0.159664;
    double M_4pi       = 12.56637;
    
    double n_H2 = n_H2fra*n_g;
    double n_He = n_Hefra*n_g;

    double mH2 = 2.0 * HYDROGEN_MASS;
    double mHe = 4.0 * HYDROGEN_MASS;
    double mHCO = 29.0*PROTON_MASS;
    // double mK = 39.098*PROTON_MASS;

    double rho_H2 = mH2*n_H2;
    double rho_He = mHe*n_He;

    // reduced mass 
    double mu_eH2 = ELECTRON_MASS * mH2/ (ELECTRON_MASS + mH2);
    double mu_eHe = ELECTRON_MASS * mHe / (ELECTRON_MASS + mHe);

    double mu_H2i = mHCO * mH2/(mHCO + mH2);
    double mu_Hei = mHCO  * mHe/(mHCO + mHe);

    // collision rate Pinto & Galli 2008
    // ELECTRON momentum transfer
    // H2
    //double momentum_eH2=sqrt(T)*(0.535+0.203*theta - 0.163 *theta*theta +0.050 *theta*theta*theta)*1.0e-9 ; // table1
    double momentum_eH2  = 3.16e-11 * std::pow(std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mu_eH2)) * 1.0e-5, 1.3); //　table2
    //He
    //double momentum_eHe=sqrt(T)* 0.428 *1.0e-9; // table1
    double momentum_eHe = 7.08e-11 * std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mu_eHe)) * 1.0e-5; // table2

    // ion momentum transfer(HCO+)
    // H2
    //double momentum_iH2 = 2.81e-9 * std::sqrt(0.804) * std::sqrt(HYDROGEN_MASS / mu_H2i); // table1
    double momentum_iH2 = 2.81e-9 * std::sqrt(0.804) * std::sqrt(HYDROGEN_MASS / mu_H2i); // table2
    // He
    //double momentum_iHe = 2.81e-9 * std::sqrt(0.207) * std::sqrt(HYDROGEN_MASS / mu_Hei); // table1
    double momentum_iHe = 2.81e-9 * std::sqrt(0.207) * std::sqrt(HYDROGEN_MASS / mu_Hei); // table2

    // dust momentum transfer
    int N_d = nd_bin;

    double momentum_dH2[N_d];
    double momentum_dHe[N_d];

    int N_th = _dust_charge_th.size();
    double momentum_dthH2[N_th];
    double momentum_dthHe[N_th];

    // sqrt(8kT/(pi*mH2)) and sqrt(8kT/(pi*mHe)) don't depend on the dust bin i,
    // so hoist them out of the loop (previously recomputed nd_bin times).
    double vth_H2 = sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mH2));
    double vth_He = sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mHe));
    for (int i=0;i<N_d;i++)
    {
        // H2
        momentum_dH2[i] = 4.0/3.0 *1.3*d_s[i] * vth_H2;
        // He
        momentum_dHe[i] = 4.0/3.0 *1.3*d_s[i] * vth_He;
    }

    // gammmaの計算
    // electron
    double gammaH2 = momentum_eH2/(mH2 + ELECTRON_MASS);
    double gammaHe = momentum_eHe/(mHe + ELECTRON_MASS);
    // ion (HCO+)
    double gamma_iH2 = momentum_iH2/(mHCO + mH2);
    double gamma_iHe = momentum_iHe/(mHCO + mHe);

    // dust (熱電離なし)
    double gamma_dH2[N_d];
    double gamma_dHe[N_d];
    for (int i=0;i<N_d;i++)
    {   
        // H2
        gamma_dH2[i] = momentum_dH2[i]/(d_mass[i] + mH2);
        //gamma_dH2[i] = momentum_dH2[i]/(m_d + mH2);
        // He
        gamma_dHe[i] = momentum_dHe[i]/(d_mass[i] + mHe);
        //gamma_dHe[i] = momentum_dHe[i]/(m_d + mHe);
    }

    // tauの計算
    // electron
    double tauH2=rho_H2*gammaH2;
    double tauHe=rho_He*gammaHe;

    //double taue = tauH2+tauHe;
    double taue = tauH2*tauHe/(tauH2+tauHe);

    // ion
    double tauiH2=rho_H2*gamma_iH2;
    double tauiHe=rho_He*gamma_iHe;

    double taui = tauiH2*tauiHe/(tauiH2+tauiHe);

    // dust (熱電離なし)
    double taudH2[N_d];
    double taudHe[N_d];

    for (int i=0;i<N_d;i++)
    {
        taudH2[i]=rho_H2*gamma_dH2[i];
        taudHe[i]=rho_He*gamma_dHe[i];
    }
    // betaの計算
    // Hall parameter electron

    double beta_e = -CHARGE_UNIT*mag / (ELECTRON_MASS*LIGHT_SPEED) * 1.0 /(tauH2+tauHe);

    // Hall parameter ion(HCO+)
    double beta_i = double(CHARGE_UNIT*mag) / (mHCO*LIGHT_SPEED) * 1.0/(tauiH2+tauiHe); 
  
    // Hall parameter dust
    // lowtauとhightauの時で効果を分けているため、分けてResistivityを計算する
    // Hall parameter low_tau
    DoubleMatrix n_dl = _dnumber_low_;
    DoubleMatrix cs_l = _dust_charge_low_;
    
    DoubleMatrix beta_dl(N_d,std :: vector<double>(2));

    for (int i=0;i<N_d;i++)
    {
        beta_dl[i][0]=double(cs_l[i][0]*mag) / (d_mass[i]*LIGHT_SPEED) * 1.0 /(taudH2[i]+taudHe[i]);
    
        beta_dl[i][1]=double(cs_l[i][2]*mag) / (d_mass[i]*LIGHT_SPEED) * 1.0 /(taudH2[i]+taudHe[i]);
    }

    // Hall parameter high_tau
    int N = _dust_charge_high_.at(0).size();
    double z,m;
    DoubleMatrix beta_dh(N_d,std :: vector<double>(N));
    //double s=0.0;
    for(int i=0;i<N_d;i++)
    {
        // Everything here except z is independent of j (m=d_mass[i] and
        // taudH2[i]+taudHe[i] only depend on the dust bin i), so hoist the
        // division out of the inner loop -- it was being redone N times per
        // dust bin (this loop is O(nd_bin*N) per grid cell).
        m = d_mass[i];
        double beta_dh_coeff = double(CHARGE_UNIT*mag) / (m*LIGHT_SPEED*(taudH2[i]+taudHe[i]));
        for(int j=0;j<N;j++)
        {
            z = _dust_charge_high_[i][j];
            beta_dh[i][j] = z*beta_dh_coeff;
        }
    }

    // conductivity
    double n_dlowm,n_dlowp,c_dlowm,c_dlowp;
    double n_dhigh,c_dhigh;
    double n_dth,c_dth;
    // sigma_O
    double sigma_O=0.0;
    // ion electron
    sigma_O = -n_e*CHARGE_UNIT*beta_e + n_i*CHARGE_UNIT*beta_i ;
    // dust
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
        
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_O += n_dlowm*c_dlowm*beta_dl[i][0] + n_dlowp*c_dlowp*beta_dl[i][1];
        
        //hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_O += n_dhigh* c_dhigh*CHARGE_UNIT * beta_dh[i][j];
        }
    }
    sigma_O = LIGHT_SPEED/mag * sigma_O;
    //sigma_P
    double sigma_H=0.0; 
    // ion electron
    sigma_H =n_e*CHARGE_UNIT* beta_e*beta_e/(1.0+beta_e*beta_e) -n_i *CHARGE_UNIT * beta_i*beta_i/(1.0+beta_i*beta_i);
    // dust
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_H -= n_dlowm*c_dlowm * (beta_dl[i][0]*beta_dl[i][0]) / (1.0+beta_dl[i][0]*beta_dl[i][0]);
        sigma_H -= n_dlowp*c_dlowp * (beta_dl[i][1]*beta_dl[i][1]) / (1.0+beta_dl[i][1]*beta_dl[i][1]);

        // hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_H -= n_dhigh*c_dhigh*CHARGE_UNIT * (beta_dh[i][j]*beta_dh[i][j]) / (1.0+beta_dh[i][j]*beta_dh[i][j]);
        }
    }
    sigma_H = LIGHT_SPEED/mag * sigma_H;

    //sigma_P
    double sigma_P=0.0;
    // ion electron
    sigma_P = -n_e*CHARGE_UNIT*beta_e / (1.0+beta_e*beta_e) + CHARGE_UNIT*n_i*beta_i / (1.0+beta_i*beta_i);
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
         
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_P += n_dlowm*c_dlowm * beta_dl[i][0]/ (1.0+beta_dl[i][0]*beta_dl[i][0]);
        sigma_P += n_dlowp*c_dlowp * beta_dl[i][1]/ (1.0+beta_dl[i][1]*beta_dl[i][1]);

        // hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_P += n_dhigh* c_dhigh*CHARGE_UNIT *beta_dh[i][j] / (1.0+beta_dh[i][j]*beta_dh[i][j]);
        }
    }    
    sigma_P = LIGHT_SPEED/mag * sigma_P;    
    
    _sigmaO_ = sigma_O;
    _sigmaH_ = sigma_H;
    _sigmaP_ = sigma_P;

    _etaO_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * 1.0/sigma_O;
    _etaH_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * sigma_H/(sigma_H*sigma_H + sigma_P*sigma_P);
    _etaA_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * sigma_P/(sigma_H*sigma_H + sigma_P*sigma_P) - _etaO_;

    res ->_eta_O_ = _etaO_;
    res ->_eta_A_ = _etaA_;

}


void CalcAlphaParameter:: CalcResistivity_ForThermal(Resistivity *res){
    // Pinto Galli 2009
    // momentum transfer
    //n_g = n_g*GAS_MOLECULAR_MASS/(1.4*HYDROGEN_MASS);
    double SQRT_1_2PI = 0.398942;
    double n_H2fra    = 0.840336;
    double n_Hefra    = 0.159664;
    double M_4pi       = 12.56637;
    
    double n_H2 = n_H2fra*n_g;
    double n_He = n_Hefra*n_g;

    double mH2 = 2.0 * HYDROGEN_MASS;
    double mHe = 4.0 * HYDROGEN_MASS;
    double mHCO = 29.0*PROTON_MASS;
    double mK = 39.098*PROTON_MASS;

    double rho_H2 = mH2*n_H2;
    double rho_He = mHe*n_He;

    // reduced mass 
    double mu_eH2 = ELECTRON_MASS * mH2/ (ELECTRON_MASS + mH2);
    double mu_eHe = ELECTRON_MASS * mHe / (ELECTRON_MASS + mHe);
    double mu_H2K = mK*mH2/(mK+mH2);
    double mu_HeK = mK*mHe/(mK+mHe);
    double mu_H2i = mHCO * mH2/(mHCO + mH2);
    double mu_Hei = mHCO  * mHe/(mHCO + mHe);

    // collision rate Pinto & Galli 2008
    // ELECTRON momentum transfer
    double momentum_eH2  = 3.16e-11 * std::pow(std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mu_eH2)) * 1.0e-5, 1.3); //　table2
    double momentum_eHe = 7.08e-11 * std::sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mu_eHe)) * 1.0e-5; // table2
    
    // ion momentum transfer(HCO+)
    double momentum_iH2 = 2.81e-9 * std::sqrt(0.804) * std::sqrt(HYDROGEN_MASS / mu_H2i); // table2
    double momentum_iHe = 2.81e-9 * std::sqrt(0.207) * std::sqrt(HYDROGEN_MASS / mu_Hei); // table2

    // ion momentum transfer(K+)
    double momentum_KH2 = 2.81e-9 * std::sqrt(0.804) * std::sqrt(HYDROGEN_MASS / mu_H2K); // table2
    double momentum_KHe = 2.81e-9 * std::sqrt(0.207) * std::sqrt(HYDROGEN_MASS / mu_HeK); // table2

    // dust momentum transfer
    int N_d = nd_bin;

    double momentum_dH2[N_d];
    double momentum_dHe[N_d];

    int N_th = _dust_charge_th.size();
    double momentum_dthH2[N_th];
    double momentum_dthHe[N_th];

    // sqrt(8kT/(pi*mH2)) and sqrt(8kT/(pi*mHe)) don't depend on the dust bin i,
    // so hoist them out of the loop (previously recomputed nd_bin times).
    double vth_H2 = sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mH2));
    double vth_He = sqrt(8.0*BOLTZMANN_CONSTANT*T/(M_PI*mHe));
    for (int i=0;i<N_d;i++)
    {
        // H2
        momentum_dH2[i] = 4.0/3.0 *1.3*d_s[i] * vth_H2;
        // He
        momentum_dHe[i] = 4.0/3.0 *1.3*d_s[i] * vth_He;
    }

    // gammmaの計算
    // electron
    double gammaH2 = momentum_eH2/(mH2 + ELECTRON_MASS);
    double gammaHe = momentum_eHe/(mHe + ELECTRON_MASS);
    // ion (HCO+)
    double gamma_iH2 = momentum_iH2/(mHCO + mH2);
    double gamma_iHe = momentum_iHe/(mHCO + mHe);
    // ion (K+)
    double gamma_KH2 = momentum_KH2/(mK + mH2);
    double gamma_KHe = momentum_KHe/(mK + mHe);
    // dust
    double gamma_dH2[N_d];
    double gamma_dHe[N_d];
    for (int i=0;i<N_d;i++)
    {   
        // H2
        gamma_dH2[i] = momentum_dH2[i]/(d_mass[i] + mH2);
        // He
        gamma_dHe[i] = momentum_dHe[i]/(d_mass[i] + mHe);
    }
    
    // tauの計算
    // electron
    double tauH2=rho_H2*gammaH2;
    double tauHe=rho_He*gammaHe;
    // ion  (HCO+)
    double tauiH2=rho_H2*gamma_iH2;
    double tauiHe=rho_He*gamma_iHe;
    // ion  (K+)
    double tauKH2=rho_H2*gamma_KH2;
    double tauKHe=rho_He*gamma_KHe;
    // dust 
    double taudH2[N_d];
    double taudHe[N_d];

    for (int i=0;i<N_d;i++)
    {
        taudH2[i]=rho_H2*gamma_dH2[i];
        taudHe[i]=rho_He*gamma_dHe[i];
    }

    // betaの計算
    // Hall parameter electron
    double beta_e = -CHARGE_UNIT*mag / (ELECTRON_MASS*LIGHT_SPEED) * 1.0 /(tauH2+tauHe);
    // Hall parameter ion(HCO+)
    double beta_i = double(CHARGE_UNIT*mag) / (mHCO*LIGHT_SPEED) * 1.0/(tauiH2+tauiHe); 
    // Hall parameter ion(K+)
    double beta_K = double(CHARGE_UNIT*mag) / (mK*LIGHT_SPEED) * 1.0/(tauKH2+tauKHe); 
    // Hall parameter dust
    // lowtauとhightauの時で効果を分けているため、分けてResistivityを計算する
    // Hall parameter low_tau
    DoubleMatrix n_dl = _dnumber_low_;
    DoubleMatrix cs_l = _dust_charge_low_;
    
    DoubleMatrix beta_dl(N_d,std :: vector<double>(2));

    for (int i=0;i<N_d;i++)
    {
        beta_dl[i][0]=double(cs_l[i][0]*mag) / (d_mass[i]*LIGHT_SPEED) * 1.0 /(taudH2[i]+taudHe[i]);
    
        beta_dl[i][1]=double(cs_l[i][2]*mag) / (d_mass[i]*LIGHT_SPEED) * 1.0 /(taudH2[i]+taudHe[i]);
    }

    // // Hall parameter high_tau
    int N = _dust_charge_high_.at(0).size();
    double z,m;
    DoubleMatrix beta_dh(N_d,std :: vector<double>(N));
    //double s=0.0;
    for(int i=0;i<N_d;i++)
    {
        // Everything here except z is independent of j (m=d_mass[i] and
        // taudH2[i]+taudHe[i] only depend on the dust bin i), so hoist the
        // division out of the inner loop -- it was being redone N times per
        // dust bin (this loop is O(nd_bin*N) per grid cell).
        m = d_mass[i];
        double beta_dh_coeff = double(CHARGE_UNIT*mag) / (m*LIGHT_SPEED*(taudH2[i]+taudHe[i]));
        for(int j=0;j<N;j++)
        {
            z = _dust_charge_high_[i][j];
            beta_dh[i][j] = z*beta_dh_coeff;
        }
    }
    
    // conductivity
    double n_dlowm,n_dlowp,c_dlowm,c_dlowp;
    double n_dhigh,c_dhigh;
    // sigma_O
    double sigma_O=0.0;
    // ion electron
    sigma_O = -n_e*CHARGE_UNIT*beta_e + n_i*CHARGE_UNIT*beta_i + nK_th*CHARGE_UNIT*beta_K;
    // dust
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
        
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_O += n_dlowm*c_dlowm*beta_dl[i][0] + n_dlowp*c_dlowp*beta_dl[i][1];
        
        //hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_O += n_dhigh* c_dhigh*CHARGE_UNIT * beta_dh[i][j];
        }
    }
    sigma_O = LIGHT_SPEED/mag * sigma_O;
    //sigma_P
    double sigma_H=0.0; 
    // ion electron
    sigma_H =n_e*CHARGE_UNIT* beta_e*beta_e/(1.0+beta_e*beta_e) - n_i *CHARGE_UNIT * beta_i*beta_i/(1.0+beta_i*beta_i) - nK_th *CHARGE_UNIT * beta_K*beta_K/(1.0+beta_K*beta_K);
    // // dust
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_H -= n_dlowm*c_dlowm * (beta_dl[i][0]*beta_dl[i][0]) / (1.0+beta_dl[i][0]*beta_dl[i][0]);
        sigma_H -= n_dlowp*c_dlowp * (beta_dl[i][1]*beta_dl[i][1]) / (1.0+beta_dl[i][1]*beta_dl[i][1]);

        // hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_H -= n_dhigh*c_dhigh*CHARGE_UNIT * (beta_dh[i][j]*beta_dh[i][j]) / (1.0+beta_dh[i][j]*beta_dh[i][j]);
        }
    }
    sigma_H = LIGHT_SPEED/mag * sigma_H;

    //sigma_P
    double sigma_P=0.0;
    sigma_P = -n_e*CHARGE_UNIT*beta_e / (1.0+beta_e*beta_e) + CHARGE_UNIT*n_i*beta_i / (1.0+beta_i*beta_i)  + CHARGE_UNIT*nK_th*beta_K / (1.0+beta_K*beta_K);
    for(int i=0;i<N_d;i++)
    {
        // lowtau dust
         
        n_dlowm =_dnumber_low_[i][0];
        n_dlowp =_dnumber_low_[i][2];

        c_dlowm =_dust_charge_low_[i][0];
        c_dlowp =_dust_charge_low_[i][2];

        sigma_P += n_dlowm*c_dlowm * beta_dl[i][0]/ (1.0+beta_dl[i][0]*beta_dl[i][0]);
        sigma_P += n_dlowp*c_dlowp * beta_dl[i][1]/ (1.0+beta_dl[i][1]*beta_dl[i][1]);

        // hightau dust
        for(int j=0 ; j< N ; j++)
        {
            n_dhigh =_dnumber_high_[i][j];
            c_dhigh =_dust_charge_high_[i][j];
            sigma_P += n_dhigh* c_dhigh*CHARGE_UNIT *beta_dh[i][j] / (1.0+beta_dh[i][j]*beta_dh[i][j]);
        }
    }    
    sigma_P = LIGHT_SPEED/mag * sigma_P;    

    _sigmaO_ = sigma_O;
    _sigmaH_ = sigma_H;
    _sigmaP_ = sigma_P;

    _etaO_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * 1.0/sigma_O;
    _etaH_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * sigma_H/(sigma_H*sigma_H + sigma_P*sigma_P);
    _etaA_ = LIGHT_SPEED*LIGHT_SPEED/M_4pi * sigma_P/(sigma_H*sigma_H + sigma_P*sigma_P) - _etaO_;

    res ->_eta_O_ = _etaO_;
    res ->_eta_A_ = _etaA_;
}

/*
    Alpha Parameter
*/

// Delage et.al 2022 (30)(31)
double CalcAlphaParameter :: GetAlphaParameter()
{
    //std :: cout << _etaO_ << " " << _etaH_ << " " << _etaA_ << std :: endl;
    V_A=mag/sqrt(4.0*M_PI*n_g*GAS_MOLECULAR_MASS); // Alfven wave

    double A   = V_A*V_A/(_etaO_*omega_k);          //Elsasser数(オーム散逸)
    double A_H = V_A*V_A/(_etaH_*omega_k);         //Elsasser数(ホール効果) 
    double A_M = V_A*V_A/(_etaA_*omega_k);         //Elsasser数(両極性拡散)
            
    double alpha_1 = 50.0/pow(A_M,1.2); 
    double alpha_2 = 8.0/pow(A_M,0.3)+1;

    double beta_min=sqrt(alpha_1*alpha_1 + alpha_2*alpha_2);//(34)
    //T.N.Delage et.al 2021(39)

    if(A>1 && p_beta>beta_min &&  p_beta>1) 
    {
        alpha=1.0/(3*p_beta);
        if (alpha < alpha_hydro) {
            alpha = alpha_hydro;
        }
    }else
    {
        alpha=alpha_hydro;
    }
    return alpha;
}