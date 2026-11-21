#ifndef _CALC_RESISTIVITY_HPP
#define _CALC_RESISTIVITY_HPP

#include "typedef.hpp"
#include "physical_constants.hpp"
#include "calculation_tool.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <set>
#include <sstream>

// Marchand 2022
// Newton法を用いて電子、イオンの存在量を計算してResistivityを計算する
// 変数としてはni,ne,phi,epsilonを計算する

class Calc_Resistivity
{
public:
    Calc_Resistivity();
    ~Calc_Resistivity();

    // ng_: gas particle number density n_g = n(H2) + n(He) [cm^-3]; converted internally to
    // n_H = NH_PER_NG * n_g (n(He) = 0.085 n_H), which the Marchand 2022 equations use (zeta n_H etc.)
    void SetGasParameter_For_Resistivity(double ionizationrate_,double temp,double ng_);
    static constexpr double NH_PER_NG = 1.0/(0.5 + 0.085);
    void SetDustParameter_For_Resistivity(DoubleVector nd_,DoubleVector m_,DoubleVector a_);
    void SetMagneticField_For_Resistivity(double magneticfield);

    // Grid location of the current call, used to annotate error/warning messages
    // from the Newton solver (i, k = grid indices; r_au, z = position; T, ng = gas state)
    void SetDebugContext(int i, int k, double r_au, double z, double T, double ng);

    void Calc_ChargeState();
    void Calc_MagneticResistivity();

    double GetEtaO()   const { return _etaO_; }
    double GetEtaA()   const { return _etaA_; }
    double GetEtaH()   const { return _etaH_; }
    double GetSigmaO() const { return _sigmaO_; }
    double GetSigmaH() const { return _sigmaH_; }
    double GetSigmaP() const { return _sigmaP_; }

    double GetNi()     const { return n_i; }
    double GetNe()     const { return n_e; }
    double GetNH_ion() const { return n_H; }
    double GetNNa()    const { return n_Na; }
    double GetNK()     const { return n_K; }
    double GetNg()     const { return n_gas; }
    double GetNH()     const { return n_g; }
    double GetZd(int j) const { return (j >= 0 && j < nd_bin) ? Z_d[j] : 0.0; }

    void SetActiveSpecies(bool doH, bool doNa, bool doK);

private:
    void Calculate_Initial_Phi0();
    void Calculate_Newton(double gosa,DoubleVector &Zd_);
    void Caluclate_F();
    void J1_Component();
    void J2_Component();
    void J3_Component();
    void J4_Component();

    double T; //temperature
    double a; //Dust mean radius

    // Jacobian
    DoubleMatrix J;
    DoubleVector F;
    DoubleVector X;

    // dust information
    DoubleVector d_n;       // dust numberdensity vector
    DoubleVector d_size;    // dust size vector
    DoubleVector d_s;       // dust cross section vector
    DoubleVector d_mass;    // dust mass vector

    // bin
    int nd_bin;     // dust number density bin

    // numberdensity
    double n_gas; // gas particle number density n_g (input)
    double n_g;   // H nuclei number density n_H = NH_PER_NG * n_gas (n_H of Marchand 2022)
    double n_d; // dust number density
    double n_H; // Hydrogen ion number density(H+)
    double n_Na;// Sodium ion number density(Na+)
    double n_K;// Pottasium ion number density(K+)
    double n_i; // ion number density
    double n_e; // electron density

    double n_s;
    double n_s0;

    double ni_ns_tot; // Total number density of ions and neutral species

    double mui_;
    double mus_;

    double mag;

    // Marchand et al. 2022 (Table1)
    //////////////////////////////////////////////////////
    double m_i = 25.0*PROTON_MASS; // ion i mass
    double m_Na = 22.99*PROTON_MASS; // Sodium mass
    double m_K = 39.09*PROTON_MASS; // Potasium mass
    double m_e = ELECTRON_MASS; // electron mass
    double m_g = GAS_MOLECULAR_MASS;  // Gas mean mass
    double m_d ; // dust mean mass

    double Na_H_rate = 3.1e-9; // n_Na/n_H (Table 1)
    double K_H_rate = 2.2e-10; // n_K/n_H
    
    double beta_Na = 1.4e-15;
    double beta_K = 6.5e-15;
    double beta_H = 2.0e-10;
    
    double T0_Na = 6.0e4;
    double T0_K = 5.1e4;
    double T0_H = 1.58e5;

    double k_i_Na = 6.2e-9; // chemical reaction rate ion to Na
    double k_Na_i = 0.0; // chemical reaction rate Na to i
    double k_i_K = 6.2e-9; // chemical reaction rate ion to K
    double k_K_i = 0.0; // chemical reaction rate K to i
    double k_i_H = 3.7e-14; // chemical reaction rate ion to H
    double k_H_i = 3.8e-9; // chemical reaction rate H to ion

    double k_i_s ; // chemical reaction rate ion to species
    double k_s_i ; // chemical reaction rate species to ion

    double sigma_v_ie    = 2.00e-7;
    double sigma_v_se_Na = 2.78e-12; // 
    double sigma_v_se_K = 2.78e-12;
    double sigma_v_se_H = 3.50e-12;
    double  sigma_v_se;

    ////////////////////////////////////////////////////

    double zeta; // ionization rate 

    double tau; // tau
    double tau_tot; // tau
    DoubleVector tau_k; // tau

    double s_i = 1.0; // collision parameter (ion)
    double s_e = 0.3; // collision parameter (electron), Umebayashi & Nakano 1990
    double s_d ; // dust cross section

    double u_i; //theral velocity (ion)
    double u_e; //thermal velocity(electron)
    
    double u_H; //theral velocity (Hydrogen)
    double u_Na; //theral velocity (Sodium)
    double u_K; //theral velocity (Pottasium)
    double u_s; //theral velocity (species)

    double f_p; // The average charge of grains
    double f_0; // The average charge of grains
    double f_m; // The average charge of grains

    DoubleVector Z_d; // // The  charge of grains
    double Z_k; // The average charge of grains
    double Z_k_tot; // total charge of grains

    double dZdEpsilon;
    double total_dZdEpsilon;

    double j_i; // cross setion (ion)
    double j_e; // cross setion (electron)
    double j_0; // neutral cross section
    double j_k; // cross setion (dust)
    double total_jk; // total cross section (dust)

    double dJdEpsilon;
    double total_dJdEpsilon;

    // charge data
    double q_e = -1.0*CHARGE_UNIT; // electron  charge
    double q_i = 1.0*CHARGE_UNIT; // ion mean charge
    double q_d; // dust mean charge
    double qis;

    // Calculation parameter
    double phi;
    double phi0 = 1.0e-2;
    double expo_phi;

    // warm start from the previous call (n_i, n_s as fractions of n_g, n_s0)
    double n_i_warmstart  = 0.0;
    double n_i_H_ws_frac  = 0.0;   // per-species HCO+ warmstart (H solve), as n_i/n_g
    double n_i_Na_ws_frac = 0.0;   // per-species HCO+ warmstart (Na solve), as n_i/n_g
    double n_i_K_ws_frac  = 0.0;   // per-species HCO+ warmstart (K solve), as n_i/n_g
    double n_H_ws_frac  = 0.0;     // n_H/n_s0(H) = n_H/n_g
    double n_Na_ws_frac = 0.0;     // n_Na/n_s0(Na)
    double n_K_ws_frac  = 0.0;     // n_K/n_s0(K)
    double epsilon_ws = 0.9999;
    double phi_ws = 0.0;
    bool   has_phi_ws = false;

    // Active species flags: which thermally ionized species to compute
    bool active_H  = true;
    bool active_Na = true;
    bool active_K  = true;

    double alpha_k;
    DoubleVector alphak_; // tau

    double epsilon;
    double delta_epsilon;
    double epsilon0;

    double Theta;
    double dndt_thermal;

    double total_cross_section;

    // Debug context (grid location) for error/warning messages
    int dbg_i_ = -1;
    int dbg_k_ = -1;
    double dbg_r_au_ = 0.0;
    double dbg_z_ = 0.0;
    double dbg_T_ = 0.0;
    double dbg_ng_ = 0.0;

    // Newton diagnostics: last kDiagN iterations, printed on failure
    struct NewtonDiagRec {
        int it;
        double F[4];
        double Xraw[4];   // Newton step before the max_phi/log step clamp
        double X[4];      // after clamp
        double gosa;
        int gosa_arg;     // which component set gosa
        double alpha;
        int ns_clamp;     // 1 if n_s was pinned to n_s0 this step
        double phi, eps, ni_ng, ns_ns0;
    };
    static const int kDiagN = 30;
    NewtonDiagRec diag_buf_[kDiagN];
    int diag_n_ = 0;
    int dbg_species_ = -1;
    double diag_init_[4] = {0.0, 0.0, 0.0, 0.0}; // phi, eps, ni/ng, ns/ns0 at entry
    int prev_dbg_i_ = -1;     // previous grid point (warm-start source)
    int prev_dbg_k_ = -1;
    bool diag_retry_ = false;
    void DumpNewtonDiag(const char *reason, DoubleVector &Zd_);

protected:
    // Resistivity 

    double _sigmaO_;
    double _sigmaP_;
    double _sigmaH_;

    double _etaO_;
    double _etaA_;
    double _etaH_;

};


#endif
