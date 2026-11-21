#ifndef ALPHA_HPP_
#define ALPHA_HPP_

#include "typedef.hpp"
#include "physical_constants.hpp"

class CalcAlphaParameter
{
public:

    CalcAlphaParameter();
    ~CalcAlphaParameter();

    void SetParameter(double ionizationrate_,double temp,double ng_, DoubleVector nd_,DoubleVector m_,DoubleVector a_);
    void SetParameterForAlpha(double plasmabeta,double magneticfield,double Omega_Kepler);
    /*
        Chemical Reaction 
    */

    // Tsukamoto et.al 2022 Section 2
    // Newton法を用いて電気的中性条件を解く
    // Calculation用のParameterの設定

    double Omega_cross();

    double Ion_cross_section();
    double Electron_cross_section();

    double total_dustcharge();
    
    double total_ion_numberdensity();
    double total_electron_numberdensity();

    double diff_cross_section_i();
    double diff_cross_section_e();

    double diff_function_f();

    double diff_total_ion_numberdensity();
    double diff_total_electron_numberdensity();
    double diff_total_dustcharge();

    void CalcThermalIonization();

    void GetChargestate(double phi_init = 1.0e-2);
    double GetPhi() const { return phi; }

    /*
        Calculate Resistivity
    */

    void CalcResistivity(Resistivity *res);
    void CalcResistivity_ForThermal(Resistivity *res);

    /*
        Calculate Alpha
    */

    double GetAlphaParameter();

private:

    // Charge State

    double T; //temperature
    double a; //Dust mean radius

    // dust information
    DoubleVector d_n;       // dust numberdensity vector
    DoubleVector d_size;    // dust size vector
    DoubleVector d_s;       // dust cross section vector
    DoubleVector d_mass;    // dust mass vector

    // bin
    int nd_bin;     // dust number density bin

    // numberdensity
    double n_g; // gas numberdensity
    double n_d; // dust numberdensity
    double n_K;// Pottasium numberdensity
    double n_i; // ion numberdensity
    double n_e; // electron density

    double p_beta;
    double mag;
    double omega_k;

    //熱電離を考慮した計算をするとき
    double nK_th;
    double ne_th;
    double ni_th;

    double m_i = 29.0*PROTON_MASS; // ion mass
    double m_K = 39.09*HYDROGEN_MASS;
    double m_e = ELECTRON_MASS; // electron mass
    double m_g = GAS_MOLECULAR_MASS;  // Gas mean mass
    double m_d ; // dust mean mass

    double zeta; // ionization rate 
    double zeta_th; // thermal ionizationrate

    double tau; // tau
    double beta; // recombination rate

    double s_i = 1.0; // collision parameter (ion)
    double s_e = 0.3; // collision parameter (electron)
    double s_d ; // dust cross section

    double u_i; //theral velocity (ion)
    double u_e; //thermal velocity(electron)

    double j_i; // cross setion (ion)
    double j_e; // cross setion (electron)
    double j_0; // neutral cross section

    // charge data
    double q_e = -1.0*CHARGE_UNIT; // electron  charge
    double q_i = 1.0*CHARGE_UNIT; // ion mean charge
    double q_d; // dust mean charge

    double qd_th; // 熱電離の時のdust charge

    // Calculation parameter
    double omega; // calculation parameter2
    double phi = 1.0e-2; // calculation parameter1
    double C1_i ;
    double C1_e; 
    double C2;
    double F;
    double func_theta;

    double diff_theta;
    double diff_f;
    double diff_j_i;
    double diff_j_e;

    double phi_th = 1.0e-2;

    // exp(phi) and exp(2*phi), refreshed once per Newton iteration in
    // GetChargestate() and reused by every helper below instead of each
    // calling exp() again with the same argument (exp(2*phi) is also just
    // exp(phi)*exp(phi), avoiding a second exp() call entirely).
    double exp_phi_ = 1.0;
    double exp_2phi_ = 1.0;

    // Resistivity
    DoubleMatrix _dnumber_low_;
    DoubleMatrix _dnumber_high_;
    DoubleMatrix _dnumber_th_;

    DoubleMatrix _dust_charge_low_;
    DoubleMatrix _dust_charge_high_;
    DoubleMatrix _dust_charge_th;

protected:
    // Resistivity 

    double _sigmaO_;
    double _sigmaP_;
    double _sigmaH_;

    double _etaO_;
    double _etaA_;
    double _etaH_;

    // Alpha Parameter
    double V_A;
    double alpha;
    double alpha_hydro = 1.0e-4; // hydrodynamic alpha
};



#endif /*ALPHA_HPP_*/