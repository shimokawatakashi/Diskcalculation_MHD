#ifndef ENERGY_HPP
#define ENERGY_HPP

// include header file
#include "physical_constants.hpp"
#include "alpha.hpp"
#include "./magnetic_field1D.hpp"
#include "typedef.hpp"
#include "utils.hpp"

// 熱平衡を考え、温度をイタレーションをして計算する
class Energy_Calc
{
public:

    Energy_Calc();
    ~Energy_Calc();

    void SetGridParameter(int ns,int ne,int n_tot);
    void SetDiskParameter(DoubleVector R_,DoubleVector Sigma_,DoubleVector T_,DoubleVector T_IRR,DoubleVector T_IRR4,DoubleVector Omega_,DoubleVector dSwdt_,DoubleVector Alpha_,DoubleVector Alpha_PZ);
    void SetMagneticField(DoubleVector Bz_,DoubleVector DBDR_,DoubleVector Eta_tot);

    void RosselandOpacity();
    void CalcViscousHeating();
    void CalcWindHeating();
    void CalcJouleheating();
    void CalcTemperature();
    void GetTemperature_Initialization(DoubleVector &T_);
    void GetTemperature(DoubleVector &T_);
    double Calculate_Dt(DoubleVector &vgr_,DoubleVector &r_bod_,DoubleVector &vol_inv);
    void GetTemperature_TimeEvolution(DoubleVector &T_,double dt);
    void SetRadDiffFactor(double f) { rad_diff_fac_ = f; }  // radial radiative diffusion strength (0 = off)
    // wind heating prefactor (2 lambda - 3)/(2 (lambda - 1)) for the magnetic lever arm lambda (lambda = 2 -> 0.5)
    void SetWindLeverArm(double lambda) { wind_heat_fac_ = (2.0*lambda - 3.0)/(2.0*(lambda - 1.0)); }

private:
    void CalculateVerticalVelocity();
    void CalculateVerticalEnergyFlux();
    int nr , nstart, nend;
    double gamma_ad = 1.4; // gamma_ad = 7.0/5.0
    double gamma_ad_1 = gamma_ad-1.0;
    double frac_1gamma_ad = 1.0/gamma_ad_1;
    double frac1_8pi = 1/(8.0*M_PI); // 1/8pi
    double rad_diff_fac_ = 1.0; // radial radiative diffusion strength (input "rad_diff_fac"; 0 = off)
    double wind_heat_fac_ = 0.5; // (2 lambda - 3)/(2 (lambda - 1)) , input "lambda_w" (default lambda = 2 -> 0.5)
    DoubleVector gosa_in,gosa_mid,gosa_out;
    DoubleVector vg_;
    DoubleVector dSigmaw_dt;
    DoubleVector inv_vol_;
    DoubleVector Ez_in,Ez_,Ez_out;
    DoubleVector T_vis4,T_irr4,T_mid4;
    DoubleVector sigmag,r_,rbod_,vz_,vz_out,vz_in;
    DoubleVector rhomid,omega,cs,alpha,alpha_pz;
    DoubleVector rhomid_in,cs_in;
    DoubleVector rhomid_out,cs_out;
    DoubleVector P_int,E_tot,E_mag,dEdt_src;
    DoubleVector Etot_flux;
    DoubleVector v_bod_; // face-centered velocity used for both Etot_flux and the -P*div(v) compression term
    DoubleVector Q_irr,Lambda_rad;
    DoubleVector T_irr,T_vis,T_mid;
    DoubleVector T_inner,T_outer,Tin_vis4,Tout_vis4,Tin_mid4,Tout_mid4;
    DoubleVector B_,dBdr_;
    DoubleVector eta_tot;
    DoubleVector F_joule,Fin_joule,Fout_joule;
    DoubleVector F_wind,Fin_wind,Fout_wind;
    DoubleVector F_rad,Fin_rad,Fout_rad;
    DoubleVector KR_,KR_in,KR_out;
    DoubleVector chi_face_; // frozen face conductivity for the implicit radial radiative-diffusion solve
    DoubleVector tauR,tauP,tauR_in,tauR_out,tauP_in,tauP_out;
};

#endif