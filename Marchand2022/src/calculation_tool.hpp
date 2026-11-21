#ifndef _CALCULATIONTOOL_
#define _CALCULATIONTOOL_

# include <iostream>
# include <iomanip>
# include <string>
# include <cmath>
# include <vector>

namespace CALC_TOOL
{
    void LU_Decomposition(const std :: vector<std :: vector <double> > A_,std :: vector<std :: vector<double> > &L_,std :: vector<std :: vector<double> > &U_); // A = LU
    void LU_Decomposition_Pivot(const std::vector<std::vector<double> > &A_,std::vector<std::vector<double> > &L_,std::vector<std::vector<double> > &U_,std::vector<int> &P_);
    std :: vector<double> Calc_LinerEquation(const std :: vector<double> &F_,const std :: vector<std :: vector<double> > &A_); // Ax = F
    std::vector<double> Calc_LinearEquation_Pivot(const std::vector<double> &F_,const std::vector<std::vector<double> > &A_);
    // std::vector<double> Calc_LinearEquation_Pivot(const std::vector<double> &F_,const std::vector<std::vector<double>> &L_,const std::vector<std::vector<double>> &U_,const std::vector<int> &P_);
}

#endif 