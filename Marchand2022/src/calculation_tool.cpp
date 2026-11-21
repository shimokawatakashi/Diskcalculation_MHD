#include "calculation_tool.hpp"

namespace CALC_TOOL
{
    void LU_Decomposition(const std::vector<std::vector<double> > A_, std::vector<std::vector<double> > &L_, std::vector<std::vector<double> > &U_)
    {
        std :: vector<double> P_;
        int n = A_.size();
        L_.resize(n, std::vector<double>(n, 0.0));
        U_.resize(n, std::vector<double>(n, 0.0));
        P_.resize(n,0.0);

        if (n != A_[0].size())
        {
            std::cerr << "We cannot decompose matrix" << std::endl;
            std::cout << "We cannot decompose matrix" << std::endl;
            exit(1);
        }

        for (int i = 0; i < n; i++)
        {
            // Uのi行目を計算
            for (int k = i; k < n; k++)
            {
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sum += L_[i][j] * U_[j][k];
                }
                U_[i][k] = A_[i][k] - sum;
            }
            if (std::abs(U_[i][i]) < 1e-12) {
                std::cerr << "Zero pivot encountered at U[" << i << "][" << i << "]." << std::endl;
                exit(1);
            }
            
            // Lのi列目を計算
            for (int k = i; k < n; k++)
            {
                if (i == k)
                {
                    L_[i][i] = 1.0; // 対角成分は1
                }
                else
                {
                    double sum = 0.0;
                    for (int j = 0; j < i; j++)
                    {
                        sum += L_[k][j] * U_[j][i];
                    }
                    L_[k][i] = (A_[k][i] - sum) / U_[i][i];
                }
            }
        }

        // if (n != A_[0].size()) {
        //     throw std::runtime_error("Matrix is not square.");
        // }

        // U_ = A_; // Copy A_ into U_ initially
        // for (int i = 0; i < n; ++i) P_[i] = i;

        // for (int i = 0; i < n; ++i) {
        //     // Pivot selection
        //     int pivot = i;
        //     double max_val = std::abs(U_[i][i]);
        //     for (int k = i + 1; k < n; ++k) {
        //         if (std::abs(U_[k][i]) > max_val) {
        //             max_val = std::abs(U_[k][i]);
        //             pivot = k;
        //         }
        //     }

        //     if (std::abs(U_[pivot][i]) < 1e-12) {
        //         throw std::runtime_error("Zero pivot encountered during decomposition.");
        //     }

        //     if (pivot != i) {
        //         std::swap(U_[i], U_[pivot]);
        //         std::swap(L_[i], L_[pivot]);
        //         std::swap(P_[i], P_[pivot]);
        //     }

        //     L_[i][i] = 1.0;

        //     for (int j = i + 1; j < n; ++j) {
        //         double factor = U_[j][i] / U_[i][i];
        //         L_[j][i] = factor;
        //         for (int k = i; k < n; ++k) {
        //             U_[j][k] -= factor * U_[i][k];
        //         }
        //     }
        // }
    }

    void LU_Decomposition_Pivot(
    const std::vector<std::vector<double> > &A_,
    std::vector<std::vector<double> > &L_,
    std::vector<std::vector<double> > &U_,
    std::vector<int> &P_
    ) {
        int n = A_.size();
        L_.assign(n, std::vector<double>(n, 0.0));
        U_ = A_;
        P_.resize(n);
        for (int i = 0; i < n; ++i) P_[i] = i;

        for (int i = 0; i < n; ++i) {

            // pivot search
            int pivot = i;
            double max_val = std::abs(U_[i][i]);
            for (int k = i+1; k < n; ++k) {
                if (std::abs(U_[k][i]) > max_val) {
                    max_val = std::abs(U_[k][i]);
                    pivot = k;
                }
            }

            // swap rows in U
            std::swap(U_[i], U_[pivot]);
            std::swap(P_[i], P_[pivot]);

            // swap rows in L "only for columns < i"
            for (int col = 0; col < i; ++col) {
                std::swap(L_[i][col], L_[pivot][col]);
            }

            L_[i][i] = 1.0;

            // elimination
            for (int row = i+1; row < n; ++row) {
                double factor = U_[row][i] / U_[i][i];
                L_[row][i] = factor;
                for (int col = i; col < n; ++col) {
                    U_[row][col] -= factor * U_[i][col];
                }
            }
        }
    }

    std::vector<double> Calc_LinerEquation(const std::vector<double> &F_, const std::vector<std::vector<double> > &A_)
    {
        // J * x = F
        // LU * x = F

        // vector x を求める
        int n = F_.size();
        std ::vector< std :: vector<double> > L_;
        std ::vector< std :: vector<double> > U_;
        std::vector<double> x_;
        std::vector<double> y_;
        L_.resize(n,std :: vector<double> (n,0.0));
        U_.resize(n,std :: vector<double> (n,0.0));
        x_.resize(n, 0.0);
        y_.resize(n, 0.0);

        // AをLU 分解
        LU_Decomposition(A_,L_,U_);

        // Ly = F（前進代入）
        for (int i = 0; i < n; i++)
        {
            y_[i] = F_[i];
            for (int j = 0; j < i; j++)
            {
                y_[i] -= L_[i][j] * y_[j];
            }
            y_[i] /= L_[i][i]; // （今回はLの対角成分が1だから省略できる場合もある）
        }

        // Ux = y（後退代入）
        for (int i = n - 1; i >= 0; i--)
        {
            x_[i] = y_[i];
            for (int j = i + 1; j < n; j++)
            {
                x_[i] -= U_[i][j] * x_[j];
            }
            x_[i] /= U_[i][i];
        }

        return x_; 
    }
    
    std::vector<double> Calc_LinearEquation_Pivot(
    const std::vector<double> &F_,
    const std::vector<std::vector<double> > &A_
    ) {
        int n = F_.size();
        std::vector<std::vector<double> > L_, U_;
        std::vector<int> P_;
        std::vector<double> y_(n), x_(n);

        LU_Decomposition_Pivot(A_, L_, U_, P_);

        // Apply permutation to F
        std::vector<double> F_perm(n);
        for (int i = 0; i < n; ++i) {
            F_perm[i] = F_[P_[i]];
        }

        // Forward substitution: L y = P F
        for (int i = 0; i < n; ++i) {
            y_[i] = F_perm[i];
            for (int j = 0; j < i; ++j) {
                y_[i] -= L_[i][j] * y_[j];
            }
        }

        // Backward substitution: U x = y
        for (int i = n - 1; i >= 0; --i) {
            x_[i] = y_[i];
            for (int j = i + 1; j < n; ++j) {
                x_[i] -= U_[i][j] * x_[j];
            }
            x_[i] /= U_[i][i];
        }

        return x_;
    }
    
}