#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <Eigen/Sparse>

// ==============================================================================
//  共役勾配法（Conjugate Gradient Method）の実装 - Template版
//  すべての精度型 (double, DD, DQ, QX) に対応
// ==============================================================================

namespace linear_algebra {

struct CGResult {
    int iter_final;
    bool is_converged;
    double time;
    std::vector<double> hist_relres_2;
    double true_relres_2;
    std::vector<double> hist_relerr_2;
    std::vector<double> hist_relerr_A;
};

// Template-based Conjugate Gradient implementation
template<typename T>
inline CGResult conjugateGradient(
    const Eigen::SparseMatrix<T>& A, 
    const Eigen::Vector<T, Eigen::Dynamic>& b, 
    Eigen::Vector<T, Eigen::Dynamic>& x, 
    const Eigen::Vector<T, Eigen::Dynamic>& x_true,
    int max_iter = 1000, 
    double tolerance = 1e-15
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    Eigen::Vector<T, Eigen::Dynamic> r = b - A * x;
    Eigen::Vector<T, Eigen::Dynamic> p = r;
    T rs_old = r.dot(r);
    
    T b_norm = sqrt(b.dot(b));
    T x_true_norm = sqrt(x_true.dot(x_true));
    
    CGResult result;
    result.hist_relres_2.reserve(max_iter + 1);
    result.hist_relerr_2.reserve(max_iter + 1);
    result.hist_relerr_A.reserve(max_iter + 1);
    
    // 初期値の履歴を記録
    T init_residual_norm = sqrt(rs_old);
    result.hist_relres_2.push_back(to_double(init_residual_norm / b_norm));
    
    Eigen::Vector<T, Eigen::Dynamic> error = x - x_true;
    T error_norm = sqrt(error.dot(error));
    result.hist_relerr_2.push_back(to_double(error_norm / x_true_norm));
    
    Eigen::Vector<T, Eigen::Dynamic> A_error = A * error;
    T A_error_norm = sqrt(A_error.dot(A_error));
    result.hist_relerr_A.push_back(to_double(A_error_norm / x_true_norm));
    
    for (int k = 0; k < max_iter; ++k) {
        Eigen::Vector<T, Eigen::Dynamic> Ap = A * p;
        T alpha = rs_old / p.dot(Ap);

        x = x + alpha * p;
        r = r - alpha * Ap;

        T rs_new = r.dot(r);
        T residual_norm = sqrt(rs_new);
        
        // 履歴記録
        result.hist_relres_2.push_back(to_double(residual_norm / b_norm));
        
        error = x - x_true;
        error_norm = sqrt(error.dot(error));
        result.hist_relerr_2.push_back(to_double(error_norm / x_true_norm));
        
        A_error = A * error;
        A_error_norm = sqrt(A_error.dot(A_error));
        result.hist_relerr_A.push_back(to_double(A_error_norm / x_true_norm));
        
        // 収束判定
        if (to_double(residual_norm) < tolerance) {
            result.iter_final = k + 1;
            result.is_converged = true;
            break;
        }

        T beta = rs_new / rs_old;
        p = r + beta * p;
        rs_old = rs_new;
        
        if (k == max_iter - 1) {
            result.iter_final = max_iter;
            result.is_converged = false;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    result.time = duration.count() / 1000.0;
    
    // 真の相対残差を計算
    Eigen::Vector<T, Eigen::Dynamic> true_residual = A * x - b;
    T true_residual_norm = sqrt(true_residual.dot(true_residual));
    result.true_relres_2 = to_double(true_residual_norm / b_norm);
    
    return result;
}

inline void print_num_results(const CGResult& results, const std::string& problem_name = "") {
    std::cout << "========================== " << std::endl;
    std::cout << "Numerical Results. " << std::endl;
    if (!problem_name.empty()) {
        std::cout << "Problem: " << problem_name << " " << std::endl;
    }
    std::cout << "========================== " << std::endl;

    if (results.is_converged) {
        std::cout << "Converged! (iter = " << results.iter_final << ")" << std::endl;
    } else {
        std::cout << "NOT converged. (max_iter = " << results.iter_final << ")" << std::endl;
    }

    std::cout << "# Iter.: " << results.iter_final << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "Time[s]: " << results.time << std::endl;
    std::cout << std::scientific << std::setprecision(2);
    std::cout << "Relres_2norm = " << results.hist_relres_2[results.iter_final] << std::endl;
    std::cout << "True_Relres_2norm = " << results.true_relres_2 << std::endl;
    std::cout << "Relerr_2norm = " << results.hist_relerr_2[results.iter_final] << std::endl;
    std::cout << "Relerr_Anorm = " << results.hist_relerr_A[results.iter_final] << std::endl;
    std::cout << "========================== " << std::endl;
    std::cout << std::endl;
}

} // namespace linear_algebra