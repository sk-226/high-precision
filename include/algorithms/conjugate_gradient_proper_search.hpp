#pragma once

#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include "algorithms/conjugate_gradient.hpp"
#include "bailey/precision_traits.hpp"

namespace algorithms::proper_search {

/// Conjugate Gradient solver variant that records lagged orthogonality
/// diagnostics. The lag value can be adjusted in this implementation for
/// experiments without changing CLI options.
template <typename T>
CGResult<T> conjugateGradient(
    const typename bailey::PrecisionTraits<T>::matrix_type& A,
    const typename bailey::PrecisionTraits<T>::vector_type& b,
    typename bailey::PrecisionTraits<T>::vector_type& x,
    const typename bailey::PrecisionTraits<T>::vector_type& x_true,
    int max_iter, double tolerance) {
    using Traits = bailey::PrecisionTraits<T>;
    using VectorType = typename Traits::vector_type;

    int lag = 1;  // ←ここを書き換えて lag を変更可能

    auto start_time = std::chrono::high_resolution_clock::now();

    // Precompute norms for relative error calculations
    T norm2_b = sqrt(b.dot(b));
    T norm2_x_true = sqrt(x_true.dot(x_true));
    T normA_x_true = sqrt(x_true.dot(A * x_true));

    CGResult<T> result;
    result.proper_search_lag = lag;
    result.hist_relres_2.reserve(max_iter + 1);
    result.hist_relerr_2.reserve(max_iter + 1);
    result.hist_relerr_A.reserve(max_iter + 1);
    result.hist_res_orthogonality.reserve(max_iter + 1);
    result.hist_search_direction_A_orthogonality.reserve(max_iter + 1);

    // Initialize residual: r = b - Ax
    VectorType r(b.size());
    r.noalias() = b - A * x;

    // buffer for the residual and A*p (size: lag+1)
    std::vector<VectorType> r_buf(lag + 1);
    std::vector<VectorType> w_buf(lag + 1);
    for (int i = 0; i <= lag; ++i) {
        r_buf[i].setZero(b.rows());
        w_buf[i].setZero(b.rows());
    }

    // Calculate initial error vector
    VectorType err = x_true - x;

    // Store initial convergence metrics
    T initial_residual_norm = sqrt(r.dot(r));
    result.initial_residual_norm = to_double(initial_residual_norm);
    result.hist_relres_2.push_back(to_double(initial_residual_norm / norm2_b));
    result.hist_relerr_2.push_back(
        to_double(sqrt(err.dot(err)) / norm2_x_true));
    result.hist_relerr_A.push_back(
        to_double(sqrt(err.dot(A * err)) / normA_x_true));

    result.hist_res_orthogonality.push_back(std::nan("NaN"));
    result.hist_search_direction_A_orthogonality.push_back(std::nan("NaN"));

    // Initialize search direction p = r
    VectorType p = r;

    // Initialize rho for beta calculation
    T rho_old = r.dot(r);

    // Main CG iteration loop
    bool is_converged = false;
    int iter_final = 0;

    // Pre-allocate temporaries used each iteration to avoid repeated
    // allocations
    VectorType w(A.rows());
    w.setZero();
    T sigma{};    // (p, A p)
    T alpha{};    // step size
    T rho_new{};  // (r_{k+1}, r_{k+1})
    T beta{};     // direction update factor

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Compute matrix-vector product (SpMV)
        w.noalias() = A * p;

        // Compute denominator for step size
        sigma = p.dot(w);

        // Compute step size α = (r,r) / (p,Ap)
        alpha = rho_old / sigma;

        // Update solution: x = x + α*p
        x = x + alpha * p;

        // Update residual: r = r - α*Ap
        r.noalias() -= alpha * w;

        // Compute new inner product immediately and reuse below
        rho_new = r.dot(r);

        // Compute current error for analysis
        err = x_true - x;

        // Record convergence metrics (reuse rho_new to avoid extra dot)
        result.hist_relres_2.push_back(to_double(sqrt(rho_new) / norm2_b));
        result.hist_relerr_2.push_back(
            to_double(sqrt(err.dot(err)) / norm2_x_true));
        result.hist_relerr_A.push_back(
            to_double(sqrt(err.dot(A * err)) / normA_x_true));

        // Update buffer: shift and push current r and w to front
        for (int j = lag; j >= 1; --j) {
            r_buf[j] = r_buf[j - 1];
            w_buf[j] = w_buf[j - 1];
        }
        r_buf[0] = r;  // current k
        w_buf[0] = w;  // current A*p_k

        // Compute orthogonality metrics
        double res_orth = std::numeric_limits<double>::quiet_NaN();
        double dirA_orth = std::numeric_limits<double>::quiet_NaN();
        if (iter >= lag) {
            const auto& r_l = r_buf[lag];  // r_{k-lag}
            const auto& w_l = w_buf[lag];  // A*p_{k-lag}

            // Residual orthogonality: |r_k^T r_{k-lag}| / (||r_k||
            // ||r_{k-lag}||)
            const double den_r = std::sqrt(to_double(r.dot(r))) *
                                 std::sqrt(to_double(r_l.dot(r_l)));
            if (den_r > 0.0) {
                const double num_r = std::abs(to_double(r.dot(r_l)));
                res_orth = num_r / den_r;
            }

            // Search direction A-orthogonality: |p_k^T A p_{k-lag}| / (||p_k||
            // ||A p_{k-lag}||)
            const double den_pAw = std::sqrt(to_double(p.dot(p))) *
                                   std::sqrt(to_double(w_l.dot(w_l)));
            if (den_pAw > 0.0) {
                const double num_pAw = std::abs(to_double(p.dot(w_l)));
                dirA_orth = num_pAw / den_pAw;
            }
        }
        result.hist_res_orthogonality.push_back(res_orth);
        result.hist_search_direction_A_orthogonality.push_back(dirA_orth);

        // Check convergence: ||r||₂ / ||b||₂ < tolerance
        if (result.hist_relres_2.back() < tolerance) {
            is_converged = true;
            iter_final = iter;
            break;
        }

        // Compute β = (r_{k+1},r_{k+1}) / (r_k,r_k)
        beta = rho_new / rho_old;

        // Update for next iteration
        rho_old = rho_new;

        // Update search direction: p = r + β*p
        p = r + beta * p;

        iter_final = iter;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);
    result.computation_time = duration.count() / 1000.0;

    // Finalize results
    result.iterations_performed = iter_final;
    result.converged = is_converged;
    result.final_residual_norm = result.hist_relres_2.back();

    // Compute true residual to check for gap with computed residual
    VectorType true_residual(b.rows());
    true_residual.noalias() = b - A * x;
    T true_residual_norm = sqrt(true_residual.dot(true_residual));
    result.true_relres_2 = to_double(true_residual_norm / norm2_b);

    return result;
}

}  // namespace algorithms::proper_search
