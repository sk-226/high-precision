#pragma once

#include <chrono>
#include <cmath>

#include "algorithms/conjugate_gradient.hpp"  // reuse CGResult and helpers
#include "bailey/precision_traits.hpp"

namespace algorithms {

// Alternative Conjugate Gradient: alpha = (p,r)/(p,Ap),
// beta = -(r_{k+1},Ap_k)/(p_k,Ap_k).
template <typename T>
CGResult<T> conjugateGradient_alt(
    const typename bailey::PrecisionTraits<T>::matrix_type& A,
    const typename bailey::PrecisionTraits<T>::vector_type& b,
    typename bailey::PrecisionTraits<T>::vector_type& x,
    const typename bailey::PrecisionTraits<T>::vector_type& x_true,
    int max_iter, double tolerance) {
    using Traits = bailey::PrecisionTraits<T>;
    using VectorType = typename Traits::vector_type;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Precompute norms for relative error calculations
    T norm2_b = sqrt(b.dot(b));
    T norm2_x_true = sqrt(x_true.dot(x_true));
    T normA_x_true = sqrt(x_true.dot(A * x_true));

    CGResult<T> result;
    result.hist_relres_2.reserve(max_iter + 1);
    result.hist_relerr_2.reserve(max_iter + 1);
    result.hist_relerr_A.reserve(max_iter + 1);

    // Initialize residual: r = b - Ax
    VectorType r(b.size());
    r.noalias() = b - A * x;

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

    // Initialize search direction p = r
    VectorType p = r;

    bool is_converged = false;
    int iter_final = 0;

    // Pre-allocate temporaries used each iteration to avoid repeated
    // allocations
    VectorType w(A.rows());
    w.setZero();
    T sigma{};  // (p, A p)
    T alpha{};  // step size
    T beta{};   // direction update factor

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Compute matrix-vector product (SpMV)
        w.noalias() = A * p;

        // Compute denominator for step size
        sigma = p.dot(w);

        // Compute step size α = (p,r) / (p,Ap)
        alpha = (p.dot(r)) / sigma;

        // Update solution: x = x + α*p
        x = x + alpha * p;

        // Update residual: r = r - α*Ap
        r.noalias() -= alpha * w;

        // Compute current error for analysis
        err = x_true - x;

        // Record convergence metrics
        result.hist_relres_2.push_back(to_double(sqrt(r.dot(r)) / norm2_b));
        result.hist_relerr_2.push_back(
            to_double(sqrt(err.dot(err)) / norm2_x_true));
        result.hist_relerr_A.push_back(
            to_double(sqrt(err.dot(A * err)) / normA_x_true));

        // Check convergence: ||r||₂ / ||b||₂ < tolerance
        if (result.hist_relres_2.back() < tolerance) {
            is_converged = true;
            iter_final = iter;
            break;
        }

        // Compute β without unary minus to support Bailey types
        beta = T(0.0);
        beta -= (r.dot(w) / sigma);

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

}  // namespace algorithms
