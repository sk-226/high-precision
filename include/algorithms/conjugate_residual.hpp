#pragma once

#include "algorithms/conjugate_gradient.hpp"

namespace algorithms {

/// Conjugate Residual solver with comprehensive convergence tracking
///
/// Solves the linear system Ax = b using the Conjugate Residual method.
/// Supports multiple precision levels through template specialization.
///
/// @param A Symmetric positive definite matrix
/// @param b Right-hand side vector
/// @param x Initial guess (modified in-place)
/// @param x_true True solution for error analysis
/// @param max_iter Maximum number of iterations
/// @param tolerance Convergence tolerance for relative residual
/// @return CGResult containing convergence history and statistics
template <typename T>
CGResult<T> conjugateResidual(
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
    VectorType Ar(b.size());
    Ar.noalias() = A * r;

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
    // Initialize q = A * p = A * r (since p = r initially)
    VectorType q = Ar;

    // Initialize rho = (r, A*r) for beta calculation
    T rho_old = r.dot(Ar);

    // Main CR iteration loop
    bool is_converged = false;
    int iter_final = 0;

    // Pre-allocate temporaries used each iteration to avoid repeated
    // allocations
    T denom{};    // (q, q) = (A*p, A*p)
    T alpha{};    // step size
    T rho_new{};  // (r_{k+1}, A*r_{k+1})
    T beta{};     // direction update factor

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Compute denominator for step size: denom = (q, q) = (A*p, A*p)
        denom = q.dot(q);

        // Compute step size α = ρ / (q, q) = (r, A*r) / (A*p, A*p)
        alpha = rho_old / denom;

        // Update solution: x = x + α*p
        x = x + alpha * p;

        // Update residual: r = r - α*q = r - α*A*p
        r.noalias() -= alpha * q;

        // Update Ar: Ar = A * r (one SpMV per iteration)
        Ar.noalias() = A * r;

        // Compute new inner product: ρ_new = (r_{k+1}, A*r_{k+1})
        rho_new = r.dot(Ar);

        // Compute current error for analysis
        err = x_true - x;

        // Record convergence metrics (use true 2-norm of r for relres_2)
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

        // Compute β = ρ_{k+1} / ρ_k = (r_{k+1}, A*r_{k+1}) / (r_k, A*r_k)
        beta = rho_new / rho_old;

        // Update for next iteration
        rho_old = rho_new;

        // Update search direction: p = r + β*p
        p = r + beta * p;

        // Update q: q = A*r + β*q
        q = Ar + beta * q;

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
