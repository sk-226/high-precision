#pragma once

#include <Eigen/Sparse>
#include <chrono>

#include "algorithms/conjugate_gradient.hpp"
#include "bailey/precision_traits.hpp"

namespace algorithms {

/// SSOR-Preconditioned Conjugate Gradient solver with mixed precision
///
/// Solves the linear system Ax = b using the Preconditioned Conjugate Gradient
/// method with SSOR (Symmetric Successive Over-Relaxation) preconditioning.
/// Uses mixed precision: the preconditioner application z = M^{-1} r is
/// computed in double precision, while the rest of the CG algorithm uses high
/// precision T.
///
/// @param A Symmetric positive definite matrix (high precision T)
/// @param b Right-hand side vector (high precision T)
/// @param x Initial guess (modified in-place, high precision T)
/// @param x_true True solution for error analysis (high precision T)
/// @param max_iter Maximum number of iterations
/// @param tolerance Convergence tolerance for relative residual
/// @param omega SSOR relaxation parameter (typically 0 < omega < 2)
/// @return CGResult containing convergence history and statistics
template <typename T>
CGResult<T> ssor_pcg_mixed_precision(
    const typename bailey::PrecisionTraits<T>::matrix_type& A,
    const typename bailey::PrecisionTraits<T>::vector_type& b,
    typename bailey::PrecisionTraits<T>::vector_type& x,
    const typename bailey::PrecisionTraits<T>::vector_type& x_true,
    int max_iter, double tolerance, double omega) {
    using Traits = bailey::PrecisionTraits<T>;
    using MatrixType = typename Traits::matrix_type;
    using VectorType = typename Traits::vector_type;

    // Double precision types for preconditioner matrices
    using DoubleMatrixType = Eigen::SparseMatrix<double>;
    using DoubleVectorType = Eigen::VectorXd;
    using DoubleTriplet = Eigen::Triplet<double>;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Precompute norms for relative error calculations
    T norm2_b = sqrt(b.dot(b));
    T norm2_x_true = sqrt(x_true.dot(x_true));
    T normA_x_true = sqrt(x_true.dot(A * x_true));

    CGResult<T> result;
    result.hist_relres_2.reserve(max_iter + 1);
    result.hist_relerr_2.reserve(max_iter + 1);
    result.hist_relerr_A.reserve(max_iter + 1);

    int n = A.rows();

    // Build double precision SSOR preconditioner matrices
    // DL_d = D + ω*StrictlyLower(A) and DU_d = D + ω*StrictlyUpper(A)
    // Extract diagonal (double precision) element-wise using to_double
    DoubleVectorType diag_d(n);
    for (int i = 0; i < n; ++i) {
        diag_d[i] = to_double(A.coeff(i, i));
    }

    // Build DL_d (lower triangular including diagonal, double precision)
    std::vector<DoubleTriplet> dl_triplets;
    dl_triplets.reserve((A.nonZeros() / 2) + n);  // Rough estimate
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double val = to_double(it.value());
            if (i == j) {
                // Diagonal: keep as is
                dl_triplets.emplace_back(i, j, val);
            } else if (i > j) {
                // Strictly lower: multiply by omega
                dl_triplets.emplace_back(i, j, omega * val);
            }
        }
    }
    DoubleMatrixType DL_d(n, n);
    DL_d.setFromTriplets(dl_triplets.begin(), dl_triplets.end());
    DL_d.makeCompressed();

    // Build DU_d (upper triangular including diagonal, double precision)
    std::vector<DoubleTriplet> du_triplets;
    du_triplets.reserve((A.nonZeros() / 2) + n);  // Rough estimate
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double val = to_double(it.value());
            if (i == j) {
                // Diagonal: keep as is
                du_triplets.emplace_back(i, j, val);
            } else if (i < j) {
                // Strictly upper: multiply by omega
                du_triplets.emplace_back(i, j, omega * val);
            }
        }
    }
    DoubleMatrixType DU_d(n, n);
    DU_d.setFromTriplets(du_triplets.begin(), du_triplets.end());
    DU_d.makeCompressed();

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

    // Main PCG iteration loop
    bool is_converged = false;
    int iter_final = 0;

    // Pre-allocate temporaries used each iteration
    VectorType w(A.rows());
    w.setZero();
    T sigma{};    // (p, A p)
    T alpha{};    // step size
    T rho_new{};  // (r_{k+1}, z_{k+1})
    T beta{};     // direction update factor

    // Double precision temporaries for mixed precision preconditioner
    DoubleVectorType r_d(n);
    DoubleVectorType t_d(n);
    DoubleVectorType u_d(n);
    DoubleVectorType z_d(n);

    // Apply initial SSOR preconditioner: z = M^{-1} r (mixed precision)
    // Cast r to double, compute preconditioner in double, cast z back to high
    // precision
    r_d = r.template cast<double>();
    t_d = DL_d.triangularView<Eigen::Lower>().solve(r_d);
    u_d = diag_d.cwiseProduct(t_d);
    z_d = DU_d.triangularView<Eigen::Upper>().solve(u_d);
    VectorType z = z_d.template cast<T>();

    // Initialize search direction p = z
    VectorType p = z;

    // Initialize rho_old = (r, z)
    T rho_old = r.dot(z);

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Compute matrix-vector product (SpMV)
        w.noalias() = A * p;

        // Compute denominator for step size
        sigma = p.dot(w);

        // Compute step size α = (r,z) / (p,Ap)
        alpha = rho_old / sigma;

        // Update solution: x = x + α*p
        x = x + alpha * p;

        // Update residual: r = r - α*Ap
        r.noalias() -= alpha * w;

        // Apply SSOR preconditioner: z = M^{-1} r (mixed precision)
        // Cast r to double, compute preconditioner in double, cast z back to
        // high precision
        r_d = r.template cast<double>();
        t_d = DL_d.triangularView<Eigen::Lower>().solve(r_d);
        u_d = diag_d.cwiseProduct(t_d);
        z_d = DU_d.triangularView<Eigen::Upper>().solve(u_d);
        z = z_d.template cast<T>();

        // Compute new inner product (r_{k+1}, z_{k+1})
        rho_new = r.dot(z);

        // Compute current error for analysis
        err = x_true - x;

        // Record convergence metrics
        T residual_norm = sqrt(r.dot(r));
        result.hist_relres_2.push_back(to_double(residual_norm / norm2_b));
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

        // Compute β = (r_{k+1},z_{k+1}) / (r_k,z_k)
        beta = rho_new / rho_old;

        // Update for next iteration
        rho_old = rho_new;

        // Update search direction: p = z + β*p
        p = z + beta * p;

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
