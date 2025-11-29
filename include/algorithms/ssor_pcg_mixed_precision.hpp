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
/// computed in precision P, while the rest of the CG algorithm uses high
/// precision T.
///
/// @tparam T Main computation precision type
/// @tparam P Preconditioner precision type
/// @param A Symmetric positive definite matrix (high precision T)
/// @param b Right-hand side vector (high precision T)
/// @param x Initial guess (modified in-place, high precision T)
/// @param x_true True solution for error analysis (high precision T)
/// @param max_iter Maximum number of iterations
/// @param tolerance Convergence tolerance for relative residual
/// @param omega SSOR relaxation parameter (typically 0 < omega < 2)
/// @return CGResult containing convergence history and statistics
template <typename T, typename P>
CGResult<T> ssor_pcg_mixed_precision(
    const typename bailey::PrecisionTraits<T>::matrix_type& A,
    const typename bailey::PrecisionTraits<T>::vector_type& b,
    typename bailey::PrecisionTraits<T>::vector_type& x,
    const typename bailey::PrecisionTraits<T>::vector_type& x_true,
    int max_iter, double tolerance, double omega) {
    using Traits = bailey::PrecisionTraits<T>;
    using MatrixType = typename Traits::matrix_type;
    using VectorType = typename Traits::vector_type;

    // Preconditioner precision types
    using PrecTraits = bailey::PrecisionTraits<P>;
    using PrecMatrixType = typename PrecTraits::matrix_type;
    using PrecVectorType = typename PrecTraits::vector_type;
    using PrecTriplet = Eigen::Triplet<P>;

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

    // Build precision P SSOR preconditioner matrices
    // DL = D + ω*StrictlyLower(A) and DU = D + ω*StrictlyUpper(A)
    // Convert from precision T to precision P element-wise
    PrecVectorType diag(n);
    for (int i = 0; i < n; ++i) {
        diag[i] = P(to_double(A.coeff(i, i)));
    }

    // Build DL (lower triangular including diagonal, precision P)
    std::vector<PrecTriplet> dl_triplets;
    dl_triplets.reserve((A.nonZeros() / 2) + n);  // Rough estimate
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            P val = P(to_double(it.value()));
            if (i == j) {
                // Diagonal: keep as is
                dl_triplets.emplace_back(i, j, val);
            } else if (i > j) {
                // Strictly lower: multiply by omega
                dl_triplets.emplace_back(i, j, P(omega) * val);
            }
        }
    }
    PrecMatrixType DL(n, n);
    DL.setFromTriplets(dl_triplets.begin(), dl_triplets.end());
    DL.makeCompressed();

    // Build DU (upper triangular including diagonal, precision P)
    std::vector<PrecTriplet> du_triplets;
    du_triplets.reserve((A.nonZeros() / 2) + n);  // Rough estimate
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            P val = P(to_double(it.value()));
            if (i == j) {
                // Diagonal: keep as is
                du_triplets.emplace_back(i, j, val);
            } else if (i < j) {
                // Strictly upper: multiply by omega
                du_triplets.emplace_back(i, j, P(omega) * val);
            }
        }
    }
    PrecMatrixType DU(n, n);
    DU.setFromTriplets(du_triplets.begin(), du_triplets.end());
    DU.makeCompressed();

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

    // Precision P temporaries for mixed precision preconditioner
    PrecVectorType r_p(n);
    PrecVectorType t_p(n);
    PrecVectorType u_p(n);
    PrecVectorType z_p(n);

    // Helper function to convert VectorType (T) to PrecVectorType (P)
    auto convert_to_prec = [](const VectorType& v) -> PrecVectorType {
        PrecVectorType result(v.size());
        for (int i = 0; i < v.size(); ++i) {
            result[i] = P(to_double(v[i]));
        }
        return result;
    };

    // Helper function to convert PrecVectorType (P) to VectorType (T)
    auto convert_from_prec = [](const PrecVectorType& v) -> VectorType {
        VectorType result(v.size());
        for (int i = 0; i < v.size(); ++i) {
            result[i] = T(to_double(v[i]));
        }
        return result;
    };

    // Apply initial SSOR preconditioner: z = M^{-1} r (mixed precision)
    // Convert r from T to P, compute preconditioner in P, convert z back to T
    r_p = convert_to_prec(r);
    t_p = DL.template triangularView<Eigen::Lower>().solve(r_p);
    u_p = diag.cwiseProduct(t_p);
    z_p = DU.template triangularView<Eigen::Upper>().solve(u_p);
    VectorType z = convert_from_prec(z_p);

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
        // Convert r from T to P, compute preconditioner in P, convert z back to T
        r_p = convert_to_prec(r);
        t_p = DL.template triangularView<Eigen::Lower>().solve(r_p);
        u_p = diag.cwiseProduct(t_p);
        z_p = DU.template triangularView<Eigen::Upper>().solve(u_p);
        z = convert_from_prec(z_p);

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
