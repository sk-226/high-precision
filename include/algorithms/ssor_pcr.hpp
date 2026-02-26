#pragma once

#include <Eigen/Sparse>
#include <chrono>

#include "algorithms/conjugate_gradient.hpp"
#include "bailey/precision_traits.hpp"

namespace algorithms {

/// SSOR-Preconditioned Conjugate Residual solver
///
/// Solves the linear system Ax = b using the Preconditioned Conjugate Residual
/// method with SSOR (Symmetric Successive Over-Relaxation) preconditioning.
/// Supports multiple precision levels through template specialization.
///
/// The SSOR preconditioner M is defined as:
///   M = (D + ω*L) * D^{-1} * (D + ω*U) / (ω*(2-ω))
/// where D is diagonal, L is strictly lower, U is strictly upper of A.
///
/// The CR method differs from CG in that it uses:
///   - rho = z' * A * z (instead of r' * z)
///   - alpha = rho / (Ap' * MAp)
///   - z is updated as z = z - alpha * MAp
///
/// @param A Symmetric positive definite matrix
/// @param b Right-hand side vector
/// @param x Initial guess (modified in-place)
/// @param x_true True solution for error analysis
/// @param max_iter Maximum number of iterations
/// @param tolerance Convergence tolerance for relative residual
/// @param omega SSOR relaxation parameter (typically 0 < omega < 2)
/// @return CGResult containing convergence history and statistics
template <typename T>
CGResult<T> ssor_pcr(
    const typename bailey::PrecisionTraits<T>::matrix_type& A,
    const typename bailey::PrecisionTraits<T>::vector_type& b,
    typename bailey::PrecisionTraits<T>::vector_type& x,
    const typename bailey::PrecisionTraits<T>::vector_type& x_true,
    int max_iter, double tolerance, double omega) {
    using Traits = bailey::PrecisionTraits<T>;
    using MatrixType = typename Traits::matrix_type;
    using VectorType = typename Traits::vector_type;
    using Triplet = Eigen::Triplet<T>;

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
    T omega_T = static_cast<T>(omega);

    // Extract diagonal
    VectorType diag = A.diagonal();

    // Build DL = D + ω*StrictlyLower(A) and DU = D + ω*StrictlyUpper(A)
    // Build DL (lower triangular including diagonal)
    std::vector<Triplet> dl_triplets;
    dl_triplets.reserve((A.nonZeros() / 2) + n);
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            T val = it.value();
            if (i == j) {
                // Diagonal: keep as is
                dl_triplets.push_back(Triplet(i, j, val));
            } else if (i > j) {
                // Strictly lower: multiply by omega
                dl_triplets.push_back(Triplet(i, j, omega_T * val));
            }
        }
    }
    MatrixType DL(n, n);
    DL.setFromTriplets(dl_triplets.begin(), dl_triplets.end());
    DL.makeCompressed();

    // Build DU (upper triangular including diagonal)
    std::vector<Triplet> du_triplets;
    du_triplets.reserve((A.nonZeros() / 2) + n);
    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            T val = it.value();
            if (i == j) {
                // Diagonal: keep as is
                du_triplets.push_back(Triplet(i, j, val));
            } else if (i < j) {
                // Strictly upper: multiply by omega
                du_triplets.push_back(Triplet(i, j, omega_T * val));
            }
        }
    }
    MatrixType DU(n, n);
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

    // Apply SSOR preconditioner: z = M^{-1} r
    // M^{-1} = (DU)^{-1} * D * (DL)^{-1}
    // Step 1: t = (DL)^{-1} * r (forward substitution)
    VectorType t = DL.template triangularView<Eigen::Lower>().solve(r);
    // Step 2: u = D * t (element-wise multiplication)
    VectorType u = diag.cwiseProduct(t);
    // Step 3: z = (DU)^{-1} * u (backward substitution)
    VectorType z = DU.template triangularView<Eigen::Upper>().solve(u);

    // Initialize search direction p = z
    VectorType p = z;

    // Ap = A * p
    VectorType Ap(n);
    Ap.noalias() = A * p;

    // Az = A * z
    VectorType Az(n);
    Az.noalias() = A * z;

    // rho_old = z' * A * z (CR uses A-inner product)
    T rho_old = z.dot(Az);

    // MAp = M^{-1} * Ap
    t = DL.template triangularView<Eigen::Lower>().solve(Ap);
    u = diag.cwiseProduct(t);
    VectorType MAp = DU.template triangularView<Eigen::Upper>().solve(u);

    // Main PCR iteration loop
    bool is_converged = false;
    int iter_final = 0;

    // Pre-allocate temporaries
    T alpha{};
    T rho_new{};
    T beta{};

    for (int iter = 1; iter <= max_iter; ++iter) {
        // alpha = rho_old / (Ap' * MAp)
        alpha = rho_old / Ap.dot(MAp);

        // x = x + alpha * p
        x = x + alpha * p;

        // r = r - alpha * Ap
        r.noalias() -= alpha * Ap;

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

        // z = z - alpha * MAp
        z.noalias() -= alpha * MAp;

        // Az = A * z
        Az.noalias() = A * z;

        // rho_new = z' * A * z
        rho_new = z.dot(Az);

        // beta = rho_new / rho_old
        beta = rho_new / rho_old;
        rho_old = rho_new;

        // p = z + beta * p
        p = z + beta * p;

        // Ap = Az + beta * Ap
        Ap = Az + beta * Ap;

        // MAp = M^{-1} * Ap
        t = DL.template triangularView<Eigen::Lower>().solve(Ap);
        u = diag.cwiseProduct(t);
        MAp = DU.template triangularView<Eigen::Upper>().solve(u);

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
