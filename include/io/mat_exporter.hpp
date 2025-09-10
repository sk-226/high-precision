#pragma once

#include "algorithms/conjugate_gradient.hpp"
#include <string>

namespace io {

/// MATLAB .mat file exporter for convergence data
/// 
/// Exports CGResult data to MATLAB-compatible .mat files using matio-cpp.
/// Writes two top-level variables for easy access in MATLAB:
///   - `metadata`: Problem info and run stats (no final metrics)
///   - `convergence`: Iteration-by-iteration histories
/// Final metrics (final residuals/errors) are intentionally omitted here and
/// should be captured in CSV for aggregation workflows.
class MatExporter {
public:
    /// Export convergence data to MATLAB .mat file
    /// 
    /// Creates a MATLAB file with the following variables:
    /// - metadata: Problem information (matrix/precision, converged, iters, time)
    /// - convergence: Iteration histories (hist_iterations, hist_*)
    /// 
    /// @param result CGResult containing convergence data
    /// @param filename Output .mat filename 
    /// @param matrix_name Name of the matrix problem
    /// @param precision_name Precision level identifier
    /// @return true if export successful, false otherwise
    template<typename T>
    static bool export_convergence_data(
        const algorithms::CGResult<T>& result,
        const std::string& filename,
        const std::string& matrix_name,
        const std::string& precision_name
    );

private:
    /// Get precision digits for metadata
    /// @param precision_name Precision level name (double, dd, dq, qx)
    /// @return Number of decimal digits for the precision level
    static int get_precision_digits(const std::string& precision_name);
};

} // namespace io

// Include implementation for template functions
#ifdef ENABLE_MAT_EXPORT
#include "io/mat_exporter_impl.hpp"
#endif
