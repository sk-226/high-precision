#pragma once

#include "algorithms/conjugate_gradient.hpp"
#include <string>

namespace io {

/// MATLAB .mat file exporter for convergence data
/// 
/// Exports CGResult data to MATLAB-compatible .mat files using matio-cpp.
/// Creates structured data with metadata and convergence history for easy
/// analysis and visualization in MATLAB.
class MatExporter {
public:
    /// Export convergence data to MATLAB .mat file
    /// 
    /// Creates a structured MATLAB file with the following hierarchy:
    /// - data.metadata: Problem information and final results
    /// - data.convergence: Iteration-by-iteration convergence history
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