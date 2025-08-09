#pragma once

// This file is only included when ENABLE_MAT_EXPORT is defined
#ifdef ENABLE_MAT_EXPORT

#include "mat_exporter.hpp"
#include <matioCpp/matioCpp.h>
#include <iostream>

namespace io {

template<typename T>
bool MatExporter::export_convergence_data(
    const algorithms::CGResult<T>& result,
    const std::string& filename,
    const std::string& matrix_name,
    const std::string& precision_name
) {
    try {
        // Create root struct for MATLAB data
        matioCpp::Struct data("data");
        
        // --- Metadata section ---
        matioCpp::Struct metadata("metadata");
        
        // Problem information
        metadata.setField(matioCpp::String("matrix_name", matrix_name));
        metadata.setField(matioCpp::String("precision_name", precision_name));
        
        // Create and set precision digits
        matioCpp::Element<double> precision_digits("precision_digits");
        precision_digits = static_cast<double>(get_precision_digits(precision_name));
        metadata.setField(precision_digits);
        
        // Convergence results
        matioCpp::Element<uint8_t> converged("converged");
        converged = result.converged ? 1 : 0;
        metadata.setField(converged);
        
        matioCpp::Element<double> iterations_performed("iterations_performed");
        iterations_performed = static_cast<double>(result.iterations_performed);
        metadata.setField(iterations_performed);
        
        matioCpp::Element<double> computation_time("computation_time");
        computation_time = result.computation_time;
        metadata.setField(computation_time);
        
        // Final convergence results
        matioCpp::Element<double> final_relres_2norm("final_relres_2norm");
        final_relres_2norm = result.final_residual_norm;
        metadata.setField(final_relres_2norm);
        
        matioCpp::Element<double> final_true_relres_2norm("final_true_relres_2norm");
        final_true_relres_2norm = result.true_relres_2;
        metadata.setField(final_true_relres_2norm);
        
        // Final error metrics (if available)
        if (!result.hist_relerr_2.empty()) {
            matioCpp::Element<double> final_relerr_2norm("final_relerr_2norm");
            final_relerr_2norm = result.hist_relerr_2.back();
            metadata.setField(final_relerr_2norm);
        }
        
        if (!result.hist_relerr_A.empty()) {
            matioCpp::Element<double> final_relerr_Anorm("final_relerr_Anorm");
            final_relerr_Anorm = result.hist_relerr_A.back();
            metadata.setField(final_relerr_Anorm);
        }
        
        data.setField(metadata);
        
        // --- Convergence history section ---
        matioCpp::Struct convergence("convergence");
        
        // Create iteration vector
        std::vector<double> iterations;
        for (size_t i = 0; i < result.hist_relres_2.size(); ++i) {
            iterations.push_back(static_cast<double>(i));
        }
        
        convergence.setField(matioCpp::Vector<double>("hist_iterations", iterations));
        convergence.setField(matioCpp::Vector<double>("hist_relres_2", result.hist_relres_2));
        convergence.setField(matioCpp::Vector<double>("hist_relerr_2", result.hist_relerr_2));
        convergence.setField(matioCpp::Vector<double>("hist_relerr_A", result.hist_relerr_A));
        
        // Add final iteration number
        matioCpp::Element<double> iter_final("iter_final");
        iter_final = static_cast<double>(result.iterations_performed);
        convergence.setField(iter_final);
        
        data.setField(convergence);
        
        // Write to file
        matioCpp::File file = matioCpp::File::Create(filename);
        if (!file.isOpen()) {
            std::cerr << "Error: Could not create file " << filename << std::endl;
            return false;
        }
        
        bool write_success = file.write(data);
        if (!write_success) {
            std::cerr << "Error: Failed to write data to " << filename << std::endl;
            return false;
        }
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error exporting to " << filename << ": " << e.what() << std::endl;
        return false;
    }
}

inline int MatExporter::get_precision_digits(const std::string& precision_name) {
    if (precision_name == "double") return 15;
    if (precision_name == "dd") return 30;
    if (precision_name == "dq") return 66;
    if (precision_name == "qx") return 33;
    return 15; // Default to double precision
}

} // namespace io

#endif // ENABLE_MAT_EXPORT