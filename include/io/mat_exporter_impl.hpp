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
        // --- Metadata section (top-level variable) ---
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

        matioCpp::Element<double> initial_residual_norm("initial_residual_norm");
        initial_residual_norm = result.initial_residual_norm;
        metadata.setField(initial_residual_norm);

        matioCpp::Element<double> final_residual_norm("final_residual_norm");
        final_residual_norm = result.final_residual_norm;
        metadata.setField(final_residual_norm);

        matioCpp::Element<double> true_relres_2("true_relres_2");
        true_relres_2 = result.true_relres_2;
        metadata.setField(true_relres_2);
        
        // Note: Final metrics are intentionally omitted from metadata.
        // They are exported to CSV separately and can be derived from histories.
        
        // --- Convergence history section (top-level variable) ---
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
        
        // Write to file (two top-level variables)
        matioCpp::File file = matioCpp::File::Create(filename);
        if (!file.isOpen()) {
            std::cerr << "Error: Could not create file " << filename << std::endl;
            return false;
        }
        bool write_ok = true;
        write_ok = write_ok && file.write(metadata);
        write_ok = write_ok && file.write(convergence);

        if (result.proper_search_lag > 0 ||
            !result.hist_res_orthogonality.empty() ||
            !result.hist_search_direction_A_orthogonality.empty()) {
            matioCpp::Struct proper_search("proper_search");

            matioCpp::Element<double> lag("lag");
            lag = static_cast<double>(result.proper_search_lag);
            proper_search.setField(lag);

            if (!result.hist_res_orthogonality.empty()) {
                proper_search.setField(matioCpp::Vector<double>(
                    "hist_res_orthogonality", result.hist_res_orthogonality));
            }

            if (!result.hist_search_direction_A_orthogonality.empty()) {
                proper_search.setField(matioCpp::Vector<double>(
                    "hist_search_direction_A_orthogonality",
                    result.hist_search_direction_A_orthogonality));
            }

            write_ok = write_ok && file.write(proper_search);
        }
        if (!write_ok) {
            std::cerr << "Error: Failed to write variables to " << filename << std::endl;
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
