#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cmath>

namespace io {

class CsvExporter {
public:
    // Append one row to CSV (create file and header if needed)
    static bool append_row(
        const std::string& filename,
        const std::string& matrix,
        int n,
        int nnz,
        const std::string& solver,
        double tolerance,
        int max_iterations,
        const std::string& rhs,
        const std::string& prec,
        const std::string& prec_kv,
        const std::string& prec_label,
        double construction_time,
        bool converged,
        int iters,
        double true_relres_2,
        double final_relres_2,
        double final_relerr_2,
        double final_relerr_A,
        double solve_time,
        const std::string& solve_status
    ) {
        try {
            bool need_header = !file_exists_and_nonempty(filename);

            std::ofstream ofs(filename, std::ios::app);
            if (!ofs.is_open()) return false;

            if (need_header) {
                ofs << header_line() << "\n";
            }

            // fixed ordering, matching example.csv
            std::ostringstream row;
            row << escape_csv(matrix) << ','
                << n << ','
                << nnz << ','
                << escape_csv(solver) << ','
                << format_double(tolerance) << ','
                << max_iterations << ','
                << escape_csv(rhs) << ','
                << escape_csv(prec) << ','
                << escape_csv(prec_kv) << ','
                << escape_csv(prec_label) << ','
                << format_double(construction_time) << ','
                << (converged ? 1 : 0) << ','
                << iters << ','
                << format_double(true_relres_2) << ','
                << format_double(final_relres_2) << ','
                << format_double(final_relerr_2) << ','
                << format_double(final_relerr_A) << ','
                << format_double(solve_time) << ','
                << escape_csv(solve_status);

            ofs << row.str() << "\n";
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    static std::string header_line() {
        return "matrix,n,nnz,solver,tolerance,max_iterations,rhs,prec,prec_kv,prec_label,construction_time,converged,iters,true_relres_2,final_relres_2,final_relerr_2,final_relerr_A,solve_time,solve_status";
    }

    static bool file_exists_and_nonempty(const std::string& path) {
        std::error_code ec;
        auto st = std::filesystem::status(path, ec);
        if (ec || !std::filesystem::exists(st)) return false;
        return std::filesystem::file_size(path, ec) > 0 && !ec;
    }

    static std::string escape_csv(const std::string& s) {
        bool need_quotes = s.find_first_of(",\"\n") != std::string::npos;
        if (!need_quotes) return s;
        std::string out = "\"";
        for (char c : s) {
            if (c == '"') out += '"';
            out += c;
        }
        out += '"';
        return out;
    }

    static std::string format_double(double v) {
        if (!std::isfinite(v)) return "NaN";
        std::ostringstream ss;
        ss.setf(std::ios::scientific);
        ss << v;
        return ss.str();
    }
};

} // namespace io
