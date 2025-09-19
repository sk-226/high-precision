#include "bailey/precision_traits.hpp"
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"
#include "algorithms/conjugate_gradient.hpp"
#include "io/matrix_market_reader.hpp"
#ifdef ENABLE_MAT_EXPORT
#include "io/mat_exporter.hpp"
#endif
#include "io/csv_exporter.hpp"

#include <charconv>
#include <system_error>
#include <stdexcept>
#include <iostream>
#include <string>
#include <string_view>
#include <variant>
#include <iomanip>
#include <filesystem>
#include <limits>

// Command line configuration
struct SolverConfig {
    std::string matrix_name;
    std::string precision_level{"qx"};  // dd, dq, qx
    double tolerance{1.0e-12};
    std::variant<int, double> max_iter{2.0};  // Default: 2*n
    std::string input_dir{"/work/inputs"};
    std::string export_mat_file;  // Empty if not specified
    std::string export_csv_file;  // Empty if not specified
};

// Command line parser
SolverConfig parseCommandLine(int argc, char** argv) {
    SolverConfig config;

    // retireve the token value next to the --option flag
    const auto next_token = [&argc, &argv](int& i) -> std::string_view {
        if (++i >= argc) {
            // i is incremented before the check, so we need to decrement it
            throw std::runtime_error(
                "Missing value after option: " + std::string(argv[i-1]) + " (try --help for usage)"
            );
        }
        std::string_view token{argv[i]}; // get the token value
        if (!token.empty() && token.front() == '-') {
            throw std::runtime_error(
                "Option " + std::string(argv[i-1]) + 
                " requires a value, but got another flag: " + std::string(token)
            );
        }
        return token;
    };
    
    const auto parse_double_strict = [](std::string_view token) -> double {
        double out{};
        const char* first = token.data();
        const char* last = first + token.size();
        auto [ptr, error_code] = std::from_chars(first, last, out);
        // check if the value is valid
        if (error_code == std::errc{} && ptr == last) {
            return out;
        }
        throw std::runtime_error(
            "Invalid floating point value: " + std::string(token)
        );
    };

    // if the token is integer, return the integer
    // if not, return the double (goto the parse_double_strict function)
    const auto parse_max_iter = [&parse_double_strict](std::string_view token) -> std::variant<int, double> {
        int fixed_max_iter{};
        const char* first = token.data();
        const char* last = first + token.size();
        auto [ptr, error_code] = std::from_chars(first, last, fixed_max_iter, 10);
        if (error_code == std::errc{} && ptr == last) {
            return fixed_max_iter;
        }
        // return the double value for coefficient (coefficient * matrix_size)
        return parse_double_strict(token);
    };

    for (int i = 1; i < argc; ++i) {
        std::string_view arg{argv[i]};
        
        if (arg == "--matrix") {
            config.matrix_name = std::string(next_token(i));
            continue;
        }

        if (arg == "--precision") {
            config.precision_level = std::string(next_token(i));
            std::string_view lvl{config.precision_level};
            if (lvl != "dd" && lvl != "dq" && lvl != "qx" && lvl != "double") {
                throw std::runtime_error("Invalid precision level. Use: dd, dq, qx, or double");
            }
            continue;
        }

        if (arg == "--tol") {
            config.tolerance = parse_double_strict(next_token(i));
            continue;
        }

        if (arg == "--max-iter") {
            config.max_iter = parse_max_iter(next_token(i));
            continue;
        }

        if (arg == "--input-dir") {
            config.input_dir = std::string(next_token(i));
            continue;
        }

        if (arg == "--export-mat") {
            config.export_mat_file = std::string(next_token(i));
            continue;
        }

        if (arg == "--export-csv") {
            config.export_csv_file = std::string(next_token(i));
            continue;
        }

        if (arg == "--help" || arg == "-h") {
            throw std::runtime_error("help");
        }

        throw std::runtime_error("Unknown argument: " + std::string(arg) + " (try --help for usage)");
    }
    
    // Matrix name MUST be specified
    if (config.matrix_name.empty()) {
        throw std::runtime_error("Matrix name is required (--matrix)");
    }
    
    return config;
}

// Helper function to resolve export file path
std::string resolveExportPath(const std::string& filename) {
    // If filename is empty, return empty
    if (filename.empty()) {
        return filename;
    }
    
    // If filename contains path separator or starts with '.', use as-is
    if (filename.find('/') != std::string::npos || filename.find('\\') != std::string::npos || 
        filename.starts_with("./") || filename.starts_with("../")) {
        return filename;
    }
    
    // Otherwise, place in outputs directory
    std::filesystem::path outputs_dir = "outputs";
    
    // Create outputs directory if it doesn't exist
    try {
        std::filesystem::create_directories(outputs_dir);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not create outputs directory: " << e.what() << '\n';
        return filename; // Fall back to current directory
    }
    
    return (outputs_dir / filename).string();
}

void printUsage(const char* program_name) {
    std::cout << "\nUsage: " << program_name << " [OPTIONS]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --matrix NAME         Matrix name (required, e.g., nos5 for nos5.mtx)\n";
    std::cout << "  --precision LEVEL     Precision level: dd, dq, qx, double (default: qx)\n";
    std::cout << "  --tol VALUE           Convergence tolerance (default: 1.0e-12)\n";
    std::cout << "  --max-iter VALUE      Maximum iterations:\n";
    std::cout << "                        - Integer: absolute number of iterations\n";
    std::cout << "                        - Float: coefficient * matrix_size (default: 2.0)\n";
    std::cout << "  --input-dir PATH      Input directory path (default: /work/inputs)\n";
    std::cout << "  --export-mat FILE     Export convergence data to MATLAB .mat file\n";
    std::cout << "  --export-csv FILE     Append final metrics to CSV (example schema)\n";
    std::cout << "  --help, -h            Show this help message\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " --matrix nos5 --precision qx --tol 1e-15\n";
    std::cout << "  " << program_name << " --matrix nos7 --precision dq --max-iter 1000\n";
    std::cout << "  " << program_name << " --matrix test --precision dd --max-iter 2.5\n";
    std::cout << "  " << program_name << " --matrix nos5 --precision double --tol 1e-10\n";
    std::cout << "  " << program_name << " --matrix nos5 --precision dq --export-mat results.mat\n";
    std::cout << "  " << program_name << " --matrix nos5 --precision dq --export-csv summary.csv\n\n";
}

// Template solver function
template<typename T>
int solveCG(const SolverConfig& config) {
    using Traits = bailey::PrecisionTraits<T>;
    using MatrixType = typename Traits::matrix_type;
    using VectorType = typename Traits::vector_type;
    
    // search matrix in the input directory (check the structure of the input directory)
    // CHECK: Download the matrix with https://github.com/sk-226/ssdownload
    std::filesystem::path matrix_dir = std::filesystem::path(config.input_dir) / config.matrix_name;
    std::filesystem::path matrix_path = matrix_dir / (matrix_dir.filename().string() + ".mtx");
    std::cout << "Loading matrix: " << matrix_path.string() << " (precision: " << Traits::name() << ")" << '\n';
    
    MatrixType A = io::loadMatrixMarket<T>(matrix_path);
    int n = A.rows();
    
    std::cout << "Matrix size: " << n << " x " << A.cols() << '\n';
    std::cout << "Non-zeros: " << A.nonZeros() << '\n';
    
    // Calculate max iterations
    int max_iterations = algorithms::resolve_max_iterations(config.max_iter, n);
    
    std::cout << "Max iterations: " << max_iterations << '\n';
    std::cout << std::scientific << std::setprecision(2) << "Tolerance: " << config.tolerance << '\n';
    
    // Set up problem: Ax = b where x_true = ones(n)
    VectorType x_true = VectorType::Ones(n);
    VectorType b = A * x_true;
    VectorType x = VectorType::Zero(n);  // Initial guess
    
    std::cout << "\nStarting CG iterations...\n";
    
    auto result = algorithms::conjugateGradient<T>(A, b, x, x_true, max_iterations, config.tolerance);
    
    // Print results
    algorithms::print_results(result, config.matrix_name + ".mtx");
    
    // Export to MATLAB .mat file if requested
    if (!config.export_mat_file.empty()) {
#ifdef ENABLE_MAT_EXPORT
        std::string export_path = resolveExportPath(config.export_mat_file);
        std::cout << "\nExporting convergence data to " << export_path << "..." << '\n';
        bool export_success = io::MatExporter::export_convergence_data(
            result, 
            export_path, 
            config.matrix_name, 
            config.precision_level
        );
        if (export_success) {
            std::cout << "Export successful." << '\n';
        } else {
            std::cerr << "Warning: Export failed." << '\n';
        }
#else
        std::cerr << "Warning: MATLAB export not available - built without matio-cpp support." << '\n';
#endif
    }
    
    // Export to CSV (final metrics) if requested
    if (!config.export_csv_file.empty()) {
        std::string csv_path = resolveExportPath(config.export_csv_file);
        bool csv_ok = io::CsvExporter::append_row(
            csv_path,
            config.matrix_name,
            n,
            A.nonZeros(),
            "cg",
            config.tolerance,
            max_iterations,
            "ones",
            config.precision_level,
            "none",
            "",
            "CG (none)",
            0.0,
            result.converged,
            result.iterations_performed,
            result.true_relres_2,
            result.final_residual_norm,
            result.hist_relerr_2.empty() ? std::numeric_limits<double>::quiet_NaN() : result.hist_relerr_2.back(),
            result.hist_relerr_A.empty() ? std::numeric_limits<double>::quiet_NaN() : result.hist_relerr_A.back(),
            result.computation_time,
            result.converged ? std::string("reached_tol") : std::string("max_iterations")
        );
        if (csv_ok) {
            std::cout << "Appended results to " << csv_path << '\n';
        } else {
            std::cerr << "Warning: Failed to append CSV results to " << csv_path << '\n';
        }
    }
    
    return result.converged ? 0 : 2;  // Exit code 2 for non-convergence (not an error)
}

// Solver dispatcher using std::variant
int runSolver(const SolverConfig& config) {
    if (config.precision_level == "dd") {
        return solveCG<bailey::DDNumber>(config);
    }
    if (config.precision_level == "dq") {
        return solveCG<bailey::DQNumber>(config);
    }
    if (config.precision_level == "qx") {
        return solveCG<bailey::QXNumber>(config);
    }
    if (config.precision_level == "double") {
        return solveCG<double>(config);
    }
    std::cerr << "Invalid precision level: " << config.precision_level << '\n';
    return 1;
}

int main(int argc, char** argv) {
    try {
        SolverConfig config = parseCommandLine(argc, argv);
        return runSolver(config);
    } catch (const std::exception& e) {
        std::string error_msg = e.what();
        if (error_msg == "help") {
            printUsage(argv[0]);
            return 0;
        }
        std::cerr << "Error: " << error_msg << '\n';
        if (error_msg.find("argument") != std::string::npos || error_msg.find("required") != std::string::npos) {
            printUsage(argv[0]);
        }
        return 1;
    }
}
