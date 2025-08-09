#include "bailey/qx_arithmetic.hpp"
#include "io/matrix_market.hpp"
#include <iostream>

int main() {
    std::cout << "=== Simple Matrix Market Test ===" << std::endl;
    
    try {
        // Test loading nos5 matrix
        std::ifstream file("/work/inputs/nos5.mtx");
        if (!file.is_open()) {
            std::cout << "Cannot open nos5.mtx" << std::endl;
            return 1;
        }
        
        std::string line;
        if (!std::getline(file, line)) {
            std::cout << "Cannot read header" << std::endl;
            return 1;
        }
        
        std::cout << "Header: " << line << std::endl;
        
        // Check if symmetric
        bool is_symmetric = line.find("symmetric") != std::string::npos;
        std::cout << "Is symmetric: " << (is_symmetric ? "yes" : "no") << std::endl;
        
        // Skip comments and read dimensions
        while (std::getline(file, line) && line[0] == '%') {}
        
        std::istringstream iss(line);
        int nrows, ncols, nnz;
        iss >> nrows >> ncols >> nnz;
        
        std::cout << "Matrix size: " << nrows << " x " << ncols << std::endl;
        std::cout << "Stored entries: " << nnz << std::endl;
        std::cout << "Expected total entries: " << (is_symmetric ? 2 * nnz - nrows : nnz) << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}