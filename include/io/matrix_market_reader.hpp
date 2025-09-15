#pragma once

#include <fast_matrix_market/app/Eigen.hpp>
#include <Eigen/Sparse>
#include <string>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include "bailey/precision_traits.hpp"

// ==============================================================================
//  Matrix Market file reader (templated)
//  Supports all precision types: double, DD, DQ, QX via PrecisionTraits
// ==============================================================================

namespace io {

template<typename T>
inline typename bailey::PrecisionTraits<T>::matrix_type loadMatrixMarket(const std::string& filename) {
    using MatrixType = typename bailey::PrecisionTraits<T>::matrix_type;

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    try {
        if constexpr (std::is_same_v<T, double>) {
            MatrixType mat;
            fast_matrix_market::read_matrix_market_eigen(file, mat);
            return mat;
        } else {
            Eigen::SparseMatrix<double> tmp;
            fast_matrix_market::read_matrix_market_eigen(file, tmp);
            MatrixType mat = tmp.template cast<T>();
            return mat;
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Matrix Market read error: ") + e.what() + " for file: " + filename);
    }
}

} // namespace io

