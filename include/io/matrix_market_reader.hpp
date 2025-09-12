#pragma once

#include <fast_matrix_market/fast_matrix_market.hpp>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include "bailey/precision_traits.hpp"

// ==============================================================================
//  Matrix Market file reader (templated)
//  Supports all precision types: double, DD, DQ, QX via PrecisionTraits
// ==============================================================================

namespace io {

template<typename T>
inline typename bailey::PrecisionTraits<T>::matrix_type loadMatrixMarket(const std::string& filename) {
    using MatrixType = typename bailey::PrecisionTraits<T>::matrix_type;
    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        fast_matrix_market::matrix_market_header header;
        std::vector<int> rows, cols;
        std::vector<double> values;

        fast_matrix_market::read_matrix_market_triplet(
            file, header, rows, cols, values
        );

        MatrixType mat(header.nrows, header.ncols);
        std::vector<Eigen::Triplet<T>> triplets;
        triplets.reserve(values.size());

        for (size_t i = 0; i < values.size(); ++i) {
            int row = rows[i];
            int col = cols[i];

            if (row < 0 || row >= static_cast<int>(header.nrows) ||
                col < 0 || col >= static_cast<int>(header.ncols)) {
                throw std::runtime_error("Index out of bounds: row=" + std::to_string(row) +
                                       ", col=" + std::to_string(col) +
                                       " for matrix " + std::to_string(header.nrows) +
                                       "x" + std::to_string(header.ncols));
            }

            triplets.emplace_back(row, col, T(values[i]));
        }

        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;

    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Matrix Market read error: ") + e.what() + " for file: " + filename);
    }
}

} // namespace io

