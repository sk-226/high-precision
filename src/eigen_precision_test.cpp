#include "bailey/precision_traits.hpp"
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

template<typename T>
double default_tolerance() {
    constexpr int digits = bailey::PrecisionTraits<T>::decimal_digits();
    const int effective_digits = std::max(6, std::min(40, digits - 4));
    return std::pow(10.0, -effective_digits);
}

template<typename T>
double to_double_safe(const T& value) {
    return to_double(value);
}

inline double to_double_safe(const double& value) {
    return value;
}

template<typename T>
bool approx_equal(const T& a, const T& b, double tol) {
    const double da = to_double_safe(a);
    const double db = to_double_safe(b);
    const double diff = std::abs(da - db);
    const double scale = std::max({1.0, std::abs(da), std::abs(db)});
    return diff <= tol * scale;
}

template<typename VectorType>
bool vector_close(const VectorType& a, const VectorType& b, double tol) {
    for (int i = 0; i < a.size(); ++i) {
        if (!approx_equal(a(i), b(i), tol)) {
            return false;
        }
    }
    return true;
}

template<typename Traits>
bool run_precision_suite(const std::string& name) {
    using T = typename Traits::scalar_type;
    using VectorType = typename Traits::vector_type;
    using MatrixType = typename Traits::matrix_type;

    const double tol = default_tolerance<T>();
    bool all_ok = true;

    auto report = [&](bool condition, const std::string& label) {
        if (condition) {
            std::cout << "[PASS] " << name << ": " << label << std::endl;
        } else {
            std::cerr << "[FAIL] " << name << ": " << label << std::endl;
            all_ok = false;
        }
    };

    VectorType v1(3);
    v1 << T(1.0), T(-1.0e-30), T(3.14159265358979323846);
    VectorType v2(3);
    v2 << T(-1.0 + 1.0e-30), T(2.0e-30), T(-2.71828182845904523536);

    // Vector addition/subtraction
    VectorType sum = v1 + v2;
    VectorType diff = v1 - v2;

    VectorType expected_sum(3);
    VectorType expected_diff(3);
    for (int i = 0; i < 3; ++i) {
        T s = v1(i);
        s = s + v2(i);
        expected_sum(i) = s;

        T d = v1(i);
        d = d - v2(i);
        expected_diff(i) = d;
    }

    report(vector_close(sum, expected_sum, tol), "Vector addition matches scalar accumulation");
    report(vector_close(diff, expected_diff, tol), "Vector subtraction matches scalar accumulation");

    // Dot product
    T dot_val = v1.dot(v2);
    T dot_expected = T(0.0);
    for (int i = 0; i < 3; ++i) {
        dot_expected = dot_expected + v1(i) * v2(i);
    }
    report(approx_equal(dot_val, dot_expected, tol), "Dot product matches manual reduction");

    // Sparse matrix-vector product
    MatrixType A(3, 3);
    std::vector<Eigen::Triplet<T>> triplets;
    triplets.emplace_back(0, 0, T(2.0));
    triplets.emplace_back(0, 1, T(-1.0e-30));
    triplets.emplace_back(0, 2, T(3.0));
    triplets.emplace_back(1, 0, T(0.5));
    triplets.emplace_back(1, 1, T(1.0));
    triplets.emplace_back(1, 2, T(-2.0));
    triplets.emplace_back(2, 0, T(1.0e30));
    triplets.emplace_back(2, 1, T(-1.0));
    triplets.emplace_back(2, 2, T(0.25));
    A.setFromTriplets(triplets.begin(), triplets.end());

    VectorType Av1 = A * v1;
    VectorType Av2 = A * v2;

    VectorType expected_Av1(3);
    VectorType expected_Av2(3);
    const std::array<std::array<double, 3>, 3> coeffs{ {
        {2.0, -1.0e-30, 3.0},
        {0.5, 1.0, -2.0},
        {1.0e30, -1.0, 0.25}
    } };

    for (int row = 0; row < 3; ++row) {
        T acc1 = T(0.0);
        T acc2 = T(0.0);
        for (int col = 0; col < 3; ++col) {
            const T coeff = T(coeffs[row][col]);
            acc1 = acc1 + coeff * v1(col);
            acc2 = acc2 + coeff * v2(col);
        }
        expected_Av1(row) = acc1;
        expected_Av2(row) = acc2;
    }

    report(vector_close(Av1, expected_Av1, tol), "A * v1 matches manual expansion");
    report(vector_close(Av2, expected_Av2, tol), "A * v2 matches manual expansion");

    VectorType Av_sum = A * (v1 + v2);
    VectorType combined = Av1 + Av2;
    report(vector_close(Av_sum, combined, tol), "Linearity: A*(v1+v2) == A*v1 + A*v2");

    return all_ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= run_precision_suite<bailey::PrecisionTraits<bailey::DDNumber>>("DD");
    ok &= run_precision_suite<bailey::PrecisionTraits<bailey::DQNumber>>("DQ");
    ok &= run_precision_suite<bailey::PrecisionTraits<bailey::QXNumber>>("QX");

    if (ok) {
        std::cout << "All Eigen high-precision checks passed." << std::endl;
        return 0;
    }

    std::cerr << "Eigen high-precision checks failed." << std::endl;
    return 1;
}

