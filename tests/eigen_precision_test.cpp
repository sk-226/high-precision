// cmake --build /work/build --config Release --target eigen_precision_test -j

#include <iomanip>
#include <iostream>

#include "bailey/precision_traits.hpp"

template <typename T>
void test_eigen_precision(void) {
    using Traits = bailey::PrecisionTraits<T>;
    using Vector = typename Traits::vector_type;
    using Matrix = typename Traits::matrix_type;
    using Scalar = typename Traits::scalar_type;
    using Triplet = Eigen::Triplet<Scalar>;

    Matrix A(2, 5);
    Scalar sqrt_2 = sqrt(Scalar(2.0));
    std::vector<Triplet> triplets = {
        {0, 0, Scalar(2.0)}, {0, 1, Scalar(-1.0 * sqrt_2)},

        {1, 2, Scalar(1.0)}, {1, 3, Scalar(1.0)},  {1, 4, Scalar(1.0)},
    };
    A.setFromTriplets(triplets.begin(), triplets.end());

    Vector x(5);
    x.setZero();
    x(0) = Scalar(1.0);
    x(1) = sqrt_2;
    x(2) = Scalar(1e-18);
    x(3) = Scalar(1.0);
    x(4) = Scalar(-1e-18);

    Vector sq(2);
    sq.setZero();
    sq(0) = Scalar(sqrt(Scalar(3.0)));

    Vector y = A * x;
    Scalar dot_product = sq.dot(sq);

    // Configure high-precision scientific output temporarily
    std::streamsize old_precision = std::cout.precision();
    std::ios_base::fmtflags old_flags = std::cout.flags();
    std::cout.setf(std::ios::scientific);
    std::cout << std::setprecision(80);

    std::cout << "================================================= \n";
    std::cout << "2 - sqrt(2)^2" << " AND " << "1e-18 + 1 - 1e-18" << " AND " << "sqrt(2)^2" << "\n\n";
    std::cout << "A * x = \n" << y << '\n';
    std::cout << "sq.dot(sq) = " << dot_product << '\n';
    std::cout << "================================================= \n\n";

    // Restore previous stream settings
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
}

int main() {
    test_eigen_precision<double>();
    test_eigen_precision<bailey::DDNumber>();
    test_eigen_precision<bailey::QXNumber>();
    test_eigen_precision<bailey::DQNumber>();

    return 0;
}
