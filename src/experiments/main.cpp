#include "bailey/qx_arithmetic.hpp"
#include "linear_algebra/conjugate_gradient.hpp"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== Bailey High-Precision Arithmetic Demo ===" << std::endl;
    
    // 基本的な算術演算のデモ
    bailey::QXNumber qa(1.0);
    bailey::QXNumber qb(3.0);
    
    std::cout << "Basic arithmetic with QXNumber:" << std::endl;
    std::cout << "qa = " << to_double(qa) << ", qb = " << to_double(qb) << std::endl;
    std::cout << "qa + qb = " << to_double(qa + qb) << std::endl;
    std::cout << "qa * qb = " << to_double(qa * qb) << std::endl;
    std::cout << "qb / qa = " << to_double(qb / qa) << std::endl;
    std::cout << "sqrt(qb) = " << std::fixed << std::setprecision(15) 
              << to_double(sqrt(qb)) << std::endl;
    
    // 小さな線形システムのテスト
    std::cout << "\n=== Small Linear System Test ===" << std::endl;
    
    Eigen::SparseMatrix<bailey::QXNumber> A(3, 3);
    std::vector<Eigen::Triplet<bailey::QXNumber>> triplets = {
        {0, 0, bailey::QXNumber(4.0)}, {0, 1, bailey::QXNumber(-1.0)}, {0, 2, bailey::QXNumber(0.0)},
        {1, 0, bailey::QXNumber(-1.0)}, {1, 1, bailey::QXNumber(4.0)}, {1, 2, bailey::QXNumber(-1.0)},
        {2, 0, bailey::QXNumber(0.0)}, {2, 1, bailey::QXNumber(-1.0)}, {2, 2, bailey::QXNumber(4.0)}
    };
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> x_true = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>::Ones(3);
    Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> b = A * x_true;
    Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> x = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>::Zero(3);
    
    auto result = linear_algebra::conjugateGradient(A, b, x, x_true, 100, 1e-15);
    
    linear_algebra::print_num_results(result, "3x3 tridiagonal");
    
    std::cout << "Solution: [" << to_double(x(0)) << ", " << to_double(x(1)) 
              << ", " << to_double(x(2)) << "]" << std::endl;
    
    return 0;
}