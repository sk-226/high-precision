#include "bailey/qx_arithmetic.hpp"
#include "linear_algebra/conjugate_gradient.hpp"
#include "matrix_io/matrix_market.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

int main() {
    std::cout << "Loading nos7.mtx matrix..." << std::endl;
    
    try {
        auto start_load = std::chrono::high_resolution_clock::now();
        auto A = matrix_io::loadMatrixMarket<bailey::QXNumber>("/work/inputs/nos7.mtx");
        auto end_load = std::chrono::high_resolution_clock::now();
        
        int n = A.rows();
        std::cout << "Matrix size: " << n << " x " << A.cols() << std::endl;
        std::cout << "Non-zeros: " << A.nonZeros() << std::endl;
        
        // 真の解とRHSベクトルの設定
        Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> x_true = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>::Ones(n);
        Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> b = A * x_true;
        Eigen::Vector<bailey::QXNumber, Eigen::Dynamic> x = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>::Zero(n);  // 初期解
        
        auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load);
        std::cout << "Matrix loading time: " << load_time.count() << " ms" << std::endl;
        
        std::cout << "\nStarting CG iterations..." << std::endl;
        
        auto result = linear_algebra::conjugateGradient(A, b, x, x_true, 1000, 1e-15);
        
        // MATLAB形式での結果出力
        linear_algebra::print_num_results(result, "nos7.mtx");
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}