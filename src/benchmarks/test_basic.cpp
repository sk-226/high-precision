#include "bailey/qx_arithmetic.hpp"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== QXNumber Basic Operations Test ===" << std::endl;
    
    bailey::QXNumber a(2.0);
    bailey::QXNumber b(6.0);
    
    std::cout << "a = " << to_double(a) << std::endl;
    std::cout << "b = " << to_double(b) << std::endl;
    
    bailey::QXNumber sum = a + b;
    bailey::QXNumber diff = b - a;
    bailey::QXNumber prod = a * b;
    bailey::QXNumber quot = b / a;
    bailey::QXNumber sq = sqrt(a);
    
    std::cout << "\nResults:" << std::endl;
    std::cout << "a + b = " << to_double(sum) << " (expected 8.0)" << std::endl;
    std::cout << "b - a = " << to_double(diff) << " (expected 4.0)" << std::endl;
    std::cout << "a * b = " << to_double(prod) << " (expected 12.0)" << std::endl;
    std::cout << "b / a = " << to_double(quot) << " (expected 3.0)" << std::endl;
    std::cout << "sqrt(a) = " << std::fixed << std::setprecision(15) << to_double(sq) 
              << " (expected ~1.414213562373095)" << std::endl;
    
    // 高精度計算の確認
    bailey::QXNumber pi_approx(3.141592653589793238462643383279502884197);
    std::cout << "\nHigh-precision value:" << std::endl;
    std::cout << "π ≈ " << std::setprecision(30) << to_double(pi_approx) << std::endl;
    
    return 0;
}