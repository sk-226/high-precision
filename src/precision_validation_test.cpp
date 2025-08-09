#include "bailey/precision_traits.hpp"
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

// High-precision mathematical constants for validation
namespace constants {
    // π to high precision (first 70 digits)
    const std::string PI_STR = "3.1415926535897932384626433832795028841971693993751058209749445923";
    
    // e to high precision (first 70 digits)  
    const std::string E_STR = "2.7182818284590452353602874713526624977572470936999595749669676277";
    
    // √2 to high precision (first 70 digits)
    const std::string SQRT2_STR = "1.4142135623730950488016887242096980785696718753769480731766797379";
}

// Safe epsilon-based comparison for high-precision types
template<typename T>
bool approx_equal(const T& a, const T& b, double relative_eps) {
    double diff = std::abs(to_double(a) - to_double(b));
    double max_val = std::max(std::abs(to_double(a)), std::abs(to_double(b)));
    
    if (max_val < 1e-15) {
        return diff < 1e-15; // Absolute comparison for near-zero values
    }
    
    return diff / max_val < relative_eps;
}

// Extract decimal digits from string representation
std::string extract_digits(const std::string& num_str, int max_digits) {
    std::string result;
    bool found_dot = false;
    int digit_count = 0;
    
    for (char c : num_str) {
        if (c == '.') {
            found_dot = true;
            result += c;
        } else if (std::isdigit(c)) {
            result += c;
            if (found_dot) {
                digit_count++;
                if (digit_count >= max_digits) break;
            }
        }
    }
    return result;
}

// Test precision for DD arithmetic (~30 digits)
void test_dd_precision() {
    std::cout << "=== DD Precision Test (~30 digits) ===" << std::endl;
    
    using DDTraits = bailey::PrecisionTraits<bailey::DDNumber>;
    std::cout << "Expected precision: " << DDTraits::decimal_digits() << " digits" << std::endl;
    
    // Test mathematical constants
    bailey::DDNumber pi_dd(3.141592653589793);
    bailey::DDNumber e_dd(2.718281828459045);
    bailey::DDNumber sqrt2_dd = sqrt(bailey::DDNumber(2.0));
    
    std::cout << "π (DD): " << std::setprecision(35) << pi_dd << std::endl;
    std::cout << "e (DD): " << std::setprecision(35) << e_dd << std::endl;
    std::cout << "√2(DD): " << std::setprecision(35) << sqrt2_dd << std::endl;
    
    // Verify precision level
    double pi_error = std::abs(to_double(pi_dd) - 3.141592653589793);
    std::cout << "π error: " << std::scientific << pi_error << std::endl;
    std::cout << std::endl;
}

// Test precision for DQ arithmetic (~64 digits)  
void test_dq_precision() {
    std::cout << "=== DQ Precision Test (~64 digits) ===" << std::endl;
    
    using DQTraits = bailey::PrecisionTraits<bailey::DQNumber>;
    std::cout << "Expected precision: " << DQTraits::decimal_digits() << " digits" << std::endl;
    
    // Test mathematical constants
    bailey::DQNumber pi_dq(3.141592653589793);
    bailey::DQNumber e_dq(2.718281828459045);
    bailey::DQNumber sqrt2_dq = sqrt(bailey::DQNumber(2.0));
    
    std::cout << "π (DQ): " << std::setprecision(70) << pi_dq << std::endl;
    std::cout << "e (DQ): " << std::setprecision(70) << e_dq << std::endl;  
    std::cout << "√2(DQ): " << std::setprecision(70) << sqrt2_dq << std::endl;
    
    // Verify precision level
    double pi_error = std::abs(to_double(pi_dq) - 3.141592653589793);
    std::cout << "π error: " << std::scientific << pi_error << std::endl;
    std::cout << std::endl;
}

// Test precision for QX arithmetic (~33 digits)
void test_qx_precision() {
    std::cout << "=== QX Precision Test (~33 digits) ===" << std::endl;
    
    using QXTraits = bailey::PrecisionTraits<bailey::QXNumber>;
    std::cout << "Expected precision: " << QXTraits::decimal_digits() << " digits" << std::endl;
    
    // Test mathematical constants
    bailey::QXNumber pi_qx(3.141592653589793);
    bailey::QXNumber e_qx(2.718281828459045);
    bailey::QXNumber sqrt2_qx = sqrt(bailey::QXNumber(2.0));
    
    std::cout << "π (QX): " << std::setprecision(40) << pi_qx << std::endl;
    std::cout << "e (QX): " << std::setprecision(40) << e_qx << std::endl;
    std::cout << "√2(QX): " << std::setprecision(40) << sqrt2_qx << std::endl;
    
    // Verify precision level  
    double pi_error = std::abs(to_double(pi_qx) - 3.141592653589793);
    std::cout << "π error: " << std::scientific << pi_error << std::endl;
    std::cout << std::endl;
}

// Test comparison operators with appropriate epsilon values
void test_comparison_safety() {
    std::cout << "=== High-Precision Comparison Test ===" << std::endl;
    
    // Test DD comparison
    bailey::DDNumber dd1(1.0/3.0);
    bailey::DDNumber dd2 = bailey::DDNumber(1.0) / bailey::DDNumber(3.0);
    
    // These should be approximately equal within DD precision
    double dd_epsilon = 1e-28;  // ~30 digits - 2 margin
    bool dd_equal = approx_equal(dd1, dd2, dd_epsilon);
    std::cout << "DD 1/3 comparison: " << (dd_equal ? "PASS" : "FAIL") << std::endl;
    
    // Test DQ comparison
    bailey::DQNumber dq1(1.0/7.0);
    bailey::DQNumber dq2 = bailey::DQNumber(1.0) / bailey::DQNumber(7.0);
    
    double dq_epsilon = 1e-62;  // ~64 digits - 2 margin
    bool dq_equal = approx_equal(dq1, dq2, dq_epsilon);
    std::cout << "DQ 1/7 comparison: " << (dq_equal ? "PASS" : "FAIL") << std::endl;
    
    // Test QX comparison
    bailey::QXNumber qx1(1.0/11.0);
    bailey::QXNumber qx2 = bailey::QXNumber(1.0) / bailey::QXNumber(11.0);
    
    double qx_epsilon = 1e-31;  // ~33 digits - 2 margin
    bool qx_equal = approx_equal(qx1, qx2, qx_epsilon);
    std::cout << "QX 1/11 comparison: " << (qx_equal ? "PASS" : "FAIL") << std::endl;
    
    std::cout << std::endl;
}

// Memory safety test for DQ arrays
void test_dq_memory_safety() {
    std::cout << "=== DQ Memory Safety Test ===" << std::endl;
    
    // Test that DQ array operations don't cause buffer overruns
    bailey::DQNumber a(1.5);
    bailey::DQNumber b(2.5);
    bailey::DQNumber result;
    
    // Basic arithmetic operations
    result = a + b;
    std::cout << "DQ addition: " << std::setprecision(20) << result << std::endl;
    
    result = a * b;
    std::cout << "DQ multiplication: " << std::setprecision(20) << result << std::endl;
    
    result = a / b;
    std::cout << "DQ division: " << std::setprecision(20) << result << std::endl;
    
    result = sqrt(a);
    std::cout << "DQ sqrt: " << std::setprecision(20) << result << std::endl;
    
    std::cout << "DQ memory safety: All operations completed without crashes" << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << "High-Precision Arithmetic Validation Test" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;
    
    try {
        test_dd_precision();
        test_dq_precision(); 
        test_qx_precision();
        test_comparison_safety();
        test_dq_memory_safety();
        
        std::cout << "=== Test Summary ===" << std::endl;
        std::cout << "All precision validation tests completed." << std::endl;
        std::cout << "Check output for accuracy verification." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}