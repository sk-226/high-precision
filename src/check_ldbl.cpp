#include <iostream>
#include <limits>
#include <cfloat>

int main() {
    using LD = long double;
    std::cout << std::boolalpha;

    std::cout << "sizeof(long double): " << sizeof(LD) << "\n";
    std::cout << "std::numeric_limits<long double>::digits: "
              << std::numeric_limits<LD>::digits << "\n";
    std::cout << "std::numeric_limits<long double>::max_exponent: "
              << std::numeric_limits<LD>::max_exponent << "\n";
    std::cout << "std::numeric_limits<long double>::is_iec559: "
              << std::numeric_limits<LD>::is_iec559 << "\n";
    std::cout << "LDBL_MANT_DIG: " << LDBL_MANT_DIG << "\n";

    bool is_ieee_quad = (sizeof(LD) == 16
                      && std::numeric_limits<LD>::digits == 113
                      && std::numeric_limits<LD>::max_exponent == 16384
                      && std::numeric_limits<LD>::is_iec559);

    std::cout << "is_ieee_quad_binary128: " << is_ieee_quad << "\n";
    return 0;
}
