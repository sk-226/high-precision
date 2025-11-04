#include "bailey/precision_traits.hpp"
#include <iostream>

int main() {
    auto pi = bailey::from_string<bailey::PrecisionTraits<bailey::DQNumber>::scalar_type>(
        "3.141592653589793238462643383279502884197169399375105820974944592307791e0" );
    std::cout << pi << '\n';
    std::cout << "parsed" << '\n';
    return 0;
}
