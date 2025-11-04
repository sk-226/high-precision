#include "bailey/precision_traits.hpp"
#include <iostream>

int main() {
    auto val = bailey::from_string<bailey::DDNumber>("3.14159265358979323846264338327950288419716939937510");
    auto result = bailey::abs(val);
    std::cout << bailey::to_string(result) << '\n';
    std::cout << result << '\n';
    return 0;
}
