#include "bailey/precision_traits.hpp"
#include <iostream>

int main() {
    using bailey::DDNumber;
    std::string s1 = "3.14159265358979323846264338327950288419716939937510";
    std::string s2 = "-6.93147180559945309417232121458176568075500134360255e-1";
    auto a = bailey::from_string<DDNumber>(s1);
    auto b = bailey::from_string<DDNumber>(s2);
    auto sum = a + b;
    std::cout << bailey::to_string(sum) << '\n';
    return 0;
}
