// Quick test to verify DQ->QX direct cast compilation
#include <iostream>

#include "bailey/precision_traits.hpp"

int main() {
    bailey::DQNumber dq_val(3.141592653589793238462643383279502884197169399375105820974944592307791e0 );
    bailey::QXNumber qx_val = bailey::precision_cast<bailey::QXNumber>(dq_val);

    std::cout << "DQ value: " << bailey::to_string(dq_val) << std::endl;
    std::cout << "QX value (from DQ): " << bailey::to_string(qx_val)
              << std::endl;

    // Verify that dq[0] was used (should match approximately)
    std::cout << "dq[0] directly: " << bailey::QXNumber(dq_val.dq[0]) << std::endl;
    std::cout << "qx value: " << bailey::QXNumber(qx_val.qx) << std::endl;
    std::cout << "to_double(DQ): " << bailey::to_double(dq_val) << std::endl;

    return 0;
}
