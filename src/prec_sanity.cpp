#include <iostream>
#include <iomanip>
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"

int main() {
    using namespace bailey;

    // Compute (sqrt(2))^2 - 2 in each precision
    // DD
    DDNumber two_dd(2.0);
    DDNumber s2_dd = sqrt(two_dd);
    DDNumber expr_dd = (s2_dd * s2_dd) - two_dd;

    // DQ
    DQNumber two_dq(2.0);
    DQNumber s2_dq = sqrt(two_dq);
    DQNumber expr_dq = (s2_dq * s2_dq) - two_dq;

    // QX
    QXNumber two_qx(2.0);
    QXNumber s2_qx = sqrt(two_qx);
    QXNumber expr_qx = (s2_qx * s2_qx) - two_qx;

    std::cout << "DD: " << expr_dd << "\n";
    std::cout << "DQ: " << expr_dq << "\n";
    std::cout << "QX: " << expr_qx << "\n";

    return 0;
}
