#include "bailey/precision_traits.hpp"
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string trim(const std::string& input) {
    const auto first = std::find_if_not(input.begin(), input.end(), [](unsigned char c) {
        return std::isspace(c) != 0;
    });
    if (first == input.end()) {
        return {};
    }
    const auto last = std::find_if_not(input.rbegin(), input.rend(), [](unsigned char c) {
        return std::isspace(c) != 0;
    }).base();
    return std::string(first, last);
}

std::vector<std::string> load_reference_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open reference file: " + path);
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

const std::string& find_line_after(const std::vector<std::string>& lines, const std::string& marker) {
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (trim(lines[i]) == marker) {
            for (std::size_t j = i + 1; j < lines.size(); ++j) {
                const std::string& candidate = lines[j];
                if (!trim(candidate).empty()) {
                    return lines[j];
                }
            }
            break;
        }
    }
    throw std::runtime_error("Failed to locate value for marker: " + marker);
}

template <typename Traits>
using Scalar = typename Traits::scalar_type;

template <typename Traits>
Scalar<Traits> parse_decimal(const std::string& source) {
    std::string cleaned = trim(source);
    if (cleaned.empty()) {
        return Scalar<Traits>{};
    }
    return bailey::from_string<Scalar<Traits>>(cleaned);
}

inline long double magnitude(const bailey::DDNumber& value) {
    const bailey::DDNumber tmp = bailey::abs(value);
    return std::abs(tmp.dd[0]) + std::abs(tmp.dd[1]);
}

inline long double magnitude(const bailey::DQNumber& value) {
    const bailey::DQNumber tmp = bailey::abs(value);
    return std::abs(tmp.dq[0]) + std::abs(tmp.dq[1]);
}

inline long double magnitude(const bailey::QXNumber& value) {
    const bailey::QXNumber tmp = bailey::abs(value);
    return std::abs(tmp.qx);
}

template <typename Traits>
long double relative_error(const Scalar<Traits>& computed, const Scalar<Traits>& reference) {
    using ScalarT = Scalar<Traits>;
    ScalarT diff = computed - reference;
    if (magnitude(reference) == 0.0L) {
        return magnitude(diff);
    }
    ScalarT ratio = diff / reference;
    return magnitude(ratio);
}

template <typename Traits>
long double tolerance();

template <> long double tolerance<bailey::PrecisionTraits<bailey::DDNumber>>() { return 1.0e-31L; }
template <> long double tolerance<bailey::PrecisionTraits<bailey::DQNumber>>() { return 1.0e-67L; }
template <> long double tolerance<bailey::PrecisionTraits<bailey::QXNumber>>() { return std::ldexp(1.0L, -110); }

std::array<bool, 3> parse_boolean_triplet(const std::string& raw_line) {
    std::array<bool, 3> results{false, false, false};
    std::istringstream iss(trim(raw_line));
    for (std::size_t idx = 0; idx < results.size(); ++idx) {
        std::string token;
        if (!(iss >> token)) {
            throw std::runtime_error("Boolean line has fewer than three entries: " + raw_line);
        }
        if (token == "T" || token == ".TRUE." || token == "TRUE") {
            results[idx] = true;
        } else if (token == "F" || token == ".FALSE." || token == "FALSE") {
            results[idx] = false;
        } else {
            throw std::runtime_error("Unexpected boolean token: " + token);
        }
    }
    return results;
}

template <typename Traits>
bool check_scalar(const std::vector<std::string>& lines,
                  const std::string& marker,
                  const Scalar<Traits>& computed,
                  long double eps,
                  const char* name) {
    const auto& line = find_line_after(lines, marker);
    const Scalar<Traits> expected = parse_decimal<Traits>(line);
    const long double rel = relative_error<Traits>(computed, expected);
    std::cout << "[scalar_precision_test] " << name << " :: " << marker
              << " -> relative error = " << std::scientific << rel
              << " (tolerance = " << eps << ")" << std::endl;
    if (rel > eps) {
        std::cout << "  expected=" << expected << " computed=" << computed << std::endl;
        return false;
    }
    return true;
}

template <typename Traits>
bool check_triplet(const std::vector<std::string>& lines,
                   const std::string& marker,
                   const std::array<bool, 3>& computed,
                   const char* name) {
    const auto& line = find_line_after(lines, marker);
    const auto expected = parse_boolean_triplet(line);
    bool ok = (computed == expected);
    std::cout << "[scalar_precision_test] " << name << " :: " << marker
              << " -> expected (" << expected[0] << ',' << expected[1] << ',' << expected[2]
              << ") computed (" << computed[0] << ',' << computed[1] << ',' << computed[2] << ')' << std::endl;
    if (!ok) {
        std::cout << "  comparison mismatch" << std::endl;
    }
    return ok;
}

template <typename Traits>
bool run_basic_suite(const std::vector<std::string>& lines, const char* name) {
    using ScalarT = Scalar<Traits>;

    const ScalarT t1 = parse_decimal<Traits>(find_line_after(lines, "t1 = pi:"));
    const ScalarT t2 = parse_decimal<Traits>(find_line_after(lines, "t2 = -log(2):"));
    const ScalarT e1 = parse_decimal<Traits>(find_line_after(lines, "e1 = 3141/8192:"));
    const ScalarT e2 = parse_decimal<Traits>(find_line_after(lines, "e2 = 6931 / 8192:"));

    const long double eps = tolerance<Traits>();
    bool ok = true;

    ok &= check_scalar<Traits>(lines, "addition: t1+t2 =", t1 + t2, eps, name);
    ok &= check_scalar<Traits>(lines, "addition: t1+e2 =", t1 + e2, eps, name);
    ok &= check_scalar<Traits>(lines, "addition: e1+t2 =", e1 + t2, eps, name);

    ok &= check_scalar<Traits>(lines, "subtraction: t1-t2 =", t1 - t2, eps, name);
    ok &= check_scalar<Traits>(lines, "subtraction: t1-e2 =", t1 - e2, eps, name);
    ok &= check_scalar<Traits>(lines, "subtraction: e1-t2 =", e1 - t2, eps, name);

    ok &= check_scalar<Traits>(lines, "multiplication: t1*t2 =", t1 * t2, eps, name);
    ok &= check_scalar<Traits>(lines, "multiplication: t1*e2 =", t1 * e2, eps, name);
    ok &= check_scalar<Traits>(lines, "multiplication: e1*t2 =", e1 * t2, eps, name);

    ok &= check_scalar<Traits>(lines, "division: t1/t2 =", t1 / t2, eps, name);
    ok &= check_scalar<Traits>(lines, "division: t1/e2 =", t1 / e2, eps, name);
    ok &= check_scalar<Traits>(lines, "division: e1/t2 =", e1 / t2, eps, name);

    const std::array<bool, 3> equal_actual{t1 == t2, e1 == t2, t1 == e2};
    ok &= check_triplet<Traits>(lines, "equal test: t1 == t2, e1 == t2, t1 == e2", equal_actual, name);

    const std::array<bool, 3> not_equal_actual{t1 != t2, e1 != t2, t1 != e2};
    ok &= check_triplet<Traits>(lines, "not-equal test: t1 /= t2, e1 /= t2, t1 =/ e2", not_equal_actual, name);

    const std::array<bool, 3> le_actual{t1 <= t2, e1 <= t2, t1 <= e2};
    ok &= check_triplet<Traits>(lines, "less-than-or-equal test: t1 <= t2, e1 <= t2, t1 <= e2", le_actual, name);

    const std::array<bool, 3> ge_actual{t1 >= t2, e1 >= t2, t1 >= e2};
    ok &= check_triplet<Traits>(lines, "greater-than-or-equal test: t1 >= t2, e1 >= t2, t1 >= e2", ge_actual, name);

    const std::array<bool, 3> lt_actual{t1 < t2, e1 < t2, t1 < e2};
    ok &= check_triplet<Traits>(lines, "less-than test: t1 < t2, e1 < t2, t1 < e2", lt_actual, name);

    const std::array<bool, 3> gt_actual{t1 > t2, e1 > t2, t1 > e2};
    ok &= check_triplet<Traits>(lines, "greater-than test: t1 > t2, e1 > t2, t1 > e2", gt_actual, name);

    ok &= check_scalar<Traits>(lines, "abs(t2) =", bailey::abs(t2), eps, name);
    ok &= check_scalar<Traits>(lines, "sqrt(t1) =", bailey::sqrt(t1), eps, name);

    return ok;
}

} // namespace

int main() {
    try {
        std::cout << "[scalar_precision_test] start" << std::endl;
        const auto dd_lines = load_reference_lines("tests/ddfun/fortran/testddfun.ref.txt");
        // for dqfun, qxfun, both use testdqfun.ref.txt
        const auto dq_lines = load_reference_lines("tests/dqfun/fortran/testdqfun.ref.txt");
        std::cout << "[scalar_precision_test] reference files loaded" << std::endl;

        bool ok = true;
        ok &= run_basic_suite<bailey::PrecisionTraits<bailey::DDNumber>>(dd_lines, "DD");
        ok &= run_basic_suite<bailey::PrecisionTraits<bailey::DQNumber>>(dq_lines, "DQ");
        ok &= run_basic_suite<bailey::PrecisionTraits<bailey::QXNumber>>(dq_lines, "QX");

        if (!ok) {
            std::cout << "[scalar_precision_test] failures detected" << std::endl;
            return 1;
        }

        std::cout << "[scalar_precision_test] all checks passed" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "scalar_precision_test aborted: " << ex.what() << std::endl;
        return 1;
    }
}
