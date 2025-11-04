#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "bailey/precision_traits.hpp"

std::string trim(const std::string& input) {
    const auto first =
        std::find_if_not(input.begin(), input.end(),
                         [](unsigned char c) { return std::isspace(c) != 0; });
    if (first == input.end()) {
        return {};
    }
    const auto last =
        std::find_if_not(input.rbegin(), input.rend(), [](unsigned char c) {
            return std::isspace(c) != 0;
        }).base();
    return std::string(first, last);
}

std::vector<std::string> load_reference_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open reference file");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

const std::string& find_line_after(const std::vector<std::string>& lines,
                                   const std::string& marker) {
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (trim(lines[i]) == marker) {
            for (std::size_t j = i + 1; j < lines.size(); ++j) {
                if (!trim(lines[j]).empty()) {
                    return lines[j];
                }
            }
        }
    }
    throw std::runtime_error("marker not found");
}

int main() {
    auto dd_lines =
        load_reference_lines("tests/ddfun/fortran/testddfun.ref.txt");
    using Traits = bailey::PrecisionTraits<bailey::DDNumber>;
    auto t1 = bailey::from_string<Traits::scalar_type>(
        find_line_after(dd_lines, "t1 = pi:"));
    std::cout << "t1 = pi:" << t1 << '\n';
    auto t2 = bailey::from_string<Traits::scalar_type>(
        find_line_after(dd_lines, "t2 = -log(2):"));
    std::cout << "t2 = -log(2):" << t2 << '\n';
    auto result = t1 + t2;
    std::cout << "t1 + t2 = " << result << '\n';
}
