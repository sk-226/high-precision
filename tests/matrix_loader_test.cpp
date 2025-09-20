#include "bailey/precision_traits.hpp"
#include "io/matrix_market_reader.hpp"

#include <filesystem>
#include <iostream>
#include <string_view>

namespace {

std::filesystem::path project_root() {
    std::filesystem::path source_path{__FILE__};
    if (!source_path.is_absolute()) {
        source_path = std::filesystem::absolute(source_path);
    }
    return source_path.parent_path().parent_path();
}

template <typename SparseMatrix>
void print_sparse_summary(const SparseMatrix& matrix, std::string_view label,
                          std::size_t max_entries = 10) {
    using Index = typename SparseMatrix::Index;

    std::cout << label << '\n';
    std::cout << "  dimensions: " << matrix.rows() << " x " << matrix.cols()
              << '\n';
    std::cout << "  nonzeros : " << matrix.nonZeros() << '\n';
    std::cout << "  sample entries (row, col, value):" << '\n';

    std::size_t count = 0;
    for (Index outer = 0; outer < matrix.outerSize() && count < max_entries; ++outer) {
        for (typename SparseMatrix::InnerIterator it(matrix, outer);
             it && count < max_entries; ++it) {
            std::cout << "    (" << it.row() << ", " << it.col() << ", "
                      << it.value() << ")" << '\n';
            ++count;
        }
    }

    const Index total_nnz = matrix.nonZeros();
    if (total_nnz > static_cast<Index>(max_entries)) {
        std::cout << "    ... " << total_nnz - static_cast<Index>(max_entries)
                  << " more entries" << '\n';
    }
    std::cout << '\n';
}

}  // namespace

int main() {
    const std::filesystem::path root = project_root();
    const std::filesystem::path matrix_dir =
        root / "inputs" / "test" / "test1";
    const std::filesystem::path matrix_path = matrix_dir / "test1.mtx";

    using DoubleMatrix = bailey::PrecisionTraits<double>::matrix_type;
    using DDMatrix = bailey::PrecisionTraits<bailey::DDNumber>::matrix_type;

    std::cout << "Loading Matrix Market file: " << matrix_path << '\n';

    DoubleMatrix double_matrix =
        io::loadMatrixMarket<double>(matrix_path.string());
    DDMatrix dd_matrix = io::loadMatrixMarket<bailey::DDNumber>(matrix_path.string());

    print_sparse_summary(double_matrix, "Double precision view");
    print_sparse_summary(dd_matrix, "DD precision view");

    return 0;
}
