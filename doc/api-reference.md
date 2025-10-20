# API Reference

This document describes the public C++ API used by the solver and examples. The project supports three precision types only:

> [!IMPORTANT]
> Supported: `dd` (double-double), `dq` (double-quad = 2×quad), `qx` (quad/extended). Not supported: `qd` (quad-double, 4×double).

## Precision Types

- `bailey::DDNumber` — DDFUN, 2×`double`, ~30 digits
- `bailey::DQNumber` — DQFUN, 2×`long double`, ~66 digits
- `bailey::QXNumber` — QXFUN, scalar `long double`, ~33 digits

Common functions are provided for all three types via free functions and operators.

### Constructors

```cpp
bailey::DDNumber xdd;                 // zero
bailey::DDNumber xdd_from_double(1.0);

bailey::DQNumber xdq;                 // zero
bailey::DQNumber xdq_from_double(1.0);

bailey::QXNumber xqx;                 // zero
bailey::QXNumber xqx_from_double(1.0);
```

### Arithmetic and Math

```cpp
// All types support +,-,*,/ and compound assignments
auto c = a + b;
a *= b;

// sqrt
auto r = sqrt(a);

// to_string / to_double utilities
std::string s = to_string(a);         // digits default: DD=32 (~30-digit accuracy), DQ=66, QX=33
double d = to_double(a);              // best-effort conversion for logging
```

> [!TIP]
> Prefer staying in one precision type during computations. Convert to `double` only for logging or comparisons against loose thresholds.

## Matrix and Vector Types

Use `PrecisionTraits<T>` to obtain the Eigen types for any precision:

```cpp
using Traits = bailey::PrecisionTraits<bailey::QXNumber>; // or DDNumber/DQNumber/double
using Matrix = Traits::matrix_type;   // Eigen::SparseMatrix<T>
using Vector = Traits::vector_type;   // Eigen::Vector<T, Eigen::Dynamic>
```

### Constructing Matrices

```cpp
const int n = 5;
Matrix A(n, n);
std::vector<Eigen::Triplet<Traits::scalar_type>> triplets;
for (int i = 0; i < n; ++i) triplets.emplace_back(i, i, Traits::scalar_type(4.0));
A.setFromTriplets(triplets.begin(), triplets.end());

Vector x = Vector::Ones(n);
Vector b = A * x;
```

## Conjugate Gradient Solver

```cpp
// include/algorithms/conjugate_gradient.hpp
namespace algorithms {
  template<typename T>
  CGResult<T> conjugateGradient(
      const bailey::PrecisionTraits<T>::matrix_type& A,
      const bailey::PrecisionTraits<T>::vector_type& b,
      bailey::PrecisionTraits<T>::vector_type& x,
      const bailey::PrecisionTraits<T>::vector_type& x_true,
      int max_iter,
      double tolerance);
}
```

Result fields include `iterations_performed`, `converged`, `hist_relres_2`, `true_relres_2`, and `computation_time`.

### Example

```cpp
using T = bailey::DQNumber;                   // choose DDNumber/DQNumber/QXNumber
using Traits = bailey::PrecisionTraits<T>;
Traits::matrix_type A = /* ... */;
Traits::vector_type b = /* ... */;
Traits::vector_type x = Traits::vector_type::Zero(A.cols());
auto result = algorithms::conjugateGradient<T>(A, b, x, Traits::vector_type::Ones(A.cols()), 1000, 1e-12);
```

## Low-level Fortran Interfaces

These are declared and used internally. Signatures (as seen in headers):

```cpp
extern "C" {
  // DD (arrays of 2 doubles)
  void ddadd_(const double*, const double*, double*);
  void ddsub_(const double*, const double*, double*);
  void ddmul_(const double*, const double*, double*);
  void dddiv_(const double*, const double*, double*);
  void dddqd_(const double* d, double* a);
  void ddsqrt_(const double* a, double* b);
  void dd_to_string(const double* a, int* n, char* c, int cl);

  // DQ (arrays of 2 long doubles)
  void dqadd_(const long double*, const long double*, long double*);
  void dqsub_(const long double*, const long double*, long double*);
  void dqmul_(const long double*, const long double*, long double*);
  void dqdiv_(const long double*, const long double*, long double*);
  void dqdqd_(const double* d, long double* a);
  void dqsqrt_(const long double* a, long double* b);
  void dq_to_string(const long double* a, int* n, char* c, int cl);

  // QX (scalar long double)
  void qxadd_(const long double*, const long double*, long double*);
  void qxsub_(const long double*, const long double*, long double*);
  void qxmul_(const long double*, const long double*, long double*);
  void qxdiv_(const long double*, const long double*, long double*);
  void qxdqd_(const double* d, long double* a);
  void qxsqrt_(const long double* a, long double* b);
  void qx_to_string(const long double* a, int* n, char* c, int cl);
}
```

> [!WARNING]
> The QX/DQ C ABI relies on `long double`. On x86_64 glibc, this is 80‑bit extended, which does not match Fortran quad. Use the provided arm64 container for DQ/QX.

## Output and Formatting

```cpp
// Stream output uses library formatting with default digits per type
std::cout << bailey::QXNumber(3.141592653589793L) << "\n";  // ~33 digits
```

To control digits, call the `*_to_string` Fortran functions (see headers) and print the returned buffer.

## Performance Notes

- Extended precision is slower than `double` by 1–2 orders of magnitude.
- Memory footprint scales with the scalar size: DD≈16 bytes, DQ≈32 bytes (arm64), QX≈16 bytes (arm64).
- Minimize precision conversions and prefer batched Eigen operations.
