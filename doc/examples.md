# Examples and Use Cases

Practical examples for DD/DQ/QX use. All snippets compile in C++17+ and use Eigen.

> [!IMPORTANT]
> Choose one precision per program section: `bailey::DDNumber`, `bailey::DQNumber`, or `bailey::QXNumber`. This project does not support `qd` (quad-double).

## Basic Arithmetic

```cpp
#include <iostream>
#include "bailey/qx_arithmetic.hpp"  // or dd_arithmetic.hpp / dq_arithmetic.hpp

int main() {
  using T = bailey::QXNumber;  // pick DDNumber/DQNumber/QXNumber
  T a(1.0), b(3.0);
  auto sum = a + b;
  auto product = a * b;
  auto quotient = b / a;
  std::cout << "1 + 3 = " << sum << "\n";
  std::cout << "1 * 3 = " << product << "\n";
  std::cout << "3 / 1 = " << quotient << "\n";
}
```

## Matrix Operations

```cpp
#include <Eigen/Sparse>
#include <vector>

template<class T>
Eigen::SparseMatrix<T> createTridiagonalMatrix(int n, T diag, T off_diag) {
  Eigen::SparseMatrix<T> A(n, n);
  std::vector<Eigen::Triplet<T>> triplets;
  for (int i = 0; i < n; ++i) {
    triplets.emplace_back(i, i, diag);
    if (i > 0)    triplets.emplace_back(i, i-1, off_diag);
    if (i < n-1)  triplets.emplace_back(i, i+1, off_diag);
  }
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

int main() {
  using T = bailey::DDNumber;
  int n = 100;
  auto A = createTridiagonalMatrix<T>(n, T(4.0), T(-1.0));
  std::cout << "Non-zeros: " << A.nonZeros() << "\n";
}
```

## Linear System Solving (CG)

```cpp
#include <iostream>
#include <Eigen/Sparse>
#include "algorithms/conjugate_gradient.hpp"

int main() {
  using T = bailey::QXNumber;  // or DDNumber/DQNumber
  int n = 100;
  Eigen::SparseMatrix<T> A(n, n);
  std::vector<Eigen::Triplet<T>> trips;
  for (int i = 0; i < n; ++i) trips.emplace_back(i, i, T(4.0));
  A.setFromTriplets(trips.begin(), trips.end());

  Eigen::Vector<T, Eigen::Dynamic> x = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
  Eigen::Vector<T, Eigen::Dynamic> x_true = Eigen::Vector<T, Eigen::Dynamic>::Ones(n);
  auto b = A * x_true;

  auto res = algorithms::conjugateGradient<T>(A, b, x, x_true, 1000, 1e-12);
  std::cout << (res.converged ? "Converged" : "Not converged")
            << ", iters=" << res.iterations_performed << "\n";
}
```

## Iterative Refinement (sketch)

```cpp
template<class T>
Eigen::Vector<T, Eigen::Dynamic>
iterativeRefinement(const Eigen::SparseMatrix<T>& A,
                    const Eigen::Vector<T,Eigen::Dynamic>& b,
                    Eigen::Vector<T,Eigen::Dynamic> x,
                    int max_refinements) {
  for (int r = 0; r < max_refinements; ++r) {
    auto residual = b - A * x;
    auto norm = sqrt(residual.dot(residual));
    if (to_double(norm) < 1e-15) break;
    auto dx = Eigen::Vector<T,Eigen::Dynamic>::Zero(A.cols());
    auto x_true = Eigen::Vector<T,Eigen::Dynamic>::Zero(A.cols()); // unknown
    algorithms::conjugateGradient<T>(A, residual, dx, x_true, 100, 1e-12);
    x += dx;
  }
  return x;
}
```

## Numerical Routines

```cpp
// Machin-like π (series)
template<class T>
T arctan_series(T x, int terms) {
  T result(0.0), xp = x, xsq = x*x;
  for (int n = 0; n < terms; ++n) {
    T term = xp / T(2*n + 1);
    result = (n % 2 == 0) ? (result + term) : (result - term);
    xp *= xsq;
  }
  return result;
}

template<class T>
T compute_pi_machin(int terms) {
  T one_fifth = T(0.2), one_239 = T(1.0/239.0);
  T pi_quarter = T(4.0) * arctan_series(one_fifth, terms) - arctan_series(one_239, terms);
  return T(4.0) * pi_quarter;
}

// Newton sqrt
template<class T>
T newton_sqrt(T x, int max_it=20) {
  T g = x;
  for (int i = 0; i < max_it; ++i) {
    T ng = g - (g*g - x) / (T(2.0) * g);
    if (to_double((ng-g)*(ng-g)) < 1e-30) return ng;
    g = ng;
  }
  return g;
}
```

## Performance Tips

```cpp
// Minimize conversions
template<class T>
double sum_as_double(const std::vector<double>& a) {
  T s(0.0); for (double v : a) s += T(v); return to_double(s);
}
```

> [!TIP]
> Convert once at the boundary (inputs/outputs). Keep all intermediate math in the chosen precision.

## Error Analysis

```cpp
template<class T>
T estimate_spectral_radius(const Eigen::SparseMatrix<T>& A, int iters=50) {
  int n = A.rows();
  Eigen::Vector<T, Eigen::Dynamic> v = Eigen::Vector<T, Eigen::Dynamic>::Random(n);
  v /= sqrt(v.dot(v));
  T eig(0.0);
  for (int i=0;i<iters;++i){ auto Av=A*v; eig=v.dot(Av); v=Av/sqrt(Av.dot(Av)); }
  return eig;
}

template<class T>
void verify_solution_accuracy() {
  int n=10;
  auto A = Eigen::SparseMatrix<T>::Identity(n,n);
  auto x_true = Eigen::Vector<T,Eigen::Dynamic>::Ones(n);
  auto b = A * x_true;
  auto x = Eigen::Vector<T,Eigen::Dynamic>::Zero(n);
  auto res = algorithms::conjugateGradient<T>(A, b, x, x_true, 100, 1e-15);
  auto err = x - x_true;
  auto err2 = sqrt(err.dot(err));
  std::cout << "Rel error: " << to_double(err2 / sqrt(x_true.dot(x_true))) << "\n";
  (void)res;
}
```

## Integrating with Double Code

```cpp
template<class T>
std::vector<double> solve_with_high_precision(
  const std::vector<std::vector<double>>& Ad,
  const std::vector<double>& bd) {
  int n = (int)Ad.size();
  Eigen::SparseMatrix<T> A(n,n);
  std::vector<Eigen::Triplet<T>> trips;
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)
    if (std::abs(Ad[i][j])>1e-15) trips.emplace_back(i,j,T(Ad[i][j]));
  A.setFromTriplets(trips.begin(), trips.end());
  Eigen::Vector<T,Eigen::Dynamic> b(n); for(int i=0;i<n;++i) b(i)=T(bd[i]);
  auto x = Eigen::Vector<T,Eigen::Dynamic>::Zero(n);
  auto x_true = Eigen::Vector<T,Eigen::Dynamic>::Zero(n);
  algorithms::conjugateGradient<T>(A,b,x,x_true,1000,1e-15);
  std::vector<double> out(n); for(int i=0;i<n;++i) out[i]=to_double(x(i));
  return out;
}
```

These examples are consistent with the current implementation and avoid `qd`/`QuadDouble` terminology.
