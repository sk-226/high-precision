# API Reference

This document provides comprehensive documentation for the high-precision arithmetic API.

## QuadDouble Class

The `QuadDouble` struct provides a C++ interface to Bailey's QX high-precision arithmetic.

### Constructor

```cpp
QuadDouble();                    // Default constructor (zero)
QuadDouble(double val);          // Convert from double precision
```

### Example
```cpp
QuadDouble zero;                 // 0.0
QuadDouble pi(3.14159265359);    // π approximation
```

## Arithmetic Operations

### Basic Operators

| Operator | Description | Example |
|----------|-------------|---------|
| `+` | Addition | `c = a + b` |
| `-` | Subtraction | `c = a - b` |
| `*` | Multiplication | `c = a * b` |
| `/` | Division | `c = a / b` |

### Assignment Operators

| Operator | Description | Example |
|----------|-------------|---------|
| `+=` | Add and assign | `a += b` |
| `-=` | Subtract and assign | `a -= b` |
| `*=` | Multiply and assign | `a *= b` |
| `/=` | Divide and assign | `a /= b` |

### Example Usage

```cpp
QuadDouble a(2.0);
QuadDouble b(3.0);

QuadDouble sum = a + b;          // 5.0
QuadDouble product = a * b;      // 6.0

a += QuadDouble(1.0);            // a becomes 3.0
b *= QuadDouble(2.0);            // b becomes 6.0
```

## Mathematical Functions

### Square Root

```cpp
QuadDouble sqrt(const QuadDouble& x);
```

Computes the square root with quad precision.

**Example:**
```cpp
QuadDouble x(2.0);
QuadDouble result = sqrt(x);     // √2 with ~128 digit precision
```

### Utility Functions

```cpp
double to_double(const QuadDouble& x);
```

Converts QuadDouble to double precision for comparisons and output.

**Example:**
```cpp
QuadDouble x(3.14159265359);
double approx = to_double(x);    // 3.14159265359 (double precision)
```

## Matrix Operations

### Type Definitions

```cpp
using SpMat_QD = Eigen::SparseMatrix<QuadDouble>;
using Vec_QD = Eigen::Vector<QuadDouble, Eigen::Dynamic>;
```

### Matrix Construction

```cpp
// Create sparse matrix
const int n = 5;
SpMat_QD A(n, n);

// Fill with triplets
std::vector<Eigen::Triplet<QuadDouble>> triplets;
for (int i = 0; i < n; ++i) {
    triplets.push_back(Eigen::Triplet<QuadDouble>(i, i, QuadDouble(4.0)));
}
A.setFromTriplets(triplets.begin(), triplets.end());
```

### Vector Operations

```cpp
// Create vector
Vec_QD x(n);
for (int i = 0; i < n; ++i) {
    x(i) = QuadDouble(1.0);
}

// Matrix-vector multiplication
Vec_QD b = A * x;

// Vector operations
QuadDouble norm_squared = x.dot(x);
QuadDouble norm = sqrt(norm_squared);
```

## Conjugate Gradient Solver

### Function Signature

```cpp
int conjugateGradient(const SpMat_QD& A, const Vec_QD& b, Vec_QD& x, 
                     int max_iter, double tolerance);
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `A` | `const SpMat_QD&` | Coefficient matrix (must be symmetric positive definite) |
| `b` | `const Vec_QD&` | Right-hand side vector |
| `x` | `Vec_QD&` | Solution vector (input: initial guess, output: solution) |
| `max_iter` | `int` | Maximum number of iterations |
| `tolerance` | `double` | Convergence tolerance |

### Return Value

Returns the number of iterations taken to converge, or `max_iter` if convergence failed.

### Example Usage

```cpp
// Create system Ax = b
SpMat_QD A = createMatrix();
Vec_QD b = createRHS();
Vec_QD x = Vec_QD::Zero(A.cols());  // Initial guess

// Solve
int iterations = conjugateGradient(A, b, x, 1000, 1e-12);

if (iterations < 1000) {
    std::cout << "Converged in " << iterations << " iterations" << std::endl;
} else {
    std::cout << "Failed to converge" << std::endl;
}
```

## Fortran Interface Functions

### Core Arithmetic (QX Functions)

These functions provide the low-level interface to Bailey's QXFUN library:

```cpp
extern "C" {
    void qxadd_(const double* a, const double* b, double* c);    // c = a + b
    void qxsub_(const double* a, const double* b, double* c);    // c = a - b  
    void qxmul_(const double* a, const double* b, double* c);    // c = a * b
    void qxdiv_(const double* a, const double* b, double* c);    // c = a / b
    void qxdqd_(const double* d, double* a);                     // a = (QuadDouble)d
    void qxsqrt_(const double* a, double* b);                    // b = sqrt(a)
    void qxtoqd_(const double* a, int* n, char* c, int cl);     // Convert to string
}
```

### Usage Notes

- All functions operate on arrays of 4 doubles (QuadDouble internal representation)
- Input parameters are typically `const double*`
- Output parameters are `double*`
- These functions are wrapped by the C++ QuadDouble interface

### Direct Usage (Advanced)

```cpp
// Direct function call (not recommended for normal use)
QuadDouble a(2.0), b(3.0), result;
qxadd_(a.qd, b.qd, result.qd);  // result = a + b
```

## Output and Formatting

### String Conversion

The `qxtoqd_` function converts QuadDouble to string representation:

```cpp
std::ostream& operator<<(std::ostream& os, const QuadDouble& q);
```

**Example:**
```cpp
QuadDouble pi(3.141592653589793);
std::cout << "π = " << pi << std::endl;
// Output: π = 3.141592653589793E+00
```

### Precision Control

The string conversion uses 15 decimal digits by default. To modify precision, you would need to call the Fortran function directly:

```cpp
QuadDouble x(3.14159265359);
char buffer[70];
int precision = 20;
qxtoqd_(x.qd, &precision, buffer, sizeof(buffer));
```

## Error Handling

### Common Issues

1. **Division by Zero**: Results in infinite or NaN values
2. **Invalid Operations**: Square root of negative numbers
3. **Convergence Failure**: CG solver may not converge for ill-conditioned matrices

### Best Practices

```cpp
// Check for valid input
if (to_double(denominator) != 0.0) {
    QuadDouble result = numerator / denominator;
}

// Verify matrix properties for CG
// (Matrix should be symmetric positive definite)

// Check convergence
int iterations = conjugateGradient(A, b, x, max_iter, tol);
if (iterations >= max_iter) {
    std::cerr << "Warning: CG failed to converge" << std::endl;
}
```

## Performance Considerations

### Computational Cost

High-precision operations are significantly more expensive:

| Operation | Relative Cost |
|-----------|---------------|
| Addition/Subtraction | ~10-50x slower than double |
| Multiplication | ~50-100x slower than double |
| Division | ~100-200x slower than double |
| Square Root | ~200-500x slower than double |

### Optimization Tips

1. **Minimize Precision Conversions**: Avoid frequent double ↔ QuadDouble conversions
2. **Batch Operations**: Use Eigen3 vector/matrix operations when possible
3. **Algorithm Choice**: Sometimes algorithmic improvements are better than higher precision
4. **Memory Access**: High-precision numbers use more memory, affecting cache performance

### Memory Usage

- `QuadDouble`: 32 bytes (4 × 8-byte doubles)
- `double`: 8 bytes
- Sparse matrices scale accordingly with chosen precision