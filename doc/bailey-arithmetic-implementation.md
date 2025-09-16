# Bailey High-Precision Arithmetic Implementation Guide

## Overview

This project implements an interface that combines C++ operator overloading with C-compatible wrappers (bind(C)) around David H. Bailey's high-precision arithmetic libraries (DDFUN, DQFUN, QXFUN) so they can be used from C++.

## Architecture

### Layer Structure

```
C++ Application Layer
     ↓ (operator overloading)
Precision types (DD/DQ/QX) + operators
     ↓ (extern "C" calls)
Fortran Wrapper Layer (ddfun_cwrap.f90 / dqfun_cwrap.f90 / qxfun_cwrap.f90)
     ↓ (bind(C) interface)
Bailey Libraries (Fortran)
     ↓ (actual computation)
High-Precision Arithmetic
```

## Precision Data Structures

```cpp
namespace bailey {
  struct DDNumber { double      dd[2];      };  // 2×double
  struct DQNumber { long double dq[2];      };  // 2×long double
  struct QXNumber { long double qx;         };  // 1×long double
}
```

**Fortran internal representation:**
- DD: `type(dd_real)` with `real(ddknd) ddr(2)`
- DQ: `type(dq_real)` with `real(dqknd) dqr(2)`
- QX: single scalar `real(qxknd)`

> [!IMPORTANT]
> QX/DQ use `long double` on the C side. On x86_64 glibc this is the 80-bit extended format and does not match the Fortran `real(16)` ABI, so run DQ/QX inside the arm64 container when possible. Use DD if reproducibility is your priority.

## DD/DQ Data Structures and Interop

### DD (Double-Double)

```fortran
type dd_real
    sequence
    real(ddknd) ddr(2)
end type
```

```cpp
namespace bailey {
struct DDNumber { double dd[2]; };
extern "C" {
    void ddadd_(const double*, const double*, double*);
    void ddsub_(const double*, const double*, double*);
    void ddmul_(const double*, const double*, double*);
    void dddiv_(const double*, const double*, double*);
    void dddqd_(const double* d, double* a);
    void ddsqrt_(const double* a, double* b);
    void dd_to_string(const double* a, int* n, char* c, int cl);
}
}
```

### DQ (Double-Quad)

```fortran
type dq_real
    sequence
    real(dqknd) dqr(2)
end type
```

```cpp
namespace bailey {
struct DQNumber { long double dq[2]; };
extern "C" {
    void dqadd_(const long double*, const long double*, long double*);
    void dqsub_(const long double*, const long double*, long double*);
    void dqmul_(const long double*, const long double*, long double*);
    void dqdiv_(const long double*, const long double*, long double*);
    void dqdqd_(const double* d, long double* a);
    void dqsqrt_(const long double* a, long double* b);
    void dq_to_string(const long double* a, int* n, char* c, int cl);
}
}
```

### Usage Notes (DD/DQ)

- On the C++ side both `DDNumber` and `DQNumber` store their two-term expansions in arrays, and the operators are thin wrappers around the Fortran routines.
- We standardized the stringification helpers to the `*_to_string` naming scheme.
- DQ relies on a `long double` ABI, so behavior can vary by platform (use the Docker image for consistency).

## Operator Overloading Implementation

### Basic Arithmetic Operators

```cpp
// Example (QX): a + b
bailey::QXNumber operator+(const bailey::QXNumber& a, const bailey::QXNumber& b) {
  bailey::QXNumber r; qxadd_(&a.qx, &b.qx, &r.qx); return r;
}
```

### Assignment Operators

```cpp
// Compound assignment: a += b
QXNumber& operator+=(QXNumber& a, const QXNumber& b) { 
    a = a + b;  // Uses operator+ defined above
    return a; 
}

// Similarly for -=, *=, /=
```

### Mathematical Functions

```cpp
// Square root
QXNumber sqrt(const QXNumber& a) { 
    QXNumber r; 
    qxsqrt_(&a.qx, &r.qx);  // Calls Bailey's qxsqrt Fortran routine
    return r; 
}
```

## Foreign Function Interface (FFI)

### Fortran-to-C Interface

The C ABI for QX uses a scalar `real(qxknd)` (received on the C++ side as a `long double*`):

```fortran
subroutine qx_add(a,b,c) bind(C,name="qxadd_")
  real(qxknd), intent(in)  :: a, b
  real(qxknd), intent(out) :: c
  c = a + b
end subroutine
```

### C++ External Declarations

```cpp
extern "C" {
    void qxadd_(const long double* a, const long double* b, long double* c);  // a + b = c
    void qxsub_(const long double* a, const long double* b, long double* c);  // a - b = c
    void qxmul_(const long double* a, const long double* b, long double* c);  // a * b = c
    void qxdiv_(const long double* a, const long double* b, long double* c);  // a / b = c
    void qxdqd_(const double* d, long double* a);                             // double → QX
    void qxsqrt_(const long double* a, long double* b);                       // sqrt(a) = b
    void qx_to_string(const long double* a, int* n, char* c, int cl);         // QX → string
}
```

## Type Conversion Functions

### Double to QXNumber

```cpp
QXNumber(double val) { 
    qx = static_cast<long double>(val);  // Cast to long double
}
```

### QXNumber to Double

```cpp
double to_double(const QXNumber& a) { 
    // Use Bailey's string conversion for accuracy
    char s[70]; 
    int d = 15; 
    qx_to_string(&a.qx, &d, s, sizeof(s));
    
    // Parse string to double
    std::string str(s, sizeof(s));
    size_t end = str.find('\0');
    if (end != std::string::npos) {
        str = str.substr(0, end);
    }
    
    try {
        return std::stod(str);
    } catch (...) {
        return static_cast<double>(a.qx);  // Fallback to most significant part
    }
}
```

### String Output

```cpp
std::ostream& operator<<(std::ostream& os, const QXNumber& q) {
  char s[128] = {0};
  int d = 33;
  qx_to_string(&q.qx, &d, s, sizeof(s));
  os << s;
  return os;
}
```

## Eigen3 Integration

### NumTraits Specialization

```cpp
namespace Eigen {
    template<> struct NumTraits<bailey::QXNumber> : GenericNumTraits<bailey::QXNumber> {
        typedef bailey::QXNumber Real; 
        typedef bailey::QXNumber NonInteger; 
        typedef bailey::QXNumber Nested;
        
        enum { 
            IsComplex = 0, 
            IsInteger = 0, 
            IsSigned = 1, 
            RequireInitialization = 1, 
            ReadCost = 1,
            AddCost = 8,
            MulCost = 16
        };
    };
}
```

This specialization lets the precision types participate in Eigen3 sparse and dense operations (see the headers for details).

## Usage Examples

### Basic Arithmetic

```cpp
QXNumber a(1.0);
QXNumber b(3.0);

QXNumber sum = a + b;        // Uses operator+, calls qxadd_
QXNumber product = a * b;    // Uses operator*, calls qxmul_
QXNumber quotient = b / a;   // Uses operator/, calls qxdiv_

std::cout << "Sum: " << sum << std::endl;  // Uses operator<<
```

### Vector Operations (with Eigen3)

```cpp
using Vec_QX = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>;

Vec_QX x = Vec_QX::Ones(100);      // Vector of QXNumber(1.0)
Vec_QX y = Vec_QX::Zero(100);      // Vector of QXNumber(0.0)

bailey::QXNumber dot_product = x.dot(y); // Eigen3 calls our operators internally
```

### Sparse Matrix Operations

```cpp
using SpMat_QX = Eigen::SparseMatrix<bailey::QXNumber>;

SpMat_QX A(100, 100);
Vec_QX x(100), b(100);

// Matrix-vector multiplication
Vec_QX result = A * x;  // Eigen3 uses our QXNumber operators
```

## Performance Considerations

### Computational Costs

1. **Operator Calls**: `+`, `-`, `*`, `/` involve C/Fortran calls plus high-precision computation.

2. **Memory Usage**: 
   - QXNumber: sizeof(long double)
   - Standard double: sizeof(double)
   - Overhead depends on platform (e.g., binary128 vs x86 extended)

3. **Precision/Speed**: QX≈33 digits, DQ≈64 digits, DD≈30 digits. Expect roughly 10–200× slower than `double`.

### Optimization Tips

1. **Minimize conversions**: Avoid frequent `to_double()` calls
2. **Vectorize operations**: Use Eigen3's vectorized operations when possible
3. **Cache results**: Store intermediate QXNumber results rather than recomputing

## Error Handling

### Common Issues

1. **Conversion failures**: String parsing in `to_double()` may fail
2. **Precision loss**: Converting QXNumber → double loses precision
3. **Memory alignment**: Fortran arrays must be properly aligned

### Debugging Tips

1. **Check string output**: Use `operator<<` to see formatted values
3. **Validate operations**: Compare with known mathematical results

## Build System Integration

### CMake Configuration

```cmake
target_link_libraries(your_target PRIVATE ddwrap ddfun dqwrap dqfun qxwrap qxfun gfortran)

# Include directories
target_include_directories(your_target PRIVATE ${QXFUN_DIR})
```

### Environment Variables

```bash
export QXFUN_DIR=/path/to/qxfun/fortran
export DQFUN_DIR=/path/to/dqfun/fortran  
export DDFUN_DIR=/path/to/ddfun/fortran
```

## Conclusion

This implementation lets you use Bailey's high-precision arithmetic libraries with natural C++ syntax. Operator overloading keeps the ergonomics close to `double` while enabling ~30 digits (DD), ~33 digits (QX), and ~64 digits (DQ).

**Key Benefits:**
- Natural C++ syntax for high-precision arithmetic
- Seamless Eigen3 integration for linear algebra
- Fortran library's computational accuracy with C++ usability
- Extensible design for additional mathematical functions
