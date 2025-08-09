# Bailey High-Precision Arithmetic Implementation Guide

## Overview

本プロジェクトでは、David H. Bailey氏のhigh-precision arithmetic libraries（DDFUN, DQFUN, QXFUN）をC++から使用するために、C++の演算子オーバーロード機能を活用してFortranライブラリへのインターフェースを実装しています。

## Architecture

### Layer Structure

```
C++ Application Layer
     ↓ (operator overloading)
QuadDouble struct + operators
     ↓ (extern "C" calls)
Fortran Wrapper Layer (qxfun_cwrap.f90)
     ↓ (bind(C) interface)
Bailey's QX Library (Fortran)
     ↓ (actual computation)
High-Precision Arithmetic
```

## QuadDouble Data Structure

### Internal Representation

```cpp
struct QuadDouble {
    double qd[4] = {0.0, 0.0, 0.0, 0.0};
    QuadDouble() = default;
    QuadDouble(double val) { qxdqd_(&val, qd); }
};
```

**Bailey QX Library Internal Format:**
- `qd[0] + qd[1] + qd[2] + qd[3]` = actual high-precision value
- `qd[0]`: Most significant part (highest order bits)
- `qd[1]`: First correction term
- `qd[2]`: Second correction term  
- `qd[3]`: Third correction term

**Precision Level:**
- QX: ~128 decimal digits (quad-precision extended)
- DQ: ~64 decimal digits (quad-double)
- DD: ~32 decimal digits (double-double)

## Operator Overloading Implementation

### Basic Arithmetic Operators

```cpp
// Addition: a + b
QuadDouble operator+(const QuadDouble& a, const QuadDouble& b) { 
    QuadDouble r; 
    qxadd_(a.qd, b.qd, r.qd);  // Calls Bailey's qxadd Fortran routine
    return r; 
}

// Subtraction: a - b  
QuadDouble operator-(const QuadDouble& a, const QuadDouble& b) { 
    QuadDouble r; 
    qxsub_(a.qd, b.qd, r.qd);  // Calls Bailey's qxsub Fortran routine
    return r; 
}

// Multiplication: a * b
QuadDouble operator*(const QuadDouble& a, const QuadDouble& b) { 
    QuadDouble r; 
    qxmul_(a.qd, b.qd, r.qd);  // Calls Bailey's qxmul Fortran routine
    return r; 
}

// Division: a / b
QuadDouble operator/(const QuadDouble& a, const QuadDouble& b) { 
    QuadDouble r; 
    qxdiv_(a.qd, b.qd, r.qd);  // Calls Bailey's qxdiv Fortran routine
    return r; 
}
```

### Assignment Operators

```cpp
// Compound assignment: a += b
QuadDouble& operator+=(QuadDouble& a, const QuadDouble& b) { 
    a = a + b;  // Uses operator+ defined above
    return a; 
}

// Similarly for -=, *=, /=
```

### Mathematical Functions

```cpp
// Square root
QuadDouble sqrt(const QuadDouble& a) { 
    QuadDouble r; 
    qxsqrt_(a.qd, r.qd);  // Calls Bailey's qxsqrt Fortran routine
    return r; 
}
```

## Foreign Function Interface (FFI)

### Fortran-to-C Interface

**File: `fortran/qxfun_cwrap.f90`**

```fortran
! Fortran wrapper with C-compatible interface
module qxfun_cwrap
    use iso_c_binding
    use qxfunmod
    implicit none

contains
    ! Addition wrapper
    subroutine qxadd_c(a, b, c) bind(c, name='qxadd_')
        real(c_double), intent(in) :: a(4), b(4)
        real(c_double), intent(out) :: c(4)
        
        type(qx_real) :: qa, qb, qc
        
        ! Convert C arrays to QX format
        qa%qxr(:) = real(a(:), qxknd)
        qb%qxr(:) = real(b(:), qxknd)
        
        ! Perform high-precision addition
        call qxadd(qa, qb, qc)
        
        ! Convert back to C format
        c(:) = real(qc%qxr(:), c_double)
    end subroutine qxadd_c
end module
```

### C++ External Declarations

```cpp
extern "C" {
    void qxadd_(const double* a, const double* b, double* c);  // a + b = c
    void qxsub_(const double* a, const double* b, double* c);  // a - b = c
    void qxmul_(const double* a, const double* b, double* c);  // a * b = c
    void qxdiv_(const double* a, const double* b, double* c);  // a / b = c
    void qxdqd_(const double* d, double* a);                   // double → QX
    void qxsqrt_(const double* a, double* b);                  // sqrt(a) = b
    void qxtoqd_(const double* a, int* n, char* c, int cl);   // QX → string
}
```

## Type Conversion Functions

### Double to QuadDouble

```cpp
QuadDouble(double val) { 
    qxdqd_(&val, qd);  // Calls Bailey's conversion routine
}
```

### QuadDouble to Double

```cpp
double to_double(const QuadDouble& a) { 
    // Use Bailey's string conversion for accuracy
    char s[70]; 
    int d = 15; 
    qxtoqd_(a.qd, &d, s, sizeof(s));
    
    // Parse string to double
    std::string str(s, sizeof(s));
    size_t end = str.find('\0');
    if (end != std::string::npos) {
        str = str.substr(0, end);
    }
    
    try {
        return std::stod(str);
    } catch (...) {
        return a.qd[0];  // Fallback to most significant part
    }
}
```

### String Output

```cpp
std::ostream& operator<<(std::ostream& os, const QuadDouble& q) {
    char s[70]; 
    int d = 15; 
    qxtoqd_(q.qd, &d, s, sizeof(s));  // Convert to formatted string
    os << std::string(s, sizeof(s)); 
    return os;
}
```

## Eigen3 Integration

### NumTraits Specialization

```cpp
namespace Eigen {
    template<> struct NumTraits<QuadDouble> : GenericNumTraits<QuadDouble> {
        typedef QuadDouble Real; 
        typedef QuadDouble NonInteger; 
        typedef QuadDouble Nested;
        
        enum { 
            IsComplex = 0, 
            IsInteger = 0, 
            IsSigned = 1, 
            RequireInitialization = 1, 
            ReadCost = 4,      // Cost of reading a value
            AddCost = 32,      // Cost of addition operation
            MulCost = 64       // Cost of multiplication operation
        };
    };
}
```

This specialization enables Eigen3 to work with QuadDouble types in:
- Sparse matrices: `Eigen::SparseMatrix<QuadDouble>`
- Dense vectors: `Eigen::Vector<QuadDouble, Eigen::Dynamic>`
- Matrix operations: dot products, matrix-vector multiplication, etc.

## Usage Examples

### Basic Arithmetic

```cpp
QuadDouble a(1.0);
QuadDouble b(3.0);

QuadDouble sum = a + b;        // Uses operator+, calls qxadd_
QuadDouble product = a * b;    // Uses operator*, calls qxmul_
QuadDouble quotient = b / a;   // Uses operator/, calls qxdiv_

std::cout << "Sum: " << sum << std::endl;  // Uses operator<<
```

### Vector Operations (with Eigen3)

```cpp
using Vec_QD = Eigen::Vector<QuadDouble, Eigen::Dynamic>;

Vec_QD x = Vec_QD::Ones(100);      // Vector of QuadDouble(1.0)
Vec_QD y = Vec_QD::Zero(100);      // Vector of QuadDouble(0.0)

QuadDouble dot_product = x.dot(y); // Eigen3 calls our operators internally
```

### Sparse Matrix Operations

```cpp
using SpMat_QD = Eigen::SparseMatrix<QuadDouble>;

SpMat_QD A(100, 100);
Vec_QD x(100), b(100);

// Matrix-vector multiplication
Vec_QD result = A * x;  // Eigen3 uses our QuadDouble operators
```

## Performance Considerations

### Computational Costs

1. **Operator Calls**: Each `+`, `-`, `*`, `/` operation involves:
   - C++ → Fortran interface call
   - Array copying (4 doubles)
   - High-precision arithmetic computation
   - Return value copying

2. **Memory Usage**: 
   - QuadDouble: 32 bytes (4 × 8 bytes)
   - Standard double: 8 bytes
   - **4x memory overhead**

3. **Speed Comparison**:
   - QuadDouble operations: ~10-100x slower than double
   - Benefit: ~128 decimal digits precision vs. ~15 digits

### Optimization Tips

1. **Minimize conversions**: Avoid frequent `to_double()` calls
2. **Vectorize operations**: Use Eigen3's vectorized operations when possible
3. **Cache results**: Store intermediate QuadDouble results rather than recomputing

## Error Handling

### Common Issues

1. **Conversion failures**: String parsing in `to_double()` may fail
2. **Precision loss**: Converting QuadDouble → double loses precision
3. **Memory alignment**: Fortran arrays must be properly aligned

### Debugging Tips

1. **Print raw values**: Access `qd[0], qd[1], qd[2], qd[3]` directly
2. **Check string output**: Use `operator<<` to see full precision
3. **Validate operations**: Compare with known mathematical results

## Build System Integration

### CMake Configuration

```cmake
# Link Bailey libraries
target_link_libraries(your_target PRIVATE qxwrap qxfun gfortran)

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

この実装により、Bailey氏の高精度算術ライブラリをC++の自然な構文で使用できるようになります。演算子オーバーロードにより、通常のdouble型を使うのと同じ感覚で128桁精度の計算が可能になります。

**Key Benefits:**
- Natural C++ syntax for high-precision arithmetic
- Seamless Eigen3 integration for linear algebra
- Fortran library's computational accuracy with C++ usability
- Extensible design for additional mathematical functions