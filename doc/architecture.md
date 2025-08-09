# Architecture

This document describes the system architecture and design decisions for the Bailey high-precision numerical analysis project.

## System Overview

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   C++ Application   │    │  Fortran Wrappers   │    │  Bailey Libraries   │
│                 │    │                  │    │                 │
│  ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│  │ QuadDouble  │ │◄──►│ │ qxfun_cwrap  │ │◄──►│ │   QXFUN     │ │
│  │   Struct    │ │    │ │              │ │    │ │ (qxmodule)  │ │
│  └─────────────┘ │    │ └──────────────┘ │    │ └─────────────┘ │
│                 │    │                  │    │                 │
│  ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│  │ Eigen3      │ │    │ │ dqfun_cwrap  │ │◄──►│ │   DQFUN     │ │
│  │ Integration │ │    │ │              │ │    │ │ (dqmodule)  │ │
│  └─────────────┘ │    │ └──────────────┘ │    │ └─────────────┘ │
│                 │    │                  │    │                 │
│  ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│  │ CG Solver   │ │    │ │ ddfun_cwrap  │ │◄──►│ │   DDFUN     │ │
│  └─────────────┘ │    │ │              │ │    │ │ (ddmodule)  │ │
└─────────────────┘    │ └──────────────┘ │    │ └─────────────┘ │
                       └──────────────────┘    └─────────────────┘
```

## Layer Architecture

### 1. Bailey Library Layer (Fortran)

**Purpose**: Provides core high-precision arithmetic operations

**Components**:
- **DDFUN**: Double-double arithmetic using pairs of double-precision numbers
- **DQFUN**: Quad-double arithmetic using pairs of quad-precision numbers  
- **QXFUN**: Extended quad arithmetic using hardware quad-precision

**Data Structures**:
```fortran
! DDFUN
type dd_real
    sequence
    real(ddknd) ddr(2)  ! ddknd = selected_real_kind(15, 307)
end type

! DQFUN  
type dq_real
    sequence
    real(dqknd) dqr(2)  ! dqknd = selected_real_kind(33, 4931)
end type

! QXFUN
! Uses raw real(qxknd) values  ! qxknd = selected_real_kind(33, 4931)
```

### 2. Fortran Wrapper Layer

**Purpose**: Provides C-compatible interface to Bailey libraries

**Design Pattern**: Foreign Function Interface (FFI) with `bind(C)`

**Example Structure**:
```fortran
subroutine qx_add(a,b,c) bind(C,name="qxadd_")
    real(qxknd), intent(in)  :: a, b
    real(qxknd), intent(out) :: c
    c = a + b
end subroutine
```

**Key Functions**:
- Arithmetic: `add`, `sub`, `mul`, `div`
- Conversion: `fromdbl`, `tostr`
- Mathematical: `sqrt`

### 3. C++ Interface Layer

**Purpose**: Provides object-oriented C++ interface with operator overloading

**Core Design**:
```cpp
struct QuadDouble {
    double qd[4] = {0.0, 0.0, 0.0, 0.0};  // QX internal representation
    
    QuadDouble() = default;
    QuadDouble(double val);                 // Conversion constructor
};
```

**Operator Overloading**:
- Arithmetic: `+`, `-`, `*`, `/`
- Assignment: `+=`, `-=`, `*=`, `/=`
- Stream: `<<` for output

### 4. Linear Algebra Integration

**Purpose**: Integrates high-precision arithmetic with Eigen3 sparse matrices

**Eigen3 NumTraits Specialization**:
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
            ReadCost = 4,     // Cost relative to double
            AddCost = 32, 
            MulCost = 64 
        };
    };
}
```

## Data Flow

### Arithmetic Operation Flow

```
C++ Code: a + b
    ↓
QuadDouble::operator+
    ↓
qxadd_(a.qd, b.qd, result.qd)
    ↓
qxfun_cwrap.f90: qx_add
    ↓
Bailey QXFUN: native quad arithmetic
    ↓
Result returned through layers
```

### Matrix-Vector Multiplication Flow

```
Eigen3: y = A * x
    ↓
Eigen3 SparseDenseProduct
    ↓
QuadDouble arithmetic operations
    ↓ (for each non-zero element)
QuadDouble::operator* and operator+=
    ↓
Bailey library calls
    ↓
High-precision computation
```

## Memory Layout

### QuadDouble Memory Structure

```
QuadDouble object (32 bytes):
┌────────────┬────────────┬────────────┬────────────┐
│  qd[0]     │  qd[1]     │  qd[2]     │  qd[3]     │
│  8 bytes   │  8 bytes   │  8 bytes   │  8 bytes   │
└────────────┴────────────┴────────────┴────────────┘
```

For QXFUN, typically only `qd[0]` and `qd[1]` are used, but the full array ensures compatibility.

### Sparse Matrix Memory Usage

```cpp
SpMat_QD A(n, n);  // n×n sparse matrix

// Memory per non-zero element:
// - Value: 32 bytes (QuadDouble)
// - Indices: 8 bytes (2 × int32)
// Total: ~40 bytes per non-zero (vs ~12 bytes for double)
```

## Build System Architecture

### CMake Configuration

```cmake
# Import Bailey libraries as static libraries
foreach(name dd dq qx)
    add_library(${name}fun STATIC IMPORTED)
    set_target_properties(${name}fun PROPERTIES
        IMPORTED_LOCATION "${${NAME_UP}FUN_DIR}/lib${name}fun.a")
    
    add_library(${name}wrap STATIC IMPORTED)  
    set_target_properties(${name}wrap PROPERTIES
        IMPORTED_LOCATION "${${NAME_UP}FUN_DIR}/lib${name}wrap.a")
endforeach()
```

### Linking Strategy

```
sample_qx executable
├── qxwrap.a        (Fortran C interface)
├── qxfun.a         (Bailey QX library)
├── gfortran        (Fortran runtime)
└── Eigen3          (Header-only, no linking)
```

## Thread Safety Design

### Bailey Library Thread Safety

- **Guarantee**: All Bailey libraries are 100% thread-safe
- **Design**: No global state, all operations on local variables
- **Implication**: Parallel algorithms can safely use QuadDouble arithmetic

### Eigen3 Integration Thread Safety

- **Matrix Operations**: Thread-safe when operating on different matrices
- **Shared Matrices**: Require external synchronization for write operations
- **Read-Only Operations**: Multiple threads can safely read the same matrix

### Implementation Considerations

```cpp
// Safe: Different threads, different matrices
#pragma omp parallel for
for (int i = 0; i < num_systems; ++i) {
    conjugateGradient(A[i], b[i], x[i], max_iter, tol);
}

// Safe: Same matrix, read-only operations
#pragma omp parallel for  
for (int i = 0; i < n; ++i) {
    Vec_QD col_i = A.col(i);  // Read-only access
    process(col_i);
}
```

## Error Handling Strategy

### Layered Error Handling

1. **Bailey Library Level**: 
   - Fortran runtime errors (overflow, underflow, invalid operations)
   - No explicit error return codes

2. **Wrapper Level**:
   - Format string errors
   - Memory access issues

3. **C++ Level**:
   - Standard C++ exceptions from Eigen3
   - Logic errors (dimension mismatches)

4. **Application Level**:
   - Convergence failures
   - Input validation

### Error Propagation

```cpp
// Bailey library errors propagate as:
// Fortran runtime error → Program termination

// Eigen3 errors propagate as:
// C++ exception → Catch and handle in application

// Application-level validation:
if (A.rows() != A.cols()) {
    throw std::invalid_argument("Matrix must be square");
}
```

## Performance Architecture

### Computational Complexity

| Operation | Double | QuadDouble | Ratio |
|-----------|--------|------------|-------|
| Add/Sub | O(1) | O(1) | ~10-50x |
| Multiply | O(1) | O(1) | ~50-100x |
| Divide | O(1) | O(1) | ~100-200x |
| Matrix×Vector | O(nnz) | O(nnz) | ~50-100x |

### Memory Hierarchy Impact

```
Cache Usage:
- Double sparse matrix: ~12 bytes/element
- QuadDouble sparse matrix: ~40 bytes/element
- Cache misses increase due to larger data structures
- Algorithm complexity remains the same
```

### Optimization Strategies

1. **Algorithmic**: Use better algorithms rather than just higher precision
2. **Vectorization**: Let Eigen3 handle SIMD where possible
3. **Memory Access**: Minimize precision conversions
4. **Sparsity**: Exploit sparse structure to reduce computation

## Extensibility Design

### Adding New Precision Types

To add support for additional Bailey libraries:

1. **Add Fortran Wrapper**: Create `newfun_cwrap.f90`
2. **Update CMake**: Add library to `foreach(name dd dq qx new)`
3. **Create C++ Interface**: Define new precision struct
4. **Specialize Eigen3**: Add NumTraits specialization

### Adding New Operations

To add new mathematical functions:

1. **Fortran Layer**: Add function to appropriate wrapper
2. **C++ Declaration**: Add `extern "C"` declaration
3. **C++ Wrapper**: Add high-level C++ interface
4. **Documentation**: Update API reference

### Integration Points

The architecture provides clean integration points for:
- Additional Bailey libraries
- Alternative linear algebra backends
- Different solver algorithms
- Custom precision types