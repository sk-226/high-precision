# Architecture

This document describes the system architecture and design decisions for the Bailey high-precision numerical analysis project.

## System Overview

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│ C++ Application │    │ Fortran Wrappers │    │ Bailey Libraries│
│                 │    │                  │    │                 │
│ ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│ │ Precision T │ │◄──►│ │ qxfun_cwrap  │ │◄──►│ │    QXFUN    │ │
│ │   Struct    │ │    │ │              │ │    │ │ (qxmodule)  │ │
│ └─────────────┘ │    │ └──────────────┘ │    │ └─────────────┘ │
│                 │    │                  │    │                 │
│ ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│ │   Eigen3    │ │    │ │ dqfun_cwrap  │ │◄──►│ │    DQFUN    │ │
│ │ Integration │ │    │ │              │ │    │ │ (dqmodule)  │ │
│ └─────────────┘ │    │ └──────────────┘ │    │ └─────────────┘ │
│                 │    │                  │    │                 │
│ ┌─────────────┐ │    │ ┌──────────────┐ │    │ ┌─────────────┐ │
│ │  CG Solver  │ │    │ │ ddfun_cwrap  │ │◄──►│ │    DDFUN    │ │
│ └─────────────┘ │    │ │              │ │    │ │ (ddmodule)  │ │
└─────────────────┘    │ └──────────────┘ │    │ └─────────────┘ │
                       └──────────────────┘    └─────────────────┘
```

## Layer Architecture

### 1. Bailey Library Layer (Fortran)

**Purpose**: Provides core high-precision arithmetic operations

**Components**:
- **DDFUN (DD)**: double-double arithmetic using pairs of double-precision numbers
- **DQFUN (DQ)**: double-quad arithmetic using pairs of quad-precision numbers (2×`real(dqknd)`)  
- **QXFUN (QX)**: quad/extended arithmetic using a single `real(qxknd)` value

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

**Core Design (actual types)**:
```cpp
namespace bailey {
  struct DDNumber { double      dd[2];      };  // DDFUN
  struct DQNumber { long double dq[2];      };  // DQFUN
  struct QXNumber { long double qx;         };  // QXFUN
}
```

**Operator Overloading**:
- Arithmetic: `+`, `-`, `*`, `/`
- Assignment: `+=`, `-=`, `*=`, `/=`
- Stream: `<<` for output

### 4. Linear Algebra Integration

**Purpose**: Integrates high-precision arithmetic with Eigen3 sparse matrices

Each type has a matching `Eigen::NumTraits<>` specialization (see headers) enabling sparse/dense operations. Algorithms use `bailey::PrecisionTraits<T>` to select types at compile time.

## Data Flow

### Arithmetic Operation Flow

```
C++ Code: a + b
    ↓
operator+(T)
    ↓
qxadd_/dqadd_/ddadd_
    ↓
Fortran C wrappers: `qxfun_cwrap`, `dqfun_cwrap`, `ddfun_cwrap`
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
High-precision arithmetic operations
    ↓ (for each non-zero element)
operator* and operator+=
    ↓
Bailey library calls
    ↓
High-precision computation
```

## Memory Layout

### Precision Memory Structures (conceptual)

- DD: 2×`double`
- DQ: 2×`long double`
- QX: 1×`long double`

### Sparse Matrix Memory Usage

```cpp
// Example (arm64): value sizes — DD≈16 B, DQ≈32 B, QX≈16 B
// Indices unchanged vs double. Sparse memory scales mainly with scalar size.
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
cg_solver executable
├── ddwrap.a dqwrap.a qxwrap.a
├── ddfun.a dqfun.a qxfun.a
├── gfortran (Fortran runtime)
└── Eigen3 (header-only)
```

## Thread Safety Design

### Bailey Library Thread Safety

- Wrappers are stateless and operate on caller-provided values.
- Use separate variables per thread; synchronize when sharing matrices/vectors.
- We have no evidence of global mutable state in these paths; however, no formal thread-safety guarantee is made by this project.

### Eigen3 Integration Thread Safety

- **Matrix Operations**: Thread-safe when operating on different matrices
- **Shared Matrices**: Require external synchronization for write operations
- **Read-Only Operations**: Multiple threads can safely read the same matrix

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

| Operation | Double | DD/DQ/QX | Ratio (rough) |
|-----------|--------|-----------|---------------|
| Add/Sub | O(1) | O(1) | 10–50× |
| Multiply | O(1) | O(1) | 50–100× |
| Divide | O(1) | O(1) | 100–200× |
| Matrix×Vector | O(nnz) | O(nnz) | 50–100× |

### Memory Hierarchy Impact

```
Cache Usage:
- Double sparse matrix: ~12 bytes/element
- High-precision sparse matrix: scales with scalar size (arm64: DD≈20, QX≈20, DQ≈36 bytes/element incl. indices (32-bit index); rough)
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
