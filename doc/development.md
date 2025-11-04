# Development Guide

This document covers development setup, workflow, and important considerations for working with the Bailey high-precision libraries.

## Development Setup

### Docker Development (Recommended)

The simplest workflow uses Make targets and the provided Dockerfile.

```bash
# Build the image and project
make build

# Run a solver from the image (note the --)
make run -- cg_solver --matrix nos5 --precision dq --tol 1e-15 --max-iter 2.0

# Quick re-run without rebuild (same binary)
make run NO_BUILD=1 -- cg_solver --matrix nos5 --precision dq --tol 1e-15 --max-iter 2.0
```

> [!NOTE]
> The container workdir is `/work`. Host `./inputs` is mounted to `/work/inputs` (read-only) and `./outputs` persists results.

> [!IMPORTANT]
> If you customize the base image or use a local toolchain, verify the floating‑point model before using DQ/QX: `cmake --build build && ./build/tests/check_ldbl`.
>
> Only `dd`, `dq`, and `qx` are supported. `qd` (quad-double) is not supported in this project.

### Local Development Setup

1. **Install Dependencies** (see [Dependencies](dependencies.md))

2. **Set Environment Variables**:
```bash
export DDFUN_DIR=/path/to/ddfun/fortran
export DQFUN_DIR=/path/to/dqfun/fortran
export QXFUN_DIR=/path/to/qxfun/fortran
```

3. **Build the Project**:
```bash
cmake -S . -B build -G Ninja \
    -DQXFUN_DIR=${QXFUN_DIR} \
    -DDQFUN_DIR=${DQFUN_DIR} \
    -DDDFUN_DIR=${DDFUN_DIR}

cmake --build build --config Release
```

## Project Structure

```
high-precision/
├── src/                    # C++ source code
│   └── cg_solver.cpp      # Main solver (CLI)
├── interfaces/bailey_wrappers/   # Fortran wrapper files
│   ├── ddfun_cwrap.f90           # DD wrapper
│   ├── dqfun_cwrap.f90           # DQ wrapper
│   └── qxfun_cwrap.f90           # QX wrapper
├── doc/                  # Documentation
├── CMakeLists.txt        # Build configuration
└── Dockerfile           # Container definition
```

## Build Process Details

### 1. Bailey Library Integration

The Dockerfile automatically:

```bash
# Download Bailey libraries
wget https://www.davidhbailey.com/dhbsoftware/{library}-{version}.tar.gz

# Extract and build
tar -xzf {library}-{version}.tar.gz
cd {library}-{version}/fortran
./gnu-complib-{prefix}.scr  # Build library
ar rcs lib{prefix}fun.a *.o  # Create static library
```

### 2. Fortran Wrapper Compilation

```bash
gfortran -c interfaces/bailey_wrappers/qxfun_cwrap.f90 -I${QXFUN_DIR}
ar rcs ${QXFUN_DIR}/libqxwrap.a qxfun_cwrap.o
```

### 3. C++ Application Build

```bash
# CMake configuration
cmake -S . -B build -G Ninja -DQXFUN_DIR=${QXFUN_DIR}

# Compilation and linking
ninja -C build
```

## Development Workflow

### Making Changes (containerized)

1. Start dev container and connect via SSH (see [Remote Development](remote_dev.md)).
2. Build inside the container: `cmake -S . -B build -G Ninja && cmake --build build -j`.
3. Run: `./build/cg_solver --matrix nos5 --precision dq --tol 1e-15`.

### Adding New Functions

To add a new high-precision function:

1. **Add Fortran Wrapper**:
```fortran
subroutine qx_newfunc(a, result) bind(C, name="qxnewfunc_")
    real(qxknd), intent(in)  :: a
    real(qxknd), intent(out) :: result
    result = newfunc(a)
end subroutine
```

2. **Declare in C++**:
```cpp
extern "C" {
    void qxnewfunc_(const long double* a, long double* result);
}
```

3. **Add C++ Wrapper**:
```cpp
bailey::QXNumber newfunc(const bailey::QXNumber& a) {
    bailey::QXNumber r;
    qxnewfunc_(&a.qx, &r.qx);
    return r;
}
```

### Code Style Guidelines

#### C++ Code

- Use `bailey::DDNumber` / `DQNumber` / `QXNumber` for high precision
- Prefer `bailey::PrecisionTraits<T>` for matrix/vector aliases
- Use `to_double()` only for logging or tolerance checks

#### Fortran Wrappers

- Use `bind(C)` for C compatibility
- Convert between Bailey types and C arrays
- Handle string formatting carefully (avoid 'E' suffix in format)
- Use appropriate precision constants (`qxknd`, `dqknd`, etc.)

## Important Considerations

### Memory Layout

- **DDFUN**: `double[2]`
- **DQFUN**: `long double[2]` (mapped from Fortran `real(dqknd)`)
- **QXFUN**: single `long double` (mapped from `real(qxknd)`)

### Precision Handling

```cpp
// Converting between precisions
double d = 1.5;
bailey::QXNumber qx(d);     // double -> QX
double back = to_double(qx); // QX -> double (for comparisons)
```

### Thread Safety

Bailey libraries are thread-safe, but:
- Each thread should have its own working variables
- Shared data structures need proper synchronization
- Eigen3 operations inherit thread safety from underlying data

### Performance Notes

- High-precision operations are significantly slower than double precision
- Use appropriate precision level for your accuracy requirements
- Consider algorithmic improvements over just increasing precision

## Testing

### Basic Functionality Test

```bash
# Run the default test
./build/sample_qx

# Expected output:
# Created test matrix. Size: 5x5
# Problem setup complete (b = A * x_true).
# Starting Conjugate Gradient solver...
# Converged after 1 iterations.
# Verification:
# Norm of (x_computed - x_true): [small number]
```

### Custom Tests

Modify `src/experiments/main.cpp` to test specific functionality:

```cpp
// Test basic arithmetic
bailey::QXNumber a(1.0);
bailey::QXNumber b(2.0);
auto c = a + b;
std::cout << "1 + 2 = " << c << std::endl;

// Test matrix operations
using T = bailey::QXNumber; // or DDNumber/DQNumber
using Traits = bailey::PrecisionTraits<T>;
Traits::matrix_type A = createTestMatrix<T>();
Traits::vector_type x = solveWithCG<T>(A, /*rhs*/);
```

## Debugging

### Common Issues

1. **Linking Errors**: Check that all Bailey libraries are built and accessible
2. **Precision Loss**: Verify using appropriate precision throughout the computation chain
3. **Format Errors**: Ensure Fortran format strings don't use 'E' exponent specification

### Debug Builds

```bash
cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug
cmake --build debug
gdb ./debug/sample_qx
```
