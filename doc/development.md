# Development Guide

This document covers development setup, workflow, and important considerations for working with the Bailey high-precision libraries.

## Development Setup

### Docker Development (Recommended)

The fastest way to get started:

```bash
# Clone and enter project directory
cd /path/to/high-precision

# Build development environment
docker build -t bailey-hp .

# Run development container
docker run -it bailey-hp /bin/bash

# Inside container - run the application
./build/sample_qx
```

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
│   └── main.cpp           # Main application
├── fortran/               # Fortran wrapper files
│   ├── ddfun_cwrap.f90   # DD wrapper
│   ├── dqfun_cwrap.f90   # DQ wrapper
│   └── qxfun_cwrap.f90   # QX wrapper
├── doc/                  # Documentation
├── CMakeLists.txt        # Build configuration
├── Dockerfile           # Container definition
└── CLAUDE.md           # AI assistant instructions
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
gfortran -c fortran/qxfun_cwrap.f90 -I${QXFUN_DIR}
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

### Making Changes

1. **Edit C++ Code**: Modify `src/main.cpp` or add new source files
2. **Update Fortran Wrappers**: Modify `fortran/*_cwrap.f90` if needed
3. **Rebuild**: `cmake --build build`
4. **Test**: `./build/sample_qx`

### Adding New Functions

To add a new high-precision function:

1. **Add Fortran Wrapper**:
```fortran
subroutine qx_newfunc(a, result) bind(C, name="qxnewfunc_")
    real(qxknd), intent(in)  :: a
    real(qxknd), intent(out) :: result
    result = newfunc(a)  ! Call Bailey library function
end subroutine
```

2. **Declare in C++**:
```cpp
extern "C" {
    void qxnewfunc_(const double* a, double* result);
}
```

3. **Add C++ Wrapper**:
```cpp
QuadDouble newfunc(const QuadDouble& a) {
    QuadDouble r;
    qxnewfunc_(a.qd, r.qd);
    return r;
}
```

### Code Style Guidelines

#### C++ Code

- Use `QuadDouble` struct for high-precision numbers
- Provide operator overloads for natural syntax
- Include Eigen3 NumTraits specialization for matrix operations
- Use `to_double()` for conversion when needed

#### Fortran Wrappers

- Use `bind(C)` for C compatibility
- Convert between Bailey types and C arrays
- Handle string formatting carefully (avoid 'E' suffix in format)
- Use appropriate precision constants (`qxknd`, `dqknd`, etc.)

## Important Considerations

### Memory Layout

- **DDFUN**: Uses `real(c_double)` arrays of size 2
- **DQFUN**: Uses `real(dqknd)` arrays of size 2  
- **QXFUN**: Uses single `real(qxknd)` values

### Precision Handling

```cpp
// Converting between precisions
double d = 1.5;
QuadDouble qd(d);           // double -> QuadDouble
double back = to_double(qd); // QuadDouble -> double (for comparisons)
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

Modify `src/main.cpp` to test specific functionality:

```cpp
// Test basic arithmetic
QuadDouble a(1.0);
QuadDouble b(2.0);
QuadDouble c = a + b;
std::cout << "1 + 2 = " << c << std::endl;

// Test matrix operations
SpMat_QD A = createTestMatrix();
Vec_QD x = solveWithCG(A, b);
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