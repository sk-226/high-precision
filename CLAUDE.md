# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a high-precision numerical analysis project that implements the Conjugate Gradient (CG) method with multiple precision levels. The project supports standard double precision and Bailey's high-precision arithmetic libraries (DD/DQ/QX), bridging Fortran-based high-precision arithmetic with C++ linear algebra using Eigen3.

## Build System

The project uses CMake with the following key commands:

```bash
# Configure build (requires environment variables for Bailey libraries)
cmake -S . -B build -G Ninja \
    -DQXFUN_DIR=${QXFUN_DIR} -DDQFUN_DIR=${DQFUN_DIR} -DDDFUN_DIR=${DDFUN_DIR}

# Build the project
cmake --build build --config Release

# Run the main CG solver
./build/cg_solver --matrix nos5 --precision qx --tol 1e-15

# Export convergence data to MATLAB
./build/cg_solver --matrix nos5 --precision dq --export-mat convergence.mat

# Run specific experiments
./build/nos5_cg
./build/nos7_cg

# Test basic arithmetic functionality
./build/tests/test_basic

# Sample programs
./build/sample_qx
```

Required environment variables:
- `QXFUN_DIR`: Path to QX function library directory
- `DQFUN_DIR`: Path to DQ function library directory  
- `DDFUN_DIR`: Path to DD function library directory

## Docker Development

The project includes a complete Docker environment:

```bash
# Build the Docker image
docker build -t bailey-hp .

# Run the container
docker run -it bailey-hp
```

The Dockerfile automatically downloads and builds Bailey's DD/DQ/QX libraries from David H. Bailey's website and sets up the complete development environment.

## Architecture

### Core Components

1. **Bailey Library Integration** (`interfaces/bailey_wrappers/`): Fortran wrapper modules that expose Bailey's high-precision arithmetic to C++
   - `ddfun_cwrap.f90`: Double-double precision wrapper (~30 digits)
   - `dqfun_cwrap.f90`: Quad-double precision wrapper (~66 digits)  
   - `qxfun_cwrap.f90`: Extended quad precision wrapper (~33 digits)

2. **Precision Type System** (`include/bailey/`): 
   - `precision_traits.hpp`: Template-based precision type system
   - `dd_arithmetic.hpp`, `dq_arithmetic.hpp`: High-precision arithmetic types
   - `quad_double.hpp`: Legacy QX precision support

3. **Core Algorithms** (`include/algorithms/`):
   - `conjugate_gradient.hpp`: Template-based CG implementation
   - Supports all precision levels through unified interface

4. **Matrix I/O and Export** (`include/io/`):
   - `matrix_market_reader.hpp`: Matrix Market reader (uses fast_matrix_market)
   - `mat_exporter.hpp`: MATLAB .mat file export for convergence data

5. **Command Line Interface** (`src/cg_solver.cpp`):
   - Unified solver supporting all precision levels
   - Performance benchmarking and result output
   - Optional MATLAB convergence data export

### Precision Levels

The project supports four precision levels:
- **double**: Standard IEEE 754 double precision (~15 decimal digits)
- **DD**: Double-double (~30 decimal digits) - Bailey's DDFUN library
- **DQ**: Quad-double (~66 decimal digits) - Bailey's DQFUN library
- **QX**: Extended quad precision (~33 decimal digits) - Bailey's QXFUN library

Usage examples:
```bash
./build/cg_solver --matrix nos5 --precision double --tol 1e-10
./build/cg_solver --matrix nos5 --precision dd --tol 1e-12  
./build/cg_solver --matrix nos5 --precision dq --tol 1e-15
./build/cg_solver --matrix nos5 --precision qx --tol 1e-20

# Export convergence data for analysis
./build/cg_solver --matrix nos5 --precision dq --export-mat convergence.mat
```

### Key Design Patterns

- **Template-Based Precision System**: Single CG implementation works with all precision types through `PrecisionTraits` specialization
- **Foreign Function Interface**: Fortran-to-C bindings with `bind(C)` for Bailey library interoperability
- **Operator Overloading**: Natural C++ arithmetic syntax for high-precision types
- **Eigen3 Integration**: NumTraits specialization enables sparse matrix operations
- **Unified Command Interface**: Single executable supports all precision levels

### Algorithm Implementation

The Conjugate Gradient implementation follows standard mathematical formulation:
- Uses relative residual convergence criterion: `||r||₂ / ||b||₂ < tolerance`
- Maintains numerical accuracy across all precision levels
- Provides comprehensive convergence history and timing measurements
- Compatible with Matrix Market sparse matrix format

## Dependencies

- CMake 3.21+
- Eigen3 3.4+
- Fortran compiler (gfortran)
- C++17/C++20 compiler
- Bailey's DD/DQ/QX libraries (auto-downloaded in Docker)
- libquadmath (for extended precision support)
- fast_matrix_market (for high-performance Matrix Market I/O)
- matio-cpp (optional, for MATLAB .mat export, auto-downloaded in Docker)

## Performance Characteristics

Relative performance compared to double precision:
- **double**: Baseline performance
- **DD**: 2-5x slower (~30 decimal digits)
- **DQ**: 10-20x slower (~66 decimal digits) 
- **QX**: 5-10x slower (~33 decimal digits)

## Test Matrices

The project includes several Matrix Market test cases:
- `nos5.mtx`: 468×468 symmetric matrix
- `nos7.mtx`: Larger test case
- `wathen120.mtx`: Additional validation matrix

## MATLAB Export

The solver can export convergence data to MATLAB .mat files:

```bash
./build/cg_solver --matrix nos5 --precision dq --export-mat results.mat
```

**Important**: Exported convergence data is stored in **double precision** format (~15 digits), regardless of the computation precision level. High-precision convergence values are converted to double during export.

Output files are automatically placed in the `outputs/` directory unless a path is explicitly specified.
