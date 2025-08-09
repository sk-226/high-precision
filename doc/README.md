# Bailey High-Precision Numerical Analysis

This project implements high-precision numerical analysis using David H. Bailey's DD/DQ/QX libraries with Eigen3 for sparse matrix computations.

## Overview

This project bridges Fortran-based high-precision arithmetic libraries with C++ linear algebra, enabling quad-precision conjugate gradient (CG) method for sparse matrices.

### Precision Levels

- **DDFUN**: Double-double precision (~32 decimal digits)
- **DQFUN**: Quad-double precision (~64 decimal digits) 
- **QXFUN**: Extended quad precision (~128+ decimal digits)

The current implementation uses QXFUN for maximum precision.

## Quick Start

### Using make + Docker (Recommended)

```bash
# Build container + project
make build

# Run the main solver (note the -- separator)
make run -- cg_solver --matrix nos5 --precision dq --tol 1e-30 --max-iter 2.0
```

Important:
- Always include `--` after `make run` so that subsequent options are passed to the solver, not to make.
- This project requires that the runtime container's `long double` is IEEE 754 binary128 (quad). Our provided Dockerfile satisfies this and we verify it at build time.

### Using raw Docker

```bash
docker build -t bailey-hp .
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision dq --tol 1e-30 --max-iter 2.0
```

### Local Development

See [Dependencies](dependencies.md) for detailed setup instructions.

## Documentation Structure

- [Dependencies](dependencies.md) - External libraries and build requirements
- [Development](development.md) - Development setup and workflow
- [API Reference](api-reference.md) - Function documentation and usage
- [Architecture](architecture.md) - System design and implementation details
- [Examples](examples.md) - Practical examples and use cases
- [Troubleshooting](troubleshooting.md) - Common issues and solutions

## Key Features

- **High-Precision Arithmetic**: Up to 128+ decimal digits of precision
- **Sparse Matrix Support**: Efficient sparse linear algebra with Eigen3
- **Cross-Language Integration**: Fortran libraries accessible from C++
- **Thread-Safe**: Bailey libraries are 100% thread-safe
- **Containerized**: Complete Docker development environment

## Example Output

```
Created test matrix. Size: 5x5
Problem setup complete (b = A * x_true).

Starting Conjugate Gradient solver...
Converged after 1 iterations.

Verification:
Norm of (x_computed - x_true): 4.677071733467427E-1
```

## License

This project builds upon the following third-party software. We are grateful to their authors and maintainers:
- David H. Bailey's DD/DQ/QX libraries (see respective DISCLAIMER.txt for licensing terms)
- Eigen3 (C++)
- fast_matrix_market (C++)
- matio / matio-cpp (optional)