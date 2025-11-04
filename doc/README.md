# Bailey High-Precision Numerical Analysis

This project implements high-precision numerical analysis using David H. Bailey's DD/DQ/QX libraries with Eigen3 for sparse matrix computations.

> [!IMPORTANT]
> Supported precisions are only: `dd` (double-double), `dq` (double-quad = 2×quad), and `qx` (quad/extended). The project does not use or provide `qd` (quad-double, 4×double). Use the exact names `dd`, `dq`, or `qx` for `--precision`.

## Overview

This project bridges Fortran-based high-precision arithmetic libraries with C++ linear algebra, enabling quad-precision conjugate gradient (CG) method for sparse matrices.

> [!IMPORTANT]
> Platform support: arm64 Linux (container) — DD/DQ/QX supported (`long double` ≈ IEEE‑754 binary128). x86_64 glibc — DQ/QX unsupported (`long double` is 80‑bit extended); use DD or run inside the provided arm64 container.

> [!TIP]
> Quick check inside the image: `make run NO_BUILD=1 -- tests/check_ldbl`.

> [!NOTE]
> Future work: plan to support DQ/QX on x86_64 via a portable REAL(16) interop (e.g., _Float128/libquadmath) rather than relying on `long double`.

### Precision Levels

- **DD (DDFUN)**: double-double (~30 decimal digits)
- **DQ (DQFUN)**: double-quad = 2×quad (~66 decimal digits)
- **QX (QXFUN)**: single quad/extended (~33 decimal digits)

The solver supports all three. Pick via `--precision dd|dq|qx`.

## Quick Start

### Using make + Docker (Recommended)

```bash
# Build container + project
make build

# Run the main solver (note the -- separator)
make run -- cg_solver --matrix nos5 --precision dq --tol 1e-30 --max-iter 2.0
```

> [!IMPORTANT]
> Always include `--` after `make run` so subsequent flags go to the solver, not to make.

> [!NOTE]
> The container workdir is `/work`. Host `./inputs` is mounted read-only to `/work/inputs`, and outputs persist to `./outputs`.

> [!TIP]
> For quick re-runs when code hasn’t changed: `make run NO_BUILD=1 -- ...` (skips image rebuild; executes the same binary as `make run`).

### Using raw Docker

```bash
docker build -t bailey-hp .
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision dq --tol 1e-30 --max-iter 2.0
```

### Local Development

See [Dependencies](dependencies.md) for detailed setup instructions, and [Remote Development over SSH](remote_dev.md) for an editor-agnostic SSH workflow.

## Documentation Structure

- [Dependencies](dependencies.md) - External libraries and build requirements
- [Development](development.md) - Development setup and workflow
- [API Reference](api-reference.md) - Function documentation and usage
- [Architecture](architecture.md) - System design and implementation details
- [Examples](examples.md) - Practical examples and use cases
- [Troubleshooting](troubleshooting.md) - Common issues and solutions

## Key Features

- **High-Precision Arithmetic**: Up to ~66 decimal digits (DQ) and ~33 digits (QX)
- **Sparse Matrix Support**: Efficient sparse linear algebra with Eigen3
- **Cross-Language Integration**: Fortran libraries accessible from C++
- **Concurrency-Friendly**: Our wrappers are stateless; using separate data per thread is safe. Coordinate access when sharing matrices/vectors.
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
