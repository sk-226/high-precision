# High-Precision Conjugate Gradient Solver

A template-based implementation of the Conjugate Gradient method supporting multiple arithmetic precision levels, from standard IEEE 754 double precision to Bailey's extended quad precision (~128 digits).

## Features

- **Multiple Precision Support**: double, DD (~32 digits), DQ (~64 digits), QX (~128 digits)
- **Template-Based Design**: Single algorithm implementation works across all precision types
- **Comprehensive Metrics**: Tracks convergence history, timing, and error analysis
- **Matrix Market Format**: Supports standard sparse matrix input files
- **Performance Benchmarking**: Fair timing comparisons across precision levels

## Usage (Must Read)

Run the solver inside Docker using make. Note the `--` separator after `make run`.

```bash
# Build container + project
make build

# Run (note the -- separator)
make run -- cg_solver --matrix nos5 --precision dq --tol 1e-30 --max-iter 2.0

# Other examples
make run -- cg_solver --matrix nos5 --precision double --tol 1e-10
make run -- cg_solver --matrix nos5 --precision dd     --tol 1e-12
make run -- cg_solver --matrix nos5 --precision dq     --tol 1e-15
make run -- cg_solver --matrix nos5 --precision qx     --tol 1e-20
```

Important:
- You must include the `--` separator so that make does not try to parse solver options as make options.
- This repository assumes that the runtime container's `long double` is IEEE 754 binary128 (quad). Our Docker image satisfies this (verified in CI); using a different base image or host toolchain where `long double` is not quad is unsupported.

### Raw Docker (alternative)

```bash
# Build the Docker environment with Bailey libraries
docker build -t bailey-hp .

# Run with different precision levels
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision double --tol 1e-10
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision dd     --tol 1e-12
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision dq     --tol 1e-15
docker run --rm bailey-hp /work/build/cg_solver --matrix nos5 --precision qx     --tol 1e-20
```

### Local Build (advanced)

```bash
# Set up Bailey library paths
export QXFUN_DIR=/path/to/qxfun/fortran
export DQFUN_DIR=/path/to/dqfun/fortran  
export DDFUN_DIR=/path/to/ddfun/fortran

# Configure and build
cmake -S . -B build -G Ninja \
    -DQXFUN_DIR=${QXFUN_DIR} -DDQFUN_DIR=${DQFUN_DIR} -DDDFUN_DIR=${DDFUN_DIR}
cmake --build build --config Release

# Run solver
./build/cg_solver --matrix nos5 --precision qx --tol 1e-15
```

## Command Line Usage

```bash
./build/cg_solver [OPTIONS]

Options:
  --matrix NAME         Matrix name (e.g., nos5 for nos5.mtx)
  --precision LEVEL     Precision: double, dd, dq, qx (default: qx)
  --tol VALUE           Convergence tolerance (default: 1.0e-12)
  --max-iter VALUE      Max iterations: integer or coefficient*size (default: 2.0)
  --input-dir PATH      Input directory (default: /work/inputs)
  --export-mat FILE     Export convergence data to MATLAB .mat file
  --help, -h            Show help message

Examples:
  ./build/cg_solver --matrix nos5 --precision double --tol 1e-10
  ./build/cg_solver --matrix nos7 --precision dq --max-iter 1000
  ./build/cg_solver --matrix test --precision qx --max-iter 2.5
  ./build/cg_solver --matrix nos5 --precision dq --export-mat convergence.mat
```

## Precision Levels

| Precision | Library | Decimal Digits | Use Case |
|-----------|---------|----------------|----------|
| `double`  | IEEE 754 | ~15 | Performance baseline |
| `dq`      | Bailey DQFUN | ~66 | High precision |
| `dd`      | Bailey DDFUN | ~30 | Extended precision |
| `qx`      | Bailey QXFUN | ~33 | Maximum precision |

## Output Format

### Console Output

```
==========================
Numerical Results.
Problem: nos5.mtx
Precision: Double (15 digits)
==========================
Converged! (iter = 459)
# Iter.: 459
Time[s]: 0.044
Relres_2norm = 9.23e-11
True_Relres_2norm = 9.23e-11
Relerr_2norm = 8.45e-09
Relerr_Anorm = 7.17e-10
==========================
```

### MATLAB Export

Temporarily omitted here. See `doc/` for the latest export guidance.

## Algorithm Details

The implementation follows the standard Conjugate Gradient formulation:

1. **Convergence Criterion**: Relative residual norm `||r||₂ / ||b||₂ < tolerance`
2. **Error Analysis**: Tracks 2-norm and A-norm relative errors against known solution
3. **Timing Measurement**: Wall-clock time excluding I/O and initialization
4. **Residual Verification**: Computes true residual `||b - Ax||₂` for numerical verification

## Project Structure

```
├── include/
│   ├── algorithms/           # CG algorithm implementation
│   ├── bailey/              # Precision traits and arithmetic types  
│   └── io/                  # Matrix Market I/O and MATLAB export
├── interfaces/
│   └── bailey_wrappers/     # Fortran-C bindings for Bailey libraries
├── src/                     # Command-line applications
├── inputs/                  # Test matrices (Matrix Market format)
├── outputs/                 # MATLAB export files (auto-created)
└── lib/                     # External dependencies
```

## Dependencies

- **CMake 3.21+**: Build system
- **Eigen3 3.4+**: Sparse linear algebra
- **C++17 Compiler**: Language support
- **Fortran Compiler**: For Bailey library wrappers
- **Bailey Libraries**: DD/DQ/QX high-precision arithmetic (auto-downloaded in Docker)
- **matio-cpp** (optional): MATLAB .mat file export (auto-downloaded in Docker)

## Acknowledgements

This project builds upon and gratefully acknowledges the following works:
- David H. Bailey — DDFUN/DQFUN/QXFUN high-precision libraries (`davidhbailey.com`)
- Eigen3 — high-quality C++ template library for linear algebra
- matio / matio-cpp — MATLAB .mat I/O libraries used for optional export support
- fast_matrix_market — performant Matrix Market I/O library bundled as a submodule

We thank the authors and maintainers of these projects for making their software available to the community.

### Licensing Notes

- Bailey libraries (DD/DQ/QX) are distributed under a Limited BSD-like license by David H. Bailey. See `third_party/licenses/bailey-LICENSE.txt` bundled in this repository. When redistributing binaries, ensure the required copyright, conditions, and disclaimer
  are preserved in accompanying documentation/materials.
- fast_matrix_market, Eigen3, matio/matio-cpp and other bundled dependencies retain their own licenses; please review their LICENSE files in their respective repositories/submodules.

## Matrix Format

Input matrices should be in Matrix Market coordinate format (`.mtx`):

```
%%MatrixMarket matrix coordinate real symmetric
468 468 5172
1 1 2.0e+00
1 2 -1.0e+00
...
```

## Performance Notes

- **Double precision**: Fastest execution, ~15 digit accuracy
- **DD precision**: ~2-5x slower than double, ~30 digit accuracy  
- **DQ precision**: ~10-20x slower than double, ~64 digit accuracy
- **QX precision**: ~5-10x slower than double, ~33 digit accuracy

Higher precision levels may converge in fewer iterations due to reduced round-off error accumulation.

## Contributing

This project implements research-grade numerical algorithms. When contributing:

1. Maintain numerical accuracy across all precision levels
2. Preserve algorithm timing measurement integrity  
3. Follow existing code documentation standards
4. Test across all supported precision types

## References

- David H. Bailey's High-Precision Software: https://www.davidhbailey.com/dhbsoftware/
- Matrix Market Exchange: https://math.nist.gov/MatrixMarket/
- Eigen3 Documentation: https://eigen.tuxfamily.org/
- https://github.com/tbeu/matio
- https://github.com/alandefreitas/matplotplusplus
