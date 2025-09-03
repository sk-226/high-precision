# Repository Guidelines

## Project Structure & Module Organization
- `src/`: CLI apps and small test executables (`cg_solver`, `experiments/`, `benchmarks/`).
- `include/`: C++ headers — `algorithms/`, `bailey/`, `io/`, `matrix_io/`, `linear_algebra/`.
- `interfaces/bailey_wrappers/`: Fortran–C bindings for Bailey DD/DQ/QX.
- `lib/fast_matrix_market/`: Submodule for Matrix Market I/O; keep vendor code untouched.
- `inputs/`: Matrix Market `.mtx` files. `outputs/`: generated `.mat` and logs.
- `doc/`, `third_party/`: documentation and licenses. Example: `run_exp.sh` batches solver runs.

## Build, Test, and Development Commands
- Docker (preferred, reproducible):
  - `make build`
  - `make run -- cg_solver --matrix nos5 --precision dq --tol 1e-15`
- Local CMake (advanced):
  - `export QXFUN_DIR=... DQFUN_DIR=... DDFUN_DIR=...`
  - `cmake -S . -B build -G Ninja -DQXFUN_DIR=$QXFUN_DIR -DDQFUN_DIR=$DQFUN_DIR -DDDFUN_DIR=$DDFUN_DIR`
  - `cmake --build build --config Release`
- Test binaries (after build): `./build/test_basic`, `./build/precision_validation_test`, `./build/check_ldbl`, `./build/prec_sanity`.

## Coding Style & Naming Conventions
- C++17/20, Fortran for wrappers; 4‑space indent, `#pragma once` in headers.
- Filenames: `snake_case.hpp/.cpp`; Types: UpperCamelCase (e.g., `CGResult`, `PrecisionTraits`).
- Namespaces: `bailey`, `algorithms`. Prefer Eigen types via `PrecisionTraits`; avoid raw `new`.
- Keep headers self‑contained; minimize changes under `lib/` and `interfaces/`.

## Testing Guidelines
- Place new micro‑tests in `src/benchmarks/` (name: `<feature>_test.cpp`).
- Use matrices from `inputs/`; don’t commit large data. Tests should return code 0 and print clear PASS/FAIL.
- Validate solver changes by running a double vs DQ/QX comparison on `nos5.mtx` with a fixed tolerance.

## Commit & Pull Request Guidelines
- Commit style: prefer Conventional Commits (`feat:`, `fix:`, `docs:`, `refactor:`). History example: `fix: digits`.
- PRs must include: rationale (what/why), commands used (copy‑pasteable), before/after logs, matrices/precision/tolerance, linked issues, and doc updates when flags or I/O change.

## Security & Configuration Tips
- Docker image enforces quad `long double`; use Docker for CI/local parity.
- Local builds require `QXFUN_DIR`, `DQFUN_DIR`, `DDFUN_DIR`; `matio-cpp` is auto‑detected (export control is optional).
- Write artifacts to `outputs/`. When using `make run`, include `--` before solver args.
