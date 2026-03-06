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
- Test binaries (after build): `./build/tests/test_basic`, `./build/tests/precision_validation_test`, `./build/tests/check_ldbl`, `./build/tests/prec_sanity`.

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

## Cursor Cloud specific instructions

### Architecture caveat (x86_64)
The Cloud VM runs on x86_64 where `long double` is 80-bit extended (not IEEE binary128). **DD precision works; DQ/QX precision produces incorrect results on x86_64.** This is a known limitation documented in the README. The `check_ldbl` test verifies this: `is_ieee_quad_binary128: false`. Test failures in QX/DQ arithmetic on x86_64 are expected and not bugs.

### Docker is required
Everything builds and runs inside Docker. Before any build/test/run, ensure Docker is running (`sudo dockerd &>/tmp/dockerd.log &` if needed, then `sudo chmod 666 /var/run/docker.sock`). The fuse-overlayfs storage driver is used in the Cloud VM (overlay2 is not supported by the kernel).

### Key commands (all from repo root)
- **Build**: `make build` (builds Docker image with all dependencies)
- **Run solver**: `make run NO_BUILD=1 -- cg_solver --matrix test1 --precision dd --tol 1e-20 --input-dir /work/inputs/test`
- **Run tests**: `make run NO_BUILD=1 -- tests/test_basic`, `make run NO_BUILD=1 -- tests/check_ldbl`, etc.
- Use `NO_BUILD=1` to skip rebuilding the Docker image when only running.
- Matrix files must follow the structure `{input_dir}/{matrix_name}/{matrix_name}.mtx`.

### Dockerfile workarounds for fuse-overlayfs
The Dockerfile requires `dnf --setopt=cachedir=/tmp/dnf-cache` for secondary `dnf install` steps because fuse-overlayfs cannot handle EPEL metadata directory renames across layers. This is already applied in the current Dockerfile.

### Build fixes for x86_64
- `libquadmath` is added to link libraries on x86_64 (CMakeLists.txt conditional).
- Fortran standalone programs (`dq_eps_and_sqrt2`, `pi_check`, `qx_io_demo`) are skipped on x86_64 due to REAL(16)/REAL(10) type mismatch.
- `cg_alt_solver` target is commented out (source file not committed to repo).
