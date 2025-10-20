# DQFUN Test Harness

Documentation for the double-quad precision regression suite in `tests/dqfun`.

## Scope
- Verifies the 66-digit arithmetic and special functions exported by `dqmodule`.
- Mirrors the double-double flow but with tighter tolerances and higher-precision data tables.

## Prerequisites
- `gfortran` reachable on `PATH`.
- Run from the repository root.

## How to run
```
bash tests/gnu-dqfun-tests.scr
```
The wrapper drops into `tests/dqfun/fortran/`, compiles the library, and executes all listed programs. Output logs (`*.txt`) stay in that directory.

## Script flow
1. **Compile library** - `tests/dqfun/fortran/gnu-complib-dq.scr` runs `gfortran -O3 -c dqfuna.f90 dqfune.f90 dqmodule.f90 second.f90`.
2. **Link and execute** - `gnu-dqfun-tests.scr` invokes `gnu-complink-dq.scr <program>`, piping stdout into `<program>.txt`, and finishes each run with `tail -1` to show the summary line.

## Core regression: `testdqfun`
- Location: `tests/dqfun/fortran/testdqfun.f90`.
- Loads `testdqfun.ref.txt` and rewinds the unit before calculations (`tests/dqfun/fortran/testdqfun.f90:43`).
- Allocates high-precision constants (pi, log2, Bernoulli table, polylog cache) and sets up test vectors (`tests/dqfun/fortran/testdqfun.f90:50`).
- Iterates through arithmetic, transcendental, and comparison checks using `dqwrite` for logging and `checkdq`/`checkdqc` for verification (`tests/dqfun/fortran/testdqfun.f90:109`). Example call:
  ```fortran
  call checkdq (nfile, 1, t1 / t2, eps, err)
  ```
- Maintains a running maximum error; the expected footer in `testdqfun.txt` is:
  ```
  Max relative error =   4.248079D-62
  ALL TESTS PASSED
  ```
- Coverage highlights:
  - Mixed-mode arithmetic against native quad constants (division, exponentiation, etc.).
  - Complex-valued operations validated through `checkdqc` (`tests/dqfun/fortran/testdqfun.f90:88`).
  - Logical comparisons (`==`, `/=`, `<=`, `>=`) tested across `dq_real` and standard double variants (`tests/dqfun/fortran/testdqfun.f90:179`).

## Auxiliary programs
- `tpslqm1dq` (`tests/dqfun/fortran/tpslqm1dq.f90:1`): PSLQM1 demonstration at 70-digit precision for harder integer relations. Reports trial counters, precision settings, CPU timings, and the recovered relation.
- `tquaddq` (`tests/dqfun/fortran/tquaddq.f90:1`): Exercises tanh-sinh, exp-sinh, and sinh-sinh quadratures using `dq_real` to confirm numerical stability with wider exponent ranges.
- `tpphixdq` (`tests/dqfun/fortran/tpphixdq.f90:1`): Recomputes the Poisson phi experiment with 70-digit arithmetic, invoking PSLQM3 for relation recovery and logging CPU costs.

Each auxiliary program writes a `<name>.txt` transcript and should terminate with the same `ALL TESTS PASSED` line contained in the corresponding `.ref.txt` file.

## Result interpretation
Expect four final summary lines from the wrapper: one for each program executed. Any deviation from `ALL TESTS PASSED` indicates either a compiler/runtime failure or a numerical mismatch with the baseline.

## Troubleshooting tips
- Link failures often stem from lingering `.o` files compiled with a different compiler version; re-run `find tests/dqfun/fortran -name '*.o' -delete` before rebuilding.
- When the relative error exceeds tolerance, re-run `./testdqfun` manually in `tests/dqfun/fortran/` to compare the full log against `testdqfun.ref.txt`.
