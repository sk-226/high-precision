# DDFUN Test Harness

Detailed notes on the vendor-provided double-double regression suite under `tests/ddfun`.

## Scope
- Validates scalar and complex arithmetic exported by `ddmodule` using `testddfun.f90`.
- Exercises high-level algorithms (`tpslqm1dd`, `tquaddd`, `tpphixdd`) to ensure the library integrates with PSLQ search, quadrature, and theta-based workflows.

## Prerequisites
- `gfortran` available on `PATH`.
- Working directory should be the repository root so relative paths resolve correctly.
- Optional: set `DDFUN_STRICT` to adjust compiler flags (`0` = `-O3 -fno-expensive-optimizations`, `1` adds strict FP guards, `2` uses `-O0`).

## How to run
```
bash tests/gnu-ddfun-tests.scr
```
This wrapper enters `tests/ddfun/fortran/`, compiles sources, and executes every regression in sequence. All generated logs (`*.txt`) appear next to their sources.

## Script flow
1. **Compile library** - `tests/ddfun/fortran/gnu-complib-dd.scr` runs `gfortran $CFLAGS -c ddfuna.f90 ddfune.f90 ddmodule.f90 second.f90` and prints the active flag block.
2. **Link and execute** - `gnu-ddfun-tests.scr` repeatedly calls `gnu-complink-dd.scr <program>` and pipes stdout into `<program>.txt`, then `tail -1` to display the pass/fail line.

## Core regression: `testddfun`
- Location: `tests/ddfun/fortran/testddfun.f90`.
- Opens the reference stream (`testddfun.ref.txt`) and rewinds it before execution (`tests/ddfun/fortran/testddfun.f90:30`).
- Seeds test data (pi, log2, Bernoulli numbers, complex constants) and writes them to the transcript (`tests/ddfun/fortran/testddfun.f90:37`).
- For each arithmetic form, prints the computed value and validates it against the reference file via `checkdd`/`checkddc` (`tests/ddfun/fortran/testddfun.f90:97`). Example call:
  ```fortran
  call checkdd (nfile, 14, t1 + t2, eps, err)
  ```
- Tracks a running maximum relative error; the final lines in `testddfun.txt` are:
  ```
  Max relative error =   9.159424D-32
  ALL TESTS PASSED
  ```
- Coverage highlights:
  - Real arithmetic: addition, subtraction, multiplication, division, exponentiation.
  - Mixed precision operations with native double values (`tests/ddfun/fortran/testddfun.f90:102`).
  - Complex arithmetic and transcendental functions through `ddwrite`/`checkddc` pairs (`tests/ddfun/fortran/testddfun.f90:76`).
  - Comparison operators (`==`, `/=`, `<=`, `>=`) and their consistency with hardware doubles (`tests/ddfun/fortran/testddfun.f90:167`).

## Auxiliary programs
- `tpslqm1dd` (`tests/ddfun/fortran/tpslqm1dd.f90:1`): builds an algebraic input vector for the PSLQM1 relation finder, reports CPU times, and prints the discovered relation. Input controls (e.g., `kq`, `kr`, `ks`, `ndp`) are defined in the parameter block and surfaced in the log header.
- `tquaddd` (`tests/ddfun/fortran/tquaddd.f90:1`): evaluates tanh-sinh, exp-sinh, and sinh-sinh quadratures for predefined integrands, verifying initialization helpers (`initqts`, `initqes`, `initqss`).
- `tpphixdd` (`tests/ddfun/fortran/tpphixdd.f90:1`): computes the Poisson phi function using theta functions and then launches PSLQM3. The output lists CPU timings (`Alpha`, `PSLQM1`, `Total`) and prints recovered relations in human-readable form.

Each auxiliary run finishes with a summary line from its `.txt` log; compare against the corresponding `.ref.txt` if discrepancies arise.

## Result interpretation
When every step succeeds, the wrapper prints four trailing `ALL TESTS PASSED` lines (one per program). Investigate any earlier compiler errors or non-zero exit codes by opening the generated log in the same folder.

## Troubleshooting tips
- Delete stale `.o` files if switching `DDFUN_STRICT` modes; the scripts do not force recompilation when the flag list changes.
- If comparisons fail, re-run `./testddfun` manually inside `tests/ddfun/fortran/` to inspect the full diff against `testddfun.ref.txt`.
