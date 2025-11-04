# QXFUN Test Harness

Notes for the quad-extended special-function suite in `tests/qxfun`.

## Scope
- Validates the QXFUN special functions (Bessel, error, gamma families, Hurwitz zeta, polylogarithms, hypergeometric, etc.).
- Confirms that the supporting Bernoulli/polylog tables and PSLQ utilities function correctly at roughly 60-digit precision.

## Prerequisites
- `gfortran` available on `PATH`.
- Run commands from the repository root.

## How to run
```
bash tests/gnu-qxfun-tests.scr
```
The wrapper changes into `tests/qxfun/fortran/`, compiles the library, then executes each regression while logging stdout to `<program>.txt`.

## Script flow
1. **Compile library** - `tests/qxfun/fortran/gnu-complib-qx.scr` runs `gfortran -O3 -c qxfune.f90 qxmodule.f90 second.f90`.
2. **Link and execute** - `gnu-qxfun-tests.scr` uses `gnu-complink-qx.scr <program>` to build binaries, records output, and shows the last log line via `tail -1`.

## Core regression: `testqxfun`
- Location: `tests/qxfun/fortran/testqxfun.f90`.
- Opens the reference file `testqxfun.ref.txt` (`tests/qxfun/fortran/testqxfun.f90:32`) and announces the run.
- Initializes Bernoulli and polylog tables, then evaluates each exported function with representative arguments.
- After every call, prints the result and verifies it through `checkqx` (`tests/qxfun/fortran/testqxfun.f90:48`). Example block:
  ```fortran
  t3 = bessel_k (pi, log2)
  call checkqx (nfile, 1, t3, err)
  ```
- Maintains the worst relative error; the reference footer in `testqxfun.txt` is:
  ```
  Max relative error =   3.729413D-34
  ALL TESTS PASSED
  ```
- Function coverage includes:
  - Modified Bessel functions `bessel_i`, `bessel_k`, plus their `n`-th order variants for small and large arguments.
  - Ordinary Bessel `bessel_j`, `bessel_y`, `bessel_jn`, `bessel_yn`.
  - Special functions `digamma`, `erf`, `erfc`, `expint`, `gamma`, `incgamma` for positive and negative domain cases.
  - Hurwitz zeta, polygamma, hypergeometric `pfq`, and polylogarithm evaluations.

## Auxiliary programs
- `tpslqm1qx` (`tests/qxfun/fortran/tpslqm1qx.f90:1`): PSLQM1 harness tuned for QXFUN's real type. Outputs precision settings, CPU timings, and the discovered relation.
- `tquadqx` (`tests/qxfun/fortran/tquadqx.f90:1`): Runs the tanh-sinh, exp-sinh, and sinh-sinh quadrature regressions with QXFUN kernels.
- `tpphixqx` (`tests/qxfun/fortran/tpphixqx.f90:1`): Re-evaluates the Poisson phi problem at QXFUN precision, printing theta ratios, PSLQM results, and CPU statistics.

Each program writes a `.txt` transcript and should end with the `ALL TESTS PASSED` line that mirrors its `.ref.txt` baseline.

## Result interpretation
A successful wrapper run prints four terminal `ALL TESTS PASSED` summaries. Any mismatch signals either runtime failure or a divergence from the vendor reference data.

## Troubleshooting tips
- Ensure no `.o` files remain from prior compiler runs when switching Fortran versions; delete and re-run the wrapper.
- For detailed diffs, execute the binary directly in `tests/qxfun/fortran/` and compare the generated `.txt` with the corresponding `.ref.txt` file.
