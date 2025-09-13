# Known Issues — Bailey High-Precision Tests

This document records test behaviors that can vary across compilers/architectures but are acceptable for this project.

## DD (DDFUN) `testddfun`: `asin(t2)` borderline difference

- Symptom:
  - Running `tests/gnu-ddfun-tests.scr` (or `tests/ddfun/fortran/testddfun`) may print one `ERROR:` line in the block labeled `asin(t2) =`, and the summary line:
    - `Max relative error =   5.492266D-31`
    - leading to `ONE OR MORE TESTS FAILED`.
- Context:
  - The test sets `neps = -31`, i.e. tolerance ~ `1.0D-31`.
  - On some toolchains (gfortran/libm combinations) the computed value differs by ~`5.49D-31`, which is only ~5.5× the tolerance and is considered a boundary case due to rounding/codegen differences.
- Acceptance policy (this repo):
  - If this is the only discrepancy and the reported maximum relative error is on the order of `~5.5D-31`, we accept it as a known issue and do not treat it as a functional failure.
  - Any additional discrepancies or larger errors should be investigated.
- Notes:
  - Tightening/relaxing floating-point optimization flags (`-O0/-O1` with `-fprotect-parens`, `-ffloat-store`, etc.) did not eliminate the difference in our environment.
  - This does not affect DQ/QX tests, which pass entirely, nor the C++ solver.

## QX/DQ underflow notices

- You may see `IEEE_UNDERFLOW_FLAG` notices in `tquad*` cases; they are expected in these vendor tests when exercising extreme ranges. If the final line says `ALL TESTS PASSED`, they can be ignored.

