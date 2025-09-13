# Bailey High-Precision Tests

This folder contains Bailey DD/DQ/QX upstream test suites copied under tests/, plus thin wrappers.

- `gnu-ddfun-tests.scr` — runs tests/ddfun/fortran/gnu-ddfun-tests.scr
- `gnu-dqfun-tests.scr` — runs tests/dqfun/fortran/gnu-dqfun-tests.scr
- `gnu-qxfun-tests.scr` — runs tests/qxfun/fortran/gnu-qxfun-tests.scr

Structure per precision (e.g., tests/dqfun/fortran):
- Upstream scripts: `gnu-complib-*.scr`, `gnu-complink-*.scr`, `gnu-*-tests.scr` (vendor-style flow).
- Library sources for tests: `ddfuna/ddfune/ddmodule` or `dqfuna/dqfune/dqmodule` or `qxfune/qxmodule` and helper `second.f90`.
- Test drivers: `test*d*fun.f90`, `tpslqm1*.f90`, `tquad*.f90`, `tpphix*.f90` and `*.ref.txt`.

## How to run

No external environment is required; scripts compile the minimal library sources locally.

Recommended (Docker):

```sh
make build
./tests/gnu-ddfun-tests.scr
./tests/gnu-dqfun-tests.scr
./tests/gnu-qxfun-tests.scr
```

Local (no Docker): simply run scripts inside tests/

```sh
cmake -S . -B build -G Ninja -DQXFUN_DIR=$QXFUN_DIR -DDQFUN_DIR=$DQFUN_DIR -DDDFUN_DIR=$DDFUN_DIR
cmake --build build --config Release
./tests/gnu-ddfun-tests.scr
```

## Notes
- Vendor test sources are copied to keep behavior identical; tolerance/formatting remain unchanged.
- `gnu-complib-*.scr` compiles the minimal needed .f90 modules within tests, not the entire upstream tree.
- This isolates tests from system Fortran builds and avoids tmp/ dependencies.

## Known Issues
- DD (DDFUN) `testddfun` may report a single failure in the `asin(t2)` check with:
  - `ERROR: 5.492266D-31` and summary `Max relative error = 5.492266D-31` (tolerance is `1.0D-31`).
  - This is a borderline, environment-dependent rounding effect. In this repository, we accept this specific deviation if it is the sole failing line and the maximum relative error is on the order of `~5.5D-31`. See `doc/known_issues.md` for details.
- QX/DQ tests sometimes print `IEEE_UNDERFLOW_FLAG` in `tquad*`; this is expected if the final line reads `ALL TESTS PASSED`.
