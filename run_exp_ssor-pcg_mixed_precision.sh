#!/usr/bin/env bash
# SSOR-PCG solver experiment script (mixed precision version)
# 
# This script uses ssor_pcg_solver_mixed_precision which supports configurable
# precision for both main computation and preconditioner.
# 
# Usage examples:
#   --precision dq --precond-precision double  (default: dq+double)
#   --precision dq --precond-precision dd      (dq+dd)
#   --precision qx --precond-precision dq     (qx+dq)
# 
# Results are labeled as "BASE+PRECOND" (e.g., "dq+dd") when different precisions
# are used, or just "BASE" when both are the same.
#
# 使いたい行列名を並べる
matrices=(
  DNVS/thread
)

precision="dq"
precond_precision="dd"
omega=(
  # 0.4
  # 0.5
  0.6
  # 0.64 # SSOR-PCG(double) で0.01刻みでomegaを変えたときの最適値
  # 0.7
  # 0.8
  # 0.9
)
# docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  for o in "${omega[@]}"; do
    make run NO_BUILD=1 -- ssor_pcg_solver_mixed_precision --matrix "$m" --precision "$precision" --precond-precision "$precond_precision" --tol 1e-12 --max-iter 2.0 --omega "$o" --export-mat "outputs/${m}_${precision}_mp_${precond_precision}_${o}.mat" --export-csv "outputs/runs.csv"
  done
done
