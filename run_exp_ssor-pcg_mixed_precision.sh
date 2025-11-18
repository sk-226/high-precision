#!/usr/bin/env bash
# SSOR-PCG solver experiment script
# Note: For mixed precision version (preconditioner in double, main computation in high precision),
#       use ssor_pcg_solver_mixed_precision instead. Results will be labeled as "precision+double"
#       (e.g., "dq+double") in CSV and MAT outputs.
#
# 使いたい行列名を並べる
matrices=(
  DNVS/thread
)

precision="dd"
omega=(
  # 0.4
  # 0.5
  0.6
  # 0.7
  # 0.8
  # 0.9
)
# docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  for o in "${omega[@]}"; do
    make run NO_BUILD=1 -- ssor_pcg_solver_mixed_precision --matrix "$m" --precision "$precision" --tol 1e-12 --max-iter 2.0 --omega "$o" --export-mat "outputs/${m}_${precision}_${o}.mat" --export-csv "outputs/runs.csv"
  done
done
