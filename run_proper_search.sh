#!/usr/bin/env bash
# 使いたい行列名を並べる
matrices=(
  Boeing/bcsstk36
  Boeing/msc23052
  Cylshell/s3rmt3m1
  Cylshell/s3rmt3m3
  DNVS/thread
  FIDAP/ex10
  FIDAP/ex10hs
  FIDAP/ex13
  FIDAP/ex9
  HB/bcsstk11
  HB/bcsstk12
  HB/bcsstk20
)

precision="dd"

# docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  make run NO_BUILD=1 -- cg_solver_proper_search --matrix "$m" --precision "$precision" --tol 1e-12 --max-iter 2.0 --export-mat "outputs/${m}_${precision}.mat" --export-csv "outputs/runs.csv"
done
