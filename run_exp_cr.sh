#!/usr/bin/env bash
# 使いたい行列名を並べる
matrices=(
  DNVS/thread
  FIDAP/ex10
)

precisions=(
  # "double"
  # "dd"
  "dq"
)

# docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  for p in "${precisions[@]}"; do
    make run NO_BUILD=1 -- cr_solver --matrix "$m" --precision "$p" --tol 1e-12 --max-iter 2.0 --export-mat "outputs/${m}_${p}.mat" --export-csv "outputs/runs.csv"
  done
done
