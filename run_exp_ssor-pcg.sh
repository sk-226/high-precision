#!/usr/bin/env bash
# 使いたい行列名を並べる
matrices=(
  DNVS/thread
)

precision="double"
omega=(
  0.4
  0.5
  # 0.6
  # 0.7
  # 0.8
  # 0.9
)
# docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  for o in "${omega[@]}"; do
    make run NO_BUILD=1 -- ssor_pcg_solver --matrix "$m" --precision "$precision" --tol 1e-12 --max-iter 2.0 --omega "$o" --export-mat "outputs/${m}_${precision}_${o}.mat" --export-csv "outputs/runs.csv"
  done
done
