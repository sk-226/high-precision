#!/usr/bin/env bash
# 使いたい行列名を並べる
matrices=(
  # diag+ssor で解けなかった行列
  # LF10000
  # LFAT5000
  # bcsstk19
  # bcsstk20
  # bloweybq
  # ex10
  # ex10hs
  # ex13
  # ex3
  # ex33
  # ex5
  # nos2
  # plat1919
  # plat362
  # thread # 計算重すぎて一旦 pass
  
  # computational fluid dynamics problem
  bcsstk13
  # ex10 # not converged
  # ex10hs # not converged
  # ex13 # not converged
  ex15
  # ex3 # not converged
  # ex33 # not converged
  # ex5 # not converged
  ex9
  Pres_Poisson
)

docker build -t bailey-hp .

for m in "${matrices[@]}"; do
  docker run --rm \
    -v "$(pwd)/outputs:/work/outputs" \
    bailey-hp \
    ./build/cg_solver --matrix "$m" --precision dd --tol 1.0e-12 --max-iter 2.0\
                      --export-mat "./outputs/${m}_dd.mat"
done
