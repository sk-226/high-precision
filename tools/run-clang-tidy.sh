#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -eq 0 ]; then
  exit 0
fi

cfg_dir=""
if [ -d build ] && [ -f build/compile_commands.json ]; then
  cfg_dir="build"
elif [ -f compile_commands.json ]; then
  cfg_dir="."
fi

for file in "$@"; do
  if [ -n "$cfg_dir" ]; then
    clang-tidy --quiet --config-file=.clang-tidy -p "$cfg_dir" "$file"
  else
    clang-tidy --quiet --config-file=.clang-tidy "$file"
  fi
done
