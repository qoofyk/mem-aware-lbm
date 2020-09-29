#!/usr/bin/env bash

# First created: 2020 Sep 11
# Last modified: 2020 Sep 11

# Author: Yuankun Fu
# email: qoofyk@gmail.com

for tile in 0 8 16 32 64 128 256 192 256; do
  for code in 2step_whole_prism_unroll; do
    for d in 672; do
      cat 0926/omp_dim_${d}_${code}_omp_${tile}.log >> omp_dim_${d}_${code}_omp_${tile}.log
    done
  done
done