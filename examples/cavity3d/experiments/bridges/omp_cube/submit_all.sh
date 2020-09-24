#!/usr/bin/env bash

# First created: 2020 Aug 09
# Last modified: 2020 Aug 09

# Author: Yuankun Fu
# email: qoofyk@gmail.com

# 112 224 336 448 560 672 784 840
for d in 112 224 336 448 560 672 784 840; do
  DIM=$d CODE=fuse bash -c 'sbatch --nodes=1 --job-name="omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'

done

# 8 16 32 64 128 192 256 320 384 448
for t in 8 16 32 64 128 192 256 320 384; do
  for d in 112 224 336 448 560 672 784 840; do
    DIM=$d CODE=fuse_prism  TILE=$t bash -c 'sbatch --nodes=1 --job-name="omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
    DIM=$d CODE=2step_whole_prism_unroll_omp TILE=$t bash -c 'sbatch --nodes=1 --job-name="omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
  done
done

# DIM=112 CODE=fuse bash -c 'sbatch --nodes=1 --job-name="omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
# DIM=7168 CODE=2step_whole_prism_unroll_omp TILE=64 bash -c 'sbatch --nodes=1 --job-name="omp-cube-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'