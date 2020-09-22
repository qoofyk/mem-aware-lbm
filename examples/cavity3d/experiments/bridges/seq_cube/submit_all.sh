#!/usr/bin/env bash

# First created: 2020 Aug 09
# Last modified: 2020 Aug 09

# Author: Yuankun Fu
# email: qoofyk@gmail.com

CODE=fuse bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}" --output="jobtime/seq-cube-${CODE}.out" mflups.sh'
CODE=2step  bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}" --output="jobtime/seq-cube-${CODE}.out" mflups.sh'

# 8 16 32 64 128 256
for t in 32 64 128 192 256 320 384 448; do
  CODE=fuse_prism  PrismSize=$t bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}${PrismSize:+_}${PrismSize}" --output="jobtime/seq-cube-${CODE}${PrismSize:+_}${PrismSize}.out" mflups.sh'
  CODE=2step_whole_prism PrismSize=$t bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}${PrismSize:+_}${PrismSize}" --output="jobtime/seq-cube-${CODE}${PrismSize:+_}${PrismSize}.out" mflups.sh'
  CODE=2step_3parts_prism PrismSize=$t bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}${PrismSize:+_}${PrismSize}" --output="jobtime/seq-cube-${CODE}${PrismSize:+_}${PrismSize}.out" mflups.sh'
  CODE=2step_whole_prism_unroll PrismSize=$t bash -c 'sbatch --nodes=1 --job-name="seq-cube-${CODE}${PrismSize:+_}${PrismSize}" --output="jobtime/seq-cube-${CODE}${PrismSize:+_}${PrismSize}.out" mflups.sh'
done