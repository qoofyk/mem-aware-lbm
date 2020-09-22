#!/usr/bin/env bash

# First created: 
# Last modified: 2018 Aug 30

# Author: Yuankun Fu
# email: qoofyk@gmail.com

Height=32
Width=32
Length=32
dim=$(( Height - 1 ))

Warm_up=0
Measure=10

PrismSize=4

# This is correct results
cd fuse
mpirun -n 1 ./cavity3d $dim $Measure $PrismSize $Warm_up $Height $Width $Length > vel_fuse_${Height}_${Measure}.dat
cd ..

# Now test and verify with different other code
# 2step 2step_3parts_prism 2step_whole_prism 2step_whole_prism_unroll
for CODE in fuse_prism 2step 2step_3parts_prism 2step_whole_prism 2step_whole_prism_unroll; do
  echo "run ${CODE}"
  cd ${CODE}
  mpirun -n 1 ./cavity3d $dim $Measure $PrismSize $Warm_up $Height $Width $Length > vel_${CODE}_${Height}_${Measure}.dat
  echo "verify ${CODE}"
  diff ../fuse/vel_fuse_${Height}_${Measure}.dat vel_${CODE}_${Height}_${Measure}.dat
  cd ..
  echo "---------------------------------------"
done

# Verify OMP
echo "Verify with 1 thread"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=1
for CODE in 2step_whole_prism_unroll_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  mpirun -n 1 ./cavity3d $dim $Measure $PrismSize $Warm_up $Height $Width $Length > vel_${CODE}_${Height}_${Measure}.dat
  echo "verify ${CODE}"
  diff ../fuse/vel_fuse_${Height}_${Measure}.dat vel_${CODE}_${Height}_${Measure}.dat
  cd ..
  echo "---------------------------------------"
done

echo "Verify with 2 thread"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=2
for CODE in 2step_whole_prism_unroll_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  mpirun -n 1 ./cavity3d $dim $Measure $PrismSize $Warm_up $Height $Width $Length > vel_${CODE}_${Height}_${Measure}.dat
  echo "verify ${CODE}"
  diff ../fuse/vel_fuse_${Height}_${Measure}.dat vel_${CODE}_${Height}_${Measure}.dat
  cd ..
  echo "---------------------------------------"
done

echo "Verify with 4 thread"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=4
for CODE in 2step_whole_prism_unroll_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  mpirun -n 1 ./cavity3d $dim $Measure $PrismSize $Warm_up $Height $Width $Length > vel_${CODE}_${Height}_${Measure}.dat
  echo "verify ${CODE}"
  diff ../fuse/vel_fuse_${Height}_${Measure}.dat vel_${CODE}_${Height}_${Measure}.dat
  cd ..
  echo "---------------------------------------"
done

rm -rf */*.dat