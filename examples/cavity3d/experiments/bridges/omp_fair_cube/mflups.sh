#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=02:30:00

export OMP_PROC_BIND=spread
total_cores=$(( $(echo $SLURM_JOB_CPUS_PER_NODE | cut -d'(' -f 1) ))
cores=$(( $total_cores))

mybin=../../../${CODE}/cavity3d

omp_square () {
  # module list
  echo $total_cores #bridges give 28
  echo $cores

  # Start Simulation
  dim=$(( $DIM - 1 ))

  # for ((k=0; k<${#threads[@]}; k++)); do
  for ((k=${#threads[@]}-1; k>=0; k--)); do

    if (( ${threads[k]} == 0 )); then
      continue
    fi

    =$DIM
    Width=$DIM
    Length=$DIM

    for repeat in 0 1 2 3 4; do
      echo "spread $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} ${Height[k]} ${Width[k]} ${Length[k]}"
      if [[ ${CODE} == "fuse"* ]]; then
        mpirun -n ${threads[k]} $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} ${Height[k]} ${Width[k]} ${Length[k]}
      else
        export OMP_NUM_THREADS=${threads[k]}
        mpirun -n 1 $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} ${Height[k]} ${Width[k]} ${Length[k]}
      fi
      echo "---------------------------------------------------------------------"
    echo
    done
  done
}


if [ $DIM == 112 ]; then
  threads=(1 2 4 8 14 16 28)
  Height=(  112 56  28  14  8 14  8)
  Width=( 112 112 224 224 224 224 224)
  Length=(  112 224 224 448 784 448 784)
  warmup_steps=(50 50 50 50 50 50 50)
  steps=(100 250 250 850 850 850 1150)

elif [ $DIM == 224 ]; then
  threads=(1 2 4 8 14 16 28)
  Height=(  224 112 56  28  16  28  16
  Width=( 224 224 448 448 448 448 448
  Length=(  224 448 448 896 1568  896 1568
  warmup_steps=(10 10 50 50 50 50 50)
  steps=(20 40 100 150 250 250 550)

elif [ $DIM == 336 ]; then
  threads=(1 2 4 6 8 12 14 16 24 28)
  Height=(  336 168 84  56  42  56  24  42  28  24)
  Width=( 336 336 672 672 672 672 672 672 672 672)
  Length=(  336 672 672 1008  1344  1008  2352  1344  2016  2352)
  warmup_steps=(4 10 10 20 20 20 20 20 20 100)
  steps=(6 10 40 80 80 140 140 140 180 500) 

elif [ $DIM == 448 ]; then
  threads=(1 2 4 8 14 16 28)
  Height=(  448 224 112 56  32  56  32)
  Width=( 448 448 896 896 896 896 896)
  Length=(  448 896 896 1792  3136  1792  3136)
  warmup_steps=(2 10 10 20 20 20 100)
  steps=(4 10 40 80 80 180 500)

elif [ $DIM == 560 ]; then
  threads=(1 2 4 8 10 14 16 20 28)
  Height=(  560 280 140 70  56  40  70  56  40)
  Width=( 560 560 1120  1120  1120  1120  1120  1120  1120)
  Length=(  560 1120  1120  2240  2800  3920  2240  2800  3920)
  warmup_steps=(2 2 10 10 10 20 20 50 50)
  steps=(4 6 10 40 90 80 80 150 150)

elif [ $DIM == 672 ]; then
  threads=(1 2 4 6 8 12 14 16 24 28)
  Height=(  672 336 168 112 84  112 48  84  56  48)
  Width=( 672 672 1344  1344  1344  1344  1344  1344  1344  1344)
  Length=(  672 1344  1344  2016  2688  2016  4704  2688  4032  4704)
  warmup_steps=(2 2 2 2 2 2 10 10 10 10)
  steps=(4 4 4 8 8 8 10 30 50 90)

elif [ $DIM == 784 ]; then
  threads=(1 2 4 8 14 16 28)
  # threads=(0 0 0 8 14 16 28)
  Height=(  784 392 196 98  56  98  56)
  Width=( 784 784 1568  1568  1568  1568  1568)
  Length=(  784 1568  1568  3136  5488  3136  5488)
  warmup_steps=(2 2 2 2 10 20 20)
  steps=(4 4 4 18 50 60 80)

elif [ $DIM == 840 ]; then
  threads=(1 2 4 6 8 10 12 14 20 24 28)
  # threads=(0 0 0 0 0 0 0 14 20 24 28)
  Height=(  840 420 210 140 105 84  140 60  84  70  60)
  Width=( 840 840 1680  1680  1680  1680  1680  1680  1680  1680  1680)
  Length=(  840 1680  1680  2520  3360  4200  2520  5880  4200  5040  5880)
  warmup_steps=(0 2 2 2 2 2 10 10 10 10 10)
  steps=(2 2 4 8 8 18 10 20 40 50 90)

  # if [[ ${CODE} == *"_tile"* ]]; then
  #   echo "contains _tile"
  #   warmup_steps=(48 48 96 96 768 768 768)
  #   steps=(96 96 192 192 600 600 600)
  # fi

else
  echo "Wrong parameters"
fi

time omp_square > omp_dim_${DIM}_${CODE}${TILE:+_}${TILE}.log