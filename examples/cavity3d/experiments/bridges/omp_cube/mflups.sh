#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=03:30:00

export OMP_PROC_BIND=spread
total_cores=$(( $(echo $SLURM_JOB_CPUS_PER_NODE | cut -d'(' -f 1) ))
cores=$(( $total_cores))

mybin=../../../${CODE}/cavity3d

omp_square () {
  # module list
  echo $total_cores #bridges give 28
  echo $cores

  # Start Simulation
  Height=$DIM
  Width=$DIM
  Length=$DIM
  dim=$(( $DIM - 1 ))

  # for ((k=0; k<${#threads[@]}; k++)); do
  for ((k=${#threads[@]}-1; k>=0; k--)); do
    

    if (( ${threads[k]} == 0 )); then
      continue
    fi

    for repeat in 0 1 2 3 4; do
      echo "spread $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} $Height $Width $Length"
      if [[ ${CODE} == "fuse"* ]]; then
        mpirun -n ${threads[k]} $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} $Height $Width $Length
      else
        export OMP_NUM_THREADS=${threads[k]}
        mpirun -n 1 $mybin $dim ${steps[k]} ${TILE:-1} ${warmup_steps[k]} $Height $Width $Length
      fi
      echo "---------------------------------------------------------------------"
    echo
    done
  done
}


if [ $DIM == 112 ]; then
  threads=(1 2 4 8 14 16 28)
  warmup_steps=(50 50 50 50 50 50 50)
  steps=(100 250 250 850 850 850 1150)

elif [ $DIM == 224 ]; then
  threads=(1 2 4 8 14 16 28)
  warmup_steps=(10 10 50 50 50 50 50)
  steps=(20 40 100 150 250 250 550)

elif [ $DIM == 336 ]; then
  threads=(1 2 4 6 8 12 14 16 24 28)
  warmup_steps=(4 10 10 10 10 10 10 20 20 20)
  steps=(6 10 40 70 70 90 90 100 180 180)

elif [ $DIM == 448 ]; then
  threads=(1 2 4 8 14 16 28)
  warmup_steps=(2 2 10 10 10 10 10)
  steps=(4 8 10 30 50 70 90)

elif [ $DIM == 560 ]; then
  threads=(1 2 4 8 10 14 16 20 28)
  warmup_steps=(2 2 10 10 10 10 10 20 20)
  steps=(4 6 6 20 20 50 50 60 80)

elif [ $DIM == 672 ]; then
  threads=(1 2 4 8 10 14 16 20 28)
  warmup_steps=(2 2 10 10 10 10 10 20 20)
  steps=(4 6 6 20 20 50 50 60 80)

elif [ $DIM == 784 ]; then
  threads=(1 2 4 8 14 16 28)
  warmup_steps=(2 2 2 2 10 10 10)
  steps=(4 4 4 8 30 50 70)

elif [ $DIM == 840 ]; then
  threads=(1 2 4 6 8 10 12 14 20 24 28)
  warmup_steps=(0 2 2 2 2 2 10 10 10 10 10)
  steps=(2 4 6 10 10 20 20 20 40 60 80)

  # if [[ ${CODE} == *"_tile"* ]]; then
  #   echo "contains _tile"
  #   warmup_steps=(48 48 96 96 768 768 768)
  #   steps=(96 96 192 192 600 600 600)
  # fi

else
  echo "Wrong parameters"
fi

time omp_square > omp_dim_${DIM}_${CODE}${TILE:+_}${TILE}.log