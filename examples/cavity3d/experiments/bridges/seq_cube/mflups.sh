#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=05:00:00

total_cores=$(( $(echo $SLURM_JOB_CPUS_PER_NODE | cut -d'(' -f 1) ))
cores=$(( $total_cores))

mybin=../../../${CODE}/cavity3d

seq_square () {
  # module list
  echo $total_cores #bridges give 28
  echo $cores

  # Start Simulation
  for ((k=0; k<${#dim[@]}; k++)); do
  # for ((k=${#dim[@]}-1; k>=0; k--)); do
    Height=${dim[k]}
    Width=${dim[k]}
    Length=${dim[k]}
    
    if (( $Height < ${PrismSize:-1} )); then
      continue
    fi

    for repeat in 0 1 2 3 4; do
      echo "$mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${PrismSize:-1}"
      mpirun -n 1 $mybin $dim ${steps[k]} ${PrismSize:-1} ${warmup_steps[k]} $Height $Width $Length
      echo "---------------------------------------------------------------------"
    echo
    done
  done
}

dim=(         64  128 256 384 512 640 768 896)
warmup_steps=(600 100 20  4   4   2   2   2)
steps=(       400 200 100 36  20  8   8   4)
time seq_square > seq_${CODE}${PrismSize:+_}${PrismSize}.log