#!/usr/bin/env bash

# First created: 
# Last modified: 2018 Aug 30

# Author: Yuankun Fu
# email: qoofyk@gmail.com

module load remora

# mybin=../../../${CODE}/cavity3d
mybin=../cavity3d

seq_square () {
  # module list
  # echo $total_cores #bridges give 28
  # echo $cores

  # Start Simulation
  for ((k=0; k<${#dim[@]}; k++)); do
  # for ((k=${#dim[@]}-1; k>=0; k--)); do
    Height=${dim[k]}
    Width=${dim[k]}
    Length=${dim[k]}

    cd ${dim[k]}
    
    if (( $Height < ${PrismSize:-1} )); then
      continue
    fi

    # 1 2 3 4
    for repeat in 0 ; do
      echo "$mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${PrismSize:-1}"
      remora mpirun -n 1 $mybin $dim ${steps[k]} ${PrismSize:-1} ${warmup_steps[k]} $Height $Width $Length
      echo "---------------------------------------------------------------------"
    echo
    done

    cd ..
  done
}

# dim=(         64  128 256 384 512 640 768 896)
# warmup_steps=(600 100 20  4   4   2   2   2)
# steps=(       400 200 100 36  20  8   8   4)

dim=(         768 896)
warmup_steps=(2   2)
steps=(       8   4)

PrismSize=32


seq_square