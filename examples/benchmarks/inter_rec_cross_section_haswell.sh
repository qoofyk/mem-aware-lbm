# idev -p skx-dev -N 1 -n 48 -m 120
# idev -p development -N 1 -n 68 -m 120

PBS_O_HOME=$HOME
PBS_O_WORKDIR=$(pwd)
export SCRATCH_DIR=${SCRATCH}/mem-aware-lbm/${SLURM_JOBID}

# bridges
groupname=$(id -Gn)
EMPTY_DIR=${SCRATCH}/empty/

BIN1=${PBS_O_WORKDIR}/cavity3d_step2_whole_omp_unroll_pyramid/cavity3d_step2
BIN2=${PBS_O_WORKDIR}/cavity3d_pymd/cavity3d

BIN3=${PBS_O_WORKDIR}/cavity3d_step2_3parts_omp_line/cavity3d_step2
BIN4=${PBS_O_WORKDIR}/cavity3d/cavity3d

max_cores=28

rec_cross_section_step2_whole_strong_perf () {
  N=$Nz

  # for ((k=0; k<${#threads[@]}; k++)); do
  for ((k=${#threads[@]}-1; k>=0; k--)); do
    export OMP_NUM_THREADS=${threads[k]}
    export OMP_PROC_BIND=spread
    # export KMP_AFFINITY=verbose

    echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
    tile=$((Nx / threads[k]))
    echo "spread $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"

    for i in `seq 3`
    do
      mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz
    done
    # export OMP_PROC_BIND=false
    # if [[ $OMP_NUM_THREADS -eq $max_cores ]]
    # then
    #   # export KMP_AFFINITY=verbose,granularity=fine,balanced
    #   export OMP_PROC_BIND=spread
    # elif [[ $OMP_NUM_THREADS -eq 16 ]]
    # then
    #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27]
    # elif [[ $OMP_NUM_THREADS -eq 24 ]]
    # then
    #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13]
    # else
    #   export OMP_PROC_BIND=spread
    # fi
    # # export KMP_AFFINITY=verbose,granularity=core,scatter
    # # export KMP_AFFINITY=verbose,granularity=core,compact
    # # export KMP_AFFINITY=none,granularity=core

    # echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
    # tile=$((Nx / threads[k]))
    # echo "explicit $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
    # mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz

    echo "---------------------------------------------------------------------"
  done
}

rec_cross_section_step2_3parts_strong_perf () {
  N=$Nz
  myNx=$((Nx + 3))
  
  # for ((k=0; k<${#threads[@]}; k++)); do
  for ((k=${#threads[@]}-1; k>=0; k--)); do
    export OMP_NUM_THREADS=${threads[k]}
    export OMP_PROC_BIND=spread
    # export KMP_AFFINITY=verbose

    ################# BIN3
    echo "spread KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
    tile=$((Nx / threads[k]))
    echo "spread $BIN3 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"

    for i in `seq 3`
    do
      mpirun -n 1 $BIN3 $N $steps $tile $warmup_steps $myNx $Ny $Nz
    done

    # export OMP_PROC_BIND=false
    # if [[ $OMP_NUM_THREADS -eq $max_cores ]]
    # then
    #   # export KMP_AFFINITY=verbose,granularity=fine,balanced
    #   export OMP_PROC_BIND=spread
    # elif [[ $OMP_NUM_THREADS -eq 16 ]]
    # then
    #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27]
    # elif [[ $OMP_NUM_THREADS -eq 24 ]]
    # then
    #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13]
    # else
    #   export OMP_PROC_BIND=spread
    # fi
    # # export KMP_AFFINITY=verbose,granularity=core,scatter
    # # export KMP_AFFINITY=verbose,granularity=core,compact
    # # export KMP_AFFINITY=none,granularity=core

    # echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
    # tile=$((Nx / threads[k]))
    # echo "explicit $BIN3 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
    # mpirun -n 1 $BIN3 $N $steps $tile $warmup_steps $Nx $Ny $Nz

    echo "---------------------------------------------------------------------"
  done
}

rec_cross_section_step1_strong_perf () {
  N=$Nz

  # for ((k=0; k<${#threads[@]}; k++)); do
  for ((k=${#threads[@]}-1; k>=0; k--)); do
      num_proc=${threads[k]}
      tile=$((Nx / threads[k]))

      # export I_MPI_DEBUG=4 #MPI runtime report process binding

      # echo "default $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      # mpirun -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "default $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      # mpirun -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "---------------------------------------------------------------------"

      for i in `seq 2`
      do
        echo "spread $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz" # spread shows same core binding with default
        mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

        echo "spread $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
        mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz
      done
      echo "---------------------------------------------------------------------"

      # mpi default's performance same as explicit
      # if [[ $num_proc -eq 16 ]]
      # then
      #   echo "explict $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      #   mpirun -env I_MPI_PIN_PROCESSOR_LIST=0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27 \
      #         -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      #   echo "explict $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      #   mpirun -env I_MPI_PIN_PROCESSOR_LIST=0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27 \
      #         -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # elif [[ $num_proc -eq 24 ]]
      # then
      #   echo "explict $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      #   mpirun -env I_MPI_PIN_PROCESSOR_LIST=0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13 \
      #         -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      #   echo "explict $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      #   mpirun -env I_MPI_PIN_PROCESSOR_LIST=0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13 \
      #         -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      #   echo "---------------------------------------------------------------------"
      # fi
  done
}

Nx=$((32 * 112))
Nz=32
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done

Nx=$((16 * 112))
Nz=64
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done

Nx=$((8 * 112))
Nz=128
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done

Nx=$((4 * 112))
Nz=256
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done

Nx=$((2 * 112))
Nz=512
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done

Nx=$((112))
Nz=1024
for Ny in {32..864..32}
do
  steps=150
  # tile=(192 96 48 32 24 16 12 8 6 4)
  warmup_steps=100
  threads=(1 2 4 8 14 16 28)
  rec_cross_section_step2_whole_strong_perf
  rec_cross_section_step2_3parts_strong_perf
  rec_cross_section_step1_strong_perf
  echo "========================================================================================="  
done