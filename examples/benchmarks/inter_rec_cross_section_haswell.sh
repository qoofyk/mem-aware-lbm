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
  N=1
  Nx=$Whole_Nx
  for i in `seq 2`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      export OMP_NUM_THREADS=${threads[k]}
      export OMP_PROC_BIND=spread
      export KMP_AFFINITY=verbose

      echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
      tile=$((Nx / threads[k]))
      echo "spread $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz

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
  done
}

rec_cross_section_step2_3parts_strong_perf () {
  N=1
  Nx=$((Whole_Nx + 3))
  for i in `seq 2`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      export OMP_NUM_THREADS=${threads[k]}
      export OMP_PROC_BIND=spread
      export KMP_AFFINITY=verbose

      ################# BIN3
      echo "spread KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
      tile=$((Nx / threads[k]))
      echo "spread $BIN3 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN3 $N $steps $tile $warmup_steps $Nx $Ny $Nz

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
  done
}

rec_cross_section_step1_strong_perf () {
  N=1
  Nx=$Whole_Nx
  for i in `seq 1`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      num_proc=${threads[k]}
      tile=$((Nx / threads[k]))

      export I_MPI_DEBUG=4 #MPI runtime report process binding

      echo "default $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      mpirun -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      echo "default $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      mpirun -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      echo "---------------------------------------------------------------------"

      # echo "spread $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz" # spread shows same core binding with default
      # mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "spread $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      # mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "---------------------------------------------------------------------"

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
  done
}

tile=8

# Whole_Nx=$((32 * 112))
# for base in 32
# do
#   for i in {10..16..2}
#   do
#     Ny=$((base * i))
#     Nz=$base
#     steps=300
#     # tile=(192 96 48 32 24 16 12 8 6 4)
#     warmup_steps=150
#     # threads=(1 2 4 6 8 12 16 24 32 48)
#     threads=(16 28)
#     rec_cross_section_step2_whole_strong_perf
#     rec_cross_section_step2_3parts_strong_perf
#     rec_cross_section_step1_strong_perf

#     echo "========================================================================================="
#   done
# done

Whole_Nx=$((16 * 112))
for base in 64
do
  for i in {1..8}
  do
    Ny=$((base * i))
    Nz=$base
    steps=150
    # tile=(192 96 48 32 24 16 12 8 6 4)
    warmup_steps=100
    # threads=(1 2 4 6 8 12 16 24 32 48)
    threads=(16 28)
    rec_cross_section_step2_whole_strong_perf
    rec_cross_section_step2_3parts_strong_perf
    rec_cross_section_step1_strong_perf

    echo "========================================================================================="
  done
done

Whole_Nx=$((8 * 112))
for base in 128
do
  for i in {1..4}
  do
    Ny=$((base * i))
    Nz=$base
    steps=150
    # tile=(192 96 48 32 24 16 12 8 6 4)
    warmup_steps=100
    # threads=(1 2 4 6 8 12 16 24 32 48)
    threads=(16 28)
    rec_cross_section_step2_whole_strong_perf
    rec_cross_section_step2_3parts_strong_perf
    rec_cross_section_step1_strong_perf

    echo "========================================================================================="
  done
done