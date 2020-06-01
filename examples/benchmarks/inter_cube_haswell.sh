# Bridges RM: interact -p RM -N 1 -n 28 -t 2:00:00
# SKX: idev -p skx-dev -N 1 -n 48 -m 120
# KNL: idev -p development -N 1 -n 68 -m 120

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

cube_step2_whole_strong_perf () {
  N=$((dim - 1))
  Nx=$dim
  Ny=$dim
  Nz=$dim
  for i in `seq 2`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      export OMP_NUM_THREADS=${threads[k]}
      export OMP_PROC_BIND=spread
      export KMP_AFFINITY=verbose

      ################# BIN1
      echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
      tile=$((Nx / threads[k]))
      echo "spread $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      export OMP_PROC_BIND=false
      if [[ $OMP_NUM_THREADS -eq $max_cores ]]
      then
        # export KMP_AFFINITY=verbose,granularity=fine,balanced
        export OMP_PROC_BIND=spread
      elif [[ $OMP_NUM_THREADS -eq 16 ]]
      then
        export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27]
      elif [[ $OMP_NUM_THREADS -eq 24 ]]
      then
        export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13]
      else
        export OMP_PROC_BIND=spread
      fi
      # export KMP_AFFINITY=verbose,granularity=core,scatter
      # export KMP_AFFINITY=verbose,granularity=core,compact
      # export KMP_AFFINITY=none,granularity=core

      echo "KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
      tile=$((Nx / threads[k]))
      echo "explicit $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      echo "---------------------------------------------------------------------"
    done
  done
}

cube_step2_3parts_strong_perf () {
  N=$((dim - 1))
  Nx=$((dim + 3))
  Ny=$dim
  Nz=$dim
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

      export OMP_PROC_BIND=false
      if [[ $OMP_NUM_THREADS -eq $max_cores ]]
      then
        # export KMP_AFFINITY=verbose,granularity=fine,balanced
        export OMP_PROC_BIND=spread
      elif [[ $OMP_NUM_THREADS -eq 16 ]]
      then
        export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,17,18,5,6,7,8,9,23,24,25,26,27]
      elif [[ $OMP_NUM_THREADS -eq 24 ]]
      then
        export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,14,1,15,2,16,17,4,18,5,19,6,7,21,8,22,9,23,24,11,25,12,26,13]
      else
        export OMP_PROC_BIND=spread
      fi
      # export KMP_AFFINITY=verbose,granularity=core,scatter
      # export KMP_AFFINITY=verbose,granularity=core,compact
      # export KMP_AFFINITY=none,granularity=core

      echo "explicit KMP_AFFINITY=${KMP_AFFINITY}, OMP_PROC_BIND=${OMP_PROC_BIND}"
      tile=$((Nx / threads[k]))
      echo "explicit $BIN3 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN3 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      echo "---------------------------------------------------------------------"
    done
  done
}

cube_step1_strong_perf () {
  N=$((dim - 1))
  Nx=$dim
  Ny=$dim
  Nz=$dim
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

      ######## MPI default binding same as spread binding
      # echo "spread $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz" #same binding with default
      # mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "spread $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      # mpirun -env I_MPI_PIN_PROCESSOR_LIST=:map=spread -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      # echo "---------------------------------------------------------------------"

      ######## MPI default's performance same as explicit
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

# dim=$((max_cores * 4)) #192
# steps=300
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=150
# # threads=(1 2 4 6 8 12 16 24 32 48)
# threads=(16 28)
# cube_step2_whole_strong_perf
# cube_step2_3parts_strong_perf
# cube_step1_strong_perf

echo "========================================================================================="

# dim=$((max_cores * 5)) #240
# steps=300
# # tile=(112 112 56 28 16 14 8)
# warmup_steps=150
# # threads=(1 2 4 6 8 10 12 16 20 24 30 40 48)
# threads=(40 48)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 6)) #288
# steps=100
# # tile=(112 84 84 42 24 21 12)
# warmup_steps=50
# # threads=(1 2 4 6 8 12 16 18 24 32 36 48)
# threads=(36 48)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 7)) #336
# steps=100
# # tile=(112 112 112 56 32 28 16)
# warmup_steps=50
# # threads=(1 2 4 6 8 12 14 16 24 28 42 48)
# threads=(42 48)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 8)) #384
# steps=100
# # tile=(112 70 70 70 40 35 20)
# warmup_steps=50
# # threads=(1 2 4 6 8 12 16 24 32 48)
# threads=(16 28)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# # dim=$((max_cores * 9)) #432
# # steps=100
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=50
# # # threads=(1 2 4 6 8 12 16 18 24 36 48)
# # threads=(36 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 10)) #480
# # steps=100
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=50
# # # threads=(1 2 4 6 8 10 12 16 20 24 30 32 40 48)
# # threads=(40 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 11)) #528
# # steps=100
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=50
# # # threads=(1 2 4 6 8 12 16 22 24 44 48)
# # threads=(44 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 12)) #576
# # steps=50
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=20
# # # threads=(1 2 4 6 8 12 16 18 24 32 36 48)
# # threads=(36 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 13)) #624
# # steps=50
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=20
# # # threads=(1 2 4 6 8 12 16 24 26 48)
# # threads=(26 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 14)) #672
# # steps=50
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=20
# # # threads=(1 2 4 6 8 12 14 16 24 28 32 42 48)
# # threads=(42 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# # dim=$((max_cores * 15)) #720
# # steps=50
# # # tile=(112 84 84 84 48 42 24)
# # warmup_steps=20
# # # threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48)
# # threads=(40 48)
# # cube_step2_whole_strong_perf
# # cube_step1_strong_perf

# # echo "========================================================================================="

# dim=$((max_cores * 16)) #720
# steps=50
# # tile=(112 84 84 84 48 42 24)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48)
# threads=(16 28)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 20)) #720
# steps=50
# # tile=(112 84 84 84 48 42 24)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48)
# threads=(20 28)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 28)) #720
# steps=50
# # tile=(112 84 84 84 48 42 24)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48)
# threads=(16 28)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="

dim=$((max_cores * 30)) #840
steps=50
# tile=(112 84 84 84 48 42 24)
warmup_steps=20
# threads=(1 2 4 6 8 10 12 14 20 24 28)
threads=(24 28)
# cube_step2_whole_strong_perf
# cube_step2_3parts_strong_perf
cube_step1_strong_perf

# echo "========================================================================================="

# dim=$((max_cores * 32)) #876 Exceed max memory
# steps=50
# # tile=(112 84 84 84 48 42 24)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 14 20 24 28)
# threads=(16 28)
# cube_step2_whole_strong_perf
# cube_step1_strong_perf

# echo "========================================================================================="
