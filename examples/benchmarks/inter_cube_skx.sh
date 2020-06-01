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

max_cores=48

step2_strong_perf () {
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
      
      # if [[ $OMP_NUM_THREADS -eq $max_cores ]]
      # then
      #   export KMP_AFFINITY=verbose,granularity=fine,balanced
      # elif [[ $OMP_NUM_THREADS -eq 26 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,1,5,9,11,7,3,13,17,21,23,19,15,25]
      # elif [[ $OMP_NUM_THREADS -eq 32 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,28,32,34,1,5,9,11,7,3,13,17,21,23,19,15,25,29,33,35]
      # elif [[ $OMP_NUM_THREADS -eq 36 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,28,32,34,30,26,1,5,9,11,7,3,13,17,21,23,19,15,25,29,33,35,31,27]
      # elif [[ $OMP_NUM_THREADS -eq 40 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,28,32,34,30,26,36,40,1,5,9,11,7,3,13,17,21,23,19,15,25,29,33,35,31,27,37,41]
      # elif [[ $OMP_NUM_THREADS -eq 42 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,28,32,34,30,26,36,40,44,1,5,9,11,7,3,13,17,21,23,19,15,25,29,33,35,31,27,37,41,45]
      # elif [[ $OMP_NUM_THREADS -eq 44 ]]
      # then
      #   export KMP_AFFINITY=verbose,explicit,granularity=fine,proclist=[0,4,8,10,6,2,12,16,20,22,18,14,24,28,32,34,30,26,36,40,44,46,1,5,9,11,7,3,13,17,21,23,19,15,25,29,33,35,31,27,37,41,45,47]
      # else
      #   export KMP_AFFINITY=verbose,granularity=core,scatter
      # fi
      # export KMP_AFFINITY=verbose,granularity=core,scatter
      # export KMP_AFFINITY=verbose,granularity=core,compact
      # export KMP_AFFINITY=none,granularity=core

      echo "KMP_AFFINITY = $KMP_AFFINITY"
      tile=$((Nx / threads[k]))
      echo " $BIN1 $N $steps $tile ${warmup_steps} $Nx $Ny $Nz"
      mpirun -n 1 $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz
    done
  done
}

step1_strong_perf () {
  N=$((dim - 1))
  Nx=$dim
  Ny=$dim
  Nz=$dim
  for i in `seq 1`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      num_proc=${threads[k]}
      tile=$((Nx / threads[k]))

      echo " $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      mpirun -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz

      echo " $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
      mpirun -n $num_proc $BIN4 $N $steps $tile $warmup_steps $Nx $Ny $Nz
    done
  done
}

dim=$((max_cores * 4)) #192
steps=300
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=150
# threads=(1 2 4 6 8 12 16 24 32 48)
threads=(32 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 5)) #240
steps=300
# tile=(112 112 56 28 16 14 8)
warmup_steps=150
# threads=(1 2 4 6 8 10 12 16 20 24 30 40 48)
threads=(40 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 6)) #288
steps=100
# tile=(112 84 84 42 24 21 12)
warmup_steps=50
# threads=(1 2 4 6 8 12 16 18 24 32 36 48)
threads=(36 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 7)) #336
steps=100
# tile=(112 112 112 56 32 28 16)
warmup_steps=50
# threads=(1 2 4 6 8 12 14 16 24 28 42 48)
threads=(42 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 8)) #384
steps=100
# tile=(112 70 70 70 40 35 20)
warmup_steps=50
# threads=(1 2 4 6 8 12 16 24 32 48)
threads=(32 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 9)) #432
steps=100
# tile=(112 84 84 84 48 42 24)
warmup_steps=50
# threads=(1 2 4 6 8 12 16 18 24 36 48)
threads=(36 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 10)) #480
steps=100
# tile=(112 84 84 84 48 42 24)
warmup_steps=50
# threads=(1 2 4 6 8 10 12 16 20 24 30 32 40 48)
threads=(40 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 11)) #528
steps=100
# tile=(112 84 84 84 48 42 24)
warmup_steps=50
# threads=(1 2 4 6 8 12 16 22 24 44 48)
threads=(44 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 12)) #576
steps=50
# tile=(112 84 84 84 48 42 24)
warmup_steps=20
# threads=(1 2 4 6 8 12 16 18 24 32 36 48)
threads=(36 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 13)) #624
steps=50
# tile=(112 84 84 84 48 42 24)
warmup_steps=20
# threads=(1 2 4 6 8 12 16 24 26 48)
threads=(26 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 14)) #672
steps=50
# tile=(112 84 84 84 48 42 24)
warmup_steps=20
# threads=(1 2 4 6 8 12 14 16 24 28 32 42 48)
threads=(42 48)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

dim=$((max_cores * 15)) #720
steps=50
# tile=(112 84 84 84 48 42 24)
warmup_steps=20
# threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48)
threads=(40 48)
step2_strong_perf
step1_strong_perf