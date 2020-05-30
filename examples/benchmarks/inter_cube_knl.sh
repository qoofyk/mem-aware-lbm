## idev -p skx-dev -N 1 -n 48 -m 120
## idev -p development -N 1 -n 68 -m 120

PBS_O_HOME=$HOME
PBS_O_WORKDIR=$(pwd)
# export SCRATCH_DIR=${SCRATCH}/mem-aware-lbm/${SLURM_JOBID}

# bridges
groupname=$(id -Gn)
EMPTY_DIR=${SCRATCH}/empty/

BIN1=${PBS_O_WORKDIR}/cavity3d_step2_whole_omp_unroll_pyramid/cavity3d_step2
BIN2=${PBS_O_WORKDIR}/cavity3d_pymd/cavity3d

BIN3=${PBS_O_WORKDIR}/cavity3d_step2_3parts_omp_line/cavity3d_step2
BIN4=${PBS_O_WORKDIR}/cavity3d/cavity3d

max_cores=68

step2_strong_perf () {
  N=$((dim - 1))
  Nx=$dim
  Ny=$dim
  Nz=$dim
  for i in `seq 1`
  do
    for ((k=0; k<${#threads[@]}; k++)); do
      export OMP_NUM_THREADS=${threads[k]}
      export OMP_PROC_BIND=spread
      export KMP_AFFINITY=verbose
      # if [[ $OMP_NUM_THREADS -eq $max_cores ]]
      # then
      #   export KMP_AFFINITY=verbose,granularity=core,compact
      # else
      #   export KMP_AFFINITY=verbose,granularity=core,scatter
      # fi
      # export KMP_AFFINITY=granularity=core,scatter
      echo "KMP_AFFINITY = $KMP_AFFINITY"

      tile=$((Nx / threads[k]))
      echo " $BIN1 $N $steps $tile $warmup_steps $Nx $Ny $Nz"
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
      mpirun -n $num_proc $BIN2 $N $steps $tile $warmup_steps $Nx $Ny $Nz
    done
  done
}

# dim=$((60 * 4))
# steps=150
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=40
# # threads=(1 2 4 6 8 10 12 16 20 24 30 40 48 60)
# threads=(48 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 4))
# steps=150
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=40
# # threads=(1 2 4 6 8 12 22 24 44 66)
# threads=(44 66)
# step2_strong_perf
# step1_strong_perf

# dim=$((68 * 4))
# steps=150
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=40
# # threads=(1 2 4 8 16 34 68)
# threads=(34 68)
# step2_strong_perf
# step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 5))
# steps=150
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=40
# # threads=(1 2 4 6 10 12 20 30 50 60)
# threads=(50 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 5))
# steps=150
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=40
# # threads=(1 2 6 10 22 30 66)
# threads=(30 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 5))
steps=150
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=40
# threads=(1 2 4 10 20 34 68)
threads=(34 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 6))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 18 20 24 30 36 40 60)
# threads=(40 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 6))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 12 18 22 36 44 66)
# threads=(44 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 6))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 6 8 12 24 34 68)
threads=(34 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 7))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 10 12 14 20 28 30 42 60)
# threads=(42 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 7))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 6 14 22 42 66)
# threads=(42 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 7))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 14 28 34 68)
threads=(34 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 8))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 16 20 24 30 32 40 48 60)
# threads=(48 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 8))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 12 16 22 24 44 48 66)
# threads=(48 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 8))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 8 16 32 34 68)
threads=(34 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 9))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 10 12 18 20 30 36 54 60)
# threads=(54 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 9))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 6 18 22 54 66)
# threads=(54 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 9))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 6 12 18 34 36 68)
threads=(36 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 10))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 20 24 30 40 50 60)
# threads=(50 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 10))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 10 12 20 22 30 44 60 66)
# threads=(60 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 10))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 8 10 20 34 40 68)
threads=(40 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 11))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 10 12 20 22 30 44 60)
# threads=(44 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 11))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 6 22 66)
# threads=(22 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 11))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 22 34 44 68)
threads=(44 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="

# dim=$((60 * 12))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 10 12 16 18 20 24 30 36 40 48 60)
# threads=(48 60)
# step2_strong_perf
# step1_strong_perf

# dim=$((66 * 12))
# steps=50
# # tile=(192 96 48 32 24 16 12 8 6 4)
# warmup_steps=20
# # threads=(1 2 4 6 8 12 18 22 24 36 44 66)
# threads=(44 66)
# step2_strong_perf
# step1_strong_perf

dim=$((68 * 12))
steps=50
# tile=(192 96 48 32 24 16 12 8 6 4)
warmup_steps=20
# threads=(1 2 4 6 8 12 16 24 34 48 68)
threads=(48 68)
step2_strong_perf
step1_strong_perf

echo "========================================================================================="