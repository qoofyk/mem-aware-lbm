#include <stdio.h>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int iam = 0, np = 1;
  int N=48;
  int thread_block = N/atoi(getenv("OMP_NUM_THREADS"));

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  #pragma omp parallel default(shared) private(iam, np)
  {
    np = omp_get_num_threads();
    iam = omp_get_thread_num();
    // printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
    //        iam, np, rank, numprocs, processor_name);

    // #pragma omp for
    #pragma omp for schedule(static, thread_block)
    for (int i=0; i<N; ++i)
      printf("i=%d, Hello from thread %d out of %d from process %d out of %d on %s\n",
           i, iam, np, rank, numprocs, processor_name);

    #pragma omp for schedule(static, thread_block)
    for (int i=0; i<N; ++i)
      printf("i=%d, Hello from thread %d out of %d from process %d out of %d on %s\n",
           i, iam, np, rank, numprocs, processor_name);
  }

  MPI_Finalize();
}
