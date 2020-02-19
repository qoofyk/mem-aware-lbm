#include <iostream>
#include <omp.h>
#include <vector>

using namespace std;

static long num_steps= 100000000;
double step;
#define    PAD      8  // assume 64 byte L1 cache line size
// #define NUM_THREADS 8

#pragma optimize ("", off)
void main () {
	int nthreads;
	double pi;
	step = 1.0/(double) num_steps;

#ifdef _OPENMP
    int NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
    omp_set_num_threads(NUM_THREADS);
#endif

  vector<vector<double>> sum(NUM_THREADS, vector<double>(PAD));

#pragma omp parallel
{
	int i, id, nthrds;
	double x;
	id = omp_get_thread_num();
	nthrds = omp_get_num_threads();
	if (id == 0)   nthreads = nthrds;

	for (i = id, sum[id][0] = 0.0; i < num_steps; i = i + nthrds) {
		x = (i + 0.5) * step;
		sum[id][0] += 4.0 / (1.0 + x * x);
	}
}

	for(int i=0, pi=0.0; i < nthreads; i++)
		pi += sum[i][0] * step;
}