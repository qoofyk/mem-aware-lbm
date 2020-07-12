/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
  * Flow in a lid-driven 3D cavity. Benchmark case
**/

#include "palabos3D.h"
#include "palabos3D.hh"   // include full template code
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "ittnotify.h"
#include <exception>
#include <cstdint>

using namespace plb;
using namespace std;

#define K 2 // current merge K = 2 steps

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

namespace plb {
    plint Nx, newNx; // Nx: Computation domain; newNx: Pillar Computation domain
    plint ykTile;  /* pillarTile. The assumed size, while "actual" memory layout on Length & Width are (ykTile + 2). E.g. Given computation domain 16*8*16, palabos actual allocated memory is (16+2) * (8+2) * (16+2) --> Let pillarTile=8, thus Equivalant pillar computation domain is 32 * 8 * 8. Actual allocated memory is (32+2) * (8+2) * (8+2). */
    plint NzTiles; // number of pillarTiles along Length (Z-direction)
    plint NyTiles; // number of pillarTiles along Width (Y-direction)
}

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                  plint Nx, plint Ny, plint Nz)
{
    // const plint nx = parameters.getNx();
    // const plint ny = parameters.getNy();
    // const plint nz = parameters.getNz();
    const plint nx = Nx;
    const plint ny = Ny;
    const plint nz = Nz;
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    // initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    // Modify by Yuankun, set init value to 0.01 to avoid 0 computation
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.01, (T)0.01, (T)0.01) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u, (T)0., u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u, 0., u) );

    lattice.initialize();
}

//OpenMP test hello
void test_omp_hello(){
    plint my_rank = 0;
    plint thread_count = 1;

#ifdef _OPENMP
    #pragma omp parallel
    {
        my_rank = omp_get_thread_num();
        thread_count = omp_get_num_threads();
        printf("Hello from thread %ld of %ld\n", my_rank, thread_count);
    }
#else
    printf("Hello from thread %ld of %ld\n", my_rank, thread_count);
#endif
}

struct MyException1 : public exception {
  const char * what () const throw () {
    return "N % OMP_NUM_THREADS != 0, Not Divisible";
  }
};

struct MyException2 : public exception {
  const char * what () const throw () {
    return "thread_block % ykBlockSize != 0, Not Divisible";
  }
};

struct MyException3 : public exception {
  const char * what () const throw () {
    return "Nz % ykTile != 0 && Ny % ykTile != 0, Not Divisible";
  }
};

int main(int argc, char* argv[]) {

    __itt_pause();
    plbInit(&argc, &argv);
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    plint N;
    plint Ny, Nz;
    plint numIter;
    plint warmUpIter;
    plint NUM_THREADS = 1;

#ifdef _OPENMP
    NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
    omp_set_num_threads(NUM_THREADS);
#endif

    // test_omp_hello();

    try {
        global::argv(1).read(N);
        global::argv(2).read(numIter);
        global::argv(3).read(ykBlockSize);
        global::argv(4).read(warmUpIter);
        global::argv(5).read(Nx);
        global::argv(6).read(Ny);
        global::argv(7).read(Nz);
        global::argv(8).read(ykTile);

        // check Nx % NUM_THREADS == 0
        if (Nx % NUM_THREADS != 0) throw MyException1();
        thread_block = Nx / NUM_THREADS;
        if (thread_block % ykBlockSize != 0) throw MyException2();
        if (Nz % ykTile != 0 && Ny % ykTile != 0)  throw MyException3();
        NzTiles = Nz / ykTile;
        NyTiles = Ny / ykTile;
        newNx = Nx * NyTiles * NzTiles;
    }
    catch (MyException1& e) {
        std::cout << e.what() << std::endl;
        exit(1);
    }
    catch (MyException2& e) {
        std::cout << e.what() << std::endl;
        exit(1);
    }
     catch (MyException3& e) {
        std::cout << e.what() << std::endl;
        exit(1);
    }
    catch (...) {
        pcout << "Wrong parameters. The syntax is " << std::endl;
        pcout << argv[0] << " N" << std::endl;
        pcout << "where N is the resolution. The benchmark cases published " << std::endl;
        pcout << "on the Palabos Wiki use N=100, N=400, N=1000, or N=4000." << std::endl;
        exit(1);
    }

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 1.,    // Re, Reynolds number
            N,         // N, resolution
            1.,        // lx
            1.,        // ly
            1.         // lz
    );

    pcout << "Starting 3D cavity benchmark with " << Nx << "x" << Ny << "x" << Nz << " grid points (Computation domain: Nx x Ny x Nz)"
          << " Estimated memory occupied " << Nx * Ny * Nz * 168 / (1024*1024) << " MB."
          << " .\n";

    // MultiBlockLattice3D<T, DESCRIPTOR> lattice (
    //         parameters.getNx(), parameters.getNy(), parameters.getNz(),
    //         new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            Nx, Ny, Nz,
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    plint numProcs = global::mpi().getSize();
    pcout << "Num of MPI Procs: " << numProcs << " Num of OpenMP threads: " << NUM_THREADS
          << "\nthread_block: " << thread_block << ", Parallelpiped ykBlockSize: " << ykBlockSize 
          << "\nPalabos Memory allocation domain: (Nx+2) x (Ny+2) x (Nz+2): " << Nx + 2 << 'x' << Ny + 2 << 'x' << Nz + 2 
          << "\npillarTile=" << ykTile << ", newNx=" << newNx << ", NzTiles=" << NzTiles << ',' << ", NyTiles=" << NyTiles 
          << "\nPillar Memory allocation domain: (newNx+2) x (pillarTile+2) x (pillarTile+2): " << newNx + 2 << 'x' << ykTile + 2 << 'x' << ykTile + 2 << '\n';
    // Current cores run approximately at 5 Mega Sus.
    T estimateSus= 5.e6 * numProcs;
    // The benchmark should run for approximately two minutes
    // (2*60 seconds).
    T wishNumSeconds = 60.;
    plint numCells = lattice.getBoundingBox().nCells();

    // Run at least three iterations.
    // plint numIter = std::max( (plint)3,
    //                           (plint)(estimateSus*wishNumSeconds/numCells+0.5));

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    global::timer("cavitySetup").start();
    // cavitySetup(lattice, parameters, *boundaryCondition);
    cavitySetup(lattice, parameters, *boundaryCondition, Nx, Ny, Nz);
    pcout << "cavitySetup time (s) = "<< global::timer("cavitySetup").getTime() << std::endl;

#if 0
    pcout << "Init: Velocity norm of the box: " << endl;
    // Box3D mybox(0, N, 0, N, 0, N);
    // pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, mybox)) << endl;
    for (plint iX=0; iX<=N; ++iX){
        for (plint iY=0; iY<=N; ++iY){
            Box3D line(iX, iX, iY, iY, 0, N);
            pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, line)) << endl;
        }
    }
#endif

    // Run the benchmark once "to warm up the machine".
    for (plint iT = 0; iT < warmUpIter; iT += K) {
    //for (plint iT=0; iT<2; iT += 2) { // use fixed value could have higher performance
        // pcout << "iT=" << iT << std::endl;
        lattice.step2collideAndStream();
    }

    // pcout << "Start bench!" << std::endl;
    // Run the benchmark for good.
    __itt_resume();
    global::timer("benchmark").start();
    // global::profiler().turnOn();
    for (plint iT = 0; iT < numIter; iT += K) {
        // pcout << "iT=" << iT << std::endl;
        lattice.step2collideAndStream();
    }

#if 0
    pcout << "After: Velocity norm of the box: " << endl;
    // pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, mybox)) << endl;
    for (plint iX = 0; iX <= N; ++iX){
        for (plint iY = 0; iY <= N; ++iY){
            Box3D line(iX, iX, iY, iY, 0, N);
            pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, line)) << endl;
        }
    }
#endif
    __itt_pause();

    pcout << "After " << numIter << " iterations: "
          << (T) (numCells*numIter) /
             global::timer("benchmark").getTime() / 1.e6
          << " Mega site updates per second.\n";
    pcout << "Running time (s) = "<< global::timer("benchmark").getTime() << "\n\n";
    // global::profiler().writeReport();

    delete boundaryCondition;
}
