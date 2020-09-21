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
#include "ittnotify.h"

using namespace plb;
using namespace std;

#include "ykglobal.h"

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                  plint Nx, plint Ny, plint Nz)
{
    // pcout << "cavitySetup\n";
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
    // pcout << "boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice)\n"; // Palabos memory outrage happens here when large num of processes are used

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    // initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    // Modify by Yuankun, set init value to 0.01 to avoid 0 computation
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.01, (T)0.01, (T)0.01) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u, (T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,0.,u) );

    lattice.initialize();
}

int main(int argc, char* argv[]) {

    __itt_pause();
    plbInit(&argc, &argv);
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    // not using Palabos pymd alg
    yk_use_pymd_flag = false;

    plint N;
    plint Nx, Ny, Nz;
    plint numIter;
    plint warmUpIter;
    try {
        global::argv(1).read(N);
        global::argv(2).read(numIter);
        global::argv(3).read(ykBlockSize);
        global::argv(4).read(warmUpIter);
        global::argv(5).read(Nx);
        global::argv(6).read(Ny);
        global::argv(7).read(Nz);
    }
    catch(...)
    {
        pcout << "Wrong parameters. The syntax is " << std::endl;
        pcout << argv[0] << " N" << std::endl;
        pcout << "where N is the resolution. The benchmark cases published " << std::endl;
        pcout << "on the Palabos Wiki use N=100, N=400, N=1000, or N=4000." << std::endl;
        exit(1);
    }

    pcout << "Nx=" << Nx << ", Ny=" << Ny << ", Nz=" << Nz << ", Memory="
          << Nx * Ny * Nz * 168 / (1024*1024) << " MB, warmUpIter=" << warmUpIter 
          << ", numIter=" << numIter << ", tile=" << ykBlockSize << '\n';

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 1.,    // Re
            N,         // N, resolution
            1.,        // lx
            1.,        // ly
            1.         // lz
    );

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            // parameters.getNx(), parameters.getNy(), parameters.getNz(),
            Nx, Ny, Nz,
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );
#if 0
    int rank; // memory outrage not happen here
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Hello from Proc %d \n", rank);
#endif

    plint numCores = global::mpi().getSize();
    pcout << "NUM_THREADS=" << numCores << '\n';

    // plint numCells = lattice.getBoundingBox().nCells();
    plint numCells = Nx * Ny * Nz;

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
    // pcout << "Warm up the machine!" << std::endl;
    // Run the benchmark once "to warm up the machine".
    for (plint iT=0; iT<warmUpIter; iT += 1) {
    //for (plint iT=0; iT<2; iT += 1) {
        // pcout << "iT=" << iT << std::endl;
        lattice.collideAndStream();
    }

    // pcout << "Start bench!" << std::endl;
    // Run the benchmark for good.
    __itt_resume();
    global::timer("benchmark").start();
    // global::profiler().turnOn();
    for (plint iT=0; iT<numIter; iT += 1) {
        // pcout << "iT=" << iT << std::endl;
        lattice.collideAndStream();
    }

#ifdef SAVE
    pcout << "After: Velocity norm of the box: " << endl;
    // pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, mybox)) << endl;
    for (plint iX=0; iX<=N; ++iX){
        for (plint iY=0; iY<=N; ++iY){
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
