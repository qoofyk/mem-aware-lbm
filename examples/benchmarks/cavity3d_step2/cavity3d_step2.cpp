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

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    // initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    // Modify by Yuankun, set init value to 0.01 to avoid 0 computation
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.01,(T)0.01,(T)0.01) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,0.,u) );

    lattice.initialize();
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    plint N;
    plint numIter;
    try {
        global::argv(1).read(N);
        global::argv(2).read(numIter);
        global::argv(3).read(ykBlockSize);
    }
    catch(...)
    {
        pcout << "Wrong parameters. The syntax is " << std::endl;
        pcout << argv[0] << " N" << std::endl;
        pcout << "where N is the resolution. The benchmark cases published " << std::endl;
        pcout << "on the Palabos Wiki use N=100, N=400, N=1000, or N=4000." << std::endl;
        exit(1);
    }

    pcout << "Starting benchmark with " << N+1 << "x" << N+1 << "x" << N+1 << " grid points "
          << "(approx. 2 minutes on modern processors)." << std::endl;


    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 1.,    // Re
            N,         // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;
    // Current cores run approximately at 5 Mega Sus.
    T estimateSus= 5.e6*numCores;
    // The benchmark should run for approximately two minutes
    // (2*60 seconds).
    T wishNumSeconds = 60.;
    plint numCells = lattice.getBoundingBox().nCells();

    // Run at least three iterations.
    // plint numIter = std::max( (plint)3,
    //                           (plint)(estimateSus*wishNumSeconds/numCells+0.5));

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

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
    // for (plint iT=0; iT<numIter; iT += 2) {
    for (plint iT=0; iT<2; iT += 2) {
        // pcout << "iT=" << iT << std::endl;
        lattice.step2collideAndStream();
    }

    // pcout << "Start bench!" << std::endl;
    // Run the benchmark for good.
    global::timer("benchmark").start();
    global::profiler().turnOn();
    for (plint iT=0; iT<numIter; iT += 2) {
        // pcout << "iT=" << iT << std::endl;
        lattice.step2collideAndStream();
    }

#if 0
    pcout << "After: Velocity norm of the box: " << endl;
    // pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, mybox)) << endl;
    for (plint iX=0; iX<=N; ++iX){
        for (plint iY=0; iY<=N; ++iY){
            Box3D line(iX, iX, iY, iY, 0, N);
            pcout << setprecision(3) << *computeVelocityNorm(*extractSubDomain(lattice, line)) << endl;
        }
    }
#endif

    pcout << "After " << numIter << " iterations: "
          << (T) (numCells*numIter) /
             global::timer("benchmark").getTime() / 1.e6
          << " Mega site updates per second." << std::endl << std::endl;
    pcout << "Running time (s) = "<< global::timer("benchmark").getTime() << std::endl;
    global::profiler().writeReport();  

    delete boundaryCondition;
}
