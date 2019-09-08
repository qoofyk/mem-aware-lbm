/** \file
  * Flow in a lid-driven 2D cavity. Benchmark case
**/

#include "palabos2D.h"
#include "palabos2D.hh"   // include full template code
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

void cavitySetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D topLid = Box2D(0, nx-1, ny-1, ny-1);
    Box2D everythingButTopLid = Box2D(0, nx-1, 0, ny-2);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,2>((T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,2>(u,(T)0.) );
    setBoundaryVelocity(lattice, topLid, Array<T,2>(u, 0.) );

    lattice.initialize();
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    plint N;
    try {
        global::argv(1).read(N);
    }
    catch(...)
    {
        pcout << "Wrong parameters. The syntax is " << std::endl;
        pcout << argv[0] << " N" << std::endl;
        pcout << "where N is the resolution. The benchmark cases published " << std::endl;
        pcout << "on the Palabos Wiki use N=100, N=400, N=1000, or N=4000." << std::endl;
        exit(1);
    }

    pcout << "Starting benchmark with " << N+1 << "x" << N+1 << " grid points "
          << "(approx. 2 minutes on modern processors)." << std::endl;


    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 1.,    // Re
            N,         // N
            1.,        // lx
            1.         // ly
    );


    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), 
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
    plint numIter = std::max( (plint)3,
                              (plint)(estimateSus*wishNumSeconds/numCells+0.5));

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    // Run the benchmark once "to warm up the machine".
    for (plint iT=0; iT<numIter; ++iT) {
        lattice.collideAndStream();
    }

    // Run the benchmark for good.
    global::timer("benchmark").start();
    global::profiler().turnOn();
    for (plint iT=0; iT<numIter; ++iT) {
        lattice.collideAndStream();
    }

    pcout << "numCells: " << numCells << " After " << numIter << " iterations: "
          << (T) (numCells*numIter) /
             global::timer("benchmark").getTime() / 1.e6
          << " Mega site updates per second." << std::endl << std::endl;

    global::profiler().writeReport();

    delete boundaryCondition;
}
