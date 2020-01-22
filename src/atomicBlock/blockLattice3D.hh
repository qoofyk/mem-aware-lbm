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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_3D_HH
#define BLOCK_LATTICE_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "core/plbTimer.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/util.h"
#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbProfiler.h"
#include <algorithm>
#include <typeinfo>
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace plb {

// Class BlockLattice3D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::BlockLattice3D (
        plint nx_, plint ny_, plint nz_,
        Dynamics<T,Descriptor>* backgroundDynamics_ )
   :  AtomicBlock3D(nx_, ny_, nz_, new BlockLatticeDataTransfer3D<T,Descriptor>()),
      backgroundDynamics(backgroundDynamics_)
{
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    // Allocate memory, and initialize dynamics.
    allocateAndInitialize();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                grid[iX][iY][iZ].attributeDynamics(backgroundDynamics);
            }
        }
    }
    // Attribute default value to the standard statistics (average uSqr,
    //   max uSqr, average rho). These have previously been subscribed
    //   in the constructor of BlockLatticeBase3D.
    std::vector<double> average, sum, max;
    std::vector<plint> intSum;
    average.push_back(Descriptor<double>::rhoBar(1.));
                            // default average rho to 1, to avoid division by
                            // zero in constRhoBGK and related models
    average.push_back(0.);  // default average uSqr to 0
    max.push_back(0.);      // default max uSqr to 0
    plint numCells = 1;     // pretend fictitious cell to evaluate statistics
    this->getInternalStatistics().evaluate (average, sum, max, intSum, numCells);
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::~BlockLattice3D()
{
    global::plbCounter("MEMORY_LATTICE").increment(-allocatedMemory());
    releaseMemory();
}

/** The whole data of the lattice is duplicated. This includes
 * both particle distribution function and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::BlockLattice3D(BlockLattice3D<T,Descriptor> const& rhs)
    : BlockLatticeBase3D<T,Descriptor>(rhs),
      AtomicBlock3D(rhs),
      backgroundDynamics(rhs.backgroundDynamics->clone())
{
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    allocateAndInitialize();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                Cell<T,Descriptor>& cell = grid[iX][iY][iZ];
                // Assign cell from rhs
                cell = rhs.grid[iX][iY][iZ];
                // Get an independent clone of the dynamics,
                //   or assign backgroundDynamics
                if (&cell.getDynamics()==rhs.backgroundDynamics) {
                    cell.attributeDynamics(backgroundDynamics);
                }
                else {
                    cell.attributeDynamics(cell.getDynamics().clone());
                }
            }
        }
    }
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated. This includes both particle distribution function
 * and external fields. 
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>& BlockLattice3D<T,Descriptor>::operator= (
        BlockLattice3D<T,Descriptor> const& rhs )
{
    BlockLattice3D<T,Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the 
 * lattice are copied, and not the lattice itself.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::swap(BlockLattice3D& rhs) {
    global::plbCounter("MEMORY_LATTICE").increment(-allocatedMemory());
    BlockLatticeBase3D<T,Descriptor>::swap(rhs);
    AtomicBlock3D::swap(rhs);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::specifyStatisticsStatus(Box3D domain, bool status) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collide(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].collide(this->getInternalStatistics());
                grid[iX][iY][iZ].revert();
            }
        }
    }
}

/** \sa collide(int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collide() {
    collide(this->getBoundingBox());
}

/** The distribution function never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should 
 * be streamed outside are simply left untouched.
 * The finalization of an iteration step is not automatically executed,
 * as it is in the method stream(). If you want it to be executed, you
 * must explicitly call the methods finalizeIteration() and
 * executeInternalProcessors().
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::stream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    bulkStream(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                     domain.y0+vicinity,domain.y1-vicinity,
                     domain.z0+vicinity,domain.z1-vicinity) );

    boundaryStream(domain, Box3D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z0,domain.z0+vicinity-1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z1-vicinity+1,domain.z1) );
}

/** At the end of this method, the methods finalizeIteration()
 * and executeInternalProcessors() are automatically invoked.
 * \sa stream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::stream()
{
    stream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
 * because memory is traversed only once instead of twice.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

    static const plint vicinity = Descriptor<T>::vicinity;

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D3Q19)
    collide(Box3D(domain.x0,domain.x0+vicinity-1,
                  domain.y0,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x1-vicinity+1,domain.x1,
                  domain.y0,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0,domain.y0+vicinity-1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y1-vicinity+1,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0+vicinity,domain.y1-vicinity,
                  domain.z0,domain.z0+vicinity-1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0+vicinity,domain.y1-vicinity,
                  domain.z1-vicinity+1,domain.z1) );

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                               domain.y0+vicinity,domain.y1-vicinity,
                               domain.z0+vicinity,domain.z1-vicinity) );

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(domain, Box3D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z0,domain.z0+vicinity-1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z1-vicinity+1,domain.z1) );
    global::profiler().stop("collStream");
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideRevertAndBoundSwapStream(Box3D bound, Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].collide(this->getInternalStatistics());
                grid[iX][iY][iZ].revert();

                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if ( nextX>=bound.x0 && nextX<=bound.x1 &&
                         nextY>=bound.y0 && nextY<=bound.y1 &&
                         nextZ>=bound.z0 && nextZ<=bound.z1 )
                    {
                        std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                                  grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideRevertAndBoundSwapStream(Box3D bound, plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    grid[iX][iY][iZ].collide(this->getInternalStatistics());
    grid[iX][iY][iZ].revert();

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                      grid[nextX][nextY][nextZ][iPop]);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::boundSwapStream(Box3D bound, plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                      grid[nextX][nextY][nextZ][iPop]);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::swapStream(plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];

        std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                  grid[nextX][nextY][nextZ][iPop]);
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_init(Box3D domain){

  static const plint vicinity = Descriptor<T>::vicinity;

  // First, do the collision on cells on bottom surface
  // (x0, x0, y0, y1, z0, z1) ==> bottom
  // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
  //   domain.x0, domain.x0, domain.y0, domain.y1, domain.z0, domain.z1);
  collideRevertAndBoundSwapStream(domain, Box3D(domain.x0, domain.x0+vicinity-1,
                domain.y0, domain.y1, domain.z0, domain.z1) );

  // compute x0+1 y-z surface
  plint iX = domain.x0+1;
  // first line of y-z surface
  plint iY = domain.y0;
  collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
              iY, iY, domain.z0, domain.z1) );
  // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
  //   iX, iX, iY, iY, domain.z0, domain.z1);

  // line [y0+1，y1-1] of each y-z surface
  for (iY = domain.y0+1; iY < domain.y1; ++iY) {

    // starting point on each line iY
    plint iZ = domain.z0;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
      // first collision on (iX, iY, iZ)
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      // second collision on (iX-1, iY-1, iZ-1)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
    }

    // ending point on each line iY
    iZ = domain.z1;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);
    // second collision on (iX-1, iY-1, iZ-1)
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
    // second collision on (iX-1, iY-1, iZ)
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
  }

  // last line of each y-z surface
  iY = domain.y1;
  collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
              iY, iY,
              domain.z0, domain.z1) );
  // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
  //   iX, iX, iY, iY, domain.z0, domain.z1);

  // second collision and stream on line y1-1 & y1 on x0
  collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
              iY-1, iY, domain.z0, domain.z1) );
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2_2nd_CollideAndStream(Box3D domain, plint iX, plint iY, plint iZ){
  if (iY == domain.y0 || iZ == domain.z0){
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
  }
  else{
    grid[iX][iY][iZ].collide(this->getInternalStatistics());
    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
  }
}

#if 0
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk(Box3D domain){

  for (plint iX = domain.x0+2; iX < domain.x1; ++iX) {

    /* 1. first line of each y-z surface */
    plint iY = domain.y0;
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                iY, iY, domain.z0, domain.z1) );
    // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
    //   iX, iX, iY, iY, domain.z0, domain.z1);

    /* 2. line [y0+1，y1-1] of each y-z surface */
    for (iY = domain.y0+1; iY <= domain.y1; ++iY) {

      // starting point on each line z
      plint iZ = domain.z0;
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
        // first collision on (iX, iY, iZ)
        if ( iY != domain.y1){
          grid[iX][iY][iZ].collide(this->getInternalStatistics());
          latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
          // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);
        }
        else{
          collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
          // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);
        }

        // second collision on (iX-1, iY-1, iZ-1)
        step2_2nd_CollideAndStream(domain, iX-1, iY-1, iZ-1);
      }

      // ending point on each line z
      iZ = domain.z1;
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      // second collision on (iX-1, iY-1, iZ-1)
      step2_2nd_CollideAndStream(domain, iX-1, iY-1, iZ-1);
      // second collision on (iX-1, iY-1, iZ)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
    }

    /* 3. last line of each y-z surface */
    // second collision and boundaryStream
    collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
                domain.y1, domain.y1, domain.z0, domain.z1) );
  }
}
#endif

#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk(Box3D domain){

  for (plint iX = domain.x0+2; iX < domain.x1; ++iX) {

    /* 1. first line of each y-z surface */
    plint iY = domain.y0;
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                iY, iY, domain.z0, domain.z1) );
    // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
    //   iX, iX, iY, iY, domain.z0, domain.z1);

    // iY = domain.y0+1;
    ++iY;
    plint iZ = domain.z0;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
      // first collision on (iX, iY, iZ)
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
      
      // second collision on (iX-1, iY-1, iZ-1)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
    }
    // iZ = domain.z1;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // second collision on last two points
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);

    /* 2. line [y0+2，y1-1] of each y-z surface */
    for (iY = domain.y0+2; iY < domain.y1; ++iY) {

      // starting point on each line z
      iZ = domain.z0;
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      //iZ = domain.z0+1;
      ++iZ;
      // first collision on (iX, iY, iZ)
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
      // second collision on (iX-1, iY-1, iZ-1)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);

      for (iZ = domain.z0+2; iZ < domain.z1; ++iZ ){
        // first collision on (iX, iY, iZ)
        grid[iX][iY][iZ].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
        
        // second collision on (iX-1, iY-1, iZ-1)
        grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
      }

      // ending point on each line z
      // iZ = domain.z1;
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      // second collision on (iX-1, iY-1, iZ-1)
      grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
      // second collision on (iX-1, iY-1, iZ)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
    }

    // iY = domain.y1
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                iY, iY, domain.z0, domain.z1) );

    /* 3. last line of each y-z surface */
    // second collision and boundaryStream on last two lines of iX-1 surface
    // can optimize on iY-1 without boundary check
    collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
                iY-1, iY, domain.z0, domain.z1) );
  }
}
#endif

plint thread_block;
#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk_omp(Box3D domain){

  // printf("domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
  //   domain.x0, domain.x0, domain.y0, domain.y1, domain.z0, domain.z1);

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  plint iX, iY, iZ, tid;

  // tid = omp_get_thread_num();
  // printf("tid: %ld, thread_block=%ld\n", tid, thread_block);

  // ---------------1. compute thread boundaries surface ---------------------
  // 1st Collide_Revert
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0+thread_block+1; iX < domain.x1; iX += thread_block){
    // printf("tid: %ld, thread_block=%ld, iX=%ld, (%ld, %ld), (%ld, %ld)\n", 
    //   tid, thread_block, iX, domain.y0, domain.y1, domain.z0, domain.z1);
    for (iY = domain.y0; iY <= domain.y1; ++iY)
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        grid[iX][iY][iZ].collide(this->getInternalStatistics());
        grid[iX][iY][iZ].revert();
      }
  }

  // ---------------2. compute bulk ---------------------
  #pragma omp for private(iX, iY, iZ) schedule(static, thread_block)
  for (iX = domain.x0+2; iX < domain.x1; ++iX) {
    int surface_id = (iX-2) % thread_block;

    // 2.1 -- x0+2, 1st surface, case-1
    if (surface_id == 1){
      /* 2.1.1 line y0 */
      iY = domain.y0;
      collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                  iY, iY, domain.z0, domain.z1) );

      /* 2.1.2 line [y0+1，y1-1]  */
      for (iY = domain.y0+1; iY < domain.y1; ++iY) {
        // starting point on each iY line
        iZ = domain.z0;
        collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

        for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
          // first collision on (iX, iY, iZ)
          grid[iX][iY][iZ].collide(this->getInternalStatistics());
          latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
        }

        // ending point on each line z
        // iZ = domain.z1;
        collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      }

      /*2.1.3 line y1*/
      // iY = domain.y1
      collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                  iY, iY, domain.z0, domain.z1) );
      
      continue;
    }

    // 2.2 -- x0+3, 2nd surface, case-2
    if (surface_id == 2){
      /* 2.2.1 line y0 */
      iY = domain.y0;
      collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                  iY, iY, domain.z0, domain.z1) );

      /* 2.2.2 line [y0+1，y1-1]  */
      for (iY = domain.y0+1; iY < domain.y1; ++iY) {
        // starting point z0 on each iY line
        iZ = domain.z0;
        collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

        for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
          // first collision on (iX, iY, iZ)
          grid[iX][iY][iZ].collide(this->getInternalStatistics());
          latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);

          // second collide_Revert, different from case-1
          grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
          grid[iX-1][iY-1][iZ-1].revert();
        }

        // ending point on each line z
        // iZ = domain.z1;
        collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

        // second collide_Revert on last two points of iX-1, iY-1
        grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
        grid[iX-1][iY-1][iZ-1].revert();
        grid[iX-1][iY-1][iZ].collide(this->getInternalStatistics());
        grid[iX-1][iY-1][iZ].revert();
      }

      /*2.2.3 line y1*/
      // iY = domain.y1
      collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                  iY, iY, domain.z0, domain.z1) );
      // second collision and revert on last two lines of iX-1 surface
      collide(Box3D(iX-1, iX-1, iY-1, iY, domain.z0, domain.z1) );

      continue;
    }

    // 2.3 last surface (thread boundaries surface), case-3
    if (surface_id == 0){
      /* 2.3.1 line y0 */ // No Collide, since we have collide at the begining
      iY = domain.y0;
      boundaryStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1));

      /* 2.3.2 line y0+1 */
      ++iY; // iY = y0+1;
      iZ = domain.z0;
      boundSwapStream(domain, iX, iY, iZ);

      for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
        // first collision on (iX, iY, iZ)
        swapStream(iX, iY, iZ);

        // second collision on (iX-1, iY-1, iZ-1)
        collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
      }
      // iZ = domain.z1;
      boundSwapStream(domain, iX, iY, iZ);
      
      // second collision on last two points of iX-1, iY-1
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1); //iY-1=y0
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
      /*------------------------ Finish 2.3.2 -----------------------------*/

      /* 2.3.3 line [y0+2，y1-1] of each y-z surface */
      for (iY = domain.y0+2; iY < domain.y1; ++iY){
        // starting point z0 on each iY line
        iZ = domain.z0;
        boundSwapStream(domain, iX, iY, iZ);

        //iZ = domain.z0+1;
        ++iZ;
        // first swap stream on (iX, iY, iZ)
        swapStream(iX, iY, iZ);
        
        // second collision on (iX-1, iY-1, iZ-1)
        collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);

        for (iZ = domain.z0+2; iZ < domain.z1; ++iZ ){
          // first swap & stream on (iX, iY, iZ)
          swapStream(iX, iY, iZ);

          // second collide, swap, Stream3D on iX-1, iY-1, iZ-1
          grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
          latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
        }

        // ending point on each line z
        // iZ = domain.z1;
        boundSwapStream(domain, iX, iY, iZ);

        // second collision on (iX-1, iY-1, iZ-1)
        grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
        // second collision on (iX-1, iY-1, iZ)
        collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
      }
      /*-------------------------- Finish 2.3.3 -------------------------*/

      /* 2.3.4 line y1 */
      // iY = domain.y1;
      boundaryStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1) );

      // second collision and boundaryStream on last two lines of iX-1 surface
      // can optimize on iY-1 without boundary check
      collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
                  iY-1, iY, domain.z0, domain.z1) );
      
      continue;
    }

    /* 2.4 bulk computation, case-4
      this case compute from surface [x0+4 + k * thread_block, x0+thread_block + k * thread_block]
      first thread block = [x0+2, x0+1thread_block]
      second = [x0+1 + thread_block + 1, x0+1 + 2 * thread_block]
      e.g., 1, 2, [3, 4, 5, 6], [7, 8, 9, 10], 11 
      this case will compute iX = 5 & 9
    */

    /* 2.4.1 first line of each y-z surface */
    iY = domain.y0; // first
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                iY, iY, domain.z0, domain.z1) );

    /* 2.4.2 second line of each y-z surface */
    // iY = domain.y0+1;
    ++iY;

    iZ = domain.z0; // 1st
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

    for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
      // 1st collision on (iX, iY, iZ)
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);

      // 2nd collision on (iX-1, iY-1, iZ-1)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
    }
    // iZ = domain.z1; // 1st
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

    // second collision on last two points
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1); //iY-1=y0
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
    /*-------------------------- Finish 2.4.2 -------------------------*/

    /* 2.4.3 line [y0+2，y1-1] of each y-z surface */
    for (iY = domain.y0+2; iY < domain.y1; ++iY) {

      // starting point on each line z
      iZ = domain.z0;
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

      //iZ = domain.z0+1;
      ++iZ;
      // first collision on (iX, iY, iZ)
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
      
      // second collision on (iX-1, iY-1, iZ-1)
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);

      for (iZ = domain.z0+2; iZ < domain.z1; ++iZ ){
        // first collision on (iX, iY, iZ)
        grid[iX][iY][iZ].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
        
        // second collision on (iX-1, iY-1, iZ-1)
        grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
      }

      // ending point on each line z
      // iZ = domain.z1; // 1st
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

      // 2nd collision on last two points of (iX-1, iY-1)
      grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
      collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
    }
    /*-------------------------- Finish 2.4.2 -------------------------*/

    /* 2.4.4 last line of each y-z surface */
    // iY = domain.y1; // 1st
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
                iY, iY, domain.z0, domain.z1) );

    // 2nd collision and boundaryStream on last two lines of iX-1 surface
    collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
                iY-1, iY, domain.z0, domain.z1) );
  }
  /*-------------------------- Finish 2 -------------------------*/

  // ---------------3. compute 2nd on rest surface---------------------
  /* 2nd collide Stream on x0+1 & thread_boundaries, case-3.1
   except the last boundary in the last thread block. 
   since we will handle its 2nd computation in the step2CollideAndStream_end()
  */
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0+1; iX < (domain.x1 - thread_block); iX += thread_block){
    iY = domain.y0; // 2nd
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1) );
    
    for (iY = domain.y0+1; iY < domain.y1; ++iY){
      iZ = domain.z0; // 2nd
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);

      for (iZ = domain.z0+1; iZ < domain.z1; ++iZ){ // 2nd
        grid[iX][iY][iZ].collide(this->getInternalStatistics());
        latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
      }

      // iZ = domain.z1; // 2nd
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    }

    // iY = domain.y1; // 2nd
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1) );
  }

  // 2nd swap Stream on 1st surface of every thread block, e.g. x0+2, x0+2+thread_block, case-3.2
  // iX = 3, 7, ... when thread_block = 4
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0+2; iX < domain.x1; iX += thread_block){
    iY = domain.y0; // 2nd
    boundaryStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1) );

    for (iY = domain.y0+1; iY < domain.y1; ++iY){
      iZ = domain.z0; // 2nd
      boundSwapStream(domain, iX, iY, iZ);

      for (iZ = domain.z0+1; iZ < domain.z1; ++iZ){ // 2nd
        swapStream(iX, iY, iZ);
      }

      // iZ = domain.z1
      boundSwapStream(domain, iX, iY, iZ);
    }

    // iY = domain.y1; // 2nd
    boundaryStream(domain, Box3D(iX, iX, iY, iY, domain.z0, domain.z1) );
  }
  /*-------------------------- Finish 3 -------------------------*/
}
#endif // end _OPENMP
}
#endif

#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk_blockwise(Box3D domain){

  // Make sure domain is contained within current lattice
  //PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();
  // printf("blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n", 
  //     blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);

  for (plint outerX=domain.x0+2; outerX<domain.x1; outerX+=blockSize) {
    // printf("outerX=%ld\n", outerX);
    for (plint outerY=domain.y0; outerY<=domain.y1+blockSize-1; outerY+=blockSize) {
      // printf("outerY=%ld\n", outerY);
      for (plint outerZ=domain.z0; outerZ<=domain.z1+2*(blockSize-1); outerZ+=blockSize) {
        // printf("outerZ=%ld\n", outerZ);
        // Inner loops.
        plint dx = 0;
        for (plint innerX=outerX; innerX <= std::min(outerX+blockSize-1, domain.x1-1);
          ++innerX, ++dx)
        {
          // Y-index is shifted in negative direction at each x-increment. to ensure
          //   that only post-collision cells are accessed during the swap-operation
          //   of the streaming.
          plint minY = outerY-dx;
          plint maxY = minY+blockSize-1;
          plint dy = 0;
          // printf("innerX=%ld, dx=%ld, dy=%ld, minY=%ld, maxY=%ld\n", 
          //   innerX, dx, dy, minY, maxY);

          for (plint innerY=std::max(minY,domain.y0);
            innerY <= std::min(maxY, domain.y1);
            ++innerY, ++dy)
          {
            // Z-index is shifted in negative direction at each x-increment. and at each
            //    y-increment, to ensure that only post-collision cells are accessed during
            //    the swap-operation of the streaming.
            plint minZ = outerZ-dx-dy;
            plint maxZ = minZ+blockSize-1;
            // printf("innerY=%ld, dx=%ld, dy=%ld, minY=%ld, maxY=%ld, minZ=%ld, maxZ=%ld\n", 
                          // innerY, dx, dy, minY, maxY, minZ, maxZ);

            for (plint innerZ=std::max(minZ,domain.z0); 
              innerZ <= std::min(maxZ, domain.z1);
              ++innerZ)
            {
              // printf("inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);

              // first on y=y0, z=z0
              if (innerY == domain.y0 || innerZ == domain.z0){
                // printf("case-1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                continue;
              }

              // first on y=y0+1
              if (innerY == domain.y0+1){
                // printf("case-2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);

                // second
                collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);

                if (innerZ == domain.z1){
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                }
                continue;
              }

              // first on y=y1
              if (innerY == domain.y1){
                // printf("case-3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                if (innerZ == domain.z0+1){
                  // printf("case-3.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  continue;
                }
                else{
                  // printf("case-3.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  // second
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX-1, innerY-1, innerZ-1);

                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2);
                }

                // z=z1, one more things to do
                if (innerZ == domain.z1){
                  // printf("case-3.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1);
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ);
                }

                continue;
              }

              if (innerZ == domain.z0+1){
                // printf("case-4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                grid[innerX][innerY][innerZ].collide (
                      this->getInternalStatistics() );
                latticeTemplates<T,Descriptor>::swapAndStream3D (
                      grid, innerX, innerY, innerZ );

                // second
                collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                continue;
              }

              if (innerZ == domain.z1){
                // printf("case-5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                // second collide
                grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                latticeTemplates<T,Descriptor>::swapAndStream3D(
                      grid, innerX-1, innerY-1, innerZ-1);
                collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                continue;
              }

              // printf("case-6 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
              // first Collide the cell.
              grid[innerX][innerY][innerZ].collide (
                      this->getInternalStatistics() );
              // Swap the populations on the cell, and then with post-collision
              //   neighboring cell, to perform the streaming step.
              latticeTemplates<T,Descriptor>::swapAndStream3D (
                      grid, innerX, innerY, innerZ );

              // second collide
              grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
              latticeTemplates<T,Descriptor>::swapAndStream3D(
                      grid, innerX-1, innerY-1, innerZ-1);
            }
          }
        }
      }
    }
  }
}
#endif

#if 0
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_end(Box3D domain){

  // last surface domain.x1
  plint iX =  domain.x1;
  /* 1. first line of each y-z surface */
  plint iY = domain.y0;
  collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
              iY, iY,
              domain.z0, domain.z1) );
  // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n", 
  //     iX, iX, iY, iY, domain.z0, domain.z1);

  /* 2. line [y0+1，y1-1] of each y-z surface */
  for (iY = domain.y0+1; iY <= domain.y1; ++iY) {

    // starting point on each line z
    plint iZ = domain.z0;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    for (iZ = domain.z0+1; iZ < domain.z1; ++iZ ){
      // first collision on (iX, iY, iZ)
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      // second collision on (iX-1, iY-1, iZ-1)
      step2_2nd_CollideAndStream(domain, iX-1, iY-1, iZ-1);
    }

    // ending point on each line z
    iZ = domain.z1;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    // second collision on (iX-1, iY-1, iZ-1)
    step2_2nd_CollideAndStream(domain, iX-1, iY-1, iZ-1);
    // second collision on (iX-1, iY-1, iZ)
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
  }

  /* 3. last line of each y-z surface */
  // second collision and boundaryStream
  collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
              domain.y1, domain.y1, domain.z0, domain.z1) );

  // second collision on top
  collideRevertAndBoundSwapStream(domain, Box3D(domain.x1, domain.x1,
                domain.y0,domain.y1, domain.z0,domain.z1) );
}
#endif

#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_end(Box3D domain){

  // last surface domain.x1
  plint iX =  domain.x1;
  /* 1. first line of each y-z surface */
  plint iY = domain.y0;
  collideRevertAndBoundSwapStream(domain, Box3D(iX, iX,
              iY, iY,
              domain.z0, domain.z1) );
  // printf("1st-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
  //     iX, iX, iY, iY, domain.z0, domain.z1);

  // iY = y0+1
  ++iY;
  plint iZ = domain.z0;
  collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
  // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

  for (iZ = domain.z0+1; iZ <= domain.z1; ++iZ ){
    // first collision on (iX, iY, iZ)
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    // second collision on (iX-1, iY-1, iZ-1)
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);
  }

  // ending point on each line z
  iZ = domain.z1;
  collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
  // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

  /* 2. line [y0+1，y1-1] of each y-z surface */
  for (iY = domain.y0+2; iY <= domain.y1; ++iY) {

    // starting point on each line z
    iZ = domain.z0;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    ++iZ;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ-1);

    for (iZ = domain.z0+2; iZ < domain.z1; ++iZ ){
      // first collision on (iX, iY, iZ)
      collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

      // second collision on (iX-1, iY-1, iZ-1)
      grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
    }

    // ending point on each line z
    // iZ = domain.z1;
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
    // printf("1st-(%ld, %ld, %ld)\n", iX, iY, iZ);

    // second collision on (iX-1, iY-1, iZ-1)
    grid[iX-1][iY-1][iZ-1].collide(this->getInternalStatistics());
    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX-1, iY-1, iZ-1);
    // second collision on (iX-1, iY-1, iZ)
    collideRevertAndBoundSwapStream(domain, iX-1, iY-1, iZ);
  }

  /* 3. last line of each y-z surface */
  // second collision and boundaryStream
  collideRevertAndBoundSwapStream(domain, Box3D(iX-1, iX-1,
              domain.y1, domain.y1, domain.z0, domain.z1) );

  // second collision on top
  collideRevertAndBoundSwapStream(domain, Box3D(domain.x1, domain.x1,
                domain.y0,domain.y1, domain.z0,domain.z1) );
}
#endif

/** This operation is 2 step collision and stream using 1 buffer
 * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
 * because memory is traversed only once instead of twice.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

    // 2 collideAndStream on x0; 1 collideAndStream on x0+1
    step2CollideAndStream_init(domain);

    // Then bulk [x0+2, x0-1]
    //if(domain.x1-domain.x0 <= 250){
      // step2CollideAndStream_bulk(domain);
      step2CollideAndStream_bulk_omp(domain);
    //}
    //else{
      //step2CollideAndStream_bulk_blockwise(domain);
    //}

    // step2CollideAndStream_bulk
    step2CollideAndStream_end(domain);

    global::profiler().stop("collStream");
}

/** At the end of this method, finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa collideAndStream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideAndStream() {
    collideAndStream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::incrementTime() {
    this->getTimeCounter().incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::allocateAndInitialize() {
    this->getInternalStatistics().subscribeAverage(); // Subscribe average rho-bar
    this->getInternalStatistics().subscribeAverage(); // Subscribe average uSqr
    this->getInternalStatistics().subscribeMax();     // Subscribe max uSqr

    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    rawData = new Cell<T,Descriptor> [nx*ny*nz];
    grid    = new Cell<T,Descriptor>** [nx];
    for (plint iX=0; iX<nx; ++iX) {
        grid[iX] = new Cell<T,Descriptor>* [ny];
        for (plint iY=0; iY<ny; ++iY) {
            grid[iX][iY] = rawData + nz*(iY+ny*iX);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::releaseMemory() {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                Dynamics<T,Descriptor>* dynamics = &grid[iX][iY][iZ].getDynamics();
                if (dynamics != backgroundDynamics) {
                    delete dynamics;
                }
            }
        }
    }
    delete backgroundDynamics;
    delete [] rawData;
    for (plint iX=0; iX<nx; ++iX) {
        delete [] grid[iX];
    }
    delete [] grid;
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::attributeDynamics (
        plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics )
{
    Dynamics<T,Descriptor>* previousDynamics = &grid[iX][iY][iZ].getDynamics();
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX][iY][iZ].attributeDynamics(dynamics);
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& BlockLattice3D<T,Descriptor>::getBackgroundDynamics() {
    return *backgroundDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& BlockLattice3D<T,Descriptor>::getBackgroundDynamics() const {
    return *backgroundDynamics;
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::resetDynamics(Dynamics<T,Descriptor> const& dynamics) {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                attributeDynamics(iX,iY,iZ, dynamics.clone());
            }
        }
    }
}

/** This method is slower than bulkStream(int,int,int,int), because it must
 * be verified which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::boundaryStream(Box3D bound, Box3D domain) {
    // Make sure bound is contained within current lattice
    PLB_PRECONDITION( contained(bound, this->getBoundingBox()) );
    // Make sure domain is contained within bound
    PLB_PRECONDITION( contained(domain, bound) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if ( nextX>=bound.x0 && nextX<=bound.x1 &&
                         nextY>=bound.y0 && nextY<=bound.y1 &&
                         nextZ>=bound.z0 && nextZ<=bound.z1 )
                    {
                        std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                                  grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

/** This method is faster than boundaryStream(int,int,int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::bulkStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                              grid[nextX][nextY][nextZ][iPop]);
                }
            }
        }
    }
}

/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::bulkCollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

     //if (Descriptor<T>::q==15 || Descriptor<T>::q==19) {
     //if (Descriptor<T>::q==19) {
         // On nearest-neighbor lattice, use the cache-efficient
         //   version of collidAndStream.
         //blockwiseBulkCollideAndStream(domain);
     //}
     //else {
        // Otherwise, use the straightforward implementation.
        //   Note that at some point, we should implement the cache-efficient
        //   version for extended lattices as well.
        linearBulkCollideAndStream(domain);
     //}
}


/** Straightforward implementation which works for all kinds of lattices,
 *  not only nearest-neighbor.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::linearBulkCollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].collide(this->getInternalStatistics());
                latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
            }
        }
    }
}


/** Sophisticated implementation which improves cache usage through block-wise
 *  loops. For now, this works only with nearest-neighbor lattices. On extended
 *  lattices, the naive version "linearBulkCollideAndStream" is used.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::blockwiseBulkCollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    // Outer loops.
    for (plint outerX=domain.x0; outerX<=domain.x1; outerX+=blockSize) {
        for (plint outerY=domain.y0; outerY<=domain.y1+blockSize-1; outerY+=blockSize) {
            for (plint outerZ=domain.z0; outerZ<=domain.z1+2*(blockSize-1); outerZ+=blockSize) {
                // Inner loops.
                plint dx = 0;
                for (plint innerX=outerX;
                     innerX <= std::min(outerX+blockSize-1, domain.x1);
                     ++innerX, ++dx)
                {
                    // Y-index is shifted in negative direction at each x-increment. to ensure
                    //   that only post-collision cells are accessed during the swap-operation
                    //   of the streaming.
                    plint minY = outerY-dx;
                    plint maxY = minY+blockSize-1;
                    plint dy = 0;
                    for (plint innerY=std::max(minY,domain.y0);
                         innerY <= std::min(maxY, domain.y1);
                         ++innerY, ++dy)
                    {
                        // Z-index is shifted in negative direction at each x-increment. and at each
                        //    y-increment, to ensure that only post-collision cells are accessed during
                        //    the swap-operation of the streaming.
                        plint minZ = outerZ-dx-dy;
                        plint maxZ = minZ+blockSize-1;
                        for (plint innerZ=std::max(minZ,domain.z0);
                             innerZ <= std::min(maxZ, domain.z1);
                             ++innerZ)
                        {
                            // Collide the cell.
                            grid[innerX][innerY][innerZ].collide (
                                    this->getInternalStatistics() );
                            // Swap the populations on the cell, and then with post-collision
                            //   neighboring cell, to perform the streaming step.
                            latticeTemplates<T,Descriptor>::swapAndStream3D (
                                    grid, innerX, innerY, innerZ );
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::implementPeriodicity() {
    static const plint vicinity = Descriptor<T>::vicinity;
    plint maxX = this->getNx()-1;
    plint maxY = this->getNy()-1;
    plint maxZ = this->getNz()-1;
    // Periodicity of planes orthogonal to x-axis.
    periodicDomain(Box3D(-vicinity,-1, 0,maxY, 0,maxZ));
    // Periodicity of planes orthogonal to y-axis.
    periodicDomain(Box3D(0,maxX, -vicinity,-1, 0,maxZ));
    // Periodicity of planes orthogonal to z-axis.
    periodicDomain(Box3D(0,maxX, 0,maxY, -vicinity,-1));

    // Periodicity of edges in y-z plane.
    // Periodicity between (-y,-z) and (+y,+z) edge.
    periodicDomain(Box3D(0,maxX, -vicinity,-1, -vicinity,-1));
    // Periodicity between (-y,+z) and (+y,-z) edge.
    periodicDomain(Box3D(0,maxX, -vicinity,-1, maxZ+1,maxZ+vicinity));

    // Periodicity of edges in x-z plane.
    // Periodicity between (-x,-z) and (+x,+z) edge.
    periodicDomain(Box3D(-vicinity,-1, 0,maxY, -vicinity,-1));
    // Periodicity between (-x,+z) and (+x,-z) edge.
    periodicDomain(Box3D(maxX+1,maxX+vicinity, 0,maxY, -vicinity,-1));

    // Periodicity of edges in x-y plane.
    // Periodicity between (-x,-y) and (+x,+y) edge.
    periodicDomain(Box3D(-vicinity,-1, -vicinity,-1, 0,maxZ));
    // Periodicity between (-x,+y) and (+x,-y) edge.
    periodicDomain(Box3D(-vicinity,-1, maxY+1,maxY+vicinity, 0,maxZ));

    // Periodicity of corners.
    // Periodicity of (+1,+1,+1) and (-1,-1,-1)
    periodicDomain(Box3D(maxX+1,maxX+vicinity, maxY+1,maxY+vicinity, maxZ+1,maxZ+vicinity));
    // Periodicity of (+1,+1,-1) and (-1,-1,+1)
    periodicDomain(Box3D(maxX+1,maxX+vicinity, maxY+1,maxY+vicinity, -vicinity,-1));
    // Periodicity of (+1,-1,+1) and (-1,+1,-1)
    periodicDomain(Box3D(maxX+1,maxX+vicinity, -vicinity,-1, maxZ+1,maxZ+vicinity));
    // Periodicity of (+1,-1,-1) and (-1,+1,+1)
    periodicDomain(Box3D(maxX+1,maxX+vicinity, -vicinity,-1, -vicinity,-1));
}

template<typename T, template<typename U> class Descriptor>
plint BlockLattice3D<T,Descriptor>::allocatedMemory() const {
    return this->getNx()*this->getNy()*this->getNz()*
           sizeof(T)* (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::periodicDomain(Box3D domain) {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                    plint prevX = iX - Descriptor<T>::c[iPop][0];
                    plint prevY = iY - Descriptor<T>::c[iPop][1];
                    plint prevZ = iZ - Descriptor<T>::c[iPop][2];

                    if ( (prevX>=0 && prevX<nx) &&
                         (prevY>=0 && prevY<ny) &&
                         (prevZ>=0 && prevZ<nz) )
                    {
                        plint nextX = (iX+nx)%nx;
                        plint nextY = (iY+ny)%ny;
                        plint nextZ = (iZ+nz)%nz;
                        std::swap (
                            grid[prevX][prevY][prevZ][indexTemplates::opposite<Descriptor<T> >(iPop)],
                            grid[nextX][nextY][nextZ][iPop] );
                    }
                }
            }
        }
    }
}

////////////////////// Class BlockLatticeDataTransfer3D /////////////////////////

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer3D<T,Descriptor>::BlockLatticeDataTransfer3D()
    : lattice(0),
      constLattice(0)
{ }

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::setBlock(AtomicBlock3D& block) {
    lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(&block);
    PLB_ASSERT(lattice);
    constLattice = lattice;
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::setConstBlock(AtomicBlock3D const& block) {
    constLattice = dynamic_cast<BlockLattice3D<T,Descriptor> const*>(&block);
    PLB_ASSERT(constLattice);
}

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer3D<T,Descriptor>* BlockLatticeDataTransfer3D<T,Descriptor>::clone() const
{
    return new BlockLatticeDataTransfer3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
plint BlockLatticeDataTransfer3D<T,Descriptor>::staticCellSize() const {
    return sizeof(T)* (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::send (
        Box3D domain, std::vector<char>& buffer, modif::ModifT kind ) const
{
    PLB_PRECONDITION( constLattice );
    PLB_PRECONDITION(contained(domain, constLattice->getBoundingBox()));
    // It's the responsibility of the functions called below to allocate
    //   the right amount of memory for the buffer.
    buffer.clear();
    switch(kind) {
        case modif::staticVariables:
            send_static(domain, buffer); break;
        case modif::dynamicVariables:
            send_dynamic(domain, buffer); break;
        // Serialization is the same no matter if the dynamics object
        //   is being regenerated or not by the recipient.
        case modif::allVariables:  
        case modif::dataStructure:
            send_all(domain,buffer); break;
        default: PLB_ASSERT(false);
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::send_static (
        Box3D domain, std::vector<char>& buffer ) const
{
    PLB_PRECONDITION( constLattice );
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells()*cellSize;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes==0) return;
    buffer.resize(numBytes);

    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                constLattice->get(iX,iY,iZ).serialize(&buffer[iData]);
                iData += cellSize;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::send_dynamic (
        Box3D domain, std::vector<char>& buffer ) const
{
    PLB_PRECONDITION( constLattice );
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // The serialize function automatically reallocates memory for buffer.
                serialize(constLattice->get(iX,iY,iZ).getDynamics(), buffer);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::send_all (
        Box3D domain, std::vector<char>& buffer ) const
{
    PLB_PRECONDITION( constLattice );
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // 1. Send dynamic info (automaic allocation of buffer memory).
                serialize(constLattice->get(iX,iY,iZ).getDynamics(), buffer);
                pluint pos = buffer.size();
                // 2. Send static info (needs manual allocation of buffer memory).
                if (staticCellSize()>0) {
                    buffer.resize(pos+staticCellSize());
                    constLattice->get(iX,iY,iZ).serialize(&buffer[pos]);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive (
        Box3D domain, std::vector<char> const& buffer,
        modif::ModifT kind, std::map<int,std::string> const& foreignIds )
{
    if (kind==modif::dataStructure && !foreignIds.empty()) {
        std::map<int,int> idIndirect;
        meta::createIdIndirection<T,Descriptor>(foreignIds, idIndirect);
        receive_regenerate(domain, buffer, idIndirect);
    }
    else {
        receive(domain, buffer, kind);
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind )
{
    PLB_PRECONDITION( lattice );
    PLB_PRECONDITION(contained(domain, lattice->getBoundingBox()));
    switch(kind) {
        case modif::staticVariables:
            receive_static(domain, buffer); break;
        case modif::dynamicVariables:
            receive_dynamic(domain, buffer); break;
        case modif::allVariables:
            receive_all(domain, buffer); break;
        case modif::dataStructure:
            receive_regenerate(domain, buffer); break;
        default:
            PLB_ASSERT( false );
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive_static (
        Box3D domain, std::vector<char> const& buffer )
{
    PLB_PRECONDITION( lattice );
    PLB_PRECONDITION( (plint) buffer.size() == domain.nCells()*staticCellSize() );
    // Avoid dereferencing uninitialized pointer.
    if (buffer.empty()) return;
    plint cellSize = staticCellSize();

    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice->get(iX,iY,iZ).unSerialize(&buffer[iData]);
                iData += cellSize;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive_dynamic (
        Box3D domain, std::vector<char> const& buffer )
{
    PLB_PRECONDITION( lattice );
    pluint serializerPos = 0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // No assert is included here, because incompatible types of
                //   dynamics are detected by asserts inside HierarchicUnserializer.
                serializerPos = 
                    unserialize (
                        lattice->get(iX,iY,iZ).getDynamics(), buffer, serializerPos );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive_all (
        Box3D domain, std::vector<char> const& buffer )
{
    PLB_PRECONDITION( lattice );
    pluint posInBuffer = 0;
    plint cellSize = staticCellSize();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // 1. Unserialize dynamic data.
                posInBuffer = 
                    unserialize (
                        lattice->get(iX,iY,iZ).getDynamics(), buffer, posInBuffer );
                // 2. Unserialize static data.
                if (staticCellSize()>0) {
                    lattice->get(iX,iY,iZ).unSerialize(&buffer[posInBuffer]);
                    posInBuffer += cellSize;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive_regenerate (
        Box3D domain, std::vector<char> const& buffer, std::map<int,int> const& idIndirect )
{
    PLB_PRECONDITION( lattice );
    pluint posInBuffer = 0;
    plint cellSize = staticCellSize();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // 1. Generate dynamics object, and unserialize dynamic data.
                std::map<int,int> const* indirectPtr = idIndirect.empty() ? 0 : &idIndirect;
                HierarchicUnserializer unserializer(buffer, posInBuffer, indirectPtr);
                Dynamics<T,Descriptor>* newDynamics =
                    meta::dynamicsRegistration<T,Descriptor>().generate(unserializer);
                posInBuffer = unserializer.getCurrentPos();
                lattice->attributeDynamics(iX,iY,iZ, newDynamics);

                // 2. Unserialize static data.
                if (staticCellSize()>0) {
                    PLB_ASSERT( !buffer.empty() );
                    PLB_ASSERT( posInBuffer+cellSize<=buffer.size() );
                    lattice->get(iX,iY,iZ).unSerialize(&buffer[posInBuffer]);
                    posInBuffer += cellSize;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind )
{
    PLB_PRECONDITION( lattice );
    PLB_PRECONDITION (typeid(from) == typeid(BlockLattice3D<T,Descriptor> const&));
    PLB_PRECONDITION(contained(toDomain, lattice->getBoundingBox()));
    BlockLattice3D<T,Descriptor> const& fromLattice = (BlockLattice3D<T,Descriptor> const&) from;
    switch(kind) {
        case modif::staticVariables:
            attribute_static(toDomain, deltaX, deltaY, deltaZ, fromLattice); break;
        case modif::dynamicVariables:
            attribute_dynamic(toDomain, deltaX, deltaY, deltaZ, fromLattice); break;
        case modif::allVariables:
            attribute_all(toDomain, deltaX, deltaY, deltaZ, fromLattice); break;
        case modif::dataStructure:
            attribute_regenerate(toDomain, deltaX, deltaY, deltaZ, fromLattice); break;
        default:
            PLB_ASSERT( false );
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute_static (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T,Descriptor> const& from )
{
    PLB_PRECONDITION( lattice );
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                lattice->get(iX,iY,iZ).attributeValues (
                        from.get(iX+deltaX,iY+deltaY,iZ+deltaZ) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute_dynamic (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T,Descriptor> const& from )
{
    PLB_PRECONDITION( lattice );
    std::vector<char> serializedData;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                serializedData.clear();
                serialize (
                    from.get(iX+deltaX,iY+deltaY,iZ+deltaZ).getDynamics(),
                    serializedData );
                unserialize (
                    lattice->get(iX,iY,iZ).getDynamics(),
                    serializedData );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute_all (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T,Descriptor> const& from )
{
    PLB_PRECONDITION( lattice );
    std::vector<char> serializedData;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                // 1. Attribute dynamic content.
                serializedData.clear();
                serialize (
                    from.get(iX+deltaX,iY+deltaY,iZ+deltaZ).getDynamics(),
                    serializedData );
                unserialize (
                    lattice->get(iX,iY,iZ).getDynamics(),
                    serializedData );

                // 2. Attribute static content.
                lattice->get(iX,iY,iZ).attributeValues (
                        from.get(iX+deltaX,iY+deltaY,iZ+deltaZ) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute_regenerate (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T,Descriptor> const& from )
{
    PLB_PRECONDITION( lattice );
    std::vector<char> serializedData;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                // 1. Generate new dynamics and attribute dynamic content.
                serializedData.clear();
                serialize (
                    from.get(iX+deltaX,iY+deltaY,iZ+deltaZ).getDynamics(),
                    serializedData );
                HierarchicUnserializer unserializer(serializedData, 0);
                Dynamics<T,Descriptor>* newDynamics =
                    meta::dynamicsRegistration<T,Descriptor>().generate(unserializer);
                lattice->attributeDynamics(iX,iY,iZ, newDynamics);

                // 2. Attribute static content.
                lattice->get(iX,iY,iZ).attributeValues (
                        from.get(iX+deltaX,iY+deltaY,iZ+deltaZ) );
            }
        }
    }
}

// add by Yuankun
plint ykBlockSize;

template<typename T, template<typename U> class Descriptor>
CachePolicy3D& BlockLattice3D<T,Descriptor>::cachePolicy() {
    static CachePolicy3D cachePolicySingleton(ykBlockSize);
    return cachePolicySingleton;
}


/////////// Free Functions //////////////////////////////

template<typename T, template<typename U> class Descriptor>
double getStoredAverageDensity(BlockLattice3D<T,Descriptor> const& blockLattice) {
    return Descriptor<T>::fullRho (
               blockLattice.getInternalStatistics().getAverage (
                  LatticeStatistics::avRhoBar ) );
}

template<typename T, template<typename U> class Descriptor>
double getStoredAverageEnergy(BlockLattice3D<T,Descriptor> const& blockLattice) {
    return 0.5 * blockLattice.getInternalStatistics().getAverage (
                        LatticeStatistics::avUSqr );
}

template<typename T, template<typename U> class Descriptor>
double getStoredMaxVelocity(BlockLattice3D<T,Descriptor> const& blockLattice) {
    return std::sqrt( blockLattice.getInternalStatistics().getMax (
                             LatticeStatistics::maxUSqr ) );
}

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_HH
