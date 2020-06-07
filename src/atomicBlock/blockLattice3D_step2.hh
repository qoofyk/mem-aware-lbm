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
#ifndef BLOCK_LATTICE_3D_STEP2_HH
#define BLOCK_LATTICE_3D_STEP2_HH

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
// add by Yuankun
plint thread_block;

// Class BlockLattice3D /////////////////////////

#ifdef STEP2_OMP
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

    // std::cout << "BlockLattice3D" << std::endl;

    #pragma omp parallel for default(shared) schedule(static)
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

    std::cout << "BlockLattice3D(rhs)" << std::endl;

    #pragma omp parallel for default(shared) schedule(static, thread_block)
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

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::specifyStatisticsStatus(Box3D domain, bool status) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
            }
        }
    }
}
#endif


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
  printf("Here! I am step2_3parts_seq_unroll_line\n");

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

#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk_omp(Box3D domain){
#if 0
  printf("domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
    domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);
#endif

  plint iX, iY, iZ;
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
   //int tid = omp_get_thread_num();
   //printf("tid: %ld, thread_block=%ld\n", tid, thread_block);
#if 1
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
    // collide(Box3D(iX, iX, domain.y0, domain.y1, domain.z0, domain.z1) );
  }
#endif
#if 1
  // ---------------2. compute bulk ---------------------
  #pragma omp for private(iX, iY, iZ) schedule(static, thread_block)
  for (iX = domain.x0+2; iX < domain.x1; ++iX) {
    int surface_id = (iX-2) % thread_block;
     //printf("tid: %ld, thread_block=%ld, iX=%ld\n", tid, thread_block, iX);

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
      first thread block = [x0+2, x0+1+thread_block]
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
#endif
#if 1
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
#endif
}
#endif // end _OPENMP
}
#endif

#if 1
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_bulk_blockwise(Box3D domain){
  // printf("Here! I am step2_3parts_seq_unroll_pyramid\n");

  // Make sure domain is contained within current lattice
  //PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();
  // printf("blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n",
  //     blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);

  for (plint outerX = domain.x0+2; outerX < domain.x1; outerX += blockSize) {
    // printf("outerX=%ld\n", outerX);
    for (plint outerY = domain.y0; outerY <= domain.y1+blockSize-1; outerY += blockSize) {
      // printf("outerY=%ld\n", outerY);
      for (plint outerZ = domain.z0; outerZ <= domain.z1+2*(blockSize-1); outerZ += blockSize) {
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
                if (innerZ == domain.z1) {
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                }
                else {
                  grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
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

  /* 2. line [y0+2，y1] of each y-z surface */
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
                domain.y0, domain.y1, domain.z0, domain.z1) );
}
#endif

#ifdef STEP2_WHOLE
//   /** This operation is 2 step collision and stream using 1 buffer
//  * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
//  * because memory is traversed only once instead of twice.
//  */
  #ifdef STEP2_OMP
    #ifdef STEP2_UNROLL
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
  // global::profiler().start("collStream");
  // global::profiler().increment("collStreamCells", domain.nCells());

  // Make sure domain is contained within current lattice
  PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();

#if 0
  printf("Here! I am step2_whole_omp_unroll_pyramid, blockSize=%ld, domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
    blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);
#endif

  plint iX, iY, iZ;
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  plint tid = omp_get_thread_num();
  // printf("tid: %ld, thread_block=%ld\n", tid, thread_block);
#if 1
  // ---------------1. compute thread boundaries surface ---------------------
  // 1st Collide_Revert
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0 + thread_block - 1; iX <= domain.x1; iX += thread_block){
    // printf("tid: %ld, thread_block=%ld, iX=%ld, (%ld, %ld), (%ld, %ld)\n",
    //   tid, thread_block, iX, domain.y0, domain.y1, domain.z0, domain.z1);
    for (iY = domain.y0; iY <= domain.y1; ++iY) {
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        grid[iX][iY][iZ].collide(this->getInternalStatistics());
        grid[iX][iY][iZ].revert();
      }
    }
  }
#endif
#if 1
  // ---------------2. compute bulk ---------------------
  #pragma omp for private(iX, iY, iZ) schedule(static, thread_block / blockSize)
  for (plint outerX = domain.x0; outerX <= domain.x1; outerX += blockSize) {
    // printf("tid: %ld, outerX=%ld\n", tid, outerX);
    for (plint outerY = domain.y0; outerY <= domain.y1+blockSize-1; outerY += blockSize) {
      // printf("outerY=%ld\n", outerY);
      for (plint outerZ = domain.z0; outerZ <= domain.z1+2*(blockSize-1); outerZ += blockSize) {
        // printf("outerZ=%ld\n", outerZ);
        // Inner loops.
        plint dx = 0;
        for (plint innerX = outerX; innerX <= std::min(outerX+blockSize-1, domain.x1);
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

          for (plint innerY=std::max(minY, domain.y0);
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

            for (plint innerZ=std::max(minZ, domain.z0);
              innerZ <= std::min(maxZ, domain.z1);
              ++innerZ)
            {
              // printf("tid%ld: inner(%ld, %ld, %ld)\n", tid, innerX, innerY, innerZ);
              int surface_id = innerX % thread_block;
              // Case-0. On x=x0, y=y0, z=z0, except last surface within a thread block
              if (surface_id != 0 && (innerX == domain.x0 || innerY == domain.y0 || innerZ == domain.z0)){
                // printf("case-0 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                continue;
              }
              
              // Case-1. On 1st surface of a thread block, innerx = domain.x0 + thread_block * n + 1
              if (surface_id == 1) {
                // 2.1 On x=x0, Or boundary y=y0,y1, z=z0, z1
                if (innerY == domain.y1 || innerZ == domain.z1){
                  // printf("case-1.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                }
                else {
                  // printf("case-1.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  grid[innerX][innerY][innerZ].collide (
                        this->getInternalStatistics() );
                  latticeTemplates<T,Descriptor>::swapAndStream3D (
                          grid, innerX, innerY, innerZ );
                }
              }
              // Case-2. On 2nd surface, innerx = domain.x0 + thread_block * n + 2
              else if (surface_id == 2) {
                // 2.1. On y=y0+1
                if (innerY == domain.y0+1){
                  // printf("case-2.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                    // second collide & revert on (innerX-1, innerY-1, innerZ)
                    grid[innerX-1][innerY-1][innerZ].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ].revert();
                  }
                  else {
                      grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);
                      
                      // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                      grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                      grid[innerX-1][innerY-1][innerZ-1].revert();
                    }
                  continue;
                }

                // 2.2. On y=y1
                if (innerY == domain.y1){
                  // first
                  // printf("case-2.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-2.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                    continue;
                  }
                  else{
                    // printf("case-2.2.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();

                    // second collide & revert on innerX-1, innerY, innerZ-2
                    grid[innerX-1][innerY][innerZ-2].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY][innerZ-2].revert();
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-2.2.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    // second collide & revert on innerX-1, innerY-1, innerZ
                    // second collide & revert on innerX-1, innerY, innerZ-1
                    // second collide & revert on innerX-1, innerY, innerZ
                    grid[innerX-1][innerY-1][innerZ].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ].revert();
                    grid[innerX-1][innerY][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY][innerZ-1].revert();
                    grid[innerX-1][innerY][innerZ].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY][innerZ].revert();
                  }

                  continue;
                }

                // 2.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  // printf("case-2.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  grid[innerX][innerY][innerZ].collide (
                        this->getInternalStatistics() );
                  latticeTemplates<T,Descriptor>::swapAndStream3D (
                        grid, innerX, innerY, innerZ );

                  // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  grid[innerX-1][innerY-1][innerZ-1].revert();
                  continue;
                }

                // 2.4. On z = z1
                if (innerZ == domain.z1){
                  // printf("case-2.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  grid[innerX-1][innerY-1][innerZ-1].revert();
                  // second collide & revert on (innerX-1, innerY-1, innerZ)
                  grid[innerX-1][innerY-1][innerZ].collide(this->getInternalStatistics());
                  grid[innerX-1][innerY-1][innerZ].revert();
                  continue;
                }

                // Other cases
                // printf("case-2.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                // first Collide & swapStream the cell.
                grid[innerX][innerY][innerZ].collide (
                      this->getInternalStatistics() );
                latticeTemplates<T,Descriptor>::swapAndStream3D (
                      grid, innerX, innerY, innerZ );
                

                // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                grid[innerX-1][innerY-1][innerZ-1].revert();
              }
              // Case-3 last surface (thread boundaries surface), case-3
              else if (surface_id == 0) {
                if (innerY == domain.y0 || innerZ == domain.z0) {
                  // printf("case-3.0 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  // first
                  boundSwapStream(domain, innerX, innerY, innerZ);
                  continue;
                }

                // 3.1 On y=y0+1
                if (innerY == domain.y0+1){
                  // printf("case-3.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerZ == domain.z1) {
                    // first
                    boundSwapStream(domain, innerX, innerY, innerZ);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }
                  else {
                    // first
                    if (innerX != domain.x1)
                      swapStream(innerX, innerY, innerZ);
                    else
                      boundSwapStream(domain, innerX, innerY, innerZ);
                    
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    }
                  continue;
                }

                // 3.2 On y=y1
                if (innerY == domain.y1){
                  // first
                  // printf("case-3.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  boundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-3.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    continue;
                  }
                  else{
                    // printf("case-3.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX-1, innerY-1, innerZ-1);

                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2);
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-3.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ);
                  }

                  continue;
                }

                // 3.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  // printf("case-3.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerX != domain.x1)
                    swapStream(innerX, innerY, innerZ);
                  else
                    boundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  continue;
                }

                // 3.4. On z = z1
                if (innerZ == domain.z1){
                  // printf("case-3.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  boundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D(
                      grid, innerX-1, innerY-1, innerZ-1);
                  
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  continue;
                }

                // Other cases
                // printf("case-3.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                // first Collide the cell.
                if (innerX != domain.x1)
                  swapStream(innerX, innerY, innerZ);
                else
                  boundSwapStream(domain, innerX, innerY, innerZ);

                // second collide
                grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                latticeTemplates<T,Descriptor>::swapAndStream3D(
                        grid, innerX-1, innerY-1, innerZ-1);
              }
              // Case 4
              else {
                // 4.1 On y=y0+1
                if (innerY == domain.y0+1){
                  // printf("case-4.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }
                  else {
                    // first
                    grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);

                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  }
                  continue;
                }

                // 4.2. On y=y1
                if (innerY == domain.y1){
                  // first
                  // printf("case-4.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-4.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    continue;
                  }
                  else{
                    // printf("case-4.2.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX-1, innerY-1, innerZ-1);

                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2);
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-4.2.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ);
                  }

                  continue;
                }

                // 4.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  // printf("case-4.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  grid[innerX][innerY][innerZ].collide (
                        this->getInternalStatistics() );
                  latticeTemplates<T,Descriptor>::swapAndStream3D (
                        grid, innerX, innerY, innerZ );
                  
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  continue;
                }

                // 2.4.4. On z = z1
                if (innerZ == domain.z1){
                  // printf("case-4.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second collide
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D(
                      grid, innerX-1, innerY-1, innerZ-1);
                  
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  continue;
                }

                // Other cases
                // printf("case-4.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
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
            } // end innerZ
          } // end innerY
        } // end innerX
      } // end outerZ
    } // end outerY
  } // end outerX
#endif // end step II
#if 1
  // ---------------3. compute 2nd on rest surface---------------------
  /* 2nd collide Stream on thread_boundaries, case 3.1 */
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0 + thread_block - 1; iX <= domain.x1; iX += thread_block){
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX, domain.y0, domain.y1, domain.z0, domain.z1) );
  }

  // 2nd swap Stream on 1st surface of every thread block, e.g. x0, x0 + thread_block, case-3.2
  // iX = 1, 5, ... when thread_block = 4
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0; iX < domain.x1; iX += thread_block){
    boundaryStream(domain, Box3D(iX, iX, domain.y0, domain.y1, domain.z0, domain.z1) );
  }
  /*-------------------------- Finish 3 -------------------------*/
#endif // end step III
}
#endif // end _OPENMP
    // global::profiler().stop("collStream");
}

#else
  // something
#endif
  #else // sequential
    #ifdef STEP2_UNROLL
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
    // printf("Here! I am step2_whole_seq_unroll_pyramid\n");

    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    // printf("blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n",
    //     blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);

    for (plint outerX = domain.x0; outerX <= domain.x1; outerX += blockSize) {
      // printf("outerX=%ld\n", outerX);
      for (plint outerY = domain.y0; outerY <= domain.y1+blockSize-1; outerY += blockSize) {
        // printf("outerY=%ld\n", outerY);
        for (plint outerZ = domain.z0; outerZ <= domain.z1+2*(blockSize-1); outerZ += blockSize) {
          // printf("outerZ=%ld\n", outerZ);
          // Inner loops.
          plint dx = 0;
          for (plint innerX = outerX; innerX <= std::min(outerX+blockSize-1, domain.x1);
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

            for (plint innerY=std::max(minY, domain.y0);
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

              for (plint innerZ=std::max(minZ, domain.z0);
                innerZ <= std::min(maxZ, domain.z1);
                ++innerZ)
              {
                // printf("inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);

                // 1. On x=x0, y=y0, z=z0
                if (innerX == domain.x0 || innerY == domain.y0 || innerZ == domain.z0){
                  // printf("case-1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  continue;
                }

                // 2. On y=y0+1
                if (innerY == domain.y0+1){
                  // printf("case-2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }
                  else {
                    // first
                    if (innerX == domain.x1){
                      collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                    }
                    else {
                      grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);
                    }

                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    }
                  continue;
                }

                // 3. On y=y1
                if (innerY == domain.y1){
                  // first
                  // printf("case-3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-3.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    continue;
                  }
                  else{
                    // printf("case-3.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    if (innerX != domain.x0 + 1){
                      grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                      latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX-1, innerY-1, innerZ-1);
                    } else {
                      collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    }

                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2); // don't forget second on line (iX-1, y1)!
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-3.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ);
                  }

                  continue;
                }

                // 4. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  // printf("case-4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  if (innerX == domain.x1){
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  }
                  else {
                    grid[innerX][innerY][innerZ].collide (
                          this->getInternalStatistics() );
                    latticeTemplates<T,Descriptor>::swapAndStream3D (
                          grid, innerX, innerY, innerZ );
                  }

                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  continue;
                }

                // 5. On z = z1
                if (innerZ == domain.z1){
                  // printf("case-5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  // second collide
                  if (innerX != domain.x0 + 1){
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(
                        grid, innerX-1, innerY-1, innerZ-1);
                  } else {
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  }
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  continue;
                }

                // Other cases
                // printf("case-6 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                // first Collide the cell.
                if (innerX == domain.x1){
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                }
                else {
                  grid[innerX][innerY][innerZ].collide (
                          this->getInternalStatistics() );
                  // Swap the populations on the cell, and then with post-collision
                  //   neighboring cell, to perform the streaming step.
                  latticeTemplates<T,Descriptor>::swapAndStream3D (
                          grid, innerX, innerY, innerZ );
                }

                // second collide
                if (innerX != domain.x0 + 1){
                  grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D(
                          grid, innerX-1, innerY-1, innerZ-1);
                }
                else{
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                }
              }
            }
          }
        }
      }
    }

    // second collision on top
    collideRevertAndBoundSwapStream(domain, Box3D(domain.x1, domain.x1,
                domain.y0, domain.y1, domain.z0, domain.z1) );

    global::profiler().stop("collStream");
}
    #else
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
    // printf("Here! I am step2_whole_seq_pyramid\n");

    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

    // Make sure domain is contained within current lattice
    //PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    // printf("blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n",
    //     blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);

    for (plint outerX = domain.x0; outerX <= domain.x1; outerX += blockSize) {
      // printf("outerX=%ld\n", outerX);
      for (plint outerY = domain.y0; outerY <= domain.y1+blockSize-1; outerY += blockSize) {
        // printf("outerY=%ld\n", outerY);
        for (plint outerZ = domain.z0; outerZ <= domain.z1+2*(blockSize-1); outerZ += blockSize) {
          // printf("outerZ=%ld\n", outerZ);
          // Inner loops.
          plint dx = 0;
          for (plint innerX = outerX; innerX <= std::min(outerX+blockSize-1, domain.x1);
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

            for (plint innerY=std::max(minY, domain.y0);
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

              for (plint innerZ = std::max(minZ, domain.z0);
                innerZ <= std::min(maxZ, domain.z1);
                ++innerZ)
              {
                // printf("inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);

                // first on y=y0, z=z0
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                // second collide & stream
                plint nextX = innerX - 1;
                plint nextY = innerY - 1;
                plint nextZ = innerZ - 1;
                if ( nextX >= domain.x0 &&
                     nextY >= domain.y0 &&
                     nextZ >= domain.z0 ) {
                  collideRevertAndBoundSwapStream(domain, nextX, nextY, nextZ);

                  if (innerZ == domain.z1){
                    collideRevertAndBoundSwapStream(domain, nextX, nextY, innerZ);
                  }

                  // innerY = y1
                  if (innerY == domain.y1 && (nextZ - 1) >= domain.z0){
                    collideRevertAndBoundSwapStream(domain, nextX, innerY, nextZ - 1);

                    // z=z1, one more things to do
                    if (innerZ == domain.z1){
                      collideRevertAndBoundSwapStream(domain, nextX, innerY, nextZ);
                      collideRevertAndBoundSwapStream(domain, nextX, innerY, innerZ);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // last surface domain.x1
    collideRevertAndBoundSwapStream(domain, Box3D(domain.x1, domain.x1,
                domain.y0, domain.y1,
                domain.z0, domain.z1) );

    global::profiler().stop("collStream");
}
    #endif
  #endif

#elif STEP2_3PARTS  // THREE_PARTS_STEP2_Implementation
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

    // step i: 2 collideAndStream on x0; 1 collideAndStream on x0+1
    step2CollideAndStream_init(domain);

    // step ii: Then bulk [x0+2, x0-1]
  #ifdef STEP2_OMP
      #ifdef STEP2_PYRAMID
      // something
      #else
      step2CollideAndStream_bulk_omp(domain);
      #endif
  #else 
      #ifdef STEP2_PYRAMID
    // if (domain.x1 - domain.x0 <= 250) {
    //   step2CollideAndStream_bulk(domain);
    // }
    // else {
      step2CollideAndStream_bulk_blockwise(domain);
    // }
      #else
        step2CollideAndStream_bulk(domain);
      #endif
  #endif

    // step iii : computing the rest
    step2CollideAndStream_end(domain);

    global::profiler().stop("collStream");
}
#else
// default
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) { }
#endif



// Add by Yuankun
plint ykPanel_len;
#if defined(STEP2_OMP)

#if defined(PANEL_MEM)
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

   // Now: use panel memory layout to compute, access grid[iX][iY + ykPanel_len * (iZ / ykPanel_len)][iZ % ykPanel_len]
    plint NewNy = ny * (nz / ykPanel_len);
    plint NewNz = ykPanel_len;
    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX = 0; iX < nx; ++iX) {
      grid[iX] = new Cell<T,Descriptor>* [NewNy];
      for (plint iY = 0; iY < NewNy; ++iY) {
        grid[iX][iY] = rawData + NewNz*(iY+NewNy*iX);
      }
    }
}
#else
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

    // Before: access grid[iX][iY][iZ]
    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX=0; iX<nx; ++iX) {
      grid[iX] = new Cell<T,Descriptor>* [ny];
      for (plint iY=0; iY<ny; ++iY) {
        grid[iX][iY] = rawData + nz*(iY+ny*iX);
      }
    }
    
}
#endif

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::releaseMemory() {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();

    // std::cout << "BlockLattice3D::releaseMemory()" << std::endl;

    #pragma omp parallel for default(shared) schedule(static)
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
#endif

#if 0
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
#endif


}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_HH
