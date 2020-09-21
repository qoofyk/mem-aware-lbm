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
#ifndef BLOCK_LATTICE_3D_STEP2_SEQ_HH
#define BLOCK_LATTICE_3D_STEP2_SEQ_HH

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
#ifdef DEBUG_PRINT  
  printf("Here! I am step2_3parts_seq_unroll_line\n");
#endif

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

#ifndef STEP2_OMP // sequential
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

  #if defined(STEP2_3PARTS) // THREE_PARTS_STEP2_Implementation
    // step i: 2 collideAndStream on x0; 1 collideAndStream on x0+1
    step2CollideAndStream_init(domain);

    // step ii: Then bulk [x0+2, x0-1]
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

    // step iii : computing the rest
    step2CollideAndStream_end(domain);

  #elif defined(STEP2_WHOLE)
    #ifdef STEP2_UNROLL
      step2CollideAndStream_seq_whole_blockwise_unroll(domain);
    #else
      step2CollideAndStream_seq_whole_blockwise(domain);
    #endif

  #endif

    global::profiler().stop("collStream");
}
#endif

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_SEQ_HH
