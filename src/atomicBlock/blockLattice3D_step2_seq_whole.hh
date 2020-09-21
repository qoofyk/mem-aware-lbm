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
#ifndef BLOCK_LATTICE_3D_STEP2_SEQ_WHOLE_HH
#define BLOCK_LATTICE_3D_STEP2_SEQ_WHOLE_HH

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
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_seq_whole_blockwise_unroll(Box3D domain) {
    // printf("Here! I am step2_whole_seq_unroll_pyramid\n");
#if 0
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    // global::profiler().start("collStream");
    // global::profiler().increment("collStreamCells", domain.nCells());
#endif

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    #if 0
    int tid = omp_get_thread_num();
    printf("Tid%d: blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n",
        tid, blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);
    #endif

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

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_seq_whole_blockwise(Box3D domain) {
    // printf("Here! I am step2_whole_seq_pyramid\n");

#if 0
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());
#endif

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    #if 0
    int tid = omp_get_thread_num();
    printf("Tid%d: blockSize=%ld, domain(%ld, %ld, %ld, %ld, %ld, %ld)\n",
        tid, blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);
    #endif

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

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_SEQ_WHOLE_HH
