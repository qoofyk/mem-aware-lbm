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
#ifndef BLOCK_LATTICE_3D_STEP2_OMP_WHOLE_HH
#define BLOCK_LATTICE_3D_STEP2_OMP_WHOLE_HH

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

#include "atomicBlock/blockLattice3D_step2_util.hh"

#ifdef _OPENMP
  #include <omp.h>
#endif

// #define DEBUG_PRINT

namespace plb {

#if 0
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_omp_whole_blockwise_unroll(Box3D domain) {

#if 0  
  global::profiler().start("collStream");
  global::profiler().increment("collStreamCells", domain.nCells());

  // Make sure domain is contained within current lattice
  PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );
#endif

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();

#if 0
  printf("Here! I am step2CollideAndStream_omp_whole_blockwise_unroll, blockSize=%ld, domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
    blockSize, domain.x0, domain.x1, domain.y0, domain.y1, domain.z0, domain.z1);
#endif

  
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  plint iX, iY, iZ;
  plint tid = omp_get_thread_num();
  plint nthreads = omp_get_num_threads();
  plint my_lx[2];
  my_lx[0] = domain.x0 + tid * thread_block;
  my_lx[1] = (tid + 1) * thread_block;

  // printf("tid: %ld, thread_block=%ld\n", tid, thread_block);
#if 1
  // ---------------1. compute thread boundaries surface ---------------------
  // 1st Collide_Revert
  iX = my_lx[1];
  for (iY = domain.y0; iY <= domain.y1; ++iY) {
    for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
      grid[iX][iY][iZ].collide(this->getInternalStatistics());
      grid[iX][iY][iZ].revert();
    }
  }
  #pragma omp barrier
#endif

#if 1
  // ---------------2. compute bulk ---------------------
  for (plint outerX = my_lx[0]; outerX <= my_lx[1]; outerX += blockSize) {  
    // printf("tid: %ld, outerX=%ld\n", tid, outerX);
    for (plint outerY = domain.y0; outerY <= domain.y1+blockSize-1; outerY += blockSize) {
      // printf("outerY=%ld\n", outerY);
      for (plint outerZ = domain.z0; outerZ <= domain.z1+2*(blockSize-1); outerZ += blockSize) {
        // printf("outerZ=%ld\n", outerZ);
        // Inner loops.
        plint dx = 0;
        plint innerX_max = std::min(outerX+blockSize-1, my_lx[1]);
        for (plint innerX = outerX; innerX <= innerX_max; ++innerX, ++dx)
        {
          // Y-index is shifted in negative direction at each x-increment. to ensure
          //   that only post-collision cells are accessed during the swap-operation
          //   of the streaming.
          plint minY = outerY-dx;
          plint maxY = minY+blockSize-1;
          plint dy = 0;
          plint innerY_start = std::max(minY, domain.y0);
          plint innerY_end = std::min(maxY, domain.y1);

          // surface_id == 1: First collide & stream
          if (innerX == my_lx[0]) {
            if (innerX == domain.x0) {
              for (plint innerY = innerY_start; innerY <= innerY_end; ++innerY, ++dy) {
                // Z-index is shifted in negative direction at each x-increment. and at each
                //    y-increment, to ensure that only post-collision cells are accessed during
                //    the swap-operation of the streaming.
                plint minZ = outerZ-dx-dy;
                plint maxZ = minZ+blockSize-1;

                plint innerZ_start = std::max(minZ, domain.z0);
                plint innerZ_end = std::min(maxZ, domain.z1);

                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                }
              }
            }
            else {
              for (plint innerY = innerY_start; innerY <= innerY_end; ++innerY, ++dy) {
                // Z-index is shifted in negative direction at each x-increment. and at each
                //    y-increment, to ensure that only post-collision cells are accessed during
                //    the swap-operation of the streaming.
                plint minZ = outerZ-dx-dy;
                plint maxZ = minZ+blockSize-1;

                plint innerZ_start = std::max(minZ, domain.z0);
                plint innerZ_end = std::min(maxZ, domain.z1);

                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  if (innerY == domain.y0 || innerY == domain.y1 || innerZ == domain.z0 || innerZ == domain.z1) {
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
                }
              }
            }
          }

          // surface_id == 2: First collide & $stream on surface_id=2; Second $collide & revert$ on surface_id=1
          else if (innerX == my_lx[0] + 1) {
            for (plint innerY = innerY_start; innerY <= innerY_end; ++innerY, ++dy) {
            
              plint minZ = outerZ-dx-dy;
              plint maxZ = minZ+blockSize-1;

              plint innerZ_start = std::max(minZ, domain.z0);
              plint innerZ_end = std::min(maxZ, domain.z1);

              // i. the first row of a surface
              if (innerY == domain.y0) {
                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                }
              }

              // ii. the second row of a surface
              else if (innerY == domain.y0 + 1) {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {

                  // 1. the first node in a row
                  if (innerZ == domain.z0) {
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  // 2. the last node in a row
                  else if  (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                    // second collide & revert on (innerX-1, innerY-1, innerZ)
                    grid[innerX-1][innerY-1][innerZ].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ].revert();
                  }
                  // 3. other nodes in a row
                  else {
                    grid[innerX][innerY][innerZ].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, innerX, innerY, innerZ);
                    
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                  }
                }
              }

              // iii. the last row of a surface
              else if (innerY == domain.y1) {
                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {

                  // first comp
                  // printf("case-2.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);

                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                  if (innerZ == domain.z0) continue;

                  // second comp
                  if (innerZ == domain.z0+1){
                    // printf("case-2.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
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
                }
              } // end iii

              // iv. other rows of a surface
              else {
                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {

                  if (innerZ == domain.z0) {
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  else if (innerZ == domain.z0+1) {
                    // first
                    // printf("case-2.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX][innerY][innerZ].collide (
                          this->getInternalStatistics() );
                    latticeTemplates<T,Descriptor>::swapAndStream3D (
                          grid, innerX, innerY, innerZ );

                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                  }

                  // 2.4. On z = z1
                  else if (innerZ == domain.z1) {
                    // printf("case-2.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ-1].revert();
                    // second collide & revert on (innerX-1, innerY-1, innerZ)
                    grid[innerX-1][innerY-1][innerZ].collide(this->getInternalStatistics());
                    grid[innerX-1][innerY-1][innerZ].revert();
                  }

                  // Other cases
                  else {
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
                }
              } // end iv
            } // end innverY
          } // end surface_id=2

          // surface_id == 0: First stream on surface_id=0; Second collide & stream on the lower surface
          else if (innerX == my_lx[1]) {

            for (plint innerY = innerY_start; innerY <= innerY_end; ++innerY, ++dy) {
          
              plint minZ = outerZ-dx-dy;
              plint maxZ = minZ+blockSize-1;

              plint innerZ_start = std::max(minZ, domain.z0);
              plint innerZ_end = std::min(maxZ, domain.z1);

              // i. the first row of a surface
              if (innerY == domain.y0) {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  boundSwapStream(domain, innerX, innerY, innerZ);
                }
              }

              // ii. the second row of a surface
              else if (innerY == domain.y0 + 1) {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  // 1. the first node in a row
                  if (innerZ == domain.z0) {
                    boundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  // 2. the last node in a row
                  else if (innerZ == domain.z1) {
                    // first
                    boundSwapStream(domain, innerX, innerY, innerZ);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }

                  // other nodes
                  else {
                    // first
                    if (innerX != domain.x1)
                      swapStream(innerX, innerY, innerZ);
                    else
                      boundSwapStream(domain, innerX, innerY, innerZ);
                    
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  }
                }
              }

              // iii. the last row of a surface
              else if (innerY == domain.y1) {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  // first
                  // printf("case-3.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  boundSwapStream(domain, innerX, innerY, innerZ);
                  if (innerZ == domain.z0) continue;

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-3.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);                  }
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
                }
              }

              // iv. other rows of a surface
              else {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  if (innerZ == domain.z0) {
                    boundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  else if (innerZ == domain.z0+1){
                    // first
                    // printf("case-3.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    if (innerX != domain.x1)
                      swapStream(innerX, innerY, innerZ);
                    else
                      boundSwapStream(domain, innerX, innerY, innerZ);

                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  }

                  else if (innerZ == domain.z1){
                    // printf("case-3.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    boundSwapStream(domain, innerX, innerY, innerZ);

                    // second
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(
                        grid, innerX-1, innerY-1, innerZ-1);
                    
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }

                  // Other cases
                  else {
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
                }
              }
            } // end innverY
          } // end surface_id=0

          // on other surfaces within each thread's sub-3D-lattice, First $collide$ \& $stream$ on (innerX, innerY, innerZ); Second $collide$ \& $stream$  on (innerX-1, innerY-1, innerZ-1)
          else {
            for (plint innerY = innerY_start; innerY <= innerY_end; ++innerY, ++dy) {
            
              plint minZ = outerZ-dx-dy;
              plint maxZ = minZ+blockSize-1;

              plint innerZ_start = std::max(minZ, domain.z0);
              plint innerZ_end = std::min(maxZ, domain.z1);

               // i. the first row of a surface
              if (innerY == domain.y0) {
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                }
              }

              // ii. the second row of a surface On y=y0+1
              else if (innerY == domain.y0+1){
                for (plint innerZ=innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  if (innerZ == domain.z0) {
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  // printf("case-4.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  else if (innerZ == domain.z1) {
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
                }
              }

              // iii. the last row of a surface
              else if (innerY == domain.y1){
                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  // first
                  // printf("case-4.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  if (innerZ == domain.z0) continue;

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-4.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
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
                }
              }

              // iv. other rows of a surface
              else {
                for (plint innerZ = innerZ_start; innerZ <= innerZ_end; ++innerZ) {
                  if (innerZ == domain.z0) {
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                  }

                  // On z = z0+1
                  else if (innerZ == domain.z0+1){
                    // first
                    // printf("case-4.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX][innerY][innerZ].collide (
                          this->getInternalStatistics() );
                    latticeTemplates<T,Descriptor>::swapAndStream3D (
                          grid, innerX, innerY, innerZ );
                    
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1);
                  }

                  // On z = z1
                  else if (innerZ == domain.z1){
                    // printf("case-4.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);

                    // second collide
                    grid[innerX-1][innerY-1][innerZ-1].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D(
                        grid, innerX-1, innerY-1, innerZ-1);
                    
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ);
                  }

                  // Other cases
                  else {
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
                }
              }
            }
          } // end other surface_id
        } // innerX
      } // outerZ
    } // outerY
  } // outerX

  #pragma omp barrier    
#endif // end step II

#if 1
  // ---------------3. compute 2nd on rest surface---------------------
#if 1
  // intersection
  iX = my_lx[1];
  if (iX == domain.x1) {
    for (iY = domain.y0; iY <= domain.y1; ++iY) {
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
      }
    }
  }
  else {
    for (iY = domain.y0; iY <= domain.y1; ++iY) {
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        if (iY == domain.y0 || iY == domain.y1 || iZ == domain.z0 || iZ == domain.z1) {
          collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
        }
        else {
          grid[iX][iY][iZ].collide (
              this->getInternalStatistics() );
          // Swap the populations on the cell, and then with post-collision
          //   neighboring cell, to perform the streaming step.
          latticeTemplates<T,Descriptor>::swapAndStream3D (
                  grid, iX, iY, iZ );
        }
      }
    }
  }
  

  #pragma omp barrier

  iX = my_lx[0];
  if (iX == domain.x0) {
    for (iY = domain.y0; iY <= domain.y1; ++iY) {
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        boundSwapStream(domain, iX, iY, iZ);
      }
    }
  }
  else {
    for (iY = domain.y0; iY <= domain.y1; ++iY) {
      for (iZ = domain.z0; iZ <= domain.z1; ++iZ ){
        if (iY == domain.y0 || iY == domain.y1 || iZ == domain.z0 || iZ == domain.z1) {
          boundSwapStream(domain, iX, iY, iZ);
        }
        else {
          swapStream(iX, iY, iZ);
        }
      }
    }
  }
#endif

#if 0
  /* 2nd collide Stream on thread_boundaries, case 3.1 */ // can be optimized
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
#endif

  /*-------------------------- Finish 3 -------------------------*/
#endif // end step III
}
#endif // end _OPENMP
    // global::profiler().stop("collStream");
}
/******************************End of STEP2_OMP + STEP2_WHOLE + STEP2_UNROLL*******************************************/
#endif

#if 1
// use STEP2_UNROLL but NOT using PILLAR_MEM
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream_omp_whole_blockwise_unroll(Box3D domain) {

#if 0  
  global::profiler().start("collStream");
  global::profiler().increment("collStreamCells", domain.nCells());

  // Make sure domain is contained within current lattice
  PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );
#endif

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();

#if 0
  printf("Here! I am step2CollideAndStream_omp_whole_blockwise_unroll, blockSize=%ld, domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
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
          plint surface_id = innerX % thread_block;
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

              // Case-0. On x=x0, y=y0, z=z0, except last surface within a thread block
              if (surface_id != 0 && (innerX == domain.x0 || innerY == domain.y0 || innerZ == domain.z0)){
                // printf("case-0 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ);
                continue;
              }
              
              // Case-1. On 1st surface of a thread block, innerX = domain.x0 + thread_block * n + 1
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
              // Case-2. On 2nd surface, innerX = domain.x0 + thread_block * n + 2
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
  /* 2nd collide Stream on thread_boundaries, case 3.1 */ // can be optimized
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
/******************************End of STEP2_OMP + STEP2_WHOLE + STEP2_UNROLL*******************************************/
#endif

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_OMP_WHOLE_HH