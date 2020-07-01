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
#ifndef BLOCK_LATTICE_3D_STEP2_OMP_HH
#define BLOCK_LATTICE_3D_STEP2_OMP_HH

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

#if 1
/* step2CollideAndStream_bulk_omp: used in step2_3parts_omp_line*/
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

#ifdef STEP2_WHOLE
//   /** This operation is 2 step collision and stream using 1 buffer
//  * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
//  * because memory is traversed only once instead of twice.
//  */
  #ifdef STEP2_OMP
    #ifdef STEP2_UNROLL
      #ifdef PILLAR_MEM

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
  // global::profiler().start("collStream");
  // global::profiler().increment("collStreamCells", domain.nCells());

  // Make sure domain is contained within current lattice
  PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

  // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
  //   the blocks, whereas the three inner loops enumerate the cells inside each block.
  const plint blockSize = cachePolicy().getBlockSize();

#ifdef DEBUG_PRINT
  printf("Here! I am step2_whole_omp_unroll_pyramid_panel_mem, blockSize=%ld, domain-(%ld, %ld) (%ld, %ld) (%ld, %ld)\n",
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
        plint iX_t = pillar_map_iX(iX, iY, iZ);
        plint iY_t = pillar_map(iY);
        plint iZ_t = pillar_map(iZ);

        grid[iX_t][iY_t][iZ_t].collide(this->getInternalStatistics());
        grid[iX_t][iY_t][iZ_t].revert();
      }
    }
  }
#endif
  #ifdef DEBUG_PRINT
  pcout << "End step I\n";
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
              plint surface_id = innerX % thread_block;

              plint innerX_t = pillar_map_iX(innerX, innerY, innerZ);
              plint innerY_t = pillar_map(innerY);
              plint innerZ_t = pillar_map(innerZ);

              // Case-0. On x=x0, y=y0, z=z0, except last surface within a thread block
              if (surface_id != 0 && (innerX == domain.x0 || innerY == domain.y0 || innerZ == domain.z0)){
                #ifdef DEBUG_PRINT
                printf("case-0 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                #endif
                collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                continue;
              }
              
              // Case-1. On 1st surface of a thread block, innerX = domain.x0 + thread_block * n + 1
              if (surface_id == 1) {
                // 2.1 On x=x0, Or boundary y=y0,y1, z=z0, z1
                if (innerY == domain.y1 || innerZ == domain.z1){
                  #ifdef DEBUG_PRINT
                  printf("case-1.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                }
                else {
                  // printf("case-1.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #ifdef CUBE_MAP
                  plint innerX_minus_1_t = pillar_map_iX(iX - 1, iY, iZ);
                  #else
                  plint innerX_minus_1_t = innerX_t - 1;
                  #endif

                  grid[innerX_t][innerY_t][innerZ_t].collide (this->getInternalStatistics());
                  latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);
                }
              }
              // Case-2. On 2nd surface, innerX = domain.x0 + thread_block * n + 2
              else if (surface_id == 2) {
                #ifdef CUBE_MAP
                plint innerX_minus_1_t = pillar_map_iX(innerX - 1, innerY, innerZ);
                #else
                plint innerX_minus_1_t = innerX_t - 1;
                #endif

                plint innerY_minus_1_t = pillar_map(innerY - 1);
                plint innerZ_minus_1_t = pillar_map(innerZ - 1);

                // 2.1. On y=y0+1
                if (innerY == domain.y0+1){
                  #ifdef DEBUG_PRINT
                  printf("case-2.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif

                  if (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
                    // second collide & revert on (innerX-1, innerY-1, innerZ)
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].revert();
                  }
                  else {
                      grid[innerX_t][innerY_t][innerZ_t].collide(this->getInternalStatistics());
                      latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);
                      
                      // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                      grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                      grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
                  }
                  continue;
                }

                // 2.2. On y=y1
                if (innerY == domain.y1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-2.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-2.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
                    continue;
                  }
                  else{
                    // printf("case-2.2.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();

                    // second collide & revert on innerX-1, innerY, innerZ-2
                    plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                    grid[innerX_minus_1_t][innerY_t][innerZ_minus_2_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_t][innerZ_minus_2_t].revert();
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-2.2.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    // second collide & revert on innerX-1, innerY-1, innerZ
                    // second collide & revert on innerX-1, innerY, innerZ-1
                    // second collide & revert on innerX-1, innerY, innerZ
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].revert();
                    grid[innerX_minus_1_t][innerY_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_t][innerZ_minus_1_t].revert();
                    grid[innerX_minus_1_t][innerY_t][innerZ_t].collide(this->getInternalStatistics());
                    grid[innerX_minus_1_t][innerY_t][innerZ_t].revert();
                  }

                  continue;
                }

                // 2.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-2.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  grid[innerX_t][innerY_t][innerZ_t].collide (
                        this->getInternalStatistics() );
                  latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);

                  // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
                  continue;
                }

                // 2.4. On z = z1
                if (innerZ == domain.z1){
                  #ifdef DEBUG_PRINT
                  printf("case-2.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
                  // second collide & revert on (innerX-1, innerY-1, innerZ)
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].collide(this->getInternalStatistics());
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_t].revert();
                  continue;
                }

                // Other cases
                #ifdef DEBUG_PRINT
                printf("case-2.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                #endif
                // first Collide & swapStream the cell.
                grid[innerX_t][innerY_t][innerZ_t].collide (
                      this->getInternalStatistics() );
                latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);

                // second collide & revert on (innerX-1, innerY-1, innerZ-1)
                grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].revert();
              }
              // Case-3 last surface (thread boundaries surface), case-3
              else if (surface_id == 0) {
                if (innerY == domain.y0 || innerZ == domain.z0) {
                  // printf("case-3.0 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  // first
                  boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                  continue;
                }

                #ifdef CUBE_MAP
                plint innerX_minus_1_t = pillar_map_iX(innerX - 1, innerY, innerZ);
                #else
                plint innerX_minus_1_t = innerX_t - 1;
                #endif

                plint innerY_minus_1_t = pillar_map(innerY - 1);
                plint innerZ_minus_1_t = pillar_map(innerZ - 1);

                // 3.1 On y=y0+1
                if (innerY == domain.y0+1){
                  #ifdef DEBUG_PRINT
                  printf("case-3.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif

                  if (innerZ == domain.z1) {
                    // first
                    boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, 
                                                    innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ,
                                                    innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                  }
                  else {
                    // first
                    if (innerX != domain.x1)
                      swapStream(innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                    else
                      boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                    
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1,
                                                    innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                  }
                  continue;
                }

                // 3.2 On y=y1
                if (innerY == domain.y1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-3.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-3.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                    continue;
                  }
                  else{
                    plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                    // printf("case-3.2.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, 
                      innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2, innerX_minus_1_t, innerY_t, innerZ_minus_2_t);
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-3.2.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ, innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1, innerX_minus_1_t, innerY_t, innerZ_minus_1_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ, innerX_minus_1_t, innerY_t, innerZ_t);
                  }

                  continue;
                }

                // 3.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-3.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif

                  if (innerX != domain.x1)
                    swapStream(innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                  else
                    boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                  continue;
                }

                // 3.4. On z = z1
                if (innerZ == domain.z1){
                  #ifdef DEBUG_PRINT
                  printf("case-3.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                  plint innerZ_minus_2_t = pillar_map(innerZ - 2); 
                  latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);
                  
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ, innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                  continue;
                }

                // Other cases
                #ifdef DEBUG_PRINT
                printf("case-3.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                #endif
                // first Collide the cell.
                if (innerX != domain.x1)
                  swapStream(innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                else
                  boundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                // second collide
                grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);
              }
              // Case 4
              else {
                #ifdef CUBE_MAP
                plint innerX_minus_1_t = pillar_map_iX(innerX - 1, innerY, innerZ);
                #else
                plint innerX_minus_1_t = innerX_t - 1;
                #endif

                plint innerY_minus_1_t = pillar_map(innerY - 1);
                plint innerZ_minus_1_t = pillar_map(innerZ - 1);

                // 4.1 On y=y0+1
                if (innerY == domain.y0+1){
                  #ifdef DEBUG_PRINT
                  printf("case-4.1 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif

                  if (innerZ == domain.z1) {
                    // first
                    collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);
                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ, innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                  }
                  else {
                    // first
                    grid[innerX_t][innerY_t][innerZ_t].collide(this->getInternalStatistics());
                    latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);

                    // second
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                  }
                  continue;
                }

                // 4.2. On y=y1
                if (innerY == domain.y1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-4.2 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second
                  if (innerZ == domain.z0+1){
                    // printf("case-4.2.1 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                    continue;
                  }
                  else{
                    // printf("case-4.2.2 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                    plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                    latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);

                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-2, innerX_minus_1_t, innerY_t, innerZ_minus_2_t);
                  }

                  // z=z1, more things to do for second
                  if (innerZ == domain.z1){
                    // printf("case-4.2.3 second inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ, innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ-1, innerX_minus_1_t, innerY_t, innerZ_minus_1_t);
                    collideRevertAndBoundSwapStream(domain, innerX-1, innerY, innerZ, innerX_minus_1_t, innerY_t, innerZ_t);
                  }

                  continue;
                }

                // 4.3. On z = z0+1
                if (innerZ == domain.z0+1){
                  // first
                  #ifdef DEBUG_PRINT
                  printf("case-4.3 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  #endif
                  grid[innerX_t][innerY_t][innerZ_t].collide (
                        this->getInternalStatistics() );
                  latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);
                  
                  // second
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t);
                  continue;
                }

                // 4.4. On z = z1
                if (innerZ == domain.z1){
                  // printf("case-4.4 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                  collideRevertAndBoundSwapStream(domain, innerX, innerY, innerZ, innerX_t, innerY_t, innerZ_t);

                  // second collide
                  grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                  plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                  latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);
                  
                  collideRevertAndBoundSwapStream(domain, innerX-1, innerY-1, innerZ, innerX_minus_1_t, innerY_minus_1_t, innerZ_t);
                  continue;
                }

                // Other cases
                #ifdef DEBUG_PRINT
                printf("case-4.5 inner(%ld, %ld, %ld)\n", innerX, innerY, innerZ);
                #endif
                // first Collide the cell.
                
                grid[innerX_t][innerY_t][innerZ_t].collide (
                        this->getInternalStatistics() );
                // Swap the populations on the cell, and then with post-collision
                //   neighboring cell, to perform the streaming step.
                latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX, innerY, innerZ, 
                                                  innerX_t, innerY_t, innerZ_t, innerX_minus_1_t);

                // second collide
                grid[innerX_minus_1_t][innerY_minus_1_t][innerZ_minus_1_t].collide(this->getInternalStatistics());
                plint innerZ_minus_2_t = pillar_map(innerZ - 2);
                latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_map(grid, innerX-1, innerY-1, innerZ-1, innerX_minus_1_t, innerY_minus_1_t, innerZ_minus_1_t, innerZ_minus_2_t);
              }
            } // end innerZ
          } // end innerY
        } // end innerX
      } // end outerZ
    } // end outerY
  } // end outerX
#endif // end step II
    // pcout << "End step II\n";
#if 1
  // ---------------3. compute 2nd on rest surface---------------------
  /* 2nd collide Stream on thread_boundaries, case 3.1 */
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0 + thread_block - 1; iX <= domain.x1; iX += thread_block){
    collideRevertAndBoundSwapStream(domain, Box3D(iX, iX, domain.y0, domain.y1, domain.z0, domain.z1) );
  }

  // pcout << "End collideRevertAndBoundSwapStream surface iX\n";

  // 2nd swap Stream on 1st surface of every thread block, e.g. x0, x0 + thread_block, case-3.2
  // iX = 1, 5, ... when thread_block = 4
  #pragma omp for private(iX, iY, iZ) schedule(static)
  for (iX = domain.x0; iX < domain.x1; iX += thread_block){
    boundaryStream(domain, Box3D(iX, iX, domain.y0, domain.y1, domain.z0, domain.z1) );
  }

  // pcout << "End boundaryStream surface iX\n";
  /*-------------------------- Finish 3 -------------------------*/
#endif // end step III
}
#endif // end _OPENMP
    // global::profiler().stop("collStream");
}
      /******************************End of STEP2_WHOLE + STEP2_UNROLL + PILLAR_MEM*******************************************/
      #else // use STEP2_UNROLL but NOT using PILLAR_MEM
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
      #endif // end PILLAR_MEM

    #else // NOT using STEP2_UNROLL
    // something
    #endif // end STEP2_UNROLL

  #endif // end STEP2_WHOLE

#elif STEP2_3PARTS  // THREE_PARTS_STEP2_Implementation
  #ifdef STEP2_OMP
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2CollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", domain.nCells());

    // step i: 2 collideAndStream on x0; 1 collideAndStream on x0+1
    step2CollideAndStream_init(domain);

    // step ii: Then bulk [x0+2, x0-1]
    #ifdef STEP2_PYRAMID
    // something
    #else
    step2CollideAndStream_bulk_omp(domain);
    #endif
  
    // step iii : computing the rest
    step2CollideAndStream_end(domain);

    global::profiler().stop("collStream");
}
  #endif
#endif

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_OMP_HH
