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
#ifndef BLOCK_LATTICE_3D_STEP2_UTIL_HH
#define BLOCK_LATTICE_3D_STEP2_UTIL_HH

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

#include "atomicBlock/blockLattice3D_pillar_mem.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

// #define DEBUG_PRINT

namespace plb {
// BlockLattice3D_step2_util /////////////////////////

// Add by Yuankun
plint thread_block;
// End add by Yuankun

inline plint pillar_map_iX(plint iX, plint iY, plint iZ) {
#ifdef CUBE_MAP
  #ifdef BIT_HACK
    return (iX & (YK_TILE - 1)) + ((iZ >> LOG_YK_TILE) << LOG_YK_TILE) + ((iY >> LOG_YK_TILE) << (log_NzTiles + LOG_YK_TILE)) + ( (iX >> LOG_YK_TILE) << (log_NzTiles + log_NyTiles + LOG_YK_TILE) );
  #else
    // return iX % ykTile + ykTile * (iZ / ykTile + (iY / ykTile) * (nz_ / ykTile) + (iX / ykTile)  * (nz_ / ykTile)  * (ny_ / ykTile));
    return iX % ykTile + ykTile * (iZ / ykTile + (iY / ykTile) * NzTiles + (iX / ykTile)  * NzTiles  * NyTiles);
  #endif
#else // PILLAR_MAP
  #ifdef BIT_HACK
    return iX + memNx * ( (iZ >> LOG_YK_TILE) + ( (iY >> LOG_YK_TILE) << log_NzTiles) );
  #else   
    return iX + memNx * (iZ / ykTile + (iY / ykTile) * NzTiles);
  #endif
#endif
}

inline plint pillar_map(plint iZ) {
  #ifdef BIT_HACK 
    return iZ & (YK_TILE - 1);
  #else
    return iZ % ykTile;
  #endif
}

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

    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX=0; iX<nx; ++iX) {
      for (plint iY=0; iY<ny; ++iY) {
        for (plint iZ=0; iZ<nz; ++iZ) {
          #ifdef PILLAR_MEM
            plint iX_t = pillar_map_iX(iX, iY, iZ);
            plint iY_t = pillar_map(iY);
            plint iZ_t = pillar_map(iZ);

            #ifdef DEBUG_PRINT
            printf("(%ld, %ld, %ld) -> (%ld, %ld, %ld)\n", iX, iY, iZ, iX_t, iY_t, iZ_t);
            #endif
            
            grid[iX_t][iY_t][iZ_t].attributeDynamics(backgroundDynamics);
            
          #else
            grid[iX][iY][iZ].attributeDynamics(backgroundDynamics);
          #endif
        }
      }
    }

    #ifdef DEBUG_PRINT
    pcout << "Pass assign attributeDynamics\n";
    #endif

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
                
              #ifdef PILLAR_MEM
                plint iX_t = pillar_map_iX(iX, iY, iZ);
                plint iY_t = pillar_map(iY);
                plint iZ_t = pillar_map(iZ);

                Cell<T,Descriptor>& cell = grid[iX_t][iY_t][iZ_t];
                // Assign cell from rhs
                cell = rhs.grid[iX_t][iY_t][iZ_t];
                
              #else
                Cell<T,Descriptor>& cell = grid[iX][iY][iZ];
                // Assign cell from rhs
                cell = rhs.grid[iX][iY][iZ];
              #endif
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
              #ifdef PILLAR_MEM
                plint iX_t = pillar_map_iX(iX, iY, iZ);
                plint iY_t = pillar_map(iY);
                plint iZ_t = pillar_map(iZ);

                grid[iX_t][iY_t][iZ_t].specifyStatisticsStatus(status);
              #else
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
              #endif
            }
        }
    }
}
#endif

#ifdef PILLAR_MEM
    
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideRevertAndBoundSwapStream(Box3D bound, Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX_t = pillar_map_iX(iX, iY, iZ);
                plint iY_t = pillar_map(iY);
                plint iZ_t = pillar_map(iZ);

                grid[iX_t][iY_t][iZ_t].collide(this->getInternalStatistics());
                grid[iX_t][iY_t][iZ_t].revert();

                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if ( nextX>=bound.x0 && nextX<=bound.x1 &&
                         nextY>=bound.y0 && nextY<=bound.y1 &&
                         nextZ>=bound.z0 && nextZ<=bound.z1 )
                    {
                        #ifdef CUBE_MAP
                        plint nextX_t = pillar_map_iX(nextX, nextY, nextZ);
                        #else
                        plint nextX_t = iX_t + Descriptor<T>::c[iPop][0];
                        #endif

                        plint nextY_t = pillar_map(nextY);
                        plint nextZ_t = pillar_map(nextZ);

                        std::swap(grid[iX_t][iY_t][iZ_t][iPop + Descriptor<T>::q/2],
                                  grid[nextX_t][nextY_t][nextZ_t][iPop]);
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

    plint iX_t = pillar_map_iX(iX, iY, iZ);
    plint iY_t = pillar_map(iY);
    plint iZ_t = pillar_map(iZ);

    grid[iX_t][iY_t][iZ_t].collide(this->getInternalStatistics());
    grid[iX_t][iY_t][iZ_t].revert();

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            #ifdef CUBE_MAP
            plint nextX_t = pillar_map_iX(nextX, nextY, nextZ);
            #else
            plint nextX_t = iX_t + Descriptor<T>::c[iPop][0];
            #endif

            plint nextY_t = pillar_map(nextY);
            plint nextZ_t = pillar_map(nextZ);

            std::swap(grid[iX_t][iY_t][iZ_t][iPop + Descriptor<T>::q/2],
                      grid[nextX_t][nextY_t][nextZ_t][iPop]);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::boundSwapStream(Box3D bound, plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    plint iX_t = pillar_map_iX(iX, iY, iZ);
    plint iY_t = pillar_map(iY);
    plint iZ_t = pillar_map(iZ);

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            #ifdef CUBE_MAP
            plint nextX_t = pillar_map_iX(nextX, nextY, nextZ);
            #else
            plint nextX_t = iX_t + Descriptor<T>::c[iPop][0];
            #endif

            plint nextY_t = pillar_map(nextY);
            plint nextZ_t = pillar_map(nextZ);

            std::swap(grid[iX_t][iY_t][iZ_t][iPop + Descriptor<T>::q/2],
                      grid[nextX_t][nextY_t][nextZ_t][iPop]);
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2_2nd_CollideAndStream(Box3D domain, plint iX, plint iY, plint iZ){
  if (iY == domain.y0 || iZ == domain.z0){
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
  }
  else{
    plint iX_t = pillar_map_iX(iX, iY, iZ);
    plint iY_t = pillar_map(iY);
    plint iZ_t = pillar_map(iZ);
    grid[iX_t][iY_t][iZ_t].collide(this->getInternalStatistics());

    latticeTemplates<T,Descriptor>::swapAndStream3D_pillar_mem(grid, iX, iY, iZ);
  }
}

    
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::swapStream(plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    plint iX_t = pillar_map_iX(iX, iY, iZ);
    plint iY_t = pillar_map(iY);
    plint iZ_t = pillar_map(iZ);

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];

        #ifdef CUBE_MAP
        plint nextX_t = pillar_map_iX(nextX, nextY, nextZ);
        #else
        plint nextX_t = iX_t + Descriptor<T>::c[iPop][0];
        #endif

        plint nextY_t = pillar_map(nextY);
        plint nextZ_t = pillar_map(nextZ);

        std::swap(grid[iX_t][iY_t][iZ_t][iPop + Descriptor<T>::q/2],
                  grid[nextX_t][nextY_t][nextZ_t][iPop]);
    }
}

#else 
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
void BlockLattice3D<T,Descriptor>::step2_2nd_CollideAndStream(Box3D domain, plint iX, plint iY, plint iZ){
  if (iY == domain.y0 || iZ == domain.z0){
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
  }
  else{
    grid[iX][iY][iZ].collide(this->getInternalStatistics());
    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY, iZ);
  }
}

#endif

#if defined(STEP2_OMP)
  #if defined(PILLAR_MEM)
  // plint ykPanel_len;
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::allocateAndInitialize() {
    this->getInternalStatistics().subscribeAverage(); // Subscribe average rho-bar
    this->getInternalStatistics().subscribeAverage(); // Subscribe average uSqr
    this->getInternalStatistics().subscribeMax();     // Subscribe max uSqr

    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    rawData = new Cell<T,Descriptor> [nx*ny*nz];

    plint newNx = nx * (ny / ykTile) * (nz / ykTile);
    grid    = new Cell<T,Descriptor>** [newNx];

    plint NewNy = ykTile;
    plint NewNz = ykTile;

    pcout << "PILLAR_MEM: allocateAndInitialize:" << nx << 'x' << ny << 'x' << nz << ' ' 
          << newNx << 'x' << NewNy << 'x' << NewNz << '\n'; 

    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX = 0; iX < newNx; ++iX) {
      grid[iX] = new Cell<T,Descriptor>* [NewNy];
      for (plint iY = 0; iY < NewNy; ++iY) {
        grid[iX][iY] = rawData + NewNz * (iY + NewNy * iX);
      }
    }

    #ifdef DEBUG_PRINT
    pcout << "PILLAR_MEM: allocateAndInitialize End\n";
    #endif
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::attributeDynamics (
        plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics )
{
    plint iX_t = pillar_map_iX(iX, iY, iZ);
    plint iY_t = pillar_map(iY);
    plint iZ_t = pillar_map(iZ);

    Dynamics<T,Descriptor>* previousDynamics = &grid[iX_t][iY_t][iZ_t].getDynamics();
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX_t][iY_t][iZ_t].attributeDynamics(dynamics);
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::releaseMemory() {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();

    // pcout << "BlockLattice3D: releaseMemory:" << nx << 'x' << ny << 'x' << nz << '\n';

    #pragma omp parallel for default(shared) schedule(static)
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                plint iX_t = pillar_map_iX(iX, iY, iZ);
                plint iY_t = pillar_map(iY);
                plint iZ_t = pillar_map(iZ);
                // printf("(%ld, %ld, %ld) -> (%ld, %ld, %ld)\n", iX, iY, iZ, 
                //         iX, iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY), iZ & (YK_PANEL_LEN - 1));
                Dynamics<T,Descriptor>* dynamics = &grid[iX_t][iY_t][iZ_t].getDynamics();

                if (dynamics != backgroundDynamics) {
                    delete dynamics;
                }
            }
        }
    }

    // pcout << "BlockLattice3D::releaseMemory() delete dynamics\n";

    delete backgroundDynamics;
    delete [] rawData;

    plint newNx = nx * (ny / ykTile) * (nz / ykTile);
    for (plint iX = 0; iX < newNx; ++iX) {
        delete [] grid[iX];
    }

    delete [] grid;
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

    pcout << "PILLAR_MEM: allocateAndInitialize:" << nx << 'x' << ny << 'x' << nz << '\n'; 

    // Before: access grid[iX][iY][iZ]
    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX=0; iX<nx; ++iX) {
      grid[iX] = new Cell<T,Descriptor>* [ny];
      for (plint iY=0; iY<ny; ++iY) {
        grid[iX][iY] = rawData + nz*(iY+ny*iX);
      }
    }
    
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
void BlockLattice3D<T,Descriptor>::releaseMemory() {
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();

    // pcout << "BlockLattice3D: releaseMemory:" << nx << 'x' << ny << 'x' << nz << '\n';

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

    // pcout << "BlockLattice3D::releaseMemory() delete dynamics\n";

    delete backgroundDynamics;
    delete [] rawData;
    for (plint iX=0; iX<nx; ++iX) {
        delete [] grid[iX];
    }
    delete [] grid;
}
  #endif

#endif

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_UTIl_HH
