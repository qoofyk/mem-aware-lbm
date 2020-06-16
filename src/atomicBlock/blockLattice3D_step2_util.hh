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

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace plb {
// BlockLattice3D_step2_util /////////////////////////

// Add by Yuankun
plint thread_block;
// End add by Yuankun

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
          #ifdef PANEL_MEM
            // printf("(%ld, %ld) -> (%ld, %ld)\n", iY, iZ, 
            //   iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY), iZ & (YK_PANEL_LEN - 1));
            grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)].attributeDynamics(backgroundDynamics);
          #else
            grid[iX][iY][iZ].attributeDynamics(backgroundDynamics);
          #endif
        }
      }
    }

    // pcout << "Pass assign\n";

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
                
              #ifdef PANEL_MEM
                Cell<T,Descriptor>& cell = grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)];
                // Assign cell from rhs
                cell = rhs.grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)];
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
              #ifdef PANEL_MEM
                grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)].specifyStatisticsStatus(status);
              #else
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
              #endif
            }
        }
    }
}
#endif

#ifdef PANEL_MEM
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideRevertAndBoundSwapStream(Box3D bound, Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
                plint iZ_p = iZ & (YK_PANEL_LEN - 1);
                grid[iX][iY_p][iZ_p].collide(this->getInternalStatistics());
                grid[iX][iY_p][iZ_p].revert();

                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if ( nextX>=bound.x0 && nextX<=bound.x1 &&
                         nextY>=bound.y0 && nextY<=bound.y1 &&
                         nextZ>=bound.z0 && nextZ<=bound.z1 )
                    {
                        plint nextY_p = nextY + (nextZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
                        plint nextZ_p = nextZ & (YK_PANEL_LEN - 1);
                        std::swap(grid[iX][iY_p][iZ_p][iPop+Descriptor<T>::q/2],
                                  grid[nextX][nextY_p][nextZ_p][iPop]);
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

    plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
    plint iZ_p = iZ & (YK_PANEL_LEN - 1);
    grid[iX][iY_p][iZ_p].collide(this->getInternalStatistics());
    grid[iX][iY_p][iZ_p].revert();

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            plint nextY_p = nextY + (nextZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
            plint nextZ_p = nextZ & (YK_PANEL_LEN - 1);
            std::swap(grid[iX][iY_p][iZ_p][iPop+Descriptor<T>::q/2],
                      grid[nextX][nextY_p][nextZ_p][iPop]);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::boundSwapStream(Box3D bound, plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
    plint iZ_p = iZ & (YK_PANEL_LEN - 1);

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if ( nextX>=bound.x0 && nextX<=bound.x1 &&
             nextY>=bound.y0 && nextY<=bound.y1 &&
             nextZ>=bound.z0 && nextZ<=bound.z1 )
        {
            plint nextY_p = nextY + (nextZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
            plint nextZ_p = nextZ & (YK_PANEL_LEN - 1);
            std::swap(grid[iX][iY_p][iZ_p][iPop+Descriptor<T>::q/2],
                      grid[nextX][nextY_p][nextZ_p][iPop]);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::swapStream(plint iX, plint iY, plint iZ) {
    // Make sure domain is contained within current lattice
    // PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
    plint iZ_p = iZ & (YK_PANEL_LEN - 1);

    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];

        plint nextY_p = nextY + (nextZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
        plint nextZ_p = nextZ & (YK_PANEL_LEN - 1);
        std::swap(grid[iX][iY_p][iZ_p][iPop+Descriptor<T>::q/2],
                  grid[nextX][nextY_p][nextZ_p][iPop]);
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::step2_2nd_CollideAndStream(Box3D domain, plint iX, plint iY, plint iZ){
  if (iY == domain.y0 || iZ == domain.z0){
    collideRevertAndBoundSwapStream(domain, iX, iY, iZ);
  }
  else{
    plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
    plint iZ_p = iZ & (YK_PANEL_LEN - 1);
    grid[iX][iY_p][iZ_p].collide(this->getInternalStatistics());
    latticeTemplates<T,Descriptor>::swapAndStream3D(grid, iX, iY_p, iZ);
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
  #if defined(PANEL_MEM)
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
    grid    = new Cell<T,Descriptor>** [nx];

    // Now: use panel memory layout to compute, access grid[iX][iY + ny * (iZ / YK_PANEL_LEN)][iZ % YK_PANEL_LEN]
    // Or use Bit hack: access by grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)]
    // plint NewNy = ny * (nz / ykPanel_len);
    plint NewNy = ny * (nz >> YK_LOG_PANEL_LEN);
    plint NewNz = YK_PANEL_LEN;

    pcout << "PANEL_MEM: allocateAndInitialize:" << nx << 'x' << ny << 'x' << nz << ' ' 
          << ' ' << NewNy << 'x' << NewNz << '\n'; 

    #pragma omp parallel for default(shared) schedule(static, thread_block)
    for (plint iX = 0; iX < nx; ++iX) {
      grid[iX] = new Cell<T,Descriptor>* [NewNy];
      for (plint iY = 0; iY < NewNy; ++iY) {
        grid[iX][iY] = rawData + NewNz * (iY + NewNy * iX);
      }
    }

    // pcout << "PANEL_MEM: allocateAndInitialize End";
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::attributeDynamics (
        plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics )
{
    plint iY_p = iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY);
    plint iZ_p = iZ & (YK_PANEL_LEN - 1);

    Dynamics<T,Descriptor>* previousDynamics = &grid[iX][iY_p][iZ_p].getDynamics();
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX][iY_p][iZ_p].attributeDynamics(dynamics);
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
              #ifdef PANEL_MEM
                // printf("(%ld, %ld, %ld) -> (%ld, %ld, %ld)\n", iX, iY, iZ, 
                //         iX, iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY), iZ & (YK_PANEL_LEN - 1));
                Dynamics<T,Descriptor>* dynamics = &grid[iX][iY + (iZ >> YK_LOG_PANEL_LEN << YK_LOG_NY)][iZ & (YK_PANEL_LEN - 1)].getDynamics();
              #else
                Dynamics<T,Descriptor>* dynamics = &grid[iX][iY][iZ].getDynamics();
              #endif
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

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_STEP2_UTIl_HH
