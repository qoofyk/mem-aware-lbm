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
 * 3D specialization of latticeTemplates functions.
 */

#ifndef LATTICE_TEMPLATES_3D_H
#define LATTICE_TEMPLATES_3D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

// Add by Yuankun
#include "atomicBlock/blockLattice3D_pillar_mem.h"
// End add by Yuankun

namespace plb {

template<typename T>
struct latticeTemplates<T, descriptors::D3Q19Descriptor> {

static void swapAndStreamCell (
      Cell<T,descriptors::D3Q19Descriptor> ***grid,
      plint iX, plint iY, plint iZ, plint nX, plint nY, plint nZ, plint iPop, T& fTmp )
{
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+9];
    grid[iX][iY][iZ][iPop+9] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
}

static void swapAndStream3D(Cell<T,descriptors::D3Q19Descriptor> ***grid,
                            plint iX, plint iY, plint iZ, plint Nx_, plint Ny_, plint Nz_)
{
    T fTmp;

    plint iX_t = cube_mem_map_iX(iX, iY, iZ, Nx_, Ny_, Nz_);
    plint iY_t = iY % ykTile;
    plint iZ_t = iZ % ykTile;

    plint iX_minus_1_t = cube_mem_map_iX(iX - 1, iY, iZ, Nx_, Ny_, Nz_);
    plint iY_minus_1_t = (iY - 1) % ykTile;
    plint iZ_minus_1_t = (iZ - 1) % ykTile;
    plint iZ_plus_1_t = (iZ + 1) % ykTile;

    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_minus_1_t, iY_t,               iZ_t,         1, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_t,         iY_minus_1_t,       iZ_t,         2, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_t,         iY_t,               iZ_minus_1_t, 3, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_minus_1_t, iY_minus_1_t,       iZ_t,         4, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_minus_1_t, (iY + 1) % ykTile,  iZ_t,         5, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_minus_1_t, iY_t,               iZ_minus_1_t, 6, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_minus_1_t, iY_t,               iZ_plus_1_t,  7, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_t,         iY_minus_1_t,       iZ_minus_1_t, 8, fTmp);
    swapAndStreamCell(grid, iX_t, iY_t, iZ_t, iX_t,         iY_minus_1_t,       iZ_plus_1_t,  9, fTmp);
}
static void swapAndStream3D(Cell<T,descriptors::D3Q19Descriptor> ***grid,
                            plint iX, plint iY, plint iZ)
{
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ,   4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ,   5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ+1, 7, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ-1, 8, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ+1, 9, fTmp);
}
};

template<typename T>
struct latticeTemplates<T, descriptors::ForcedD3Q19Descriptor> {

static void swapAndStreamCell (
      Cell<T,descriptors::ForcedD3Q19Descriptor> ***grid,
      plint iX, plint iY, plint iZ, plint nX, plint nY, plint nZ, plint iPop, T& fTmp )
{
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+9];
    grid[iX][iY][iZ][iPop+9] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
}

static void swapAndStream3D(Cell<T,descriptors::ForcedD3Q19Descriptor> ***grid,
                            plint iX, plint iY, plint iZ)
{
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ,   4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ,   5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ+1, 7, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ-1, 8, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ+1, 9, fTmp);
}

};

template<typename T>
struct latticeTemplates<T, descriptors::D3Q15Descriptor> {

static void swapAndStreamCell (
      Cell<T,descriptors::D3Q15Descriptor> ***grid,
      plint iX, plint iY, plint iZ, plint nX, plint nY, plint nZ, plint iPop, T& fTmp )
{
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+7];
    grid[iX][iY][iZ][iPop+7] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
}

static void swapAndStream3D(Cell<T,descriptors::D3Q15Descriptor> ***grid,
                            plint iX, plint iY, plint iZ)
{
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);

    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 7, fTmp);
}

};

template<typename T>
struct latticeTemplates<T, descriptors::ForcedD3Q15Descriptor> {

static void swapAndStreamCell (
      Cell<T,descriptors::ForcedD3Q15Descriptor> ***grid,
      plint iX, plint iY, plint iZ, plint nX, plint nY, plint nZ, plint iPop, T& fTmp )
{
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+7];
    grid[iX][iY][iZ][iPop+7] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
}

static void swapAndStream3D(Cell<T,descriptors::ForcedD3Q15Descriptor> ***grid,
                            plint iX, plint iY, plint iZ)
{
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);

    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 7, fTmp);
}

};

}  // namespace plb

#endif

