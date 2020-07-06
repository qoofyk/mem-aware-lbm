/*
 * Description:
 *
 * First created: 
 * Last modified: 
 *
 * Author: Yuankun Fu
 * e-mail: qoofyk@gmail.com
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#define DEBUG_PRINT

#ifdef DEBUG_PRINT
#define debug cout
#else
#define debug 0 && cout
#endif

using namespace std;

// #define ykTile 32

typedef int plint;
plint*** palabos_grid;
plint* palabos_rawData;

// Using Pillar memory layout
plint*** pillar_grid;
plint* pillar_rawData;
plint newNx;
plint ykTile; /* pillarTile. The assumed size, while "actual" memory layout on Length & Width are (ykTile + 2). E.g. Given computation domain 16*8*16, palabos actual allocated memory is (16+2) * (8+2) * (16+2) --> Let pillarTile=8, thus Equivalant pillar computation domain is 32 * 8 * 8. Actual allocated memory is (32+2) * (8+2) * (8+2). */
plint NzTiles; // number of pillarTiles along Length (Z-direction)
plint NyTiles; // number of pillarTiles along Width (Y-direction)

// Using parallelepiped block Traverse
plint blockSize; // Length of parallelogram. 

// #define block_map_iX(iX, iY, iZ) (iX % ykTile + ykTile * (iZ / ykTile + (iY / ykTile) * NzTiles + (iX / ykTile) * NzTiles * NyTiles))

plint myLog2 (plint x) {
    int targetlevel = 0;
    while (x >>= 1) ++targetlevel;
    return targetlevel;
}

static inline uint32_t myAsmLog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

#if 1
plint*** allocateAndInitialize(plint Nx, plint Ny, plint Nz) {
    cout << "Palabos allocateAndInitialize\n";
    plint allocNx = Nx + 2;
    plint allocNy = Ny + 2;
    plint allocNz = Nz + 2;

    palabos_rawData = new plint [allocNx * allocNy * allocNz];
    plint* rawData = palabos_rawData;

    plint*** grid  = new plint** [allocNx];
    
    int cnt = 0;   
    for (plint iX = 0; iX < allocNx; ++iX) {
      grid[iX] = new int* [allocNy];
      for (plint iY = 0; iY < allocNy; ++iY) {
          grid[iX][iY] = rawData + allocNz * (iY + allocNy * iX);
          for (plint iZ = 0; iZ < allocNz; ++ iZ) {
              grid[iX][iY][iZ] = cnt++;
              // cout << grid[iX][iY][iZ] << '\n';
          }
      }        
    }
    cout << "Palabos allocateAndInitialize Done\n";
    return grid;
}
#endif

#if 1
plint*** pillar_memory_allocateAndInitialize(plint Nx, plint Ny, plint Nz) {
    cout << "pillar_memory_allocateAndInitialize\n";

    plint allocNx = newNx + 2;
    plint allocNy = ykTile + 2;
    plint allocNz = ykTile + 2;
    cout << newNx << ' ' << allocNy << ' ' << allocNz << '\n';

    pillar_rawData = new plint [allocNx * allocNy * allocNz];
    plint* rawData = pillar_rawData;

    plint*** grid  = new plint** [allocNx];
    
    int cnt = 0;   
    for (plint iX = 0; iX < allocNx; ++iX) {
      grid[iX] = new int* [allocNy];
      for (plint iY = 0; iY < allocNy; ++iY) {
          grid[iX][iY] = rawData + allocNz * (iY + allocNy * iX);
          for (plint iZ = 0; iZ < allocNz; ++ iZ) {
            grid[iX][iY][iZ] = cnt++;
          }
      }        
    }

    return grid;
}
#endif

#if 1
void orginal_access_pattern(plint*** grid, plint Nx, plint Ny, plint Nz) {
    cout << "orginal_access_pattern\n";

    for (plint iX = 1; iX <= Nx; ++iX) {

      for (plint iY = 1; iY <= Ny; ++iY) {
        
        for (plint iZ = 1; iZ <= Nz; ++iZ) {

          // Previous access patten grid[iX][iY][iZ]
          // debug << "Org:(" << iX << ',' << iY << ',' << iZ << ")= ";

          // New memory map access pattern
          // 1. pillar memory layout
          plint iX_t = iX + Nx * ((iZ - 1) / ykTile + ((iY - 1) / ykTile) * NzTiles);
          plint iY_t = (iY - 1) % ykTile + 1;
          plint iZ_t = (iZ - 1) % ykTile + 1;

          // debug << "Org:(" << iX << ',' << iY << ',' << iZ
          //     << ")--> New:(" << iX_t << ',' << iY_t << ',' << iZ_t << ")= ";

          int point_val = grid[iX_t][iY_t][iZ_t];
          debug << point_val << ' ';
          
        }
        debug << '\n';
    }
    debug << "----------------------------------------\n";
  }
}
#endif

// Without using %, *, /
#if 1
void pillar_access_pattern(plint*** grid, plint Nx, plint Ny, plint Nz) {
    cout << "pillar_access_pattern\n";

    plint pillarRowTiles = Nx * NzTiles;

    for (plint iX = 1; iX <= Nx; ++iX) {
      
      plint out_mem_iX = iX;
      plint iY = 1;
      for (plint tY = 0; tY < NyTiles; ++tY, out_mem_iX += pillarRowTiles) {

        for (plint iiY = 1; iiY <= ykTile; ++iiY, ++iY) {
          plint mem_iX = out_mem_iX;

          plint iZ = 1;
          for (plint tZ = 0; tZ < NzTiles; ++tZ, mem_iX += Nx) {                
                 
            for (plint iiZ = 1; iiZ <= ykTile; ++iiZ, ++iZ) {
              // Previous access patten grid[iX][iY][iZ]

              // Now 
              // cout << '(' << iX << ',' << iY << ',' << iZ << ") -> "
              //      << '[' << mem_iX << ',' << iiY << ',' << iiZ << "]; ";
              int point_val = grid[mem_iX][iiY][iiZ];
              debug << point_val << ' ';
              // for inner points' neighbors grid[mem_iX+-1][iiY+-1][iiZ+-1];
              // for outer memory points' neighbors need extra computation. 
            }
          }
          debug << '\n';
        }
      }
      debug << "----------------------------------------\n";
    }
}
#endif

void original_parallelepiped_access_pattern(plint*** grid, plint Nx, plint Ny, plint Nz) {
  cout << "original_parallelepiped_access_pattern\n";

  plint XBlockID = 1;
  for (plint outerX = 1; outerX <= Nx; outerX += blockSize, XBlockID += (Nz / blockSize + 2) * (Nz / blockSize + 1)) 
  {
    // printf("tid: %ld, outerX=%ld\n", tid, outerX);
    plint YBlockId = XBlockID;
    for (plint outerY = 1; outerY <= Ny+blockSize-1; outerY += blockSize, YBlockId += Nz / blockSize + 2) 
    {
      // printf("outerY=%ld\n", outerY);
      plint innerBlockId = YBlockId;
      for (plint outerZ = 1; outerZ <= Nz+2*(blockSize-1); outerZ += blockSize, ++innerBlockId) 
      {
        // printf("outerZ=%ld\n", outerZ);
        cout << "****** innerBlockId=" << innerBlockId << " outerZ=" << outerZ << "\n=========================================\n";
        // Inner loops.
        plint dx = 0;
        for (plint innerX = outerX; innerX <= std::min(outerX+blockSize-1, Nx);
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

          for (plint innerY=std::max(minY, 1); innerY <= std::min(maxY, Ny);
            ++innerY, ++dy)
          {
            // Z-index is shifted in negative direction at each x-increment. and at each
            //    y-increment, to ensure that only post-collision cells are accessed during
            //    the swap-operation of the streaming.
            plint minZ = outerZ-dx-dy;
            plint maxZ = minZ+blockSize-1;
            // printf("innerY=%ld, dx=%ld, dy=%ld, minY=%ld, maxY=%ld, minZ=%ld, maxZ=%ld\n",
                          // innerY, dx, dy, minY, maxY, minZ, maxZ);

            for (plint innerZ=std::max(minZ, 1); innerZ <= std::min(maxZ, Nz);
              ++innerZ)
            {
              // Previous access patten grid[innerX][innerY][innerZ];

              // New memory map access pattern
              // 1. pillar memory layout
              plint iX_t = innerX + Nx * ((innerZ - 1) / ykTile + ((innerY - 1) / ykTile) * NzTiles);
              plint iY_t = (innerY - 1) % ykTile + 1;
              plint iZ_t = (innerZ - 1) % ykTile + 1;

              debug << "Org:(" << innerX << ',' << innerY << ',' << innerZ << "); ";
                   // << ")--> New:(" << iX_t << ',' << iY_t << ',' << iZ_t << ")=";

              int point_val = grid[iX_t][iY_t][iZ_t];

              // debug << point_val << ' ';
            }
            debug << '\n';
          }
          debug << "----------------------------------------\n";
        }
      }
    }
  }
}

void pillar_parallelepiped_access_pattern(plint*** grid, plint Nx, plint Ny, plint Nz) {
  cout << "pillar_parallelepiped_access_pattern\n";

  for (plint outerX = 1; outerX <= Nx; outerX += blockSize) 
  {
    // printf("tid: %ld, outerX=%ld\n", tid, outerX);
    for (plint outerY = 1; outerY <= Ny+blockSize-1; outerY += blockSize) 
    {
      // printf("outerY=%ld\n", outerY);
      for (plint outerZ = 1; outerZ <= Nz+2*(blockSize-1); outerZ += blockSize) 
      {
        // printf("outerZ=%ld\n", outerZ);
        // Inner loops.
        plint dx = 0;
        for (plint innerX = outerX; innerX <= outerX+blockSize-1;
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

          for (plint innerY=std::max(minY, 1); innerY <= std::min(maxY, Ny);
            ++innerY, ++dy)
          {
            // Z-index is shifted in negative direction at each x-increment. and at each
            //    y-increment, to ensure that only post-collision cells are accessed during
            //    the swap-operation of the streaming.
            plint minZ = outerZ-dx-dy;
            plint maxZ = minZ+blockSize-1;
            // printf("innerY=%ld, dx=%ld, dy=%ld, minY=%ld, maxY=%ld, minZ=%ld, maxZ=%ld\n",
                          // innerY, dx, dy, minY, maxY, minZ, maxZ);

            for (plint innerZ=std::max(minZ, 1); innerZ <= std::min(maxZ, Nz);
              ++innerZ)
            {
              // New memory map access pattern
              // 1. pillar memory layout
              plint iX_t = innerX + Nx * ((innerZ - 1) / ykTile + ((innerY - 1) / ykTile) * NzTiles);
              plint iY_t = (innerY - 1) % ykTile + 1;
              plint iZ_t = (innerZ - 1) % ykTile + 1;

              debug << "Org:(" << innerX << ',' << innerY << ',' << innerZ << "); ";
                   // << ")--> New:(" << iX_t << ',' << iY_t << ',' << iZ_t << ")=";

              int point_val = grid[iX_t][iY_t][iZ_t];

              // debug << point_val << ' ';
            }
            debug << '\n';
          }
          debug << "----------------------------------------\n";
        }
      }
    }
  }
}

#if 0
void access_one_line(plint out_mem_iX) {
  plint mem_iX = out_mem_iX;

  tZ = 0;
  for (iiZ = 1; iiZ < ykTile; ++iiZ) {
    // cout << mem_iX << ',' << iiY << ',' << iiZ << "=\n";

    int point_val = grid[mem_iX][iiY][iiZ];
    debug << point_val << ' ';
  }
  mem_iX += newNx;
  
  // tz = [1, NzTiles - 2]
  for (tZ = 1; tZ < NzTiles - 1; ++tZ) {                
    for (iiZ = 0; iiZ < ykTile; ++iiZ) {
      // cout << mem_iX << ',' << iiY << ',' << iiZ << "=\n";

      int point_val = grid[mem_iX][iiY][iiZ];
      debug << point_val << ' ';
    }

    mem_iX += newNx;
  }

  tZ = NzTiles - 1;
  for (iiZ = 1; iiZ < ykTile; ++iiZ) {
      // cout << mem_iX << ',' << iiY << ',' << iiZ << "=\n";

      int point_val = grid[mem_iX][iiY][iiZ];
      debug << point_val << ' ';
  }
  mem_iX += newNx;

}

// Without using %, *, /
void pillar_access_pattern2(plint Nx, plint Ny, plint Nz) {
    cout << "pillar_access_pattern\n";
    plint iX, iY, iZ;
    plint tY, iiY, tZ, iiZ;

    for (plint iX = 1; iX <= Nx; ++iX) {
      
      plint out_mem_iX = iX;

      //------------------------------------------
      tY = 0;
      for (iiY = 1; iiY < ykTile; ++iiY) {
        access_one_line(out_mem_iX);
        debug << '\n';
      }
      out_mem_iX += newNx * NzTiles;
      //------------------------------------------


      // tY = [1, NyTiles - 2]
      for (tY = 1; tY < NyTiles - 1; ++tY) {

        for (iiY = 0; iiY < ykTile; ++iiY) {
          
          plint mem_iX = out_mem_iX;
          for (tZ = 0; tZ < NzTiles; ++tZ, mem_iX += Nx) {                
                 
            for (iiZ = 0; iiZ < ykTile; ++iiZ) {
                // cout << mem_iX << ',' << iiY << ',' << iiZ << "=\n";

                int point_val = grid[mem_iX][iiY][iiZ];
                debug << point_val << ' ';
            }
          }
          debug << '\n';
        }

        out_mem_iX += Nx * NzTiles;
      }
      debug << "----------------------------------------\n";
    }
}
#endif

void releaseMemory(plint*** grid, plint* rawData, plint Nx) {
    plint allocNx = Nx + 2;
    delete [] rawData;
    for (plint iX = 0; iX < allocNx; ++iX) {
      delete [] grid[iX];
    }
    delete [] grid;
}


int main(int argc, char** argv)
{
    cout<<"Test memory map\n";

    if (argc < 6) {
      cout << "Input wrong, please recheck! Correct: ./pillar_test Nx Ny Nz ykTile blockSize\n";
      exit(1);
    }

    int Nx = atoi(argv[1]);
    int Ny = atoi(argv[2]);
    int Nz = atoi(argv[3]);
    ykTile = atoi(argv[4]); // pillarTile;
    blockSize = atoi(argv[5]);

    if (Ny % ykTile != 0 || Nz % ykTile != 0 ) {
      cout << "Ny % ykTile != 0 || Nz % ykTile != 0\n";
      exit(1);
    }
    NyTiles = Ny / ykTile;
    NzTiles = Nz / ykTile;

    // newNx = Nx * (Ny / ykTile) * (Nz / ykTile);
    newNx = Nx * NyTiles * NzTiles;

    cout << "Computation domain: Nx x Ny x Nz: " << Nx << 'x' << Ny << 'x' << Nz
         << "\nPalabos Memory allocation domain: (Nx+2) x (Ny+2) x (Nz+2): " << Nx + 2 << 'x' << Ny + 2 << 'x' << Nz + 2 
         << "\npillarTile=" << ykTile << ", newNx=" << newNx << ", NzTiles=" << NzTiles << ',' << ", NyTiles=" << NyTiles 
         << "\nPillar Memory allocation domain: (newNx+2) x (pillarTile+2) x (pillarTile+2): " << newNx + 2 << 'x' << ykTile + 2 << 'x' << ykTile + 2 << '\n';

    cout << (plint) log2(NzTiles) << ' ' << myLog2(NzTiles) << ' ' << myAsmLog2(NzTiles) << '\n';
    
    // allocateAndInitialize(palabos_grid, palabos_rawData, Nx, Ny, Nz);
    palabos_grid = allocateAndInitialize(Nx, Ny, Nz);
    pillar_grid  = pillar_memory_allocateAndInitialize(Nx, Ny, Nz);

    orginal_access_pattern(pillar_grid, Nx, Ny, Nz);
    pillar_access_pattern(pillar_grid, Nx, Ny, Nz);

    original_parallelepiped_access_pattern(pillar_grid, Nx, Ny, Nz);

    #ifdef PILLAR_PYMD
    pillar_access_pymd_pattern(Nx, Ny, Nz);
    #endif

    releaseMemory(palabos_grid, palabos_rawData, Nx);
    releaseMemory(pillar_grid, pillar_rawData, newNx);
    
    cout<<"Done!\n";

    return 0;
}