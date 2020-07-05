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
plint*** grid;
plint* rawData;
plint newNx;
plint ykTile; // actual memory lay out on Long & Width
plint ykBlockH; // actual memory lay out on Height
plint NzTiles;
plint NyTiles;
plint allocNx;
plint allocNy;
plint allocNz;
plint blockSize;

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
void allocateAndInitialize(plint Nx, plint Ny, plint Nz) {

    rawData = new plint [allocNx * allocNy * allocNz];
    
    // newNx = allocNx * (allocNy / ykTile) * (allocNz / ykTile);
    newNx = allocNx * NyTiles * NzTiles;

    grid  = new plint** [newNx];

    plint NewNy = ykTile;
    plint NewNz = ykTile;
    
    cout << newNx << ' ' << NewNy << ' ' << NewNz << '\n';
    
    int cnt = 0;
    
    for (plint iX = 0; iX < newNx; ++iX) {
      grid[iX] = new int* [NewNy];
      for (plint iY = 0; iY < NewNy; ++iY) {
          grid[iX][iY] = rawData + NewNz * (iY + NewNy * iX);
          for (plint iZ = 0; iZ < NewNz; ++ iZ) {
              grid[iX][iY][iZ] = cnt++;
          }
      }        
    }
}
#endif

#if 0
void allocateAndInitialize(plint Nx, plint Ny, plint Nz) {

    rawData = new plint [allocNx * allocNy * allocNz];

    grid  = new plint** [newNx];
    
    int cnt = 0;   
    for (plint iX = 0; iX < allocNx; ++iX) {
      grid[iX] = new int* [allocNy];
      for (plint iY = 0; iY < allocNz; ++iY) {
          grid[iX][iY] = rawData + allocNz * (iY + allocNx * iX);
          for (plint iZ = 0; iZ < allocNz; ++ iZ) {
              grid[iX][iY][iZ] = cnt++;
          }
      }        
    }
}
#endif

#if 0
void orginal_access_pattern(plint Nx, plint Ny, plint Nz) {
    cout << "orginal_access_pattern\n";

    for (plint iX = 1; iX <= Nx; ++iX) {
      for (plint iY = 1; iY <= Ny; ++iY) {
        for (plint iZ = 1; iZ <= Nz; ++iZ) {

          // Previous access patten grid[iX][iY][iZ]

          // New memory map access pattern
          // 1. pillar memory layout
          // plint iX_t = iX + allocNx * (iZ / ykTile + (iY / ykTile) * (Nz / ykTile));
          plint iX_t = iX + allocNx * (iZ / ykTile + (iY / ykTile) * NzTiles);

          // 2. cube memory layout
          // plint iX_t = iX % ykBlockH + ykTile * (iZ / ykTile + (iY / ykTile) * (Nz / ykTile) + (iX / ykBlockH)  * (Nz / ykTile)  * (ny / ykTile));
          // plint iX_t = iX % ykBlockH + ykTile * (iZ / ykTile + (iY / ykTile) * NzTiles + (iX / ykBlockH) * NzTiles * NyTiles);
          // plint iX_t = block_map_iX(iX, iY, iZ);

          plint iY_t = iY % ykTile;
          plint iZ_t = iZ % ykTile;

          // debug << "Org:(" << iX << ',' << iY << ',' << iZ 
          //      << ")--> New:(" << iX_t << ',' << iY_t << ',' << iZ_t << ")=";

          int point_val = grid[iX_t][iY_t][iZ_t];

          debug << point_val << ' ';
        }
        debug << '\n';
    }
    debug << "----------------------------------------\n";
  }
}
#endif

void orginal_wedge_access_pattern(plint Nx, plint Ny, plint Nz) {
  cout << "orginal_wedge_access_pattern\n";
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
              // plint iX_t = innerX + allocNx * (innerZ / ykTile + (innerY / ykTile) * (Nz / ykTile));
              plint iX_t = innerX + allocNx * (innerZ / ykTile + (innerY / ykTile) * NzTiles);
              plint iY_t = innerY % ykTile;
              plint iZ_t = innerZ % ykTile;

              debug << "Org:(" << innerX << ',' << innerY << ',' << innerZ << ") ";
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

// Without using %, *, /
void pillar_access_pattern(plint Nx, plint Ny, plint Nz) {
    cout << "pillar_access_pattern\n";

    plint pillarRowTiles = allocNx * NzTiles;
   
    for (plint iX = 1; iX <= Nx; ++iX) {
      
      plint out_mem_iX = iX;
      plint iY = 0;
      for (plint tY = 0; tY < NyTiles; ++tY, out_mem_iX += pillarRowTiles) {

        for (plint iiY = 0; iiY < ykTile; ++iiY, ++iY) {
          if (iY == 0 || iY == Ny + 1) continue;

          plint mem_iX = out_mem_iX;

          plint iZ = 0;
          for (plint tZ = 0; tZ < NzTiles; ++tZ, mem_iX += allocNx) {                
                 
            for (plint iiZ = 0; iiZ < ykTile; ++iiZ, ++iZ) {
                if (iZ == 0 || iZ == Nz + 1) continue;
                // Previous access patten grid[iX][iY][iZ]

                // Now 
                // cout << mem_iX << ',' << iiY << ',' << iiZ << "=\n";
                int point_val = grid[mem_iX][iiY][iiZ];
                debug << point_val << ' ';
            }
          }
          debug << '\n';
        }
      }
      debug << "----------------------------------------\n";
    }
}

void pillar_wedge_access_pattern(plint Nx, plint Ny, plint Nz) {
  cout << "orginal_access_pymd_pattern\n";

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
              // Previous access patten grid[innerX][innerY][innerZ];
              // New memory map access pattern
              // 1. pillar memory layout
              // plint iX_t = innerX + allocNx * (innerZ / ykTile + (innerY / ykTile) * (Nz / ykTile));
              plint iX_t = innerX + allocNx * (innerZ / ykTile + (innerY / ykTile) * NzTiles);
              plint iY_t = innerY % ykTile;
              plint iZ_t = innerZ % ykTile;

              // debug << "Org:(" << innerX << ',' << innerY << ',' << innerZ 
              //      << ")--> New:(" << iX_t << ',' << iY_t << ',' << iZ_t << ")=";

              int point_val = grid[iX_t][iY_t][iZ_t];

              debug << point_val << ' ';
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

void releaseMemory(plint Nx, plint ny, plint Nz) {

    delete [] rawData;
    for (plint iX = 0; iX < newNx; ++iX) {
        delete [] grid[iX];
    }
    delete [] grid;
}


int main(int argc, char** argv)
{
    cout<<"Test memory map\n";

    int Nx = atoi(argv[1]);
    int Ny = atoi(argv[2]);
    int Nz = atoi(argv[3]);
    ykTile = atoi(argv[4]);
    blockSize = atoi(argv[5]);

    allocNx = Nx + 2;
    allocNy = Ny + 2;
    allocNz = Nz + 2;

    if (allocNy % ykTile || allocNz % ykTile ) {
      cout << "(Ny + 2) % ykTile != 0 || (Nz + 2) % ykTile != 0\n";
      exit(1);
    }
    NyTiles = allocNy / ykTile;
    NzTiles = allocNz / ykTile;

    cout << "Computation domain: Nx x Ny x Nz: " << Nx << 'x' << Ny << 'x' << Nz
         << "\nMemory allocation domain: (Nx+2) x (Ny+2) x (Nz+2): " << allocNx << 'x' << allocNy << 'x' << allocNz 
         << "\nNzTiles=" << NzTiles << ',' << "NyTiles=" << NyTiles << '\n';

    cout << (plint) log2(NzTiles) << ' ' << myLog2(NzTiles) << ' ' << myAsmLog2(NzTiles) << '\n';
    
    allocateAndInitialize(Nx, Ny, Nz);
    
    #ifdef ORG    
    orginal_access_pattern(Nx, Ny, Nz);
    #endif

    // #ifdef ORG_PYMD
    orginal_wedge_access_pattern(Nx, Ny, Nz);
    // #endif

    #ifdef PILLAR
    pillar_access_pattern(Nx, Ny, Nz);
    #endif

    #ifdef PILLAR_PYMD
    pillar_access_pymd_pattern(Nx, Ny, Nz);
    #endif

    releaseMemory(Nx, Ny, Nz);
    
    cout<<"Done!\n";

    return 0;
}