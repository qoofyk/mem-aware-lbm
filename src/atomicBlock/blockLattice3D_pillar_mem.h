#ifndef BLOCK_LATTICE_3D_PILLAR_MEM_H
#define BLOCK_LATTICE_3D_PILLAR_MEM_H
namespace plb {

// Add by Yuankun
// #define BIT_HACK
// #define YK_TILE 32 // 2^5=5
// #define LOG_YK_TILE 5 // 2^5=5
// #define LOG_NZ_TILES 4// 2^9 512
// #define YK_PANEL_LEN_MASK (~0 << YK_LOG_PANEL_LEN) //0xffffffe0
// End add by Yuankun

extern plint ykTile;
extern plint NzTiles;
extern plint NyTiles;
extern plint Nx, newNx;

// #define cube_mem_map_iX(iX, iY, iZ) ((iX) % ykTile + ykTile * ((iZ) / ykTile + ((iY) / ykTile) * NzTiles + ((iX) / ykTile) * NzTiles * NyTiles))

extern inline plint pillar_map_iX (plint iX, plint iY, plint iZ);
// #define pillar_map_iX(iX, iY, iZ) ((iX) + memNx * ( ((iZ) >> LOG_YK_TILE) + ( ((iY) >> LOG_YK_TILE) << log_NzTiles) ))
// #define pillar_map_iX(iX, iY, iZ) ((iX) + memNx * ( ((iZ) >> LOG_YK_TILE) + ( ((iY) >> LOG_YK_TILE) << LOG_NZ_TILES) ))

extern inline plint pillar_map (plint iZ);
// #define pillar_map(iZ) ((iZ) & (YK_TILE - 1))

}

#endif  // BLOCK_LATTICE_3D_PILLAR_MEM_H